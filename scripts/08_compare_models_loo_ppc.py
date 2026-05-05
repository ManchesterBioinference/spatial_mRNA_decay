#!/usr/bin/env python
"""
Compare degradation models using LOO-CV and posterior predictive checks.

Compares Delayed (Gaussian Random Walk age-dependent) and Biphasic (mechanistic
poly-A) models against the Null (constant-D) and Spatial (AP-binned constant-D)
baselines using Leave-One-Out cross-validation (LOO-CV) via ArviZ.

Outputs:
    - ELPD comparison table (CSV)
    - LOO comparison bar chart
    - Posterior predictive check figures for each model
    - PPC overlay figure comparing all models

Usage:
    conda run -n spatial-mrna-decay python scripts/08_compare_models_loo_ppc.py \\
        --delayed-dir results_1200/stripe3/e8_9um_age \\
        --biphasic-dir results_1200/stripe3/e8_9um_biphasic \\
        --null-dir results_1200/stripe3/e8_9um_null \\
        --spatial-dir results_1200/stripe3/e8_9um \\
        --transcription data/processed_transcription_data/transcription_traces_no_ids_stripe3_1200.csv \\
        --mrna results_1200/data/processed_mRNA_data_stripe3/e8_9um_sass_formodel.csv \\
        --output-dir results_1200/comparison/stripe3/e8_9um
"""
import argparse
import json
import os
import re
from dataclasses import dataclass
from typing import Literal
import numpy as np
import pandas as pd
import arviz as az
import matplotlib.pyplot as plt


MODEL_A_NAME = "Delayed"
MODEL_B_NAME = "Biphasic"


@dataclass
class ModelSpec:
    """Specification for a model to compare."""
    name: str
    model_type: Literal["age", "null_constant", "spatial_ap"]
    results_dir: str


def load_matrix(path: str) -> np.ndarray:
    return pd.read_csv(path, header=None).values


def load_vector(path: str) -> np.ndarray:
    return pd.read_csv(path, header=None).values.flatten()


def get_d_columns(df: pd.DataFrame) -> list[str]:
    pattern = re.compile(r"^D\[(\d+)\]$")
    indexed = []
    for col in df.columns:
        match = pattern.match(col)
        if match:
            indexed.append((int(match.group(1)), col))
    indexed.sort(key=lambda item: item[0])
    return [col for _, col in indexed]


def load_chain_df(chain_csv: str) -> pd.DataFrame:
    """Load chain CSV into DataFrame."""
    return pd.read_csv(chain_csv)


def parse_chain_params(
    df: pd.DataFrame,
    model_type: str,
    n_time: int,
    n_ap_bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse chain parameters for different model types.
    
    Args:
        df: Chain DataFrame from CSV
        model_type: "age", "null_constant", or "spatial_ap"
        n_time: Number of timepoints (for age model)
        n_ap_bins: Number of AP bins (for spatial model)
    
    Returns:
        D_samples: Degradation parameters with shape:
            - age: [n_samples, n_time]
            - null_constant: [n_samples] (scalar per sample)
            - spatial_ap: [n_samples, n_ap_bins]
        gamma: [n_samples]
        sigma: [n_samples]
    """
    if "gamma" not in df.columns or "sigma" not in df.columns:
        raise ValueError(f"Expected columns 'gamma' and 'sigma' in chain CSV.")
    
    gamma = df["gamma"].to_numpy(dtype=np.float64)
    sigma = df["sigma"].to_numpy(dtype=np.float64)
    
    d_cols = get_d_columns(df)
    if len(d_cols) == 0:
        raise ValueError(f"No D[i] columns found in chain CSV.")
    
    if model_type == "age":
        if len(d_cols) != n_time:
            raise ValueError(
                f"D columns ({len(d_cols)}) do not match transcription timepoints ({n_time}) for age model."
            )
        D_samples = df[d_cols].to_numpy(dtype=np.float64)
        
    elif model_type == "null_constant":
        if len(d_cols) != 1:
            raise ValueError(
                f"Null constant model expects exactly 1 D column, found {len(d_cols)}."
            )
        D_samples = df[d_cols[0]].to_numpy(dtype=np.float64)
        
    elif model_type == "spatial_ap":
        if len(d_cols) != n_ap_bins:
            raise ValueError(
                f"D columns ({len(d_cols)}) do not match n_ap_bins ({n_ap_bins}) for spatial model."
            )
        D_samples = df[d_cols].to_numpy(dtype=np.float64)
        
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    return D_samples, gamma, sigma


def build_expected_mrna_for_model(
    model_type: str,
    D_samples: np.ndarray,
    gamma: np.ndarray,
    transcription: np.ndarray,
    dt: float,
    n_ap_bins: int,
    n_dv_bins: int,
) -> np.ndarray:
    """
    Build expected mRNA for different model types.
    
    Args:
        model_type: "age", "null_constant", or "spatial_ap"
        D_samples: Degradation parameters
        gamma: Scaling factors [n_samples]
        transcription: Transcription data [n_obs, n_time]
        dt: Time step in minutes
        n_ap_bins: Number of AP bins
        n_dv_bins: Number of DV bins
    
    Returns:
        mu: Expected mRNA [n_samples, n_obs]
    """
    n_obs, n_time = transcription.shape
    
    # Common precompute: trapezoid weights
    weights = np.ones(n_time, dtype=np.float64)
    weights[0] = 0.5
    weights[-1] = 0.5
    weights_scaled = weights * dt
    weighted_transcription = transcription * weights_scaled[None, :]
    
    if model_type == "age":
        # Age-dependent: D_samples shape [n_samples, n_time]
        # Keep existing implementation
        cumulative_hazard = np.cumsum(D_samples, axis=1) * dt
        survival = np.exp(-cumulative_hazard)
        survival_reversed = survival[:, ::-1]
        integrals = survival_reversed @ weighted_transcription.T
        mu = gamma[:, None] * integrals
        
    elif model_type == "null_constant":
        # Null constant: D_samples is scalar per sample [n_samples]
        # Expand to time axis and compute kernel directly
        # kernel(t) = exp(-D0 * (T - t))
        n_samples = len(D_samples)
        time_from_end = (n_time - 1 - np.arange(n_time)) * dt
        
        # Broadcast: [n_samples, 1] * [1, n_time] = [n_samples, n_time]
        kernel = np.exp(-D_samples[:, None] * time_from_end[None, :])
        
        # Integrate: [n_samples, n_time] @ [n_obs, n_time].T = [n_samples, n_obs]
        integrals = kernel @ weighted_transcription.T
        mu = gamma[:, None] * integrals
        
    elif model_type == "spatial_ap":
        # Spatial AP: D_samples shape [n_samples, n_ap_bins]
        # Each observation k maps to AP bin b = k // n_dv_bins
        n_samples = D_samples.shape[0]
        mu = np.zeros((n_samples, n_obs), dtype=np.float64)
        
        time_from_end = (n_time - 1 - np.arange(n_time)) * dt
        
        for b in range(n_ap_bins):
            # Indices for this AP bin
            idx_b = [k for k in range(n_obs) if k // n_dv_bins == b]
            if len(idx_b) == 0:
                continue
            
            # Kernel for this bin: exp(-D[b] * time_from_end)
            # [n_samples, 1] * [1, n_time] = [n_samples, n_time]
            kernel_b = np.exp(-D_samples[:, b][:, None] * time_from_end[None, :])
            
            # Integrate for traces in this bin
            # [n_samples, n_time] @ [n_traces_in_bin, n_time].T = [n_samples, n_traces_in_bin]
            integrals_b = kernel_b @ weighted_transcription[idx_b].T
            
            # Apply gamma scaling
            mu[:, idx_b] = gamma[:, None] * integrals_b
    
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    return mu


def infer_chain_draw_shape(n_rows: int, n_chains: int, n_draws: int | None) -> tuple[int, int]:
    if n_draws is None:
        if n_rows % n_chains != 0:
            raise ValueError(f"Rows ({n_rows}) not divisible by n_chains ({n_chains}).")
        return n_chains, n_rows // n_chains
    if n_chains * n_draws != n_rows:
        raise ValueError(
            f"n_chains * n_draws ({n_chains} * {n_draws} = {n_chains * n_draws}) does not match rows ({n_rows})."
        )
    return n_chains, n_draws


def build_log_likelihood(mu: np.ndarray, sigma: np.ndarray, observed: np.ndarray) -> np.ndarray:
    sigma2 = np.square(sigma)[:, None]
    residual2 = np.square(observed[None, :] - mu)
    return -0.5 * np.log(2.0 * np.pi * sigma2) - 0.5 * residual2 / sigma2


def get_pareto_k(loo_result: az.ELPDData) -> np.ndarray:
    """Extract Pareto k-hat values from a LOO result as a plain numpy array."""
    k = loo_result.pareto_k
    if hasattr(k, "values"):
        return k.values
    return np.asarray(k)


def build_idata_from_chain(
    chain_csv: str,
    model_name: str,
    model_type: str,
    transcription: np.ndarray,
    observed: np.ndarray,
    dt: float,
    n_chains: int,
    n_draws: int | None,
    n_ap_bins: int,
    n_dv_bins: int,
) -> tuple[az.InferenceData, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build InferenceData from chain CSV for any model type.
    
    Args:
        chain_csv: Path to chain CSV file
        model_name: Display name for model
        model_type: "age", "null_constant", or "spatial_ap"
        transcription: Transcription matrix [n_obs, n_time]
        observed: Observed mRNA [n_obs]
        dt: Time step in minutes
        n_chains: Number of MCMC chains
        n_draws: Draws per chain (None to infer)
        n_ap_bins: Number of AP bins
        n_dv_bins: Number of DV bins
    
    Returns:
        idata: ArviZ InferenceData
        mu: Expected mRNA [n_samples, n_obs]
        gamma: Gamma samples [n_samples]
        sigma: Sigma samples [n_samples]
    """
    n_time = transcription.shape[1]
    
    # Load and parse chain
    df = load_chain_df(chain_csv)
    D_samples, gamma, sigma = parse_chain_params(df, model_type, n_time=n_time, n_ap_bins=n_ap_bins)
    
    # Build expected mRNA
    mu = build_expected_mrna_for_model(
        model_type, D_samples, gamma, transcription, dt, n_ap_bins=n_ap_bins, n_dv_bins=n_dv_bins
    )
    
    # Build log likelihood
    log_lik = build_log_likelihood(mu, sigma, observed)
    
    # Reshape for ArviZ
    chains, draws = infer_chain_draw_shape(len(df), n_chains, n_draws)
    
    # Prepare posterior dict based on model type
    if model_type == "age":
        posterior = {
            "D": D_samples.reshape(chains, draws, n_time),
            "gamma": gamma.reshape(chains, draws),
            "sigma": sigma.reshape(chains, draws),
        }
    elif model_type == "null_constant":
        # D0 should be 2D [chains, draws] with no extra dimension
        posterior = {
            "D0": D_samples.reshape(chains, draws),
            "gamma": gamma.reshape(chains, draws),
            "sigma": sigma.reshape(chains, draws),
        }
    elif model_type == "spatial_ap":
        posterior = {
            "D": D_samples.reshape(chains, draws, n_ap_bins),
            "gamma": gamma.reshape(chains, draws),
            "sigma": sigma.reshape(chains, draws),
        }
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    log_lik_reshaped = log_lik.reshape(chains, draws, log_lik.shape[1])
    
    idata = az.from_dict(
        posterior=posterior,
        log_likelihood={"m_obs": log_lik_reshaped},
        observed_data={"m_obs": observed},
    )
    idata.posterior.attrs["model_name"] = model_name
    
    return idata, mu, gamma, sigma


# ---------------------------------------------------------------------------
# LOO refit wrapper  (used by az.reloo for high-Pareto-k observations)
# ---------------------------------------------------------------------------

class MRNADecaySamplingWrapper(az.SamplingWrapper):
    """
    SamplingWrapper for ``az.reloo`` on spatial mRNA decay models.

    Implements exact LOO via random-walk Metropolis-Hastings (RWMH) for the
    small number of observations whose PSIS Pareto k-hat exceeds the threshold
    (typically 1-2 per model).  All sampling uses pure NumPy — no additional
    packages required.  Chains are seeded from the full-data posterior mean so
    burnin converges quickly.

    Parameters
    ----------
    idata : az.InferenceData
        Full-data InferenceData.  Used both to seed RWMH and as the base
        for the initial PSIS-LOO call inside ``az.reloo``.
    model_type : str
        One of ``"age"``, ``"null_constant"``, ``"spatial_ap"``.
    transcription : ndarray  (n_obs, n_time)
    observed : ndarray  (n_obs,)
    dt : float — time step in minutes
    n_ap_bins, n_dv_bins : int
    n_burnin : int — discarded steps per chain
    n_production : int or None
        Production steps per chain.  ``None`` → use original draws-per-chain,
        with one RWMH chain per original MCMC chain.
    seed : int
    """

    def __init__(
        self,
        idata: az.InferenceData,
        model_type: str,
        transcription: np.ndarray,
        observed: np.ndarray,
        dt: float,
        n_ap_bins: int,
        n_dv_bins: int,
        n_burnin: int = 500,
        n_production: int | None = None,
        seed: int = 42,
    ) -> None:
        super().__init__(model=None, idata_orig=idata)
        self.model_type    = model_type
        self.transcription = transcription
        self.observed      = observed
        self.dt            = dt
        self.n_ap_bins     = n_ap_bins
        self.n_dv_bins     = n_dv_bins
        self.n_burnin      = n_burnin
        self.seed          = seed
        self._holdout_idx: int | None = None

        # ── parameter layout by model type ──────────────────────────────────
        n_time = transcription.shape[1]
        if model_type == "age":
            self._d_dim = n_time
            d_flat = idata.posterior["D"].values.reshape(-1, n_time)
        elif model_type == "null_constant":
            self._d_dim = 1
            d_flat = idata.posterior["D0"].values.reshape(-1, 1)
        elif model_type == "spatial_ap":
            self._d_dim = n_ap_bins
            d_flat = idata.posterior["D"].values.reshape(-1, n_ap_bins)
        else:
            raise ValueError(f"Unknown model_type: {model_type}")

        gamma_flat = idata.posterior["gamma"].values.reshape(-1)
        sigma_flat = idata.posterior["sigma"].values.reshape(-1)
        self._n_params = self._d_dim + 2

        # Full-posterior statistics for seeding the RWMH chains
        self._pmean = np.concatenate([
            d_flat.mean(0), [gamma_flat.mean()], [sigma_flat.mean()]
        ])
        self._pstd = np.maximum(
            np.concatenate([d_flat.std(0), [gamma_flat.std()], [sigma_flat.std()]]),
            1e-8,
        )

        # Number of independent RWMH chains (matches original sampling)
        if n_production is None:
            n_chains_orig = idata.posterior.dims["chain"]
            n_draws_orig  = idata.posterior.dims["draw"]
            self._n_chains_reloo   = n_chains_orig
            self._n_production     = n_draws_orig
        else:
            self._n_chains_reloo   = 1
            self._n_production     = n_production

    # ── abstract method implementations ────────────────────────────────────

    def sel_observations(self, idx) -> tuple:
        """Store hold-out index; return (data_without_i, excluded_obs_index)."""
        holdout = int(np.asarray(idx).ravel()[0])
        self._holdout_idx = holdout
        # new_obs: dict telling sample() which obs to hold out
        # excluded_obs: int index passed to log_likelihood__i
        return {"holdout_idx": holdout}, holdout

    def _log_prob(self, theta: np.ndarray, holdout_idx: int) -> float:
        """Log-posterior with observation ``holdout_idx`` excluded."""
        d_vals = theta[:self._d_dim]
        gamma  = theta[self._d_dim]
        sigma  = theta[self._d_dim + 1]

        if np.any(d_vals <= 0) or gamma <= 0 or sigma <= 0:
            return -np.inf

        # Weakly informative half-normal prior (scale = 5 × full posterior std)
        prior_scale = 5.0 * self._pstd
        log_prior = float(np.sum(-0.5 * (theta / prior_scale) ** 2))

        # Single-sample forward pass
        if self.model_type == "null_constant":
            D_batch = np.array([d_vals[0]])   # shape (1,) scalar
        else:
            D_batch = d_vals[np.newaxis, :]   # shape (1, n_d)

        mu = build_expected_mrna_for_model(
            self.model_type,
            D_batch,
            np.array([gamma]),
            self.transcription,
            self.dt,
            self.n_ap_bins,
            self.n_dv_bins,
        )  # (1, n_obs)
        log_lik = build_log_likelihood(mu, np.array([sigma]), self.observed)  # (1, n_obs)

        mask = np.ones(len(self.observed), dtype=bool)
        mask[holdout_idx] = False
        return float(log_lik[0, mask].sum()) + log_prior

    @staticmethod
    def _rwmh_chain(
        log_prob,
        start: np.ndarray,
        step: np.ndarray,
        n_burnin: int,
        n_production: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """
        Single random-walk Metropolis-Hastings chain.

        Uses adaptive step-size during burnin (targets ~35 % acceptance).
        Returns production samples of shape ``(n_production, n_params)``.
        """
        n_params = len(start)
        current   = start.copy()
        current_lp = log_prob(current)
        step = step.copy()

        # Burnin with step-size adaptation every 100 steps
        accepted = 0
        for i in range(n_burnin):
            proposal = np.abs(current + step * rng.standard_normal(n_params)) + 1e-9
            prop_lp  = log_prob(proposal)
            if np.log(max(rng.random(), 1e-300)) < prop_lp - current_lp:
                current    = proposal
                current_lp = prop_lp
                accepted  += 1
            if (i + 1) % 100 == 0:
                rate = accepted / 100
                step *= 1.1 if rate > 0.44 else (0.9 if rate < 0.23 else 1.0)
                accepted = 0

        # Production
        chain = np.empty((n_production, n_params))
        for i in range(n_production):
            proposal = np.abs(current + step * rng.standard_normal(n_params)) + 1e-9
            prop_lp  = log_prob(proposal)
            if np.log(max(rng.random(), 1e-300)) < prop_lp - current_lp:
                current    = proposal
                current_lp = prop_lp
            chain[i] = current
        return chain

    def sample(self, modified_observed_data: dict) -> np.ndarray:
        """
        Refit with independent RWMH chains holding out
        ``modified_observed_data["holdout_idx"]``.

        Returns array of shape ``(n_chains * n_production, n_params)``
        with all chains concatenated (flat chain).
        """
        holdout_idx = int(modified_observed_data["holdout_idx"])
        rng = np.random.default_rng(self.seed + holdout_idx)

        # Initial step size ≈ 10 % of posterior std
        step = 0.1 * self._pstd

        flat_chains = []
        for c in range(self._n_chains_reloo):
            start = np.abs(
                self._pmean + 1e-4 * self._pstd * rng.standard_normal(self._n_params)
            ) + 1e-9
            print(
                f"    [reloo obs {holdout_idx}] chain {c+1}/{self._n_chains_reloo}  "
                f"burnin={self.n_burnin}  production={self._n_production}"
            )
            chain = self._rwmh_chain(
                lambda theta, hid=holdout_idx: self._log_prob(theta, hid),
                start, step, self.n_burnin, self._n_production, rng,
            )
            flat_chains.append(chain)

        return np.concatenate(flat_chains, axis=0)  # (n_chains*n_production, n_params)

    def get_inference_data(self, flat_chain: np.ndarray) -> az.InferenceData:
        """
        Build InferenceData from the RWMH flat chain.

        Log-likelihood is computed for **all** observations (including the
        held-out one) so that ``log_likelihood__i`` can extract the exact LOO
        log p(y_i | θ_{-i}).
        """
        D_all     = flat_chain[:, :self._d_dim]
        gamma_all = flat_chain[:, self._d_dim]
        sigma_all = flat_chain[:, self._d_dim + 1]

        D_samples = D_all[:, 0] if self.model_type == "null_constant" else D_all

        mu      = build_expected_mrna_for_model(
            self.model_type, D_samples, gamma_all,
            self.transcription, self.dt, self.n_ap_bins, self.n_dv_bins,
        )  # (n_samples, n_obs)
        log_lik = build_log_likelihood(mu, sigma_all, self.observed)  # (n_samples, n_obs)

        # ArviZ expects (chains, draws, n_obs) — treat all RWMH samples as 1 chain
        ll_3d = log_lik[np.newaxis, :, :]

        if self.model_type == "age":
            posterior = {
                "D":     D_samples[np.newaxis, :, :],
                "gamma": gamma_all[np.newaxis, :],
                "sigma": sigma_all[np.newaxis, :],
            }
        elif self.model_type == "null_constant":
            posterior = {
                "D0":    D_samples[np.newaxis, :],
                "gamma": gamma_all[np.newaxis, :],
                "sigma": sigma_all[np.newaxis, :],
            }
        else:  # spatial_ap
            posterior = {
                "D":     D_samples[np.newaxis, :, :],
                "gamma": gamma_all[np.newaxis, :],
                "sigma": sigma_all[np.newaxis, :],
            }

        return az.from_dict(
            posterior=posterior,
            log_likelihood={"m_obs": ll_3d},
            observed_data={"m_obs": self.observed},
        )

    def log_likelihood__i(
        self,
        excluded_obs: int,
        idata__i: az.InferenceData,
    ):
        """
        Return log p(y_i | theta_{-i}) as a DataArray of shape (chain, draw).

        ``excluded_obs`` is the integer observation index returned by
        ``sel_observations``; ``idata__i`` is the LOO InferenceData from
        ``get_inference_data`` whose log_likelihood includes all n_obs obs.
        """
        # log_likelihood["m_obs"] has dims (chain, draw, m_obs_dim_0)
        return idata__i.log_likelihood["m_obs"].isel(m_obs_dim_0=excluded_obs)


# ---------------------------------------------------------------------------
# Convenience wrapper: run LOO and escalate to reloo if k > threshold
# ---------------------------------------------------------------------------

def maybe_reloo(
    idata: az.InferenceData,
    wrapper: MRNADecaySamplingWrapper,
    k_thresh: float = 0.7,
    scale: str = "deviance",
    verbose: bool = True,
) -> tuple[az.ELPDData, list[int]]:
    """
    Compute PSIS-LOO; run exact LOO (``az.reloo``) for any observation with
    Pareto k > k_thresh.

    Returns
    -------
    loo_result : ELPDData — corrected where reloo was applied
    reloo_indices : list[int] — observation indices that were reloo-corrected
    """
    loo_result = az.loo(idata, pointwise=True, scale=scale)
    pareto_k   = get_pareto_k(loo_result)
    bad        = list(map(int, np.where(pareto_k > k_thresh)[0]))

    if not bad:
        return loo_result, []

    if verbose:
        print(
            f"  Pareto k > {k_thresh} for {len(bad)} obs {bad} — "
            "running az.reloo (exact LOO via RWMH) ..."
        )
    loo_result = az.reloo(wrapper, loo_orig=loo_result, k_thresh=k_thresh, verbose=verbose)
    return loo_result, bad


def load_models_from_config(
    config_path: str,
    chain_relative_override: str | None = None,
) -> list[ModelSpec]:
    """
    Load model specifications from JSON or YAML configuration file.
    
    Args:
        config_path: Path to JSON or YAML config file
        chain_relative_override: Optional override for chain_relative_path from file
    
    Returns:
        List of ModelSpec objects
    
    Config file format (JSON or YAML):
    {
        "chain_relative_path": "chains/degradation_chain.csv",  # Optional
        "models": [
            {"name": "ModelName", "type": "age|null_constant|spatial_ap", "dir": "path/to/results"},
            ...
        ]
    }
    """
    # Detect file type by extension
    ext = os.path.splitext(config_path)[1].lower()
    
    with open(config_path, 'r') as f:
        if ext in ['.json']:
            config = json.load(f)
        elif ext in ['.yaml', '.yml']:
            try:
                import yaml
                config = yaml.safe_load(f)
            except ImportError:
                raise ImportError(
                    "PyYAML is required to load YAML files. "
                    "Install with: pip install pyyaml\n"
                    "Alternatively, use JSON format instead."
                )
        else:
            raise ValueError(f"Unsupported config file extension: {ext}. Use .json, .yaml, or .yml")
    
    if "models" not in config:
        raise ValueError("Config file must contain 'models' key")
    
    models = []
    for model_dict in config["models"]:
        if "name" not in model_dict or "type" not in model_dict or "dir" not in model_dict:
            raise ValueError(
                f"Each model must have 'name', 'type', and 'dir' fields. Got: {model_dict}"
            )
        
        models.append(ModelSpec(
            name=model_dict["name"],
            model_type=model_dict["type"],
            results_dir=model_dict["dir"],
        ))
    
    return models


def update_or_add_model(
    models: list[ModelSpec],
    name: str,
    model_type: str,
    results_dir: str,
) -> list[ModelSpec]:
    """
    Add a new model or update an existing one in the model list.
    
    Args:
        models: Existing list of ModelSpec objects
        name: Model name
        model_type: Model type ("age", "null_constant", "spatial_ap")
        results_dir: Path to model results directory
    
    Returns:
        Updated list of ModelSpec objects
    """
    # Check if model with this name already exists
    for i, model in enumerate(models):
        if model.name == name:
            # Update existing model
            models[i] = ModelSpec(name=name, model_type=model_type, results_dir=results_dir)
            return models
    
    # Add new model
    models.append(ModelSpec(name=name, model_type=model_type, results_dir=results_dir))
    return models


def save_ppc_plot(
    model_name: str,
    mu: np.ndarray,
    sigma: np.ndarray,
    observed: np.ndarray,
    output_path: str,
    n_pp_samples: int,
    random_seed: int,
) -> None:
    rng = np.random.default_rng(random_seed)
    n_total = mu.shape[0]
    n_take = min(n_pp_samples, n_total)
    idx = rng.choice(n_total, size=n_take, replace=False)
    y_rep = rng.normal(loc=mu[idx], scale=sigma[idx, None])

    ppc_idata = az.from_dict(
        posterior_predictive={"m_obs": y_rep[None, :, :]},
        observed_data={"m_obs": observed},
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    az.plot_ppc(ppc_idata, num_pp_samples=min(100, n_take), ax=ax)
    ax.set_title(f"Posterior Predictive Check: {model_name}")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare multiple models with LOO and PPC.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # All four models via CLI
  python %(prog)s --delayed-dir results_age --biphasic-dir results_biphasic \\
      --null-dir results_null --spatial-dir results_spatial \\
      --transcription data.csv --mrna mrna.csv --output-dir output

  # Using JSON config (can be overridden by CLI args)
  python %(prog)s --models-config models.json \\
      --transcription data.csv --mrna mrna.csv --output-dir output
"""
    )
    
    # Model directory arguments
    parser.add_argument("--delayed-dir", help="Path to Delayed (GRW) model results (age model)")
    parser.add_argument("--biphasic-dir", help="Path to Biphasic (poly-A) model results (age model)")
    
    # New model directory arguments (Phase 3)
    parser.add_argument("--null-dir", help="Path to null constant model results")
    parser.add_argument("--spatial-dir", help="Path to spatial AP model results")
    
    # Configuration file argument (Phase 3)
    parser.add_argument(
        "--models-config",
        help="Path to JSON/YAML config file defining models (CLI args override this)"
    )
    parser.add_argument(
        "--chain-relative-path",
        default="chains/degradation_chain.csv",
        help="Relative path to chain CSV within each model directory (default: chains/degradation_chain.csv)"
    )
    
    # Required data arguments
    parser.add_argument("--transcription", required=True, help="Path to transcription CSV")
    parser.add_argument("--mrna", required=True, help="Path to observed mRNA CSV")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    
    # Optional parameters
    parser.add_argument("--n-chains", type=int, default=4)
    parser.add_argument("--n-draws", type=int, default=None)
    parser.add_argument("--n-pp-samples", type=int, default=500)
    parser.add_argument("--seed", type=int, default=14)
    parser.add_argument("--n-ap-bins", type=int, default=5, help="Number of AP bins (for spatial models)")
    parser.add_argument("--n-dv-bins", type=int, default=5, help="Number of DV bins (for spatial models)")

    # reloo options
    parser.add_argument(
        "--reloo-threshold", type=float, default=0.7,
        help="Pareto k threshold above which exact LOO via RWMH sampling is run (default: 0.7)",
    )
    parser.add_argument(
        "--no-reloo", action="store_true",
        help="Disable reloo; use plain PSIS-LOO even when k > threshold",
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Build model list with priority: CLI > config file > defaults
    models = []
    
    # Step 1: Load from config file if provided
    if args.models_config:
        models = load_models_from_config(args.models_config, args.chain_relative_path)
        print(f"Loaded {len(models)} models from config: {args.models_config}")
    
    # Step 2: Apply CLI arguments (add or override)
    if args.delayed_dir:
        models = update_or_add_model(models, "Delayed", "age", args.delayed_dir)
    if args.biphasic_dir:
        models = update_or_add_model(models, "Biphasic", "age", args.biphasic_dir)
    if args.null_dir:
        models = update_or_add_model(models, "Constant", "null_constant", args.null_dir)
    if args.spatial_dir:
        models = update_or_add_model(models, "Spatial", "spatial_ap", args.spatial_dir)
    
    # Step 3: Validate model list
    if len(models) == 0:
        parser.error(
            "No models specified. Provide either:\n"
            "  - Model directories via CLI args (--convolution-dir, --biopolya-dir, etc.)\n"
            "  - A config file via --models-config"
        )
    
    if len(models) < 2:
        parser.error(
            f"Model comparison requires at least 2 models, but only {len(models)} was specified.\n"
            f"Specified model: {models[0].name} ({models[0].model_type})\n"
            f"Add at least one more model via CLI args or config file."
        )
    
    print(f"\nComparing {len(models)} models:")
    for model in models:
        print(f"  - {model.name} ({model.model_type}): {model.results_dir}")
    print()

    # Load data
    transcription = load_matrix(args.transcription)
    observed = load_vector(args.mrna)

    if transcription.shape[0] != observed.shape[0]:
        raise ValueError(
            f"Transcription rows ({transcription.shape[0]}) must equal observed mRNA length ({observed.shape[0]})."
        )

    dt = 20.0 / 60.0

    # Build InferenceData for each model
    idata_dict = {}
    mu_dict = {}
    sigma_dict = {}
    
    for model in models:
        chain_csv = os.path.join(model.results_dir, args.chain_relative_path)
        
        if not os.path.exists(chain_csv):
            raise FileNotFoundError(
                f"Chain file not found for model '{model.name}': {chain_csv}\n"
                f"Check that:\n"
                f"  1. Model directory exists: {model.results_dir}\n"
                f"  2. Chain path is correct: {args.chain_relative_path}"
            )
        
        print(f"Loading {model.name} from {chain_csv}...")
        
        idata, mu, gamma, sigma = build_idata_from_chain(
            chain_csv,
            model.name,
            model.model_type,
            transcription,
            observed,
            dt,
            args.n_chains,
            args.n_draws,
            args.n_ap_bins,
            args.n_dv_bins,
        )
        
        idata_dict[model.name] = idata
        mu_dict[model.name] = mu
        sigma_dict[model.name] = sigma

    # Perform LOO comparison (with optional reloo for high-k observations)
    print("\nPerforming LOO comparison...")
    loo_dict: dict[str, az.ELPDData] = {}
    reloo_info: dict[str, list[int]] = {}  # model_name -> reloo-corrected obs indices

    for model in models:
        idata = idata_dict[model.name]
        if args.no_reloo:
            loo_result = az.loo(idata, pointwise=True, scale="deviance")
            reloo_indices: list[int] = []
        else:
            wrapper = MRNADecaySamplingWrapper(
                idata,
                model.model_type,
                transcription,
                observed,
                dt,
                args.n_ap_bins,
                args.n_dv_bins,
            )
            print(f"  LOO for {model.name} ...")
            loo_result, reloo_indices = maybe_reloo(
                idata, wrapper,
                k_thresh=args.reloo_threshold,
                verbose=True,
            )
        loo_dict[model.name] = loo_result
        reloo_info[model.name] = reloo_indices

    comparison = az.compare(loo_dict)

    # Present models in a fixed, manuscript-friendly order.
    preferred_order = ["Constant", "Spatial", "Delayed", "Biphasic"]
    ordered = [name for name in preferred_order if name in comparison.index]
    remaining = [name for name in comparison.index if name not in ordered]
    comparison = comparison.loc[ordered + remaining]
    print(comparison)

    comparison_csv = os.path.join(args.output_dir, "loo_comparison.csv")
    comparison.to_csv(comparison_csv)

    best_elpd = float(comparison.loc[comparison["rank"] == 0, "elpd_loo"].iloc[0])
    MM = 1/25.4
    fig, ax = plt.subplots(figsize=(140*MM, 50*MM))
    az.plot_compare(
        comparison, ax=ax, order_by_rank=False,
        plot_kwargs={"color_ls_min_ic": "none"},  # suppress arviz's line (always at iloc[0])
    )
    ax.set_title("LOO model comparison (lower is better)", fontsize=10)
    ax.set_ylabel(None)
    ax.axvline(best_elpd, ls="--", color="grey", lw=1)  # draw at best-ranked model
    fig.tight_layout()
    compare_plot = os.path.join(args.output_dir, "loo_compare.pdf")
    fig.savefig(compare_plot, dpi=200)
    plt.close(fig)

    # Generate PPC plots for each model
    ppc_paths = {}
    for model in models:
        # Sanitize filename
        safe_name = model.name.lower().replace(" ", "_")
        ppc_path = os.path.join(args.output_dir, f"ppc_{safe_name}.pdf")
        save_ppc_plot(
            model.name,
            mu_dict[model.name],
            sigma_dict[model.name],
            observed,
            ppc_path,
            args.n_pp_samples,
            args.seed
        )
        ppc_paths[model.name] = ppc_path

    # Generate summary using model rank (not display row order).
    ranked_comparison = comparison.sort_values("rank")
    top_model = ranked_comparison.index[0]
    top_deviance = float(ranked_comparison.iloc[0]["elpd_loo"])
    second_deviance = float(ranked_comparison.iloc[1]["elpd_loo"])
    d_loo = abs(second_deviance - top_deviance)
    d_se = float(ranked_comparison.iloc[1]["dse"])
    evidence = "clear" if d_loo > 2.0 * d_se else "weak_or_indistinguishable"

    summary_txt = os.path.join(args.output_dir, "comparison_summary.txt")
    with open(summary_txt, "w", encoding="utf-8") as handle:
        handle.write("LOO Model Comparison\n")
        handle.write("====================\n")
        handle.write(comparison.to_string())
        handle.write("\n\n")
        handle.write(f"Winner (rank 0): {top_model}\n")
        handle.write(f"d_loo (deviance gap runner-up vs winner): {d_loo:.4f}\n")
        handle.write(f"SE of difference (runner-up dse): {d_se:.4f}\n")
        handle.write(f"Evidence strength (d_loo > 2*se): {evidence}\n")
        handle.write(f"\nreloo threshold: {args.reloo_threshold}"
                     f"  (disabled: {args.no_reloo})\n")
        handle.write("Reloo-corrected observations per model:\n")
        for model_name, ridx in reloo_info.items():
            if ridx:
                handle.write(f"  {model_name}: obs {ridx}\n")
            else:
                handle.write(f"  {model_name}: none (all k <= {args.reloo_threshold})\n")
        handle.write(f"\nCompare plot: {compare_plot}\n")
        handle.write("\nPPC plots:\n")
        for model_name, ppc_path in ppc_paths.items():
            handle.write(f"  {model_name}: {ppc_path}\n")

    print("\n" + "="*60)
    print(comparison)
    print()
    print(f"Winner: {top_model}")
    print(f"d_loo (deviance gap): {d_loo:.4f} | se: {d_se:.4f} | evidence: {evidence}")
    for model_name, ridx in reloo_info.items():
        if ridx:
            print(f"  reloo applied to {model_name}: obs {ridx}")
    print(f"\nSaved: {comparison_csv}")
    print(f"Saved: {compare_plot}")
    for model_name, ppc_path in ppc_paths.items():
        print(f"Saved: {ppc_path}")
    print(f"Saved: {summary_txt}")
    print("="*60)



if __name__ == "__main__":
    main()
