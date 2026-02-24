#!/usr/bin/env python
"""
LOO diagnostics for Bayesian model comparison.

For each model (specified via --models-config or CLI):
  1. Rebuild ArviZ InferenceData from chain CSV.
  2. Run az.loo(idata, pointwise=True) to get per-observation Pareto k-hat.
  3. Compute per-observation log-likelihood variance across posterior samples.
  4. Map observation indices to AP bins (obs k -> AP bin k // n_dv_bins).
  5. Produce diagnostic plots and a summary text file.

Usage:
    micromamba run -p .snakemake/conda/... python scripts/09_loo_diagnostics.py \\
        --models-config config/models_stripe3.json \\
        --transcription data/processed_transcription_data/transcription_traces_no_ids_stripe3_1200.csv \\
        --mrna results_1200/data/processed_mRNA_data_stripe3/e8_9um_sass_formodel.csv \\
        --output-dir comparison_output/diagnostics \\
        --n-ap-bins 5 --n-dv-bins 5
"""

import argparse
import importlib.util
import os
import sys

import numpy as np
import arviz as az
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Import from 08_compare_models_loo_ppc.py
# A leading digit in the filename prevents standard import, so use importlib.
# ---------------------------------------------------------------------------

def _import_compare_module():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, "08_compare_models_loo_ppc.py")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cannot find {path}")
    spec = importlib.util.spec_from_file_location("compare_models", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_cmp = _import_compare_module()

build_idata_from_chain = _cmp.build_idata_from_chain
load_models_from_config = _cmp.load_models_from_config
update_or_add_model = _cmp.update_or_add_model
load_matrix = _cmp.load_matrix
load_vector = _cmp.load_vector


# ---------------------------------------------------------------------------
# Colour palette for AP bins (up to 8)
# ---------------------------------------------------------------------------

_AP_COLORS = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2",
    "#59a14f", "#edc948", "#b07aa1", "#ff9da7",
]


# ---------------------------------------------------------------------------
# Diagnostic helpers
# ---------------------------------------------------------------------------

def map_obs_to_ap_bin(n_obs: int, n_dv_bins: int) -> np.ndarray:
    """Return AP bin index for each observation index (shape n_obs)."""
    return np.array([k // n_dv_bins for k in range(n_obs)])


def map_obs_to_dv_bin(n_obs: int, n_dv_bins: int) -> np.ndarray:
    """Return DV bin index for each observation index (shape n_obs)."""
    return np.array([k % n_dv_bins for k in range(n_obs)])


def compute_loglik_variance(idata: az.InferenceData) -> np.ndarray:
    """
    Per-observation variance of log-likelihood across all posterior samples.

    Extracts idata.log_likelihood["m_obs"] (chains, draws, n_obs), flattens
    the first two axes, and returns Var_s[log p(y_i | theta_s)] for each obs i.
    Shape: (n_obs,).
    """
    ll = idata.log_likelihood["m_obs"].values  # (chains, draws, n_obs)
    ll_flat = ll.reshape(-1, ll.shape[-1])      # (n_samples, n_obs)
    return ll_flat.var(axis=0)


def run_pointwise_loo(idata: az.InferenceData) -> az.ELPDData:
    """Run PSIS-LOO with pointwise=True and return ELPDData."""
    return az.loo(idata, pointwise=True)


def get_pareto_k(loo_result: az.ELPDData) -> np.ndarray:
    """Extract Pareto k-hat values as a plain numpy array."""
    k = loo_result.pareto_k
    if hasattr(k, "values"):   # xarray DataArray in newer ArviZ
        return k.values
    return np.asarray(k)


def flag_observations_vs_pp(
    mu_samples: np.ndarray,
    sigma_samples: np.ndarray,
    observed: np.ndarray,
    pareto_k: np.ndarray,
    ap_bins: np.ndarray,
    dv_bins: np.ndarray,
    threshold: float = 0.7,
    seed: int = 42,
) -> list[str]:
    """
    For observations with Pareto k > threshold, compare the observed value
    against the Monte-Carlo posterior predictive interval to diagnose whether
    high k is driven by outlier-like observations.

    Parameters
    ----------
    mu_samples    : (n_samples, n_obs) posterior predicted means per sample.
    sigma_samples : (n_samples,) posterior noise sigma per sample.
    observed      : (n_obs,) observed mRNA values.
    pareto_k      : (n_obs,) Pareto k-hat values from PSIS-LOO.
    ap_bins       : (n_obs,) AP bin index per observation.
    dv_bins       : (n_obs,) DV bin index per observation.
    threshold     : Flag observations with k above this value (default 0.7).
    seed          : RNG seed for reproducibility.

    Returns
    -------
    List of formatted text lines suitable for the summary file.
    """
    rng = np.random.default_rng(seed)
    n_samples, n_obs = mu_samples.shape

    # Draw one PP sample per posterior sample × observation
    noise = rng.standard_normal((n_samples, n_obs))
    pp = mu_samples + sigma_samples[:, np.newaxis] * noise   # (n_samples, n_obs)

    pp_lo  = np.percentile(pp, 2.5,  axis=0)
    pp_med = np.median(pp, axis=0)
    pp_hi  = np.percentile(pp, 97.5, axis=0)

    bad = np.where(pareto_k > threshold)[0]

    header = (
        f"  {'obs_i':>5}  {'AP':>3}  {'DV':>3}  {'k_hat':>6}  "
        f"{'y_obs':>8}  {'pp_2.5':>8}  {'pp_50':>8}  {'pp_97.5':>8}  "
        f"{'in_CI':>6}  {'pp_rank%':>8}  {'pattern'}"
    )
    lines = [
        f"  High-k observations vs posterior predictive (k > {threshold}):",
        header,
        "  " + "-" * 95,
    ]

    if len(bad) == 0:
        lines.append(f"  None — all observations have k < {threshold}.")
        return lines

    n_outside = 0
    for i in bad:
        in_ci   = bool(pp_lo[i] <= observed[i] <= pp_hi[i])
        pp_rank = float(np.mean(pp[:, i] < observed[i])) * 100

        if not in_ci:
            pattern = "OUTLIER"
            n_outside += 1
        elif pp_rank > 90 or pp_rank < 10:
            pattern = "LEVERAGE (tail)"    # plausible but at PP boundary
        else:
            pattern = "LEVERAGE (central)" # posterior pulled toward this point

        lines.append(
            f"  {i:>5d}  {int(ap_bins[i]):>3d}  {int(dv_bins[i]):>3d}  "
            f"{pareto_k[i]:>6.3f}  {float(observed[i]):>8.4f}  "
            f"{pp_lo[i]:>8.4f}  {pp_med[i]:>8.4f}  {pp_hi[i]:>8.4f}  "
            f"{'yes' if in_ci else 'NO':>6}  {pp_rank:>8.1f}  {pattern}"
        )

    lines.append(
        f"  {n_outside}/{len(bad)} high-k observations fall outside the 95% PP CI."
    )
    return lines


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_pp_check(
    mu_samples: np.ndarray,
    sigma_samples: np.ndarray,
    observed: np.ndarray,
    pareto_k: np.ndarray,
    ap_bins: np.ndarray,
    dv_bins: np.ndarray,
    model_name: str,
    output_path: str,
    threshold: float = 0.7,
    seed: int = 42,
) -> None:
    """
    For observations with Pareto k > threshold, plot observed value vs the
    95 % posterior predictive interval (error bar) and the PP median.

    Points that fall outside the 95 % PP CI are highlighted in red; those
    inside are shown in grey.  The observation index, AP bin, and DV bin
    are annotated on each point.
    """
    rng = np.random.default_rng(seed)
    n_samples, n_obs = mu_samples.shape

    noise = rng.standard_normal((n_samples, n_obs))
    pp = mu_samples + sigma_samples[:, np.newaxis] * noise

    pp_lo  = np.percentile(pp, 2.5,  axis=0)
    pp_med = np.median(pp, axis=0)
    pp_hi  = np.percentile(pp, 97.5, axis=0)

    bad = np.where(pareto_k > threshold)[0]
    if len(bad) == 0:
        return  # nothing to plot

    in_ci  = np.array([pp_lo[i] <= observed[i] <= pp_hi[i] for i in bad])
    colors = np.where(in_ci, "#888888", "#e15759")

    fig, ax = plt.subplots(figsize=(max(6, len(bad) * 0.9), 5))

    for idx, (i, c) in enumerate(zip(bad, colors)):
        # Error bar spans the 95% PP CI; centre is PP median
        ax.errorbar(
            idx, pp_med[i],
            yerr=[[pp_med[i] - pp_lo[i]], [pp_hi[i] - pp_med[i]]],
            fmt="none", color="#cccccc", capsize=5, linewidth=2, zorder=2,
        )
        # PP median marker
        ax.scatter(idx, pp_med[i], color="#cccccc", s=60, zorder=3,
                   label="PP median" if idx == 0 else None)
        # Observed value
        ax.scatter(idx, observed[i], color=c, s=80, zorder=4, marker="D",
                   label=("observed (in CI)" if (c == "#888888" and idx == 0)
                          else ("observed (outside CI)" if (c == "#e15759" and idx == 0)
                                else None)))
        label = f"i={i}\nAP{int(ap_bins[i])} DV{int(dv_bins[i])}\nk={pareto_k[i]:.2f}"
        ax.annotate(
            label, xy=(idx, observed[i]),
            xytext=(0, 10), textcoords="offset points",
            ha="center", va="bottom", fontsize=7,
        )

    ax.set_xticks(range(len(bad)))
    ax.set_xticklabels([f"obs {i}" for i in bad], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("mRNA value")
    ax.set_title(
        f"High-k observations vs 95% PP interval (k > {threshold}) — {model_name}"
    )

    # De-duplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=8, loc="best")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_pareto_k(
    pareto_k: np.ndarray,
    ap_bins: np.ndarray,
    model_name: str,
    n_ap_bins: int,
    output_path: str,
    dv_bins: np.ndarray | None = None,
) -> None:
    """
    Scatter plot of per-observation Pareto k-hat.

    AP bins shown as background bands; threshold lines at k=0.5 and k=0.7.
    Each point is annotated with its DV bin index if dv_bins is provided.
    """
    n_obs = len(pareto_k)
    n_dv = n_obs // n_ap_bins

    fig, ax = plt.subplots(figsize=(10, 4))

    # Background shading per AP bin
    for b in range(n_ap_bins):
        x0 = b * n_dv - 0.5
        x1 = (b + 1) * n_dv - 0.5
        ax.axvspan(x0, x1, alpha=0.08, color=_AP_COLORS[b % len(_AP_COLORS)], linewidth=0)

    # Scatter: one colour per AP bin
    for b in range(n_ap_bins):
        mask = ap_bins == b
        obs_idx = np.where(mask)[0]
        ax.scatter(
            obs_idx, pareto_k[mask],
            color=_AP_COLORS[b % len(_AP_COLORS)],
            label=f"AP bin {b}",
            s=60, zorder=3,
        )

    # Annotate each point with its DV bin index
    if dv_bins is not None:
        for i in range(n_obs):
            ax.annotate(
                str(dv_bins[i]),
                xy=(i, pareto_k[i]),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center", va="bottom",
                fontsize=7, color="#333333",
            )

    # Threshold lines
    ax.axhline(0.5, color="orange", linestyle="--", linewidth=1.2, label="k=0.5 (marginal)")
    ax.axhline(0.7, color="red",    linestyle="--", linewidth=1.2, label="k=0.7 (unreliable)")

    ax.set_xlabel("Observation index")
    ax.set_ylabel("Pareto k\u0302")
    ax.set_title(f"PSIS Pareto k\u0302 per observation \u2014 {model_name}")
    ax.set_xlim(-0.5, n_obs - 0.5)
    ax.legend(fontsize=8, ncol=min(n_ap_bins + 2, 8), loc="upper left")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_loglik_variance(
    ll_var: np.ndarray,
    ap_bins: np.ndarray,
    model_name: str,
    n_ap_bins: int,
    output_path: str,
    dv_bins: np.ndarray | None = None,
) -> None:
    """Per-observation log-likelihood variance, same layout as plot_pareto_k."""
    n_obs = len(ll_var)
    n_dv = n_obs // n_ap_bins

    fig, ax = plt.subplots(figsize=(10, 4))

    for b in range(n_ap_bins):
        x0 = b * n_dv - 0.5
        x1 = (b + 1) * n_dv - 0.5
        ax.axvspan(x0, x1, alpha=0.08, color=_AP_COLORS[b % len(_AP_COLORS)], linewidth=0)

    for b in range(n_ap_bins):
        mask = ap_bins == b
        obs_idx = np.where(mask)[0]
        ax.scatter(
            obs_idx, ll_var[mask],
            color=_AP_COLORS[b % len(_AP_COLORS)],
            label=f"AP bin {b}",
            s=60, zorder=3,
        )

    # Annotate each point with its DV bin index
    if dv_bins is not None:
        for i in range(n_obs):
            ax.annotate(
                str(dv_bins[i]),
                xy=(i, ll_var[i]),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center", va="bottom",
                fontsize=7, color="#333333",
            )

    ax.set_xlabel("Observation index")
    ax.set_ylabel("Var[log p(y | \u03b8)]")
    ax.set_title(f"Per-observation log-likelihood variance \u2014 {model_name}")
    ax.set_xlim(-0.5, n_obs - 0.5)
    ax.legend(fontsize=8, ncol=min(n_ap_bins, 8), loc="upper left")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_k_by_ap_bin(
    mean_k_by_bin: dict,
    n_ap_bins: int,
    output_path: str,
) -> None:
    """
    Grouped bar chart: mean Pareto k per AP bin, one group of bars per model.
    """
    models = list(mean_k_by_bin.keys())
    n_models = len(models)
    x = np.arange(n_ap_bins)
    width = 0.8 / n_models
    model_colors = plt.cm.tab10(np.linspace(0, 0.9, n_models))

    fig, ax = plt.subplots(figsize=(8, 4))

    for i, name in enumerate(models):
        offsets = x + (i - (n_models - 1) / 2.0) * width
        ax.bar(offsets, mean_k_by_bin[name], width=width * 0.9,
               color=model_colors[i], label=name, alpha=0.85)

    ax.axhline(0.5, color="orange", linestyle="--", linewidth=1.2, label="k=0.5")
    ax.axhline(0.7, color="red",    linestyle="--", linewidth=1.2, label="k=0.7")

    ax.set_xlabel("AP bin")
    ax.set_ylabel("Mean Pareto k\u0302")
    ax.set_title("Mean PSIS Pareto k\u0302 per AP bin by model")
    ax.set_xticks(x)
    ax.set_xticklabels([f"bin {b}" for b in range(n_ap_bins)])
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_posterior_D_by_ap_bin(
    idata: az.InferenceData,
    model_name: str,
    n_ap_bins: int,
    output_path: str,
) -> None:
    """
    Posterior degradation rate D and implied mRNA half-life per AP bin.

    Left panel: D (degradation rate, min^-1) median +/- 95% CI.
    Right panel: half-life = ln(2) / D median +/- 95% CI.

    Note: CI bounds are inverted for half-life because larger D means shorter half-life.
    """
    if "D" not in idata.posterior:
        print(f"  [skip] Posterior 'D' not found in {model_name} — no D-by-AP-bin plot produced")
        return

    D_vals = idata.posterior["D"].values     # (chains, draws, n_ap_bins)
    D_flat = D_vals.reshape(-1, n_ap_bins)   # (n_samples, n_ap_bins)

    D_median = np.median(D_flat, axis=0)
    D_lo     = np.percentile(D_flat, 2.5,  axis=0)
    D_hi     = np.percentile(D_flat, 97.5, axis=0)

    # Half-life = ln(2) / D; inverting the CI bounds correctly
    HL_median = np.log(2) / D_median
    HL_lo     = np.log(2) / D_hi   # larger D -> shorter half-life
    HL_hi     = np.log(2) / D_lo

    bins = np.arange(n_ap_bins)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    # Left: D per AP bin
    ax = axes[0]
    ax.errorbar(
        bins, D_median,
        yerr=[D_median - D_lo, D_hi - D_median],
        fmt="o-", color="#4e79a7", capsize=5, linewidth=1.5,
    )
    ax.set_xlabel("AP bin  (0 = most anterior)")
    ax.set_ylabel("D (min\u207b\u00b9)")
    ax.set_title(f"Posterior degradation rate \u2014 {model_name}")
    ax.set_xticks(bins)

    # Right: half-life per AP bin
    ax = axes[1]
    ax.errorbar(
        bins, HL_median,
        yerr=[HL_median - HL_lo, HL_hi - HL_median],
        fmt="o-", color="#e15759", capsize=5, linewidth=1.5,
    )
    ax.set_xlabel("AP bin  (0 = most anterior)")
    ax.set_ylabel("Half-life (min)")
    ax.set_title(f"Implied mRNA half-life \u2014 {model_name}")
    ax.set_xticks(bins)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="LOO diagnostics: per-observation Pareto k and log-lik variance.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--models-config",       help="Path to JSON/YAML model config file")
    parser.add_argument("--convolution-dir",     help="Convolution model results dir (age)")
    parser.add_argument("--biopolya-dir",        help="BioPolyA model results dir (age)")
    parser.add_argument("--null-dir",            help="NullConstant model results dir")
    parser.add_argument("--spatial-dir",         help="Spatial model results dir")
    parser.add_argument(
        "--chain-relative-path",
        default="chains/degradation_chain.csv",
        help="Relative path to chain CSV within each model directory",
    )
    parser.add_argument("--transcription", required=True, help="Transcription matrix CSV")
    parser.add_argument("--mrna",          required=True, help="Observed mRNA CSV")
    parser.add_argument("--output-dir",    required=True, help="Output directory for plots/summary")
    parser.add_argument("--n-chains",  type=int, default=4)
    parser.add_argument("--n-draws",   type=int, default=None,
                        help="Draws per chain (inferred from chain file if omitted)")
    parser.add_argument("--n-ap-bins", type=int, default=5)
    parser.add_argument("--n-dv-bins", type=int, default=5)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Build model list (config file first, then CLI overrides)
    models = []
    if args.models_config:
        models = load_models_from_config(args.models_config, args.chain_relative_path)
        print(f"Loaded {len(models)} models from {args.models_config}")
    if args.convolution_dir:
        models = update_or_add_model(models, "Convolution",   "age",          args.convolution_dir)
    if args.biopolya_dir:
        models = update_or_add_model(models, "BioPolyA",      "age",          args.biopolya_dir)
    if args.null_dir:
        models = update_or_add_model(models, "NullConstant",  "null_constant", args.null_dir)
    if args.spatial_dir:
        models = update_or_add_model(models, "Spatial",       "spatial_ap",    args.spatial_dir)

    if not models:
        parser.error(
            "No models specified. Provide --models-config or individual --*-dir flags."
        )

    print(f"\nDiagnosing {len(models)} models:")
    for m in models:
        print(f"  - {m.name} ({m.model_type}): {m.results_dir}")
    print()

    # Load shared data
    transcription = load_matrix(args.transcription)
    observed      = load_vector(args.mrna)
    dt = 20.0 / 60.0  # minutes per time step

    n_obs   = observed.shape[0]
    ap_bins = map_obs_to_ap_bin(n_obs, args.n_dv_bins)   # (n_obs,)
    dv_bins = map_obs_to_dv_bin(n_obs, args.n_dv_bins)   # (n_obs,)

    # Containers for cross-model summary plot
    mean_k_by_bin: dict[str, np.ndarray] = {}

    summary_lines: list[str] = [
        "LOO Diagnostics Summary",
        "=" * 60,
        f"n_obs={n_obs}  n_ap_bins={args.n_ap_bins}  n_dv_bins={args.n_dv_bins}",
        "",
    ]

    for model in models:
        chain_csv = os.path.join(model.results_dir, args.chain_relative_path)

        if not os.path.exists(chain_csv):
            print(f"[WARNING] Chain not found for {model.name}: {chain_csv} — skipping")
            summary_lines.append(f"Model: {model.name}  SKIPPED (chain not found: {chain_csv})")
            summary_lines.append("")
            continue

        print(f"Processing {model.name} ({model.model_type})")
        print(f"  Chain: {chain_csv}")

        # Rebuild InferenceData
        idata, mu_samples, _gamma, sigma_samples = build_idata_from_chain(
            chain_csv, model.name, model.model_type,
            transcription, observed, dt,
            args.n_chains, args.n_draws,
            args.n_ap_bins, args.n_dv_bins,
        )

        # Pointwise LOO
        print("  Running az.loo(pointwise=True) ...")
        loo_result = run_pointwise_loo(idata)
        pareto_k   = get_pareto_k(loo_result)
        ll_var     = compute_loglik_variance(idata)

        # AP-bin aggregates
        k_by_bin = np.array([
            pareto_k[ap_bins == b].mean() for b in range(args.n_ap_bins)
        ])
        mean_k_by_bin[model.name] = k_by_bin

        n_bad_05 = int((pareto_k > 0.5).sum())
        n_bad_07 = int((pareto_k > 0.7).sum())

        print(f"  p_loo={loo_result.p_loo:.3f}  elpd_loo={loo_result.elpd_loo:.3f}")
        print(f"  Pareto k: max={pareto_k.max():.3f}  mean={pareto_k.mean():.3f}")
        print(f"  n k>0.5: {n_bad_05}/{n_obs}   n k>0.7: {n_bad_07}/{n_obs}")

        # Posterior predictive check for high-k observations
        pp_flag_lines = flag_observations_vs_pp(
            mu_samples, sigma_samples, observed, pareto_k,
            ap_bins, dv_bins, threshold=0.7,
        )
        for line in pp_flag_lines:
            print(line)

        summary_lines += [
            f"Model: {model.name}  ({model.model_type})",
            f"  p_loo    = {loo_result.p_loo:.3f}",
            f"  elpd_loo = {loo_result.elpd_loo:.3f}",
            f"  Pareto k: max={pareto_k.max():.3f}  mean={pareto_k.mean():.3f}",
            f"  n k>0.5: {n_bad_05}/{n_obs}",
            f"  n k>0.7: {n_bad_07}/{n_obs}",
            "  Mean Pareto k per AP bin:",
        ]
        for b, kv in enumerate(k_by_bin):
            summary_lines.append(f"    AP bin {b}: {kv:.3f}")
        summary_lines += [""] + pp_flag_lines + [""]

        safe_name = model.name.lower().replace(" ", "_")

        # Plot: Pareto k per observation
        out_k = os.path.join(args.output_dir, f"pareto_k_{safe_name}.png")
        plot_pareto_k(pareto_k, ap_bins, model.name, args.n_ap_bins, out_k, dv_bins=dv_bins)
        print(f"  Saved: {out_k}")

        # Plot: log-likelihood variance per observation
        out_var = os.path.join(args.output_dir, f"loglik_variance_{safe_name}.png")
        plot_loglik_variance(ll_var, ap_bins, model.name, args.n_ap_bins, out_var, dv_bins=dv_bins)
        print(f"  Saved: {out_var}")

        # Plot: observed vs PP interval for high-k observations
        if n_bad_07 > 0:
            out_pp = os.path.join(args.output_dir, f"pp_check_high_k_{safe_name}.png")
            plot_pp_check(
                mu_samples, sigma_samples, observed, pareto_k,
                ap_bins, dv_bins, model.name, out_pp, threshold=0.7,
            )
            print(f"  Saved: {out_pp}")

        # Plot: posterior D and half-life by AP bin (spatial model only)
        if model.model_type == "spatial_ap":
            out_D = os.path.join(args.output_dir, f"posterior_D_{safe_name}.png")
            plot_posterior_D_by_ap_bin(idata, model.name, args.n_ap_bins, out_D)
            print(f"  Saved: {out_D}")

    # Cross-model summary plot (mean k per AP bin)
    if mean_k_by_bin:
        out_summary_k = os.path.join(args.output_dir, "pareto_k_by_ap_bin_summary.png")
        plot_k_by_ap_bin(mean_k_by_bin, args.n_ap_bins, out_summary_k)
        print(f"\nSaved: {out_summary_k}")

    # Write summary text
    summary_path = os.path.join(args.output_dir, "loo_diagnostics_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(summary_lines))
    print(f"Saved: {summary_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
