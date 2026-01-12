#!/usr/bin/env python
"""
Infer spatially-varying mRNA degradation rates using Bayesian inference.

This script implements a Bayesian model to infer degradation rates (D) across spatial bins
by fitting transcription input functions (F) to observed mRNA counts (m) using an ODE model:
    dm/dt = γ*F(t) - D*m

The model is fit using PyMC for MCMC sampling.

Prior choices (matching Julia implementation in infer_D_across_stripe2.jl):
- D (degradation rates): TruncatedNormal(mu=0.0, sigma=1.0, lower=0.0) - unbounded above
- gamma (transcription scaling): InverseGamma(alpha=2, beta=3)
- sigma (observation noise): InverseGamma(alpha=2, beta=3)

Based on infer_D_across_stripe2.jl
"""

import os
import sys
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from scipy.integrate import odeint, trapezoid
from scipy.interpolate import interp1d
from scipy.stats import invgamma
import pytensor
import pytensor.tensor as pt
from pytensor.graph import Op, Apply
import argparse


def smfish_ode(m, t, F_interp, D, gamma):
    """
    ODE model for mRNA dynamics with transcription input and degradation.
    
    dm/dt = γ*F(t) - D*m
    
    Args:
        m: mRNA concentration
        t: time
        F_interp: interpolated transcription function
        D: degradation rate
        gamma: transcription scaling factor
    
    Returns:
        dm/dt
    """
    F_t = F_interp(t)
    dmdt = gamma * F_t - D * m
    return dmdt


def integrate_ode_solution_numerical(D, gamma, F_interp, t_array):
    """
    Solve ODE numerically and return final mRNA concentration.
    
    Args:
        D: degradation rate (scalar)
        gamma: transcription scaling factor (scalar)
        F_interp: interpolated transcription function
        t_array: time points
    
    Returns:
        mRNA concentration at final time point
    """
    m0 = 0.0  # initial condition
    sol = odeint(smfish_ode, m0, t_array, args=(F_interp, D, gamma))
    return sol[-1, 0]


def solve_ode_analytical(D, gamma, F_values, t_array):
    """
    Analytical solution to dm/dt = γ*F(t) - D*m with m(0) = 0.
    
    Using integrating factor method with numerical stability improvements:
    m(t) = γ * exp(-D*t) * ∫[0 to t] F(s) * exp(D*s) ds
    
    To avoid overflow, we compute: exp(-D*T) * ∫ F(s) * exp(D*s) ds
    by working with exp(D*(s-T)) = exp(-D*(T-s)) which stays bounded.
    
    Args:
        D: degradation rate (scalar or array)
        gamma: transcription scaling factor (scalar)
        F_values: transcription values at time points (array)
        t_array: time points (array)
    
    Returns:
        mRNA concentration at final time point
    """
    # Ensure inputs are numpy arrays for numerical computation
    D = np.asarray(D, dtype=np.float64)
    gamma = np.asarray(gamma, dtype=np.float64)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)
    
    # Clip D to reasonable range to avoid numerical issues
    # Lower bound prevents division by zero, upper bound prevents exp() overflow
    # With unbounded prior, sampler may explore large D values, so clip conservatively
    D = np.clip(D, 1e-6, 10.0)
    
    t_final = t_array[-1]
    
    # Compute integrand using exp(D*(t - t_final)) = exp(-D*(t_final - t))
    # This keeps the exponential bounded since t <= t_final
    exp_arg = D * (t_array - t_final)
    exp_term = np.exp(np.clip(exp_arg, -700, 700))  # Clip to prevent overflow/underflow
    integrand = F_values * exp_term
    
    # Integrate using trapezoidal rule
    integral = trapezoid(integrand, t_array)
    
    # Final solution: m(T) = γ * integral
    # (the exp(-D*T) and exp(D*T) terms cancel in the reformulation)
    m_final = gamma * integral
    
    # Ensure output is valid and within reasonable bounds
    if not np.isfinite(m_final):
        m_final = 0.0
    m_final = np.clip(m_final, 0, 1e6)  # Reasonable upper bound for mRNA counts
    
    return float(m_final)


class ODESolverOp(Op):
    """
    Custom PyTensor Op for solving the mRNA ODE analytically.
    
    This allows PyMC to use the ODE solver within the computational graph.
    """
    
    def __init__(self, F_values, t_array):
        """
        Args:
            F_values: Transcription values at time points (numpy array)
            t_array: Time points (numpy array)
        """
        self.F_values = np.asarray(F_values, dtype=np.float64)
        self.t_array = np.asarray(t_array, dtype=np.float64)
    
    def make_node(self, D, gamma):
        # Convert inputs to tensor variables
        D = pt.as_tensor_variable(D)
        gamma = pt.as_tensor_variable(gamma)
        # Output is a scalar
        output = pt.dscalar()
        return Apply(self, [D, gamma], [output])
    
    def perform(self, node, inputs, outputs):
        D, gamma = inputs
        # Compute the ODE solution
        result = solve_ode_analytical(D, gamma, self.F_values, self.t_array)
        outputs[0][0] = np.array(result, dtype=np.float64)
    
    def grad(self, inputs, output_grads):
        # Numerical gradients (PyMC will handle this automatically)
        # Return None to use automatic differentiation
        return [pytensor.gradient.grad_undefined(self, i, inp) 
                for i, inp in enumerate(inputs)]


def calculate_expected_mRNA(D_array, gamma, F_data_arrays, t_array):
    """
    Calculate expected mRNA for all spatial bins and traces given parameters.
    
    Matches Julia implementation: For each of 5 bins with degradation rate D[i],
    calculate expected mRNA for all 5 transcription traces in that bin.
    This produces 25 expected values total (5 bins × 5 traces per bin).
    
    Args:
        D_array: Array of degradation rates for each bin [n_bins]
        gamma: transcription scaling factor
        F_data_arrays: List of transcription data arrays for each bin [n_bins][n_traces_per_bin, n_timepoints]
        t_array: time points
    
    Returns:
        Array of expected mRNA concentrations [n_bins * n_traces_per_bin]
    """
    n_bins = len(D_array)
    expected_m_all = []
    
    for i in range(n_bins):
        # For bin i with degradation rate D[i], calculate expected mRNA
        # for EACH transcription trace (not averaged)
        n_traces = F_data_arrays[i].shape[0]
        for j in range(n_traces):
            F_trace = F_data_arrays[i][j, :]
            m_expected = solve_ode_analytical(D_array[i], gamma, F_trace, t_array)
            expected_m_all.append(m_expected)
    
    return np.array(expected_m_all)


def build_pymc_model(F_data_arrays, m_data, t_array, n_bins=5):
    """
    Build PyMC Bayesian model for inferring degradation rates.
    
    Matches Julia Turing model structure:
    - For each of 5 spatial bins, sample one degradation rate D[i]
    - For each bin i, calculate expected mRNA for all 5 transcription traces using D[i]
    - This produces 25 expected values compared to 25 observed mRNA values
    
    Prior distributions (matching Julia):
        D ~ TruncatedNormal(0, 1) for each spatial bin
        γ ~ InverseGamma(2, 3) - transcription scaling
        σ ~ InverseGamma(2, 3) - observation noise
    
    Likelihood:
        m_obs ~ MvNormal(expected_m_all, σ²*I)
    
    Args:
        F_data_arrays: List of transcription data arrays [n_bins][n_traces_per_bin, n_timepoints]
        m_data: Observed mRNA data [n_bins * n_traces_per_bin] = 25 values
        t_array: time points
        n_bins: number of spatial bins
    
    Returns:
        PyMC model
    """
    print("Building PyMC model...")
    
    # Get structure info
    n_traces_per_bin = F_data_arrays[0].shape[0]
    total_observations = n_bins * n_traces_per_bin
    
    print(f"  {n_bins} spatial bins")
    print(f"  {n_traces_per_bin} traces per bin")
    print(f"  {total_observations} total expected mRNA values")
    
    with pm.Model() as model:
        # Priors for degradation rates (one per spatial bin)
        # 
        # Matching Julia implementation: D ~ filldist(truncated(Normal(0, 1)), 5)
        # - truncated(Normal(0, 1)) in Julia = lower bound at 0, no upper bound
        # - This allows D to range from 0 to ∞, exploring full posterior
        # - Less restrictive than biological constraints but numerically more stable
        #   with undefined gradients (flatter tails → better-conditioned mass matrix)
        #
        D = pm.TruncatedNormal(
            'D', 
            mu=0.0,
            sigma=1.0, 
            lower=0.0,
            shape=n_bins
        )
        
        # Prior for transcription scaling factor
        # Matching Julia implementation: γ ~ InverseGamma(2, 3)
        # InverseGamma(alpha=2, beta=3): Mean = 3, variance = 9
        gamma = pm.InverseGamma('gamma', alpha=2, beta=3)
        
        # Prior for observation noise
        # Matching Julia implementation: σ ~ InverseGamma(2, 3)
        # InverseGamma(alpha=2, beta=3): Mean = 3, variance = 9
        sigma = pm.InverseGamma('sigma', alpha=2, beta=3)
        
        # Expected mRNA concentrations for all bins and traces
        # For each bin i (with D[i]), compute expected mRNA for all traces in that bin
        ode_ops_all = []
        expected_m_list = []
        
        for i in range(n_bins):
            # For bin i, create ODE solver for each of its transcription traces
            for j in range(n_traces_per_bin):
                F_trace = F_data_arrays[i][j, :]
                ode_op = ODESolverOp(F_trace, t_array)
                m_ij = ode_op(D[i], gamma)  # Use D[i] for all traces in bin i
                expected_m_list.append(m_ij)
        
        expected_m = pm.Deterministic(
            'expected_m',
            pt.stack(expected_m_list)
        )
        
        # Likelihood: MATCHES JULIA m ~ MvNormal(integrand_all, σ^2*I)
        # Using independent Normal distributions (equivalent to MvNormal with diagonal covariance)
        m_obs = pm.Normal(
            'm_obs', 
            mu=expected_m, 
            sigma=sigma, 
            observed=m_data
        )
    
    return model


def run_mcmc_inference(model, n_samples=10000, n_chains=4, target_accept=0.8):
    """
    Run MCMC sampling using NUTS algorithm.
    
    Args:
        model: PyMC model
        n_samples: number of samples per chain
        n_chains: number of parallel chains
        target_accept: target acceptance rate for NUTS
    
    Returns:
        InferenceData object with samples
    """
    print(f"Running MCMC inference: {n_samples} samples × {n_chains} chains...")
    
    with model:
        # Sample using NUTS (No-U-Turn Sampler)
        trace = pm.sample(
            draws=n_samples,
            chains=n_chains,
            cores=n_chains,
            target_accept=target_accept,
            init='adapt_diag',  # Better initialization for numerical stability
            return_inferencedata=True,
            random_seed=14
        )
    
    return trace


def save_chain_results(trace, output_path):
    """
    Save MCMC chain samples to CSV.
    
    Args:
        trace: InferenceData object
        output_path: path to save CSV
    """
    print(f"Saving chain results to: {output_path}")
    
    # Extract samples and convert to DataFrame
    # Use group='posterior' and flatten=True to get all dimensions as separate columns
    posterior = trace.posterior
    
    # Stack chains and draws
    samples_dict = {}
    for var_name in ['D', 'gamma', 'sigma']:
        if var_name in posterior:
            var_data = posterior[var_name].values
            
            # Flatten chain and draw dimensions
            if var_name == 'D':
                # D has shape (chains, draws, n_bins)
                n_chains, n_draws, n_bins = var_data.shape
                flattened = var_data.reshape(-1, n_bins)
                for i in range(n_bins):
                    samples_dict[f'D[{i}]'] = flattened[:, i]
            else:
                # gamma and sigma are scalars
                samples_dict[var_name] = var_data.flatten()
    
    df = pd.DataFrame(samples_dict)
    
    # Save to CSV
    df.to_csv(output_path, index=False)
    
    print(f"Saved {len(df)} samples")


def plot_trace(trace, output_path):
    """Plot MCMC trace plots for diagnostics."""
    print(f"Generating trace plot: {output_path}")
    
    az.plot_trace(trace, var_names=['D', 'gamma', 'sigma'])
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def print_summary_statistics(trace):
    """Print summary statistics of posterior distributions."""
    print("\n=== Posterior Summary Statistics ===")
    summary = az.summary(trace, var_names=['D', 'gamma', 'sigma'])
    print(summary)
    
    # Print half-lives
    D_means = summary.loc[['D[0]', 'D[1]', 'D[2]', 'D[3]', 'D[4]'], 'mean'].values
    print("\n=== mRNA Half-lives by Spatial Bin ===")
    print("Each bin has one degradation rate applied to 5 transcription traces:")
    for i, D_mean in enumerate(D_means):
        halflife = np.log(2) / D_mean
        print(f"  Bin {i+1}: t_1/2 = {halflife:.2f} min (D = {D_mean:.3f} min⁻¹)")


def load_data(transcription_path, mrna_path):
    """
    Load transcription and mRNA data.
    
    Args:
        transcription_path: Path to transcription data CSV (no headers)
        mrna_path: Path to mRNA data CSV (no headers)
    
    Returns:
        F_data: transcription data array [n_bins * n_traces_per_bin, n_timepoints]
        m_data: mRNA data array [n_bins]
    """
    print(f"Loading transcription data: {transcription_path}")
    F_data = pd.read_csv(transcription_path, header=None).values
    
    print(f"Loading mRNA data: {mrna_path}")
    m_data = pd.read_csv(mrna_path, header=None).values.flatten()
    
    print(f"Transcription data shape: {F_data.shape}")
    print(f"mRNA data shape: {m_data.shape}")
    
    # Data scaling for numerical stability
    print("Scaling data for numerical stability...")
    F_mean, F_std = F_data.mean(), F_data.std()
    m_mean, m_std = m_data.mean(), m_data.std()
    print(f"  F_data: mean={F_mean:.3f}, std={F_std:.3f}")
    print(f"  m_data: mean={m_mean:.3f}, std={m_std:.3f}")
    
    # Z-score normalization
    F_data = (F_data - F_mean) / F_std
    m_data = (m_data - m_mean) / m_std
    
    # Store scaling factors for potential denormalization (not used in model, but for reference)
    print("  Data normalized using z-score")
    
    return F_data, m_data


def organize_transcription_data(F_data, n_ap_bins=5, n_dv_bins=5):
    """
    Organize transcription data by AP position.
    
    The preprocessing creates a 5×5 grid where rows are ordered by spatial position.
    We group by AP position to infer AP-specific degradation rates.
    
    Assuming row order from preprocessing is:
    - Rows 0-4: AP position 1, DV positions 1-5
    - Rows 5-9: AP position 2, DV positions 1-5
    - Rows 10-14: AP position 3, DV positions 1-5
    - Rows 15-19: AP position 4, DV positions 1-5
    - Rows 20-24: AP position 5, DV positions 1-5
    
    Args:
        F_data: Transcription data [25 spatial bins, n_timepoints]
        n_ap_bins: Number of AP bins (default: 5)
        n_dv_bins: Number of DV bins per AP position (default: 5)
    
    Returns:
        List of arrays [n_ap_bins][n_dv_bins, n_timepoints]
        Each array contains the 5 DV traces for one AP position
    """
    total_bins = n_ap_bins * n_dv_bins
    
    if len(F_data) != total_bins:
        raise ValueError(
            f"Expected {total_bins} bins ({n_ap_bins} AP × {n_dv_bins} DV), "
            f"got {len(F_data)} rows"
        )
    
    F_data_arrays = []
    for ap_idx in range(n_ap_bins):
        start_idx = ap_idx * n_dv_bins
        end_idx = (ap_idx + 1) * n_dv_bins
        F_data_arrays.append(F_data[start_idx:end_idx, :])
    
    print(f"Organized transcription data:")
    print(f"  {n_ap_bins} AP positions")
    print(f"  {n_dv_bins} DV bins per AP position")
    print(f"  Total: {total_bins} spatial bins")
    print(f"  Each bin has one averaged transcription trace")
    
    return F_data_arrays


def organize_mrna_data(m_data, n_ap_bins=5, n_dv_bins=5):
    """
    Validate mRNA data matches expected structure from preprocessing.
    
    The smFISH preprocessing produces exactly 25 values:
    - 5 AP bins × 5 DV bins = 25 spatial bins
    - Each value = average mRNA/nucleus in that spatial bin
    - Ordered same as transcription data (by AP then DV)
    
    Args:
        m_data: mRNA data [25] - one value per spatial bin
        n_ap_bins: Number of AP bins (default: 5)
        n_dv_bins: Number of DV bins per AP position (default: 5)
    
    Returns:
        mRNA data unchanged [25]
    
    Raises:
        ValueError if data size doesn't match 25 bins
    """
    total_bins = n_ap_bins * n_dv_bins
    
    if len(m_data) != total_bins:
        raise ValueError(
            f"Expected {total_bins} mRNA values ({n_ap_bins} AP × {n_dv_bins} DV), "
            f"got {len(m_data)}. Check that smFISH preprocessing uses same 5×5 binning "
            f"as transcription preprocessing."
        )
    
    print(f"mRNA data validated:")
    print(f"  {len(m_data)} values (one per spatial bin)")
    print(f"  {n_ap_bins} AP positions × {n_dv_bins} DV bins")
    print(f"  Each value = average mRNA/nucleus in that bin")
    
    return m_data


def main():
    parser = argparse.ArgumentParser(
        description='Infer spatially-varying mRNA degradation rates'
    )
    parser.add_argument( '--transcription', required=True, help='Path to transcription data CSV (no headers)')
    parser.add_argument( '--mrna', required=True, help='Path to mRNA count data CSV (no headers)')
    parser.add_argument( '--output-chain', required=True, help='Output path for MCMC chain CSV')
    parser.add_argument( '--output-trace', required=True, help='Output path for trace plot')
    parser.add_argument( '--n-samples', type=int, default=10000, help='Number of MCMC samples per chain (default: 10000)')
    parser.add_argument( '--n-chains', type=int, default=4, help='Number of parallel MCMC chains (default: 4)')
    parser.add_argument( '--n-ap-bins', type=int, default=5, help='Number of anterior-posterior bins (default: 5)')
    parser.add_argument( '--n-dv-bins', type=int, default=5, help='Number of dorsal-ventral bins (default: 5)')

    args = parser.parse_args()
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Load data
    F_data, m_data_raw = load_data(args.transcription, args.mrna)
    
    # Organize data by AP position
    # Structure: 5 AP bins, each containing 5 DV transcription traces and 5 DV mRNA values
    # Model infers 5 degradation rates (one per AP position)
    # Each D[i] applies to all 5 DV bins at that AP position
    F_data_arrays = organize_transcription_data(F_data, args.n_ap_bins, args.n_dv_bins)
    m_data = organize_mrna_data(m_data_raw, args.n_ap_bins, args.n_dv_bins)
    
    # Time array (0 to 1200 seconds, 20 second intervals)
    t_array = np.arange(0, 1201, 20)
    
    # Build Bayesian model
    model = build_pymc_model(F_data_arrays, m_data, t_array, args.n_ap_bins)
    
    # Run MCMC inference
    trace = run_mcmc_inference(model, args.n_samples, args.n_chains)
    
    # Save results
    save_chain_results(trace, args.output_chain)
    plot_trace(trace, args.output_trace)
    
    # Print summary
    print_summary_statistics(trace)
    
    # Convergence diagnostics
    print("\n=== Convergence Diagnostics ===")
    rhat = az.rhat(trace, var_names=['D', 'gamma', 'sigma'])
    print("R-hat values (should be < 1.01 for convergence):")
    print(rhat)
    
    ess = az.ess(trace, var_names=['D', 'gamma', 'sigma'])
    print("\nEffective sample sizes:")
    print(ess)
    
    print("\n=== Inference complete ===")


if __name__ == '__main__':
    main()
