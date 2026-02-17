#!/usr/bin/env python
"""
Infer constant (null model) mRNA degradation rate using Bayesian inference.

This script implements a NULL MODEL where the degradation rate is CONSTANT across
all molecular ages and all spatial positions. This serves as a baseline for comparison
against more complex models (age-dependent and spatial models).

Mathematical formulation:
    - D(age) = D₀ (constant for all ages)
    - S(age) = exp(-D₀ · age) = exponential survival
    - m(T) = γ∫[0 to T] F(t)·exp(-D₀·(T-t))dt = convolution with exponential kernel

Key insight: This is the simplest possible model - all mRNA molecules decay at the
same constant rate regardless of when they were transcribed or where they are located.

The model is fit using PyMC for MCMC sampling with NUTS.

Reference: plans/compare-age-null-spatial-models-plan.md (Phase 1)
"""

import os
import sys
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
import pytensor.tensor as pt
import argparse


def solve_convolution_constant_D(D0, gamma, F_values, t_array):
    """
    Analytical solution for constant degradation using convolution.
    
    Calculates m(T) = gamma * Integral( F(t) * exp(-D₀*(T-t)) dt )
    
    Where the survival probability is simply: S(age) = exp(-D₀ * age)
    This is equivalent to exponential decay with constant rate D₀.
    
    Args:
        D0: Constant degradation rate (scalar)
        gamma: Transcription scaling factor (scalar)
        F_values: Transcription values at time points (array)
        t_array: Time points (array)
    
    Returns:
        mRNA concentration at final time point
    """
    # Ensure inputs are numpy arrays for numerical computation
    D0 = np.asarray(D0, dtype=np.float64)
    gamma = np.asarray(gamma, dtype=np.float64)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)
    
    # Clip D to reasonable range to avoid numerical issues
    D0 = np.clip(D0, 1e-6, 1.0)
    
    # Calculate survival curve: S(age) = exp(-D₀ * age)
    t_final = t_array[-1]
    # Age of each mRNA transcribed at t when measured at T:
    # mRNA from t=0 has age T, mRNA from t=T has age 0
    ages = t_final - t_array
    survival_profile = np.exp(-D0 * ages)
    
    # Align with transcription: reverse to match F(t) ordering
    survival_profile_reversed = survival_profile[::-1]
    
    # Convolve (integrate product)
    integrand = F_values * survival_profile_reversed
    integral = trapezoid(integrand, t_array)
    
    # Final solution: m(T) = γ * integral
    m_final = gamma * integral
    
    # Ensure output is valid
    if not np.isfinite(m_final):
        m_final = 0.0
    
    return float(m_final)


def calculate_expected_mRNA(D0, gamma, F_data_arrays, t_array):
    """
    Calculate expected mRNA for all spatial bins using CONSTANT degradation.
    
    Unlike spatial or age-dependent models, here D₀ is a single scalar value
    shared by ALL observations (all spatial positions and all ages).
    
    Args:
        D0: Single constant degradation rate (scalar)
        gamma: Transcription scaling factor (scalar)
        F_data_arrays: List of transcription data arrays for each bin [n_bins][n_traces_per_bin, n_timepoints]
        t_array: Time points array
    
    Returns:
        Array of expected mRNA concentrations [n_bins * n_traces_per_bin]
    """
    n_bins = len(F_data_arrays)
    expected_m_all = []
    
    for i in range(n_bins):
        # For each transcription trace, calculate expected mRNA
        # using the same constant degradation rate
        n_traces = F_data_arrays[i].shape[0]
        for j in range(n_traces):
            F_trace = F_data_arrays[i][j, :]
            m_expected = solve_convolution_constant_D(D0, gamma, F_trace, t_array)
            expected_m_all.append(m_expected)
    
    return np.array(expected_m_all)


def build_pymc_model(F_data_arrays, m_data, t_array, n_bins=5):
    """
    Build PyMC model for NULL CONSTANT degradation.
    
    D₀ is a single scalar shared by all observations. This is the simplest
    possible model and serves as a baseline for comparison.
    
    Mathematical formulation:
        - D₀ = constant degradation rate for all ages and positions
        - S(age) = exp(-D₀ · age) = exponential survival
        - m(T) = γ∫[0 to T] F(t)·exp(-D₀·(T-t))dt
    
    Prior distributions (matching spatial model conventions):
        D0 ~ HalfNormal(1.0) - constant degradation rate
        γ ~ InverseGamma(2, 3) - transcription scaling
        σ ~ InverseGamma(2, 3) - observation noise
    
    Likelihood:
        m_obs ~ Normal(expected_m_all, σ²)
    
    Args:
        F_data_arrays: List of transcription data arrays [n_bins][n_traces_per_bin, n_timepoints]
        m_data: Observed mRNA data [n_bins * n_traces_per_bin] = 25 values
        t_array: Time points
        n_bins: Number of spatial bins
    
    Returns:
        PyMC model
    """
    print("Building PyMC model with Constant Degradation (Null Model)...")
    
    # Get structure info
    n_traces_per_bin = F_data_arrays[0].shape[0]
    n_timepoints = len(t_array)
    total_observations = n_bins * n_traces_per_bin
    
    print(f"  {n_bins} spatial bins")
    print(f"  {n_traces_per_bin} traces per bin")
    print(f"  {total_observations} total observations")
    print(f"  Model: Single constant D₀ shared by all observations")
    
    # Pre-calculate trapezoidal weights for integration
    dt = t_array[1] - t_array[0]
    weights = np.ones_like(t_array)
    weights[0] = 0.5
    weights[-1] = 0.5
    weights_scaled = pt.as_tensor_variable(weights * dt)
    
    with pm.Model() as model:
        # --- Priors ---
        
        # D0: Single constant degradation rate
        # Using HalfNormal(1.0) as suggested in plan
        D0 = pm.HalfNormal('D0', sigma=1.0)
        
        # Store as D[0] for compatibility with comparison scripts
        pm.Deterministic('D', pt.stack([D0]))
        
        # Scalar parameters matching spatial model conventions
        gamma = pm.InverseGamma('gamma', alpha=2, beta=3)
        sigma = pm.InverseGamma('sigma', alpha=2, beta=3)
        
        # --- Calculate Expected mRNA using Vectorized Convolution ---
        
        # Survival probability: S(age) = exp(-D₀ * age)
        t_final = t_array[-1]
        ages = t_final - t_array  # Age of mRNA transcribed at each timepoint
        survival_prob = pm.math.exp(-D0 * ages)
        survival_profile_reversed = survival_prob[::-1]
        
        # Expected mRNA for each observation
        expected_m_list = []
        
        for i in range(n_bins):
            n_traces = F_data_arrays[i].shape[0]
            for j in range(n_traces):
                F_trace = F_data_arrays[i][j, :]
                # Convolution: integrate F(t) * S(T-t)
                integrand = F_trace * survival_profile_reversed
                integral = pm.math.sum(integrand * weights_scaled)
                m_ij = gamma * integral
                expected_m_list.append(m_ij)
        
        expected_m = pm.Deterministic('expected_m', pt.stack(expected_m_list))
        
        # --- Likelihood ---
        m_obs = pm.Normal('m_obs', mu=expected_m, sigma=sigma, observed=m_data)
    
    return model


def run_mcmc_inference(model, n_samples=2000, n_chains=4, target_accept=0.9):
    """
    Run MCMC sampling using NUTS algorithm.
    
    Args:
        model: PyMC model
        n_samples: Number of samples per chain
        n_chains: Number of parallel chains
        target_accept: Target acceptance rate for NUTS
    
    Returns:
        InferenceData object with samples
    """
    print(f"\nRunning MCMC inference: {n_samples} samples × {n_chains} chains...")
    print(f"Target accept: {target_accept}")
    
    with model:
        # Sample using NUTS (No-U-Turn Sampler)
        trace = pm.sample(
            draws=n_samples,
            chains=n_chains,
            cores=n_chains,
            target_accept=target_accept,
            return_inferencedata=True,
            random_seed=42
        )
    
    return trace


def save_chain_results(trace, output_path):
    """
    Save MCMC chain samples to CSV.
    
    Output format for null constant model:
        - D[0]: The single constant degradation rate
        - gamma: Transcription scaling factor
        - sigma: Observation noise
    
    Args:
        trace: InferenceData object
        output_path: Path to save CSV (will create parent directories)
    """
    print(f"\nSaving chain results to: {output_path}")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    posterior = trace.posterior
    samples_dict = {}
    
    # Extract D (stored as Deterministic with shape [1])
    # We want to save it as D[0] for compatibility with comparison scripts
    if 'D' in posterior:
        var_data = posterior['D'].values  # Shape: (chains, draws, 1)
        n_chains, n_draws, _ = var_data.shape
        flattened = var_data.reshape(-1, 1)
        samples_dict['D[0]'] = flattened[:, 0]
    elif 'D0' in posterior:
        # Fallback: extract D0 directly
        var_data = posterior['D0'].values  # Shape: (chains, draws)
        samples_dict['D[0]'] = var_data.flatten()
    
    # Extract scalar parameters
    for var_name in ['gamma', 'sigma']:
        if var_name in posterior:
            var_data = posterior[var_name].values
            samples_dict[var_name] = var_data.flatten()
    
    df = pd.DataFrame(samples_dict)
    
    # Save to CSV
    df.to_csv(output_path, index=False)
    
    print(f"Saved {len(df)} samples")
    print(f"Columns: {', '.join(df.columns.tolist())}")


def plot_trace(trace, output_path):
    """
    Plot MCMC trace plots for diagnostics.
    
    Args:
        trace: InferenceData object
        output_path: Path to save plot
    """
    print(f"\nGenerating trace plot: {output_path}")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Plot D0 (the underlying parameter) instead of D (deterministic)
    var_names = ['D0', 'gamma', 'sigma']
    
    az.plot_trace(trace, var_names=var_names)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Trace plot saved")


def print_summary_statistics(trace):
    """
    Print summary statistics of posterior distributions.
    
    Args:
        trace: InferenceData object
    """
    print("\n" + "="*60)
    print("POSTERIOR SUMMARY STATISTICS")
    print("="*60)
    
    summary = az.summary(trace, var_names=['D0', 'gamma', 'sigma'])
    print(summary)
    
    # Print interpretation
    posterior = trace.posterior
    D0_samples = posterior['D0'].values.flatten()
    D0_mean = D0_samples.mean()
    D0_std = D0_samples.std()
    
    # Convert to half-life for interpretation
    halflife_mean = np.log(2) / D0_mean if D0_mean > 0 else np.inf
    
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)
    print(f"Constant Degradation Rate D₀:")
    print(f"  Mean: {D0_mean:.4f} min⁻¹ (± {D0_std:.4f})")
    print(f"  Half-life: {halflife_mean:.2f} minutes")
    print("\nThis null model assumes all mRNA molecules decay at the same")
    print("constant rate regardless of age or spatial position.")
    print("="*60)


def load_data(transcription_path, mrna_path):
    """
    Load transcription and mRNA data.
    
    Args:
        transcription_path: Path to transcription data CSV (no headers)
        mrna_path: Path to mRNA data CSV (no headers)
    
    Returns:
        F_data: Transcription data array [n_bins * n_traces_per_bin, n_timepoints]
        m_data: mRNA data array [n_bins * n_traces_per_bin]
    """
    print(f"Loading transcription data: {transcription_path}")
    F_data = pd.read_csv(transcription_path, header=None).values
    
    print(f"Loading mRNA data: {mrna_path}")
    m_data = pd.read_csv(mrna_path, header=None).values.flatten()
    
    print(f"Transcription data shape: {F_data.shape}")
    print(f"mRNA data shape: {m_data.shape}")
    
    return F_data, m_data


def organize_transcription_data(F_data, n_ap_bins=5, n_dv_bins=5):
    """
    Organize transcription data by spatial bins (for compatibility with other models).
    
    Even though the null model doesn't use spatial structure, we maintain the same
    data organization as other models for consistency.
    
    Args:
        F_data: Transcription data [25 spatial bins, n_timepoints]
        n_ap_bins: Number of AP bins (default: 5)
        n_dv_bins: Number of DV bins per AP position (default: 5)
    
    Returns:
        List of arrays [n_ap_bins][n_dv_bins, n_timepoints]
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
    
    print(f"\nOrganized transcription data:")
    print(f"  {n_ap_bins} AP positions × {n_dv_bins} DV bins = {total_bins} total bins")
    print(f"  (Null model uses single D₀ for all bins)")
    
    return F_data_arrays


def organize_mrna_data(m_data, n_ap_bins=5, n_dv_bins=5):
    """
    Validate mRNA data matches expected structure.
    
    Args:
        m_data: mRNA data [25] - one value per spatial bin
        n_ap_bins: Number of AP bins (default: 5)
        n_dv_bins: Number of DV bins per AP position (default: 5)
    
    Returns:
        mRNA data unchanged
    """
    total_bins = n_ap_bins * n_dv_bins
    
    if len(m_data) != total_bins:
        raise ValueError(
            f"Expected {total_bins} mRNA values ({n_ap_bins} AP × {n_dv_bins} DV), "
            f"got {len(m_data)}"
        )
    
    print(f"\nmRNA data validated: {len(m_data)} observations")
    
    return m_data


def main():
    parser = argparse.ArgumentParser(
        description='Infer constant (null model) mRNA degradation rate',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default settings
  python %(prog)s \\
    --transcription data/trans.csv \\
    --mrna data/mrna.csv \\
    --output-chain results/chains/degradation_chain.csv \\
    --output-trace results/trace.png

  # Quick test run (fewer samples)
  python %(prog)s \\
    --transcription data/trans.csv \\
    --mrna data/mrna.csv \\
    --output-chain results/chains/degradation_chain.csv \\
    --output-trace results/trace.png \\
    --n-samples 500 \\
    --n-chains 2
"""
    )
    
    # Required arguments
    parser.add_argument(
        '--transcription', 
        required=True, 
        help='Path to transcription data CSV (no headers)'
    )
    parser.add_argument(
        '--mrna', 
        required=True, 
        help='Path to mRNA count data CSV (no headers)'
    )
    parser.add_argument(
        '--output-chain', 
        required=True, 
        help='Output path for MCMC chain CSV (will create directory)'
    )
    parser.add_argument(
        '--output-trace', 
        required=True, 
        help='Output path for trace plot PNG'
    )
    
    # Sampling arguments
    parser.add_argument(
        '--n-samples', 
        type=int, 
        default=2000, 
        help='Number of MCMC samples per chain (default: 2000)'
    )
    parser.add_argument(
        '--n-chains', 
        type=int, 
        default=4, 
        help='Number of parallel MCMC chains (default: 4)'
    )
    parser.add_argument(
        '--target-accept', 
        type=float, 
        default=0.9, 
        help='Target acceptance rate for NUTS (default: 0.9)'
    )
    
    # Binning arguments (for data organization compatibility)
    parser.add_argument(
        '--n-ap-bins', 
        type=int, 
        default=5, 
        help='Number of anterior-posterior bins (default: 5)'
    )
    parser.add_argument(
        '--n-dv-bins', 
        type=int, 
        default=5, 
        help='Number of dorsal-ventral bins (default: 5)'
    )

    args = parser.parse_args()
    
    print("="*60)
    print("NULL CONSTANT DEGRADATION MODEL INFERENCE")
    print("="*60)
    print("Model: D(age) = D₀ (constant for all ages and positions)")
    print("="*60)
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Load data
    F_data, m_data_raw = load_data(args.transcription, args.mrna)
    
    # Organize data (maintains compatibility with other models)
    F_data_arrays = organize_transcription_data(F_data, args.n_ap_bins, args.n_dv_bins)
    m_data = organize_mrna_data(m_data_raw, args.n_ap_bins, args.n_dv_bins)
    
    # Time array in minutes (assuming 20-second intervals, 0-20 minutes)
    n_timepoints = F_data.shape[1]
    t_array = np.arange(0, n_timepoints * 20, 20) / 60.0
    print(f"\nTime array: {n_timepoints} points from {t_array[0]:.2f} to {t_array[-1]:.2f} minutes")
    
    # Build Bayesian model
    model = build_pymc_model(F_data_arrays, m_data, t_array, args.n_ap_bins)
    
    # Run MCMC inference
    trace = run_mcmc_inference(model, args.n_samples, args.n_chains, args.target_accept)
    
    # Save results
    save_chain_results(trace, args.output_chain)
    plot_trace(trace, args.output_trace)
    
    # Print summary
    print_summary_statistics(trace)
    
    print("\n" + "="*60)
    print("INFERENCE COMPLETE")
    print("="*60)
    print(f"Chain CSV: {args.output_chain}")
    print(f"Trace plot: {args.output_trace}")
    print("="*60)


if __name__ == '__main__':
    main()
