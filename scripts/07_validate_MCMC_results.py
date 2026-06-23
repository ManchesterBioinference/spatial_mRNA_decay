#!/usr/bin/env python3
"""
Posterior predictive check for MCMC degradation rate inference.

Validates inferred degradation rates by comparing model-predicted mRNA counts
against observed smFISH data. Supports both the spatial (constant-D per AP bin)
and age-dependent (D(age)) model variants.

The predicted mRNA m̂ is computed by re-running the analytical ODE solution using
posterior mean D and gamma, then compared visually to observed mRNA counts.

Usage:
    python 07_validate_MCMC_results.py \\
        --chain results/stripe2/e6/chain.csv \\
        --transcription data/processed_transcription_data/transcription_traces_no_ids_stripe2_1200.csv \\
        --mrna results/data/processed_mRNA_data_stripe2/e6_sass_formodel.csv \\
        --n-ap-bins 5 --n-dv-bins 5 \\
        --output results/stripe2/e6/posterior_predictive_check.pdf
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.integrate import cumulative_trapezoid, trapezoid

def solve_spatial(D_scalar, gamma, F_values, t_array):
    """
    Analytical solution for spatial (constant-D) model: dm/dt = γ*F(t) - D*m.

    Matches solve_ode_analytical in 02_infer_degradation_rates_spatial.py:
        m(T) = γ * ∫ F(s) * exp(D*(s - T)) ds

    Args:
        D_scalar: constant degradation rate for this AP bin (scalar)
        gamma: transcription scaling factor (scalar)
        F_values: transcription values at time points (array, length n_timepoints)
        t_array: time points (array, length n_timepoints)

    Returns:
        mRNA concentration at final time point (scalar)
    """
    D_scalar = float(np.clip(D_scalar, 1e-6, 10.0))
    gamma = float(gamma)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)

    t_final = t_array[-1]
    exp_term = np.exp(D_scalar * (t_array - t_final))
    integral = trapezoid(F_values * exp_term, t_array)
    return gamma * integral


def solve_analytical(D_age, gamma, F_values, t_array):
    """
    Re-implementation of age-dependent degradation model for validation.
    
    MUST match the numerical integration in the inference script exactly:
    - Uses cumsum(D_age) * dt for cumulative hazard (rectangular rule)
    - Uses trapezoidal rule for final convolution integral
    
    m(T) = gamma * Integral( F(t) * S(T-t) dt )
    where S(age) = exp(-Integral D(tau) dtau)
    """
    # Ensure arrays
    D_age = np.asarray(D_age, dtype=np.float64)
    gamma = float(gamma)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)
    
    # Calculate survival curve - MUST match inference script
    dt = t_array[1] - t_array[0]
    cumulative_hazard = np.cumsum(D_age) * dt  # Rectangular rule, same as PyMC
    survival_prob = np.exp(-cumulative_hazard)
    
    # Reverse to align with transcription times
    # mRNA made at t=0 has age=T, mRNA made at t=T has age=0
    survival_profile_reversed = survival_prob[::-1]
    
    # Convolve using trapezoidal rule (matches weights in inference script)
    integrand = F_values * survival_profile_reversed
    integral = trapezoid(integrand, t_array)
    
    return gamma * integral

def main():
    parser = argparse.ArgumentParser(
        description='Validate age-dependent degradation MCMC results with posterior predictive check'
    )
    parser.add_argument('--chain', required=True, help='MCMC chain CSV')
    parser.add_argument('--transcription', required=True, help='Original transcription CSV')
    parser.add_argument('--mrna', required=True, help='Original mRNA CSV')
    parser.add_argument('--n-ap-bins', type=int, required=True, help='Number of AP bins')
    parser.add_argument('--n-dv-bins', type=int, required=True, help='Number of DV bins')
    parser.add_argument('--output', default='posterior_predictive_check.pdf')
    args = parser.parse_args()

    n_ap_bins = args.n_ap_bins
    n_dv_bins = args.n_dv_bins

    # 1. Load Data
    df_samples = pd.read_csv(args.chain)
    F_data = pd.read_csv(args.transcription, header=None).values
    m_obs = pd.read_csv(args.mrna, header=None).values.flatten()
    
    # Derive time array from transcription data shape (must match inference script)
    n_timepoints = F_data.shape[1]
    t_array = np.arange(0, n_timepoints * 20, 20) / 60.0

    # 2. Extract D_age posterior means
    # For age-dependent model, D columns represent degradation at different molecular ages
    D_columns = [col for col in df_samples.columns if col.startswith('D[') and col.endswith(']')]
    D_columns_sorted = sorted(D_columns, key=lambda x: int(x.split('[')[1].split(']')[0]))
    
    n_ages = len(D_columns_sorted)
    print(f"Found {n_ages} age bins in chain")
    
    # Get posterior median for each age
    D_age_medians_raw = np.array([np.median(df_samples[col]) for col in D_columns_sorted])

    # Detect model type from number of D columns
    is_spatial_model = False
    if n_ages == 1:
        # Null constant model: single D → expand to constant array over time
        print("Null constant model detected (single D[0]): expanding to constant D(age) array")
        D_age_medians = np.full(n_timepoints, D_age_medians_raw[0])
    elif n_ages == n_ap_bins:
        # Spatial model: one constant D per AP bin
        print(f"Spatial model detected ({n_ages} D values = {n_ap_bins} AP bins): using constant D per bin")
        is_spatial_model = True
        D_spatial_medians = D_age_medians_raw  # shape (n_ap_bins,)
        D_age_medians = D_age_medians_raw  # kept for printing below
    else:
        D_age_medians = D_age_medians_raw

    gamma_median = float(np.median(df_samples['gamma']))
    
    print(f"Median gamma: {gamma_median:.2f}")
    print(f"Median D(age=0): {D_age_medians[0]:.3f} min⁻¹")
    print(f"Median D(age=max): {D_age_medians[-1]:.3f} min⁻¹")
    
    # Check if this is the mechanistic poly-A protection model
    polya_params = ['NA_0', 'deadenylation_rate', 'beta', 'D_protected', 'D_unprotected']
    is_polya_model = all(param in df_samples.columns for param in polya_params)
    
    if is_polya_model:
        print("\n" + "="*60)
        print("MECHANISTIC POLY-A PROTECTION MODEL DETECTED")
        print("="*60)
        print("\nInferred Biological Parameters:")
        print(f"  Initial poly-A tail length (NA_0):")
        print(f"    Mean: {df_samples['NA_0'].mean():.1f} adenosines")
        print(f"    Std:  {df_samples['NA_0'].std():.1f} adenosines")
        print(f"    95% CI: [{np.percentile(df_samples['NA_0'], 2.5):.1f}, {np.percentile(df_samples['NA_0'], 97.5):.1f}]")
        
        print(f"\n  Deadenylation rate:")
        print(f"    Mean: {df_samples['deadenylation_rate'].mean():.3f} As/min")
        print(f"    Std:  {df_samples['deadenylation_rate'].std():.3f} As/min")
        print(f"    95% CI: [{np.percentile(df_samples['deadenylation_rate'], 2.5):.3f}, {np.percentile(df_samples['deadenylation_rate'], 97.5):.3f}]")
        
        print(f"\n  Protection sharpness (β):")
        print(f"    Mean: {df_samples['beta'].mean():.4f}")
        print(f"    Std:  {df_samples['beta'].std():.4f}")
        print(f"    95% CI: [{np.percentile(df_samples['beta'], 2.5):.4f}, {np.percentile(df_samples['beta'], 97.5):.4f}]")
        print(f"    (Literature value: 0.096)")
        
        print(f"\n  Protected degradation rate (D_protected):")
        print(f"    Mean: {df_samples['D_protected'].mean():.4f} min⁻¹")
        print(f"    Std:  {df_samples['D_protected'].std():.4f} min⁻¹")
        print(f"    Half-life: {np.log(2) / df_samples['D_protected'].mean():.1f} min")
        
        print(f"\n  Unprotected degradation rate (D_unprotected):")
        print(f"    Mean: {df_samples['D_unprotected'].mean():.4f} min⁻¹")
        print(f"    Std:  {df_samples['D_unprotected'].std():.4f} min⁻¹")
        print(f"    Half-life: {np.log(2) / df_samples['D_unprotected'].mean():.1f} min")
        
        # Calculate time when poly-A tail reaches critical length (~20 As)
        na0_mean = df_samples['NA_0'].mean()
        deadenyl_mean = df_samples['deadenylation_rate'].mean()
        time_to_critical = (na0_mean - 20) / deadenyl_mean if deadenyl_mean > 0 else np.inf
        print(f"\n  Time to critical poly-A length (~20 As): {time_to_critical:.1f} min")
        print("="*60)
    else:
        print("\n(Delayed age-dependent model - no mechanistic poly-A parameters)")

    # 3. Predict mRNA for every trace
    m_pred = []
    for i in range(n_ap_bins):  # AP spatial bins
        for j in range(n_dv_bins):  # DV traces per bin
            idx = i * n_dv_bins + j
            F_trace = F_data[idx, :]
            if is_spatial_model:
                # Spatial model: use constant D[i] for this AP bin
                val = solve_spatial(D_spatial_medians[i], gamma_median, F_trace, t_array)
            else:
                # Age-dependent or null model: D(age) curve shared across all bins
                val = solve_analytical(D_age_medians, gamma_median, F_trace, t_array)
            m_pred.append(val)
    
    m_pred = np.array(m_pred)

    # 4. Calculate fit statistics
    residuals = m_obs - m_pred
    rmse = np.sqrt(np.mean(residuals**2))
    r_squared = 1 - (np.sum(residuals**2) / np.sum((m_obs - np.mean(m_obs))**2))
    
    print(f"\nModel fit statistics:")
    print(f"  RMSE: {rmse:.2f}")
    print(f"  R²: {r_squared:.3f}")

    # 5. Save figure data CSV
    # Derive figure_data/ directory from output path (e.g. results_1200/stripe2/e8_9um/... → results_1200/stripe2/e8_9um/figure_data/)
    output_dir = os.path.dirname(args.output)
    figure_data_dir = os.path.join(output_dir, 'figure_data')
    os.makedirs(figure_data_dir, exist_ok=True)

    # Build dataframe with per-observation data
    bin_indices = np.repeat(np.arange(n_ap_bins), n_dv_bins)
    df_figure_data = pd.DataFrame({
        'observed_mRNA': m_obs,
        'predicted_mRNA': m_pred,
        'residuals': residuals,
        'spatial_bin': bin_indices + 1,
    })

    figure_data_path = os.path.join(figure_data_dir, 'posterior_predictive_check.csv')
    df_figure_data.to_csv(figure_data_path, index=False)
    print(f"Figure data saved to {figure_data_path}")

    # 6. Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left panel: Posterior Predictive Check
    max_val = max(max(m_obs), max(m_pred))
    ax1.plot([0, max_val], [0, max_val], color='red', linestyle='--', alpha=0.5, 
             linewidth=2, label='Perfect Fit')
    
    # Plot all bins (they share the same D_age, so color by spatial bin for reference)
    # Generate colors dynamically for any number of AP bins
    import matplotlib.cm as cm
    colors = cm.tab10(np.linspace(0, 1, max(n_ap_bins, 10)))
    
    for i in range(n_ap_bins):
        start, end = i * n_dv_bins, (i + 1) * n_dv_bins
        ax1.scatter(m_obs[start:end], m_pred[start:end], 
                    color=colors[i], label=f'Spatial Bin {i+1}', 
                    s=100, edgecolors='k', alpha=0.7)

    ax1.set_xlabel('Observed mRNA (Data)', fontsize=12)
    ax1.set_ylabel('Predicted mRNA (Model)', fontsize=12)
    ax1.set_title(f'Posterior Predictive Check\nR² = {r_squared:.3f}, RMSE = {rmse:.2f}', 
                  fontsize=13)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Right panel: Residuals
    bin_indices = np.repeat(np.arange(n_ap_bins), n_dv_bins)
    ax2.scatter(m_pred, residuals, c=bin_indices, cmap='viridis', 
                s=100, edgecolors='k', alpha=0.7)
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=2)
    ax2.set_xlabel('Predicted mRNA', fontsize=12)
    ax2.set_ylabel('Residuals (Observed - Predicted)', fontsize=12)
    ax2.set_title('Residual Plot', fontsize=13)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nValidation plot saved to {args.output}")

if __name__ == "__main__":
    main()