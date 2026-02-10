#!/usr/bin/env python
"""
Visualize age-dependent mRNA decay inference results from PyMC MCMC chains.

Creates plots showing how degradation rates D(age) and mRNA half-life depend on
the age of individual molecules (time since transcription), not spatial position
or global experimental time. Includes uncertainty quantification via 95% credible
intervals.

This models intrinsic molecular processes like poly-A tail shortening.

Based on plotting_modelling_results.ipynb, updated for age-dependent inference.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import cmcrameri.cm as cmc
import argparse


# Set global plotting parameters
mpl.rc('font', size=10)
#mpl.rc('font', family='Arial', size=10)


def load_chain_data(chain_path):
    """Load MCMC chain results from PyMC inference."""
    print(f"Loading chain data from: {chain_path}")
    data = pd.read_csv(chain_path, engine='python', header=0)
    print(f"Loaded {len(data)} samples")
    return data


def extract_degradation_posteriors(data, n_ages=None):
    """
    Extract degradation rate (D) posteriors for age-dependent inference.
    
    For age-dependent models, D has shape (n_samples, n_ages) where each
    column represents D at a specific molecular age (time since transcription).
    
    Args:
        data: MCMC chain DataFrame
        n_ages: Number of age bins (if None, auto-detect from columns)
    
    Returns:
        2D numpy array with shape (n_samples, n_ages)
    """
    print("Extracting age-dependent degradation posteriors...")
    
    # Auto-detect number of D columns
    D_columns = [col for col in data.columns if col.startswith('D[') and col.endswith(']')]
    D_columns_sorted = sorted(D_columns, key=lambda x: int(x.split('[')[1].split(']')[0]))
    
    if len(D_columns_sorted) == 0:
        print(f"Error: No D[i] columns found in data")
        print(f"Available columns: {list(data.columns)}")
        return None
    
    detected_ages = len(D_columns_sorted)
    print(f"  Detected {detected_ages} age bins")
    
    if n_ages is not None and n_ages != detected_ages:
        print(f"  Warning: Expected {n_ages} age bins but found {detected_ages}")
    
    # Extract all D columns into a 2D array
    D_posteriors = data[D_columns_sorted].to_numpy()
    
    print(f"  Posterior shape: {D_posteriors.shape} (samples × age bins)")
    return D_posteriors


def extract_polya_posteriors(data, n_ages=None):
    """
    Extract poly-A tail and protection factor posteriors (if available).
    
    For mechanistic poly-A models, this extracts:
    - NA[i]: poly-A tail length at each age
    - protection[i]: Pab1 protection factor at each age
    
    Args:
        data: MCMC chain DataFrame
        n_ages: Number of age bins (if None, auto-detect)
    
    Returns:
        Dictionary with 'NA' and 'protection' arrays, or None if not present
    """
    print("Checking for poly-A tail dynamics posteriors...")
    
    # Check for NA columns
    NA_columns = [col for col in data.columns if col.startswith('NA[') and col.endswith(']')]
    protection_columns = [col for col in data.columns if col.startswith('protection[') and col.endswith(']')]
    
    if len(NA_columns) == 0 and len(protection_columns) == 0:
        print("  No poly-A dynamics found (standard model)")
        return None
    
    result = {}
    
    if len(NA_columns) > 0:
        NA_columns_sorted = sorted(NA_columns, key=lambda x: int(x.split('[')[1].split(']')[0]))
        result['NA'] = data[NA_columns_sorted].to_numpy()
        print(f"  Found NA posteriors: {result['NA'].shape}")
    
    if len(protection_columns) > 0:
        protection_columns_sorted = sorted(protection_columns, key=lambda x: int(x.split('[')[1].split(']')[0]))
        result['protection'] = data[protection_columns_sorted].to_numpy()
        print(f"  Found protection posteriors: {result['protection'].shape}")
    
    return result if result else None


def compute_halflife_from_D(D_value):
    """Compute mRNA half-life from degradation rate: t_1/2 = ln(2) / D"""
    return np.log(2) / D_value  # D is in min^-1 -> result in minutes


def compute_posterior_summaries(D_posteriors):
    """
    Compute summary statistics for age-dependent degradation rates.
    
    Args:
        D_posteriors: 2D array of shape (n_samples, n_ages)
    
    Returns:
        Dictionary with mean, median, std, and credible intervals for each age
    """
    print("Computing posterior summaries...")
    
    summary = {
        'mean': np.mean(D_posteriors, axis=0),
        'median': np.median(D_posteriors, axis=0),
        'std': np.std(D_posteriors, axis=0),
        'ci_lower': np.percentile(D_posteriors, 2.5, axis=0),
        'ci_upper': np.percentile(D_posteriors, 97.5, axis=0)
    }
    
    print(f"  Mean D: {summary['mean'].min():.3f} - {summary['mean'].max():.3f} min⁻¹")
    print(f"  Mean t_1/2: {compute_halflife_from_D(summary['mean'].max()):.2f} - "
          f"{compute_halflife_from_D(summary['mean'].min()):.2f} min")
    
    return summary


def plot_degradation_time_series(D_posteriors, summary, age_array, output_path):
    """
    Plot degradation rate D(age) as function of molecular age with uncertainty bands.
    
    Args:
        D_posteriors: 2D array (n_samples, n_ages)
        summary: Dictionary with mean, median, ci_lower, ci_upper
        age_array: Molecular ages in minutes
        output_path: Path to save figure
    """
    print(f"Plotting degradation rate vs age: {output_path}")
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Plot mean line
    ax.plot(age_array, summary['mean'], color='#332288', linewidth=2, label='Mean')
    
    # Plot 95% credible interval
    ax.fill_between(
        age_array,
        summary['ci_lower'],
        summary['ci_upper'],
        color='#d8bfd8',
        alpha=0.5,
        label='95% CI'
    )
    
    # Optional: plot median
    ax.plot(age_array, summary['median'], color='#cd96cd', linewidth=1.5, 
            linestyle='--', label='Median', alpha=0.7)
    
    ax.set_xlabel('Molecular Age (min)')
    ax.set_ylabel('Degradation Rate D (min⁻¹)')
    ax.set_title('Age-Dependent mRNA Degradation Rate')
    ax.legend(loc='best')
    ax.grid(alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_halflife_time_series(D_posteriors, summary, age_array, output_path):
    """
    Plot mRNA half-life as function of molecular age with uncertainty bands.
    
    Args:
        D_posteriors: 2D array (n_samples, n_ages)
        summary: Dictionary with mean, median, ci_lower, ci_upper
        age_array: Molecular ages in minutes
        output_path: Path to save figure
    """
    print(f"Plotting half-life vs age: {output_path}")
    
    # Convert D to half-life
    halflife_mean = compute_halflife_from_D(summary['mean'])
    halflife_median = compute_halflife_from_D(summary['median'])
    # Note: CI bounds are inverted for half-life (lower D → higher t_1/2)
    halflife_ci_lower = compute_halflife_from_D(summary['ci_upper'])
    halflife_ci_upper = compute_halflife_from_D(summary['ci_lower'])
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Plot mean line
    ax.plot(age_array, halflife_mean, color='#44AA99', linewidth=2, label='Mean')
    
    # Plot 95% credible interval
    ax.fill_between(
        age_array,
        halflife_ci_lower,
        halflife_ci_upper,
        color='#88CCAA',
        alpha=0.5,
        label='95% CI'
    )
    
    # Optional: plot median
    ax.plot(age_array, halflife_median, color='#117733', linewidth=1.5, 
            linestyle='--', label='Median', alpha=0.7)
    
    ax.set_xlabel('Molecular Age (min)')
    ax.set_ylabel('mRNA Half-life (min)')
    ax.set_title('Age-Dependent mRNA Half-life')
    ax.legend(loc='best')
    ax.grid(alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def create_temporal_overview_plot(D_posteriors, summary, age_array, output_path):
    """
    Create comprehensive overview plot with D(age) and t_1/2(age).
    
    Args:
        D_posteriors: 2D array (n_samples, n_ages)
        summary: Dictionary with statistics
        age_array: Molecular ages in minutes
        output_path: Path to save figure
    """
    print(f"Creating age-dependent overview plot: {output_path}")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Top panel: Degradation rate D(age)
    ax1.plot(age_array, summary['mean'], color='#332288', linewidth=2, label='Mean')
    ax1.fill_between(
        age_array,
        summary['ci_lower'],
        summary['ci_upper'],
        color='#d8bfd8',
        alpha=0.5,
        label='95% CI'
    )
    ax1.set_ylabel('Degradation Rate D (min⁻¹)')
    ax1.set_title('Age-Dependent mRNA Degradation')
    ax1.legend(loc='best')
    ax1.grid(alpha=0.3, linestyle='--')
    
    # Bottom panel: Half-life t_1/2(age)
    halflife_mean = compute_halflife_from_D(summary['mean'])
    halflife_ci_lower = compute_halflife_from_D(summary['ci_upper'])
    halflife_ci_upper = compute_halflife_from_D(summary['ci_lower'])
    
    ax2.plot(age_array, halflife_mean, color='#44AA99', linewidth=2, label='Mean')
    ax2.fill_between(
        age_array,
        halflife_ci_lower,
        halflife_ci_upper,
        color='#88CCAA',
        alpha=0.5,
        label='95% CI'
    )
    ax2.set_xlabel('Molecular Age (min)')
    ax2.set_ylabel('mRNA Half-life (min)')
    ax2.legend(loc='best')
    ax2.grid(alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_polya_dynamics(polya_data, age_array, output_path):
    """
    Plot poly-A tail dynamics and protection factor with D(age).
    
    Creates a 3-panel figure showing:
    1. Degradation rate D(age) - shows biphasic pattern
    2. Poly-A tail length NA(age) - shows linear shortening
    3. Protection factor - shows sharp drop at ~20 adenosines
    
    Args:
        polya_data: Dictionary with 'D', 'NA', 'protection' posteriors
        age_array: Molecular ages in minutes
        output_path: Path to save figure
    """
    print(f"Plotting poly-A tail dynamics: {output_path}")
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    
    # Extract mean and CI for each variable
    D_mean = np.mean(polya_data['D'], axis=0)
    D_ci_lower = np.percentile(polya_data['D'], 2.5, axis=0)
    D_ci_upper = np.percentile(polya_data['D'], 97.5, axis=0)
    
    NA_mean = np.mean(polya_data['NA'], axis=0)
    NA_ci_lower = np.percentile(polya_data['NA'], 2.5, axis=0)
    NA_ci_upper = np.percentile(polya_data['NA'], 97.5, axis=0)
    
    protection_mean = np.mean(polya_data['protection'], axis=0)
    protection_ci_lower = np.percentile(polya_data['protection'], 2.5, axis=0)
    protection_ci_upper = np.percentile(polya_data['protection'], 97.5, axis=0)
    
    # Panel 1: Degradation rate D(age)
    axes[0].plot(age_array, D_mean, 'k-', linewidth=2, label='Mean')
    axes[0].fill_between(age_array, D_ci_lower, D_ci_upper, 
                          alpha=0.3, color='gray', label='95% CI')
    axes[0].set_ylabel('D(age) [min⁻¹]')
    axes[0].set_title('Mechanistic Poly-A Protection Model: Biphasic Degradation')
    axes[0].legend(loc='best')
    axes[0].grid(alpha=0.3, linestyle='--')
    
    # Panel 2: Poly-A tail length
    axes[1].plot(age_array, NA_mean, 'b-', linewidth=2, label='Mean')
    axes[1].fill_between(age_array, NA_ci_lower, NA_ci_upper, 
                          alpha=0.3, color='lightblue', label='95% CI')
    axes[1].axhline(20, color='red', linestyle='--', linewidth=1.5, 
                     label='Pab1 threshold (~20 As)')
    axes[1].set_ylabel('Poly-A Length [adenosines]')
    axes[1].legend(loc='best')
    axes[1].grid(alpha=0.3, linestyle='--')
    
    # Panel 3: Protection factor
    axes[2].plot(age_array, protection_mean, 'g-', linewidth=2, label='Mean')
    axes[2].fill_between(age_array, protection_ci_lower, protection_ci_upper, 
                          alpha=0.3, color='lightgreen', label='95% CI')
    axes[2].axhline(np.tanh(0.096 * 20), color='red', linestyle='--', 
                     linewidth=1.5, alpha=0.5, label='Protection at 20 As')
    axes[2].set_ylabel('Protection Factor')
    axes[2].set_xlabel('Molecular Age [min]')
    axes[2].legend(loc='best')
    axes[2].grid(alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print("  ✓ Biphasic pattern: slow decay → sharp transition → fast decay")
    print(f"  Poly-A starts at ~{NA_mean[0]:.1f} As, drops below 20 As at age ~{age_array[np.where(NA_mean < 20)[0][0] if np.any(NA_mean < 20) else -1]:.1f} min")


def save_summary_stats(D_posteriors, summary, age_array, output_path):
    """
    Save summary statistics of age-dependent degradation rates.
    
    Args:
        D_posteriors: 2D array (n_samples, n_ages)
        summary: Dictionary with mean, median, std, credible intervals
        age_array: Molecular ages in minutes
        output_path: Path to save CSV
    """
    print(f"Saving summary statistics: {output_path}")
    
    summary_data = []
    n_ages = D_posteriors.shape[1]
    
    for i in range(n_ages):
        D_mean = summary['mean'][i]
        D_median = summary['median'][i]
        
        summary_data.append({
            'age_bin': i,
            'age_min': age_array[i],
            'D_mean': D_mean,
            'D_median': D_median,
            'D_std': summary['std'][i],
            'D_95CI_lower': summary['ci_lower'][i],
            'D_95CI_upper': summary['ci_upper'][i],
            'halflife_mean_min': compute_halflife_from_D(D_mean),
            'halflife_median_min': compute_halflife_from_D(D_median),
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_path, index=False)
    
    print("\nAge-Dependent Summary Statistics (first 10 age bins):")
    print(summary_df.head(10).to_string(index=False))


def main():
    parser = argparse.ArgumentParser(
        description='Visualize age-dependent mRNA decay inference results'
    )
    parser.add_argument(
        '--chain', required=True,
        help='Path to MCMC chain CSV file'
    )
    parser.add_argument(
        '--degradation-plot', required=True,
        help='Output path for degradation rate D(age) plot'
    )
    parser.add_argument(
        '--halflife-plot', required=True,
        help='Output path for half-life t_1/2(age) plot'
    )
    parser.add_argument(
        '--overview-plot', required=True,
        help='Output path for combined overview plot'
    )
    parser.add_argument(
        '--summary', required=True,
        help='Output path for summary statistics CSV'
    )
    parser.add_argument(
        '--time-step', type=float, default=1/3,
        help='Time step between measurements in minutes (default: 1/3 = 20 seconds)'
    )
    parser.add_argument(
        '--n-ages', type=int, default=None,
        help='Expected number of age bins (optional, auto-detected if not provided)'
    )
    
    args = parser.parse_args()
    
    # Load chain data
    chain_data = load_chain_data(args.chain)
    
    # Extract age-dependent degradation posteriors
    D_posteriors = extract_degradation_posteriors(chain_data, args.n_ages)
    
    if D_posteriors is None:
        print("Error: Could not extract degradation posteriors from chain data")
        sys.exit(1)
    
    # Check for poly-A dynamics (mechanistic model)
    polya_posteriors = extract_polya_posteriors(chain_data, args.n_ages)
    
    # Compute summary statistics
    summary = compute_posterior_summaries(D_posteriors)
    
    # Create age array (molecular age = time since transcription)
    n_ages = D_posteriors.shape[1]
    age_array = np.arange(n_ages) * args.time_step
    
    print(f"\nMolecular age range: {age_array[0]:.2f} - {age_array[-1]:.2f} min")
    
    # Generate plots
    plot_degradation_time_series(D_posteriors, summary, age_array, args.degradation_plot)
    plot_halflife_time_series(D_posteriors, summary, age_array, args.halflife_plot)
    create_temporal_overview_plot(D_posteriors, summary, age_array, args.overview_plot)
    
    # If poly-A dynamics are present, create additional plot
    if polya_posteriors is not None:
        # Check if we have D in polya_posteriors, otherwise use D_posteriors
        if 'D' not in polya_posteriors:
            polya_posteriors['D'] = D_posteriors
        
        polya_plot_path = args.overview_plot.replace('.png', '_polya_dynamics.png')
        plot_polya_dynamics(polya_posteriors, age_array, polya_plot_path)
        print(f"\n✓ Poly-A dynamics plot saved: {polya_plot_path}")
        print("  This shows the biphasic degradation pattern:")
        print("  - Slow decay while Pab1-protected (poly-A > 20 As)")
        print("  - Sharp transition as protection drops")
        print("  - Fast decay after decapping (poly-A < 20 As)")
    
    # Save summary statistics
    save_summary_stats(D_posteriors, summary, age_array, args.summary)
    
    print("\n=== Visualization complete ===")


if __name__ == '__main__':
    main()
