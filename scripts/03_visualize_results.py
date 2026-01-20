#!/usr/bin/env python
"""
Visualize mRNA decay inference results from Julia MCMC chains.

Based on plotting_modelling_results.ipynb
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
    """Load MCMC chain results from Julia inference."""
    print(f"Loading chain data from: {chain_path}")
    data = pd.read_csv(chain_path, engine='python', header=0)
    print(f"Loaded {len(data)} samples")
    return data


def extract_degradation_posteriors(data, n_bins=5):
    """
    Extract degradation rate (D) posteriors for each spatial bin.
    
    Args:
        data: MCMC chain DataFrame
        n_bins: Number of spatial bins (default 5)
    
    Returns:
        List of arrays containing D posteriors for each bin
    """
    print(f"Extracting degradation posteriors for {n_bins} bins...")
    
    D_posteriors = []
    
    # Try different column naming conventions
    for i in range(n_bins):
        col_name = None
        
        # Try PyMC 0-based bracket format: D[0], D[1], ... (our new format)
        if f'D[{i}]' in data.columns:
            col_name = f'D[{i}]'
        # Try Julia 1-based format: D[1], D[2], ...
        elif f'D[{i+1}]' in data.columns:
            col_name = f'D[{i+1}]'
        # Try ArviZ format: D_dim_0, D_dim_1, ...
        elif f'D_dim_{i}' in data.columns:
            col_name = f'D_dim_{i}'
        
        if col_name:
            D_posteriors.append(data[col_name].to_numpy())
        else:
            print(f"Warning: Could not find column for bin {i}")
    
    if len(D_posteriors) == 0:
        print(f"Available columns: {list(data.columns)}")
    
    return D_posteriors


def compute_halflife_from_D(D_value):
    """Compute mRNA half-life from degradation rate: t_1/2 = ln(2) / D"""
    return np.log(2) / D_value  # D is in min^-1 -> result in minutes


def compute_posterior_modes(D_posteriors):
    """
    Compute mode of each degradation rate posterior.
    
    Uses kernel density estimation to find mode.
    """
    print("Computing posterior modes...")
    
    modes = []
    for i, D_posterior in enumerate(D_posteriors):
        # Use histogram to approximate mode
        hist, bin_edges = np.histogram(D_posterior, bins=50)
        mode_idx = np.argmax(hist)
        mode_value = (bin_edges[mode_idx] + bin_edges[mode_idx + 1]) / 2
        modes.append(mode_value)
        
        print(f"  Bin {i+1}: D mode = {mode_value:.3f}, "
              f"half-life = {compute_halflife_from_D(mode_value):.2f} min")
    
    return modes


def plot_degradation_posteriors(D_posteriors, modes, output_path):
    """
    Plot histograms of degradation rate posteriors with mode lines.
    
    Args:
        D_posteriors: List of D posterior arrays
        modes: List of mode values
        output_path: Path to save figure
    """
    print(f"Plotting degradation posteriors: {output_path}")
    
    n_bins = len(D_posteriors)
    fig, axes = plt.subplots(
        nrows=n_bins, ncols=1, 
        figsize=(6, 8), 
        sharex=True
    )
    
    if n_bins == 1:
        axes = [axes]
    
    for i, (D_posterior, mode) in enumerate(zip(D_posteriors, modes)):
        # D_posterior and mode are already in min⁻¹
        ax = axes[i]
        ax.hist(
            D_posterior, 
            bins=50, 
            edgecolor='#cd96cd', 
            color='#d8bfd8',
            alpha=0.7
        )
        ax.axvline(
            mode, 
            color='#332288', 
            linestyle='dashed', 
            linewidth=1.3,
            label=f'Mode: {mode:.2f}'
        )
        ax.set_ylabel('Frequency')
        ax.set_title(f'Bin {i+1}')
        ax.legend(loc='upper right')
    # Ensure x-axis range consistent across all subplots
    for ax in axes:
        ax.set_xlim(-0.25, 5.25)

    axes[-1].set_xlabel('Degradation Rate (D, min⁻¹)')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def create_halflife_heatmap(modes, n_ap_bins, n_dv_bins, output_path):
    """
    Create spatial heatmap of mRNA half-lives.
    
    Args:
        modes: List of D mode values
        n_ap_bins: Number of AP bins
        n_dv_bins: Number of DV bins
        output_path: Path to save heatmap
    """
    print(f"Creating half-life heatmap: {output_path}")
    
    # Compute half-lives from degradation rates
    halflives = [compute_halflife_from_D(D) for D in modes]
    
    # Reshape into 2D array (AP varies fastest). If only AP modes are provided,
    # tile across DV bins for visualization.
    halflives_arr = np.array(halflives)
    expected_total = n_ap_bins * n_dv_bins
    if len(halflives_arr) == n_ap_bins and n_dv_bins > 1:
        halflife_array = np.tile(halflives_arr, (n_dv_bins, 1))
    elif len(halflives_arr) == expected_total:
        halflife_array = halflives_arr.reshape(n_dv_bins, n_ap_bins)
    else:
        raise ValueError(
            f"Cannot reshape {len(halflives_arr)} halflife values into {n_dv_bins}×{n_ap_bins}."
        )
    
    # Create heatmap
    plt.figure(figsize=(7, 6))
    sns.heatmap(
        halflife_array,
        cmap=cmc.acton_r,
        annot=True,
        fmt='.2f',
        cbar_kws={'label': 'mRNA Half-life (min)'}
    )
    plt.title('Spatial Pattern of mRNA Stability')
    plt.xlabel('AP Bin')
    plt.ylabel('DV Bin')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def save_summary_stats(D_posteriors, modes, output_path):
    """Save summary statistics of inferred parameters."""
    print(f"Saving summary statistics: {output_path}")
    
    summary_data = []
    for i, (D_posterior, mode) in enumerate(zip(D_posteriors, modes)):
        summary_data.append({
            'bin': i + 1,
            'D_mode': mode,
            'D_mean': np.mean(D_posterior),
            'D_median': np.median(D_posterior),
            'D_std': np.std(D_posterior),
            'halflife_mode_min': compute_halflife_from_D(mode),
            'halflife_mean_min': compute_halflife_from_D(np.mean(D_posterior)),
            'D_95CI_lower': np.percentile(D_posterior, 2.5),
            'D_95CI_upper': np.percentile(D_posterior, 97.5)
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_path, index=False)
    
    print("\nSummary Statistics:")
    print(summary_df.to_string(index=False))


def main():
    parser = argparse.ArgumentParser( description='Visualize mRNA decay inference results')
    parser.add_argument( '--chain', required=True, help='Path to Julia MCMC chain CSV file')
    parser.add_argument( '--posteriors-plot', required=True, help='Output path for degradation posteriors plot')
    parser.add_argument( '--heatmap-plot', required=True, help='Output path for half-life heatmap')
    parser.add_argument( '--summary', required=True, help='Output path for summary statistics CSV')
    parser.add_argument( '--n-ap-bins', type=int, default=5, help='Number of AP bins for heatmap display (default: 5)')
    parser.add_argument( '--n-dv-bins', type=int, default=1, help='Number of DV bins for heatmap display (default: 1)')
    
    args = parser.parse_args()
    
    # Load and process chain data
    chain_data = load_chain_data(args.chain)
    
    # Extract degradation posteriors
    D_posteriors = extract_degradation_posteriors(chain_data, args.n_ap_bins)
    
    if len(D_posteriors) == 0:
        print("Error: No degradation posteriors found in chain data")
        sys.exit(1)
    
    # Compute posterior modes
    modes = compute_posterior_modes(D_posteriors)
    
    # Generate plots
    plot_degradation_posteriors(D_posteriors, modes, args.posteriors_plot)
    create_halflife_heatmap(modes, args.n_ap_bins, args.n_dv_bins, args.heatmap_plot)
    
    # Save summary statistics
    save_summary_stats(D_posteriors, modes, args.summary)
    
    print("\n=== Visualization complete ===")


if __name__ == '__main__':
    main()
