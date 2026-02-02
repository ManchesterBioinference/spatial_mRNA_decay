#!/usr/bin/env python3
"""
Aggregate and validate mRNA counts from SASS position_data-intense.txt.

This script:
1. Bins mRNA spot data into an AP×DV spatial grid matching transcription data binning
2. Computes k-nearest neighbor metrics for nuclei density validation
3. Computes bin nuclei counts and validates against Berrocal_2020 distribution
4. Generates QC visualizations and validation report

Usage:
    python 04_aggregate_and_validate_mrna.py <position_data_file> <stripe> <embryo_id> <config_yaml> <output_file> <validation_report>

Output:
    - Binned mRNA counts CSV file (n_ap_bins × n_dv_bins values, no header) for inference
    - Validation report (text file with PASS/FAIL status)
    - Heatmap visualization of spatial binning
    - X-Y scatter plot showing spatial distribution of spots and nuclei
    - Expression density plot (spot-level) along X-axis
    - Ridge plot comparing bin nuclei count distributions (transcription vs mRNA)
"""

import os
import subprocess
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from ruamel.yaml import YAML


def load_config(config_path: str) -> dict:
    """Load config.yaml and return as dict."""
    yaml = YAML()
    with open(config_path, 'r') as f:
        return yaml.load(f)


def load_position_data(filepath: str) -> pd.DataFrame:
    """
    Load SASS position_data-intense.txt file.
    
    Expected format: tab-separated with columns including:
    - nucx, nucy, nucz: nuclear positions
    - spotx, spoty, spotz: spot positions
    - nuc: nuclear ID
    - num_spots: number of spots per nucleus
    - time: time point
    """
    df = pd.read_csv(filepath, sep='\t')
    return df


def bin_mrna_data(df: pd.DataFrame, n_ap_bins: int = 5, n_dv_bins: int = 5, no_groupByNuclei: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Bin mRNA data into spatial grid using the data's own spatial extent.
    
    Strategy:
    1. Aggregate spots by nucleus
    2. Normalize AP and DV coordinates to [0, 1] using actual data range
    3. Create n_ap_bins × n_dv_bins spatial grid
    4. Compute average mRNA count per nucleus in each bin
    
    Args:
        df: DataFrame with nuclear positions and spot counts
        n_ap_bins: Number of bins along AP axis (default: 5)
        n_dv_bins: Number of bins along DV axis (default: 5)
    
    Returns:
        Tuple of (binned_data DataFrame, nuc_data_with_bins DataFrame, spatial_ranges dict)
        - binned_data: DataFrame with columns: apBin, yBin, avg_mrna_count
        - nuc_data_with_bins: DataFrame with one row per nucleus including bin assignments
        - spatial_ranges: dict with keys 'ap_min', 'ap_max', 'dv_min', 'dv_max'
    """
    # Normalize coordinates to [0, 1] using actual data extent
    ap_min_data = df['nucx'].min()
    ap_max_data = df['nucx'].max()
    ap_range = ap_max_data - ap_min_data
    dv_min = df['nucy'].min()
    dv_max = df['nucy'].max()
    dv_range = dv_max - dv_min
    
    # Aggregate to one row per nucleus with total spot count
        
    if no_groupByNuclei:
        nuc_data_filtered = df[['spot','spotx','spoty','nuc', 'nucx', 'nucy',]].copy()
        nuc_data_filtered = nuc_data_filtered[(nuc_data_filtered['spotx'] >= ap_min_data) & (nuc_data_filtered['spotx'] <= ap_max_data) & (nuc_data_filtered['spoty'] >= dv_min) & (nuc_data_filtered['spoty'] <= dv_max)]
        nuc_data_filtered = nuc_data_filtered.drop_duplicates(subset=['spot'])
        nuc_data_filtered['ap_norm'] = (nuc_data_filtered['spotx'] - ap_min_data) / ap_range
        nuc_data_filtered['dv_norm'] = (nuc_data_filtered['spoty'] - dv_min) / dv_range
    else:
        # Get positions at t_min
        nuc_data_filtered = df[['nuc', 'nucx', 'nucy','num_spots']].copy()
        nuc_data_filtered = nuc_data_filtered.drop_duplicates(subset=['nuc'])
        nuc_data_filtered['ap_norm'] = (nuc_data_filtered['nucx'] - ap_min_data) / ap_range
        nuc_data_filtered['dv_norm'] = (nuc_data_filtered['nucy'] - dv_min) / dv_range
    
    # Assign spatial bins
    nuc_data_filtered['apBin'] = pd.cut(
        nuc_data_filtered['ap_norm'],
        bins=n_ap_bins,
        labels=False,
        include_lowest=True
    )
    nuc_data_filtered['yBin'] = pd.cut(
        nuc_data_filtered['dv_norm'],
        bins=n_dv_bins,
        labels=False,
        include_lowest=True
    )
    
    if no_groupByNuclei:
        # Count total spots in each bin
        binned_data = nuc_data_filtered.groupby(['apBin', 'yBin']).agg({
            'spot': 'count'
        }).reset_index()
    else:
        # Compute average mRNA count per nucleus in each bin
        binned_data = nuc_data_filtered.groupby(['apBin', 'yBin']).agg({
            'num_spots': 'mean'
        }).reset_index()
    
    binned_data.columns = ['apBin', 'yBin', 'avg_mrna_count']
    
    # Ensure all bins are represented (fill missing with 0)
    all_bins = pd.DataFrame([
        (ap, y) for ap in range(n_ap_bins) for y in range(n_dv_bins)
    ], columns=['apBin', 'yBin'])
    
    binned_data = all_bins.merge(binned_data, on=['apBin', 'yBin'], how='left')
    binned_data['avg_mrna_count'] = binned_data['avg_mrna_count'].fillna(0)
    
    # Sort by apBin, then yBin for consistent ordering
    binned_data = binned_data.sort_values(['apBin', 'yBin']).reset_index(drop=True)
    
    # Return spatial ranges used for binning
    spatial_ranges = {
        'ap_min': ap_min_data,
        'ap_max': ap_max_data,
        'dv_min': dv_min,
        'dv_max': dv_max
    }

    # Ensure we have nuc level data in nuc_data_filtered even if we don't group by nuclei
    nuc_data_filtered = nuc_data_filtered.drop_duplicates(subset=['nuc'])
    
    return binned_data, nuc_data_filtered, spatial_ranges


def trim_ap_axis(df: pd.DataFrame, trim_fraction: float = 0.05) -> tuple[pd.DataFrame, float, float]:
    """
    Trim the AP axis evenly from both ends to reduce nuclei density.
    
    This removes data symmetrically from the left and right edges of the AP axis
    to keep the data centered while reducing the number of nuclei.
    
    Args:
        df: DataFrame with nucx column
        trim_fraction: Fraction to trim from EACH end (total trimming = 2 * trim_fraction)
    
    Returns:
        Tuple of (trimmed_df, new_ap_min, new_ap_max)
    """
    ap_min = df['nucx'].min()
    ap_max = df['nucx'].max()
    ap_range = ap_max - ap_min
    
    # Calculate new AP bounds (trim evenly from both ends)
    trim_amount = ap_range * trim_fraction
    new_ap_min = ap_min + trim_amount
    new_ap_max = ap_max - trim_amount
    
    # Filter data to keep only nuclei within new bounds
    trimmed_df = df[(df['nucx'] >= new_ap_min) & (df['nucx'] <= new_ap_max)].copy()
    
    return trimmed_df, new_ap_min, new_ap_max


def attempt_trimming_for_density(
    df: pd.DataFrame,
    thresholds: dict,
    n_ap_bins: int,
    n_dv_bins: int,
    no_groupByNuclei: bool,
    k: int = 4,
    max_iterations: int = 10,
    trim_fraction: float = 0.05
) -> tuple[pd.DataFrame, dict, dict]:
    """
    Attempt to correct nuclei density by iteratively trimming the AP axis.
    
    Only applies when nuclei are too dense (median_nn_distance < min_nn_distance).
    Trims evenly from both ends to keep data centered.
    
    Args:
        df: Original DataFrame
        thresholds: Validation thresholds
        n_ap_bins: Number of AP bins
        n_dv_bins: Number of DV bins
        no_groupByNuclei: Whether to skip grouping by nuclei
        k: Number of nearest neighbors for density calculation
        max_iterations: Maximum number of trimming attempts
        trim_fraction: Fraction to trim from each end per iteration
    
    Returns:
        Tuple of (trimmed_df or original_df, trim_info, final_metrics)
        - If trimming not needed or fails: returns original df, empty trim_info, original metrics
        - If trimming succeeds: returns trimmed df, trim_info with details, new metrics
    """
    # First, compute metrics on original data
    _, nuc_data_original, _ = bin_mrna_data(df, n_ap_bins, n_dv_bins, no_groupByNuclei)
    original_metrics = compute_nuclei_density_metrics(nuc_data_original, k=k)
    nn_passed, _ = validate_nn_metrics(original_metrics, thresholds)
    
    # Check if nuclei are too dense (median NN distance below minimum)
    min_nn = thresholds['min_nn_distance']
    if original_metrics['median_nn_distance'] >= min_nn:
        # Density is acceptable or too sparse - no trimming needed
        return df, {}, original_metrics
    
    # Nuclei are too dense - attempt trimming
    print(f"  Nuclei too dense: median NN distance {original_metrics['median_nn_distance']:.4f} < {min_nn:.4f}")
    print(f"  Attempting AP axis trimming...")
    
    current_df = df.copy()
    original_ap_min = df['nucx'].min()
    original_ap_max = df['nucx'].max()
    
    for iteration in range(1, max_iterations + 1):
        # Trim AP axis evenly from both ends
        trimmed_df, new_ap_min, new_ap_max = trim_ap_axis(current_df, trim_fraction)
        
        # Check if we still have enough data
        if len(trimmed_df) == 0:
            print(f"  WARNING: Trimming resulted in no data. Reverting to previous iteration.")
            break
        
        # Compute new metrics
        _, nuc_data_trimmed, _ = bin_mrna_data(trimmed_df, n_ap_bins, n_dv_bins, no_groupByNuclei)
        
        if len(nuc_data_trimmed) == 0:
            print(f"  WARNING: No nuclei after trimming. Reverting to previous iteration.")
            break
        
        trimmed_metrics = compute_nuclei_density_metrics(nuc_data_trimmed, k=k)
        nn_passed, _ = validate_nn_metrics(trimmed_metrics, thresholds)
        
        total_trim_percent = (1 - (new_ap_max - new_ap_min) / (original_ap_max - original_ap_min)) * 100
        
        print(f"    Iteration {iteration}: median NN = {trimmed_metrics['median_nn_distance']:.4f}, "
              f"n_nuclei = {trimmed_metrics['n_nuclei']}, trimmed {total_trim_percent:.1f}%")
        
        if nn_passed:
            # Success! Density is now acceptable
            trim_info = {
                'applied': True,
                'iterations': iteration,
                'original_ap_range': (original_ap_min, original_ap_max),
                'trimmed_ap_range': (new_ap_min, new_ap_max),
                'trim_fraction_per_iteration': trim_fraction,
                'total_trim_percent': total_trim_percent,
                'original_n_nuclei': original_metrics['n_nuclei'],
                'final_n_nuclei': trimmed_metrics['n_nuclei'],
                'original_median_nn': original_metrics['median_nn_distance'],
                'final_median_nn': trimmed_metrics['median_nn_distance']
            }
            print(f"  ✓ Trimming successful after {iteration} iteration(s)")
            return trimmed_df, trim_info, trimmed_metrics
        
        # Update for next iteration
        current_df = trimmed_df
    
    # Max iterations reached without success
    print(f"  ✗ Failed to correct density after {max_iterations} iterations")
    return df, {}, original_metrics


def compute_nuclei_density_metrics(nuc_data: pd.DataFrame, k: int = 4) -> dict:
    """
    Compute nuclei density metrics using k-nearest neighbors.
    
    Args:
        nuc_data: DataFrame with one row per nucleus (nucx, nucy columns)
        k: Number of nearest neighbors to compute
    
    Returns:
        Dictionary with validation metrics
    """
    n_nuclei = len(nuc_data)
    
    if n_nuclei == 0:
        raise ValueError("No nuclei found in position data")
    
    # Normalize coordinates to [0, 1]
    x_min = nuc_data['nucx'].min()
    x_max = nuc_data['nucx'].max()
    y_min = nuc_data['nucy'].min()
    y_max = nuc_data['nucy'].max()
    
    nucx_norm = (nuc_data['nucx'] - x_min) / (x_max - x_min)
    nucy_norm = (nuc_data['nucy'] - y_min) / (y_max - y_min)
    
    # Compute k-nearest neighbors in normalized coordinates
    coords = np.column_stack([nucx_norm, nucy_norm])
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(coords)
    distances, _ = nbrs.kneighbors(coords)
    
    # Remove first column (distance to self = 0)
    distances = distances[:, 1:]
    
    # Compute metrics
    avg_distance = np.mean(distances)
    std_distance = np.std(distances)
    median_distance = np.median(distances)
    
    return {
        'n_nuclei': n_nuclei,
        'avg_nn_distance': avg_distance,
        'std_nn_distance': std_distance,
        'median_nn_distance': median_distance,
        'min_nn_distance': np.min(distances),
        'max_nn_distance': np.max(distances)
    }


def compute_bin_nuclei_counts(nuc_data: pd.DataFrame, n_ap_bins: int = 5, n_dv_bins: int = 5) -> list[int]:
    """
    Compute the number of nuclei in each spatial bin.
    
    Args:
        nuc_data: DataFrame with one row per nucleus including apBin, yBin columns
        n_ap_bins: Number of bins along AP axis (default: 5)
        n_dv_bins: Number of bins along DV axis (default: 5)
    
    Returns:
        List of nuclei counts (one per bin), ordered by [apBin, yBin]
    """
    # Count nuclei per bin
    bin_counts = nuc_data.groupby(['apBin', 'yBin']).size().reset_index(name='nuclei_count')
    
    # Ensure all bins are represented
    all_bins = pd.DataFrame([
        (ap, y) for ap in range(n_ap_bins) for y in range(n_dv_bins)
    ], columns=['apBin', 'yBin'])
    
    bin_counts = all_bins.merge(bin_counts, on=['apBin', 'yBin'], how='left')
    bin_counts['nuclei_count'] = bin_counts['nuclei_count'].fillna(0).astype(int)
    
    # Sort by apBin, then yBin for consistent ordering
    bin_counts = bin_counts.sort_values(['apBin', 'yBin']).reset_index(drop=True)
    
    return bin_counts['nuclei_count'].tolist()


def load_validation_thresholds(config: dict, stripe: str) -> dict:
    """
    Load validation thresholds for a specific stripe from config.yaml.
    
    Args:
        config_path: Path to config.yaml
        stripe: Stripe name (e.g., 'stripe2', 'stripe3')
    
    Returns:
        Dictionary with validation thresholds
    """
    validation_thresholds = config.get('validation_thresholds', {})
    if stripe not in validation_thresholds:
        raise ValueError(
            f"No validation thresholds found for {stripe} in config.yaml. "
            f"Run 06_compute_validation_thresholds.py first."
        )
    
    return validation_thresholds[stripe]


def validate_nn_metrics(metrics: dict, thresholds: dict) -> tuple[bool, list[str]]:
    """
    Validate k-NN metrics against stripe-specific Berrocal_2020 benchmarks.
    
    Args:
        metrics: Computed metrics from position data
        thresholds: Expected thresholds from config.yaml
    
    Returns:
        (passed, warnings): Boolean indicating if validation passed, list of warning messages
    """
    warnings = []
    passed = True
    
    min_nn = thresholds['min_nn_distance']
    max_nn = thresholds['max_nn_distance']
    expected_median = thresholds['expected_nn_distance_median']
    expected_count = thresholds['expected_nuclei_count_mean']
    
    # Hard gate: median NN distance must be within observed range
    if metrics['median_nn_distance'] < min_nn:
        warnings.append(
            f"FAIL: Median NN distance {metrics['median_nn_distance']:.4f} < "
            f"minimum observed {min_nn:.4f} (nuclei too densely packed)"
        )
        passed = False
    
    if metrics['median_nn_distance'] > max_nn:
        warnings.append(
            f"FAIL: Median NN distance {metrics['median_nn_distance']:.4f} > "
            f"maximum observed {max_nn:.4f} (nuclei too sparsely packed)"
        )
        passed = False
    
    # Soft warnings for expected values
    if abs(metrics['median_nn_distance'] - expected_median) > 0.02:
        warnings.append(
            f"WARNING: Median NN distance {metrics['median_nn_distance']:.4f} differs from "
            f"expected {expected_median:.4f} by more than 0.02"
        )
    
    if abs(metrics['n_nuclei'] - expected_count) > 50:
        warnings.append(
            f"WARNING: Nuclei count {metrics['n_nuclei']} differs from "
            f"expected {expected_count:.1f} by more than 50"
        )
    
    return passed, warnings


def validate_bin_counts(
    mrna_bin_counts: list[int],
    thresholds: dict,
    n_ap_bins: int,
    n_dv_bins: int
) -> tuple[bool, list[str]]:
    """
    Validate mRNA bin nuclei counts against transcription data distribution.
    
    Uses 2-sigma range: expected ± 2*std for each bin.
    
    Args:
        mrna_bin_counts: List of nuclei counts per bin for mRNA data
        thresholds: Validation thresholds including bin_count_distribution
        n_ap_bins: Number of bins along AP axis
        n_dv_bins: Number of bins along DV axis
    
    Returns:
        (passed, warnings): Boolean indicating if validation passed, list of warning messages
    """
    warnings = []
    passed = True
    
    if 'bin_count_distribution' not in thresholds:
        warnings.append("WARNING: No bin count distribution found in validation thresholds. Skipping bin count validation.")
        return True, warnings
    
    bin_dist = thresholds['bin_count_distribution']
    stats = bin_dist.get('stats', {})
    
    if not stats:
        warnings.append("WARNING: No bin count statistics found. Skipping bin count validation.")
        return True, warnings
    
    mean = np.array(stats['mean'])
    std = np.array(stats['std'])
    expected_bins = n_ap_bins * n_dv_bins

    if len(mean) != expected_bins or len(std) != expected_bins:
        warnings.append(
            "FAIL: Bin count statistics length does not match current binning "
            f"({len(mean)} stats vs {expected_bins} bins). "
            "Re-run 06_compute_validation_thresholds.py with the updated config."
        )
        return False, warnings
    
    # Compute 2-sigma bounds
    lower_bound = mean - 2 * std
    upper_bound = mean + 2 * std
    
    # Count how many bins fall outside 2-sigma range
    mrna_counts = np.array(mrna_bin_counts)
    outside_range = (mrna_counts < lower_bound) | (mrna_counts > upper_bound)
    n_outside = outside_range.sum()
    
    if n_outside > max(1, int(0.2 * expected_bins)):
        warnings.append(
            f"FAIL: {n_outside}/{expected_bins} bins have nuclei counts outside 2-sigma range "
            f"(>{20}% threshold)"
        )
        passed = False
    elif n_outside > 0:
        warnings.append(
            f"WARNING: {n_outside}/{expected_bins} bins have nuclei counts outside 2-sigma range"
        )
    
    return passed, warnings


def plot_heatmap(binned_data: pd.DataFrame, output_path: str, n_ap_bins: int, n_dv_bins: int, 
                 binned_data_original: pd.DataFrame = None) -> None:
    """
    Create and save a heatmap visualization of the spatial binning.
    
    Args:
        binned_data: DataFrame with columns apBin, yBin, avg_mrna_count (trimmed or final)
        output_path: Path to save the heatmap PNG file
        n_ap_bins: Number of bins along AP axis
        n_dv_bins: Number of bins along DV axis
        binned_data_original: Optional original data for comparison (if trimming was applied)
    """
    # Determine if we need dual-panel comparison
    is_comparison = binned_data_original is not None
    
    if is_comparison:
        # Create dual-panel figure
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # Original data (left panel)
        heatmap_data_orig = binned_data_original.pivot(index='yBin', columns='apBin', values='avg_mrna_count')
        sns.heatmap(
            heatmap_data_orig,
            annot=True,
            fmt='.2f',
            cmap='viridis',
            cbar_kws={'label': 'Avg mRNA count/nucleus'},
            ax=axes[0],
            linewidths=0.5,
            linecolor='white'
        )
        axes[0].set_xlabel('AP bin', fontsize=12)
        axes[0].set_ylabel('DV bin', fontsize=12)
        axes[0].set_title(f'Original Data\nmRNA abundance: {n_ap_bins}×{n_dv_bins} spatial binning', fontsize=14)
        
        # Trimmed data (right panel)
        heatmap_data = binned_data.pivot(index='yBin', columns='apBin', values='avg_mrna_count')
        sns.heatmap(
            heatmap_data,
            annot=True,
            fmt='.2f',
            cmap='viridis',
            cbar_kws={'label': 'Avg mRNA count/nucleus'},
            ax=axes[1],
            linewidths=0.5,
            linecolor='white'
        )
        axes[1].set_xlabel('AP bin', fontsize=12)
        axes[1].set_ylabel('DV bin', fontsize=12)
        axes[1].set_title(f'After Trimming\nmRNA abundance: {n_ap_bins}×{n_dv_bins} spatial binning', fontsize=14)
        
    else:
        # Single panel figure
        fig, ax = plt.subplots(figsize=(8, 7))
        heatmap_data = binned_data.pivot(index='yBin', columns='apBin', values='avg_mrna_count')
        
        sns.heatmap(
            heatmap_data,
            annot=True,
            fmt='.2f',
            cmap='viridis',
            cbar_kws={'label': 'Avg mRNA count/nucleus'},
            ax=ax,
            linewidths=0.5,
            linecolor='white'
        )
        
        ax.set_xlabel('AP bin', fontsize=12)
        ax.set_ylabel('DV bin', fontsize=12)
        ax.set_title(f'mRNA abundance: {n_ap_bins}×{n_dv_bins} spatial binning', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    heatmap_data.to_csv(output_path.replace('.png', '_data.csv'))


def plot_spatial_distribution(df: pd.DataFrame, output_path: str, spatial_ranges: dict,
                              n_ap_bins: int = 5, n_dv_bins: int = 5,
                              df_original: pd.DataFrame = None, 
                              spatial_ranges_original: dict = None) -> None:
    """
    Create X-Y scatter plot showing spatial distribution of spots and nuclei with grid overlay.
    
    Args:
        df: DataFrame with spot and nuclear positions (trimmed or final)
        output_path: Path to save the scatter plot PNG file
        spatial_ranges: Dict with 'ap_min', 'ap_max', 'dv_min', 'dv_max' from binning
        n_ap_bins: Number of bins along AP axis (default: 5)
        n_dv_bins: Number of bins along DV axis (default: 5)
        df_original: Optional original data for comparison (if trimming was applied)
        spatial_ranges_original: Optional original spatial ranges for comparison
    """
    # Determine if we need dual-panel comparison
    is_comparison = df_original is not None
    
    def plot_single_distribution(data, ranges, ax, title_prefix=""):
        """Helper function to plot a single spatial distribution."""
        sample_t = data['time'].min()
        sample_data = data[data['time'] == sample_t].copy()
        
        ap_data_min = sample_data['nucx'].min()
        ap_data_max = sample_data['nucx'].max()
        dv_min = sample_data['nucy'].min()
        dv_max = sample_data['nucy'].max()
        
        unique_nucs = sample_data['nuc'].unique()
        shuffled_nucs = np.random.permutation(unique_nucs)
        color_map = {nuc: i for i, nuc in enumerate(shuffled_nucs)}
        colors = sample_data['nuc'].map(color_map)
        
        ax.scatter(sample_data['spotx'], sample_data['spoty'], 
                   alpha=0.5, s=10, c=colors, cmap='tab20', label='Spots')
        ax.scatter(sample_data['nucx'], sample_data['nucy'], 
                   alpha=0.7, s=50, c='steelblue', marker='s', label='Nuclei')
        
        ap_range = ap_data_max - ap_data_min
        dv_range = dv_max - dv_min
        
        for i in range(n_ap_bins + 1):
            x = ap_data_min + i * (ap_range / n_ap_bins)
            ax.axvline(x, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        
        for i in range(n_dv_bins + 1):
            y = dv_min + i * (dv_range / n_dv_bins)
            ax.axhline(y, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        
        ax.set_xlabel('X Position (AP axis)', fontsize=12)
        ax.set_ylabel('Y Position (DV axis)', fontsize=12)
        ax.set_title(
            f'{title_prefix}Spatial Distribution (X-Y) with {n_ap_bins}×{n_dv_bins} Grid\nt={sample_t:.1f}',
            fontsize=14
        )
        ax.set_aspect('equal')
        ax.legend()
    
    if is_comparison:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        plot_single_distribution(df_original, spatial_ranges_original, axes[0], "Original Data\n")
        plot_single_distribution(df, spatial_ranges, axes[1], "After Trimming\n")
    else:
        fig, ax = plt.subplots(figsize=(10, 8))
        plot_single_distribution(df, spatial_ranges, ax)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_expression_density(df: pd.DataFrame, output_path: str, bin_size: float = 0.75,
                           df_original: pd.DataFrame = None) -> None:
    """
    Create spot-level expression density plot across X-axis.
    
    Args:
        df: DataFrame with spot positions (trimmed or final)
        output_path: Path to save the density plot PNG file
        bin_size: Size of bins along X-axis (default: 0.75)
        df_original: Optional original data for comparison (if trimming was applied)
    """
    import scipy.signal as signal
    
    # Determine if we need dual-panel comparison
    is_comparison = df_original is not None
    
    def compute_density(data):
        """Helper function to compute density for a dataset."""
        min_x = data['spotx'].min()
        max_x = data['spotx'].max()
        bins = np.arange(min_x, max_x + bin_size, bin_size)
        df_copy = data.copy()
        df_copy['x_bin'] = pd.cut(df_copy['spotx'], bins)
        
        spot_density = df_copy.groupby('x_bin').size().reset_index(name='spot_count')
        spot_density['x_bin_center'] = spot_density['x_bin'].apply(lambda x: x.mid)
        spot_density = spot_density.dropna()
        
        window = 3
        spot_density['spot_count_smooth'] = np.convolve(
            spot_density['spot_count'], 
            np.ones(window)/window, 
            mode='same'
        )
        
        peaks, _ = signal.find_peaks(spot_density['spot_count_smooth'].values, prominence=200)
        return spot_density, peaks
    
    def plot_single_density(data, ax, title_prefix=""):
        """Helper function to plot density on an axis."""
        spot_density, peaks = compute_density(data)
        
        ax.plot(spot_density['x_bin_center'].values, 
                spot_density['spot_count_smooth'].values, 
                'o-', linewidth=2, markersize=5, color='darkorange', 
                label='Spot count per bin')
        ax.fill_between(spot_density['x_bin_center'].values, 
                         (spot_density['spot_count'] - spot_density['spot_count'].std()).values, 
                         (spot_density['spot_count'] + spot_density['spot_count'].std()).values, 
                         alpha=0.2, color='darkorange')
        ax.plot(spot_density['x_bin_center'].iloc[peaks].values, 
                spot_density['spot_count_smooth'].iloc[peaks].values, 
                'ro', markersize=8, label='Peaks')
        
        ax.set_xlabel('X Position', fontsize=12)
        ax.set_ylabel('Number of Spots', fontsize=12)
        ax.set_title(f'{title_prefix}Expression Density - Spot Level (Spot Count)', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    if is_comparison:
        fig, axes = plt.subplots(2, 1, figsize=(12, 12))
        plot_single_density(df_original, axes[0], "Original Data\n")
        plot_single_density(df, axes[1], "After Trimming\n")
    else:
        fig, ax = plt.subplots(figsize=(12, 6))
        plot_single_density(df, ax)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_nuclei_expression_density(df: pd.DataFrame, output_path: str, bin_size: float = 0.75,
                                   df_original: pd.DataFrame = None) -> None:
    """
    Create nucleus-level expression density plot across X-axis.

    Spots are summed per nucleus (counts per `nuc`) and bins average across
    nuclei in each X-bin (mean spots per nucleus), rather than averaging
    across individual spots.

    Args:
        df: DataFrame with spot rows including `nuc` and `nucx` columns (trimmed or final)
        output_path: Path to save the density plot PNG file
        bin_size: Size of bins along X-axis (default: 0.75)
        df_original: Optional original data for comparison (if trimming was applied)
    """
    import scipy.signal as signal

    # Determine if we need dual-panel comparison
    is_comparison = df_original is not None
    
    def compute_nuclei_density(data):
        """Helper function to compute nuclei-level density."""
        nuc_counts = (
            data.groupby('nuc', dropna=True)
              .agg(nucx=('nucx', 'first'), spots_per_nuc=('spotx', 'count'))
              .reset_index()
        )

        if nuc_counts.empty:
            return None, None

        min_x = nuc_counts['nucx'].min()
        max_x = nuc_counts['nucx'].max()
        bins = np.arange(min_x, max_x + bin_size, bin_size)
        nuc_counts['x_bin'] = pd.cut(nuc_counts['nucx'], bins)

        bin_stats = nuc_counts.groupby('x_bin').agg(
            mean_spots_per_nuc=('spots_per_nuc', 'mean'),
            std_spots_per_nuc=('spots_per_nuc', 'std'),
            n_nuclei=('spots_per_nuc', 'size')
        ).reset_index()
        bin_stats['x_bin_center'] = bin_stats['x_bin'].apply(lambda x: x.mid)
        bin_stats = bin_stats.dropna()

        if bin_stats.empty:
            return None, None

        window = 3
        bin_stats['mean_smooth'] = np.convolve(
            bin_stats['mean_spots_per_nuc'].values, np.ones(window) / window, mode='same'
        )

        prom = max(1.0, float(bin_stats['mean_smooth'].std() * 1.5))
        peaks, _ = signal.find_peaks(bin_stats['mean_smooth'].values, prominence=prom)
        
        return bin_stats, peaks
    
    def plot_single_nuclei_density(data, ax, title_prefix=""):
        """Helper function to plot nuclei density on an axis."""
        bin_stats, peaks = compute_nuclei_density(data)
        
        if bin_stats is None:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
            return
        
        ax.plot(
            bin_stats['x_bin_center'].values,
            bin_stats['mean_smooth'].values,
            'o-', linewidth=2, markersize=5, color='seagreen', label='Mean spots per nucleus (binned)'
        )

        std_vals = bin_stats['std_spots_per_nuc'].fillna(0).values
        ax.fill_between(
            bin_stats['x_bin_center'].values,
            (bin_stats['mean_spots_per_nuc'] - std_vals).values,
            (bin_stats['mean_spots_per_nuc'] + std_vals).values,
            alpha=0.2, color='seagreen'
        )

        if len(peaks) > 0:
            ax.plot(
                bin_stats['x_bin_center'].iloc[peaks].values,
                bin_stats['mean_smooth'].iloc[peaks].values,
                'ro', markersize=8, label='Peaks'
            )

        ax.set_xlabel('X Position', fontsize=12)
        ax.set_ylabel('Mean spots per nucleus', fontsize=12)
        ax.set_title(f'{title_prefix}Expression Density - Nucleus Level (Mean spots per nucleus)', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    if is_comparison:
        fig, axes = plt.subplots(2, 1, figsize=(6, 12))
        plot_single_nuclei_density(df_original, axes[0], "Original Data\n")
        plot_single_nuclei_density(df, axes[1], "After Trimming\n")
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_single_nuclei_density(df, ax)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_bin_count_ridge(mrna_bin_counts: list[int], thresholds: dict, 
                         embryo_id: str, output_path: str,
                         mrna_bin_counts_original: list[int] = None) -> None:
    """
    Create ridge plot comparing bin nuclei count distributions.
    
    Shows transcription embryos (top rows), all transcription combined (middle),
    and mRNA embryo (bottom, highlighted). If trimming was applied, shows both
    original and trimmed mRNA distributions.
    
    Args:
        mrna_bin_counts: List of nuclei counts per bin for mRNA data (trimmed or final)
        thresholds: Validation thresholds including bin_count_distribution
        embryo_id: ID of mRNA embryo for labeling
        output_path: Path to save the ridge plot PNG file
        mrna_bin_counts_original: Optional original bin counts (if trimming was applied)
    """
    if 'bin_count_distribution' not in thresholds:
        print("  WARNING: No bin count distribution found. Skipping ridge plot.")
        return
    
    bin_dist = thresholds['bin_count_distribution']
    
    # Prepare data for plotting
    plot_data = []
    
    # Add transcription embryo distributions (each as a separate row)
    embryo_ids = sorted([k for k in bin_dist.keys() if k not in ['all_combined', 'stats']])
    for emb_id in embryo_ids:
        counts = bin_dist[emb_id]
        for count in counts:
            plot_data.append({'embryo': f'Transcription {emb_id}', 'count': count, 'type': 'transcription'})
    
    # Add all combined transcription
    if 'all_combined' in bin_dist:
        for count in bin_dist['all_combined']:
            plot_data.append({'embryo': 'Transcription (all)', 'count': count, 'type': 'transcription_all'})
    
    # Add mRNA embryo distribution(s)
    if mrna_bin_counts_original is not None:
        # Show both original and trimmed
        for count in mrna_bin_counts_original:
            plot_data.append({'embryo': f'mRNA {embryo_id} (original)', 'count': count, 'type': 'mrna_original'})
        for count in mrna_bin_counts:
            plot_data.append({'embryo': f'mRNA {embryo_id} (trimmed)', 'count': count, 'type': 'mrna'})
    else:
        # Show only final data
        for count in mrna_bin_counts:
            plot_data.append({'embryo': f'mRNA {embryo_id}', 'count': count, 'type': 'mrna'})
    
    df = pd.DataFrame(plot_data)
    
    # Calculate global min/max for x-axis
    global_min = df['count'].min()
    global_max = df['count'].max()
    
    # Create figure
    categories = df['embryo'].unique()
    n_categories = len(categories)
    
    fig, axes = plt.subplots(n_categories, 1, figsize=(6, 1 * n_categories), sharex=True)
    if n_categories == 1:
        axes = [axes]
    
    for idx, (ax, category) in enumerate(zip(axes, categories)):
        data = df[df['embryo'] == category]['count']
        data_type = df[df['embryo'] == category]['type'].iloc[0]
        
        # Determine color
        if data_type == 'mrna':
            color = 'red'
            alpha = 0.7
        elif data_type == 'mrna_original':
            color = 'orange'
            alpha = 0.6
        elif data_type == 'transcription_all':
            color = 'blue'
            alpha = 0.5
        else:
            color = 'gray'
            alpha = 0.4
        
        # Plot KDE
        ax.fill_between(
            np.linspace(data.min() - 2, data.max() + 2, 100),
            0,
            [0] * 100,  # Will be replaced by KDE
            alpha=alpha,
            color=color
        )
        
        # Use seaborn kdeplot for smooth density
        sns.kdeplot(data=data, ax=ax, fill=True, alpha=alpha, color=color, linewidth=2)
        
        # Add label
        ax.text(0.02, 0.5, category, transform=ax.transAxes, fontsize=10, 
                verticalalignment='center')
        
        # Remove y-axis labels and ticks
        ax.set_ylabel('')
        ax.set_yticks([])
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # Only show x-axis on bottom plot
        if idx < n_categories - 1:
            ax.set_xlabel('')
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])
    
    # Set x-axis limits and ticks on bottom plot
    axes[-1].set_xlim(global_min - 1, global_max + 1)
    tick_step = max(1, (global_max - global_min) // 5)
    axes[-1].set_xticks(range(global_min, global_max + 1, tick_step))
    
    # Set x-axis label
    fig.text(0.5, 0.001, 'Nuclei count per bin', ha='center', fontsize=12)
    axes[-1].tick_params(axis='x', labelbottom=True)
    
    fig.suptitle('Bin Nuclei Count Distributions', fontsize=14, y=0.995)
    fig.subplots_adjust(bottom=0.75)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    if len(sys.argv) != 9:
        print("Usage: python 04_aggregate_and_validate_mrna.py <position_data_file> <stripe> <embryo_id> <config_yaml> <output_file> <validation_report> <figure_directory>")
        sys.exit(1)
    
    position_data_file = sys.argv[1]
    stripe = sys.argv[2]
    embryo_id = sys.argv[3]
    config_yaml = sys.argv[4]
    output_file = sys.argv[5]
    validation_report = sys.argv[6]
    fig_dir = Path(sys.argv[7]).parent
    no_groupByNuclei = sys.argv[8] == 'True'  # Set to True to skip grouping by nuclei for binning

    # Print full command to reproduce the run
    print("\n=== Command to reproduce this script run ===")
    print(" ".join(sys.argv))
    
    print(f"Processing {position_data_file}")
    print(f"  Stripe: {stripe}")
    print(f"  Embryo: {embryo_id}")
    
    config = load_config(config_yaml)
    n_ap_bins = int(config.get('n_ap_bins', 5))
    n_dv_bins = int(config.get('n_dv_bins', 5))

    # if stripe == 'stripe2':
    #     if n_ap_bins == 5:
    #         subprocess.run(['ln','-s','data/processed_mRNA_data_stripe2',os.path.join(output_file.split('/')[0],'data')])
    #     else:
    #         sys.exit("ERROR: stripe2 only supports n_ap_bins=5")

    # Load validation thresholds
    try:
        thresholds = load_validation_thresholds(config, stripe)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # Load and process data
    df = load_position_data(position_data_file)
    print(f"  Loaded {len(df)} spot records")
    
    # Attempt AP axis trimming if nuclei are too dense
    k = thresholds.get('k', 4)
    df_processed, trim_info, initial_metrics = attempt_trimming_for_density(
        df,
        thresholds=thresholds,
        n_ap_bins=n_ap_bins,
        n_dv_bins=n_dv_bins,
        no_groupByNuclei=no_groupByNuclei,
        k=k,
        max_iterations=10,
        trim_fraction=0.05
    )
    
    # Bin mRNA data (using potentially trimmed data)
    binned_data, nuc_data_with_bins, spatial_ranges = bin_mrna_data(
        df_processed,
        n_ap_bins=n_ap_bins,
        n_dv_bins=n_dv_bins,
        no_groupByNuclei=no_groupByNuclei
    )
    print(f"  Binned into {len(binned_data)} spatial bins ({n_ap_bins}×{n_dv_bins})")
    print(f"  Spatial extent - AP: [{spatial_ranges['ap_min']:.2f}, {spatial_ranges['ap_max']:.2f}], DV: [{spatial_ranges['dv_min']:.2f}, {spatial_ranges['dv_max']:.2f}]")
    
    # Compute validation metrics (will match initial_metrics if trimming was applied)
    nn_metrics = compute_nuclei_density_metrics(nuc_data_with_bins, k=k)
    print(f"  Computed k={k} NN metrics: median distance = {nn_metrics['median_nn_distance']:.4f}")
    
    # Compute bin nuclei counts
    bin_counts = compute_bin_nuclei_counts(nuc_data_with_bins, n_ap_bins=n_ap_bins, n_dv_bins=n_dv_bins)
    print(f"  Computed bin nuclei counts: mean = {np.mean(bin_counts):.1f}, range = [{min(bin_counts)}, {max(bin_counts)}]")
    
    # Validate metrics
    nn_passed, nn_warnings = validate_nn_metrics(nn_metrics, thresholds)
    bin_passed, bin_warnings = validate_bin_counts(bin_counts, thresholds, n_ap_bins, n_dv_bins)
    
    overall_passed = nn_passed and bin_passed
    all_warnings = nn_warnings + bin_warnings
    
    # Write binned mRNA output (for inference)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    binned_data['avg_mrna_count'].to_csv(
        output_file,
        index=False,
        header=False
    )
    print(f"  Wrote output to {output_file}")
    
    # Generate visualizations
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    # If trimming was applied, compute original data for comparison
    binned_data_original = None
    spatial_ranges_original = None
    bin_counts_original = None
    
    if trim_info.get('applied', False):
        print("  Computing original data for comparison plots...")
        binned_data_original, nuc_data_original, spatial_ranges_original = bin_mrna_data(
            df,  # Original untrimmed data
            n_ap_bins=n_ap_bins,
            n_dv_bins=n_dv_bins,
            no_groupByNuclei=no_groupByNuclei
        )
        bin_counts_original = compute_bin_nuclei_counts(nuc_data_original, n_ap_bins=n_ap_bins, n_dv_bins=n_dv_bins)
    
    # 1. Heatmap
    heatmap_path = fig_dir / f"{embryo_id}_sass_formodel_heatmap.png"
    plot_heatmap(binned_data, str(heatmap_path), n_ap_bins, n_dv_bins, binned_data_original)
    print(f"  Wrote heatmap to {heatmap_path}")
    
    # 2. X-Y spatial distribution with grid overlay
    scatter_path = fig_dir / f"{embryo_id}_sass_formodel_spatial_xy.png"
    plot_spatial_distribution(df_processed, str(scatter_path), spatial_ranges, n_ap_bins, n_dv_bins,
                            df_original=df if trim_info.get('applied', False) else None,
                            spatial_ranges_original=spatial_ranges_original)
    print(f"  Wrote X-Y scatter plot to {scatter_path}")
    
    # 3. Expression density (spot-level)
    density_path = fig_dir / f"{embryo_id}_sass_formodel_expression_density.png"
    plot_expression_density(df_processed, str(density_path), 
                          df_original=df if trim_info.get('applied', False) else None)
    print(f"  Wrote expression density plot to {density_path}")

    # 3b. Expression density (nuclei-level)
    nuclei_density_path = fig_dir / f"{embryo_id}_sass_formodel_nuclei_expression_density.png"
    plot_nuclei_expression_density(df_processed, str(nuclei_density_path),
                                  df_original=df if trim_info.get('applied', False) else None)
    print(f"  Wrote nuclei-level expression density plot to {nuclei_density_path}")
    
    # 4. Ridge plot (bin count distributions)
    ridge_path = fig_dir / f"{embryo_id}_bin_count_ridge.png"
    plot_bin_count_ridge(bin_counts, thresholds, embryo_id, str(ridge_path), bin_counts_original)
    print(f"  Wrote ridge plot to {ridge_path}")
    
    # Write validation report
    report_path = Path(validation_report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    
    status = "PASS" if overall_passed else "FAIL"
    
    with open(validation_report, 'w') as f:
        f.write(f"mRNA Data Validation Report\n")
        f.write(f"=" * 80 + "\n\n")
        f.write(f"Input file: {position_data_file}\n")
        f.write(f"Stripe: {stripe}\n")
        f.write(f"Embryo: {embryo_id}\n")
        f.write(f"Status: {status}\n\n")
        
        # Report AP axis trimming if applied
        if trim_info.get('applied', False):
            f.write(f"AP Axis Trimming Applied:\n")
            f.write(f"  Reason: Nuclei too dense (median NN distance below threshold)\n")
            f.write(f"  Iterations: {trim_info['iterations']}\n")
            f.write(f"  Original AP range: [{trim_info['original_ap_range'][0]:.2f}, {trim_info['original_ap_range'][1]:.2f}]\n")
            f.write(f"  Trimmed AP range: [{trim_info['trimmed_ap_range'][0]:.2f}, {trim_info['trimmed_ap_range'][1]:.2f}]\n")
            f.write(f"  Total trimmed: {trim_info['total_trim_percent']:.1f}% of AP axis\n")
            f.write(f"  Original nuclei count: {trim_info['original_n_nuclei']}\n")
            f.write(f"  Final nuclei count: {trim_info['final_n_nuclei']}\n")
            f.write(f"  Original median NN distance: {trim_info['original_median_nn']:.4f}\n")
            f.write(f"  Final median NN distance: {trim_info['final_median_nn']:.4f}\n\n")
        
        f.write(f"Nuclei Density Metrics (k={k} NN):\n")
        f.write(f"  Number of nuclei: {nn_metrics['n_nuclei']}\n")
        f.write(f"  Average k={k} NN distance: {nn_metrics['avg_nn_distance']:.4f}\n")
        f.write(f"  Median k={k} NN distance: {nn_metrics['median_nn_distance']:.4f}\n")
        f.write(f"  Std k={k} NN distance: {nn_metrics['std_nn_distance']:.4f}\n")
        f.write(f"  Range: {nn_metrics['min_nn_distance']:.4f} - {nn_metrics['max_nn_distance']:.4f}\n\n")
        
        f.write(f"Bin Nuclei Count Statistics:\n")
        f.write(f"  Mean: {np.mean(bin_counts):.1f}\n")
        f.write(f"  Std: {np.std(bin_counts):.1f}\n")
        f.write(f"  Range: {min(bin_counts)} - {max(bin_counts)}\n\n")
        
        f.write(f"Validation Criteria ({stripe} - Berrocal_2020):\n")
        f.write(f"  Expected nuclei count: ~{thresholds['expected_nuclei_count_mean']:.1f}\n")
        f.write(f"  Expected median NN distance: ~{thresholds['expected_nn_distance_median']:.4f}\n")
        f.write(f"  Valid NN range: {thresholds['min_nn_distance']:.4f} - {thresholds['max_nn_distance']:.4f}\n")
        f.write(f"  Based on {thresholds['n_embryos_used']} transcription embryos\n\n")
        
        if all_warnings:
            f.write(f"Validation Messages:\n")
            for warning in all_warnings:
                f.write(f"  {warning}\n")
    
    print(f"\n{'✓' if overall_passed else '✗'} Validation {'PASSED' if overall_passed else 'FAILED'}")
    print(f"  Report written to {validation_report}")
    
    # Exit with error code if validation failed
    sys.exit(0 if overall_passed else 1)


if __name__ == '__main__':
    main()
