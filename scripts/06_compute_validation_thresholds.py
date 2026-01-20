#!/usr/bin/env python3
"""
Compute nuclei density validation thresholds from Berrocal_2020 data.

This script analyzes the Berrocal_2020 transcription data to compute expected
nuclei density metrics for each stripe. These metrics are used to validate
processed mRNA data quality.

For each stripe, computes:
- Expected nuclei count (mean across embryos)
- k=4 nearest neighbor distance statistics (mean, median)
- Valid range (min/max observed across embryos)

Results are saved to config.yaml under 'validation_thresholds' for each stripe.

Usage:
    python 06_compute_validation_thresholds.py --input <berrocal_data.csv> --config <config.yaml>
    
References:
    Based on analysis in notebooks/afterHandoff/exploreTranscriptionData.ipynb
    (cells 9-10)
"""

import argparse
import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedSeq

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_berrocal_data(input_path: str) -> pd.DataFrame:
    """Load and prepare Berrocal_2020 data."""
    logger.info(f"Loading Berrocal_2020 data from {input_path}")
    df = pd.read_csv(input_path)
    logger.info(f"  Loaded {len(df)} records")
    return df


def prepare_position_data(df: pd.DataFrame, max_time_seconds: int = 1200) -> pd.DataFrame:
    """
    Aggregate fluorescence and position data per nucleus.
    
    Args:
        df: Raw Berrocal data with time-series
        max_time_seconds: Maximum time to include in sum (default: 1200s)
    
    Returns:
        DataFrame with one row per nucleus containing position and summed fluorescence
    """
    # Pivot fluorescence traces
    fluo_traces = df.pivot(index='nucleus_id', columns='time', values='fluo')
    fluo_traces = fluo_traces.fillna(0)
    fluo_traces.reset_index(inplace=True)
    
    # Sum fluorescence up to max_time (time columns are 0, 20, 40, ...)
    time_cols = [col for col in fluo_traces.columns if isinstance(col, (int, float))]
    time_cols_in_range = [col for col in time_cols if col <= max_time_seconds]
    fluo_traces['sum_fluo_b41200'] = fluo_traces[time_cols_in_range].sum(axis=1)
    
    # Get position data (median position over time window)
    pos_data = df.drop(columns=['particle_id', 'set_id', 'ap_raw', 'burst_stripe_id', 
                                 'time', 'fluo', 'inference_stripe_id', 'inference_fluo_id', 
                                 'v_state', 'v_fluo'])
    pos_data = pos_data.groupby('nucleus_id').median().reset_index()
    
    # Merge fluorescence and position
    pos_sum_df = fluo_traces[['nucleus_id', 'sum_fluo_b41200']].merge(pos_data, on='nucleus_id')
    
    return pos_sum_df


def compute_bin_nuclei_counts(nuclei_data: pd.DataFrame, n_ap_bins: int = 5, n_dv_bins: int = 5) -> list[int]:
    """
    Compute the number of nuclei in each spatial bin.
    
    Follows identical binning logic to 01_preprocess_eve_data.py to ensure
    consistency between transcription and mRNA data binning.
    
    Args:
        nuclei_data: DataFrame with ap_registered, yPos columns (already filtered to stripe)
        n_ap_bins: Number of bins along AP axis (default: 5)
        n_dv_bins: Number of bins along DV axis (default: 5)
    
    Returns:
        List of 25 nuclei counts (one per bin), ordered by [apBin, yBin]
    """
    # Normalize AP and DV positions to [0, 1]
    minap = nuclei_data['ap_registered'].min()
    maxap = nuclei_data['ap_registered'].max()
    nuclei_data = nuclei_data.copy()
    nuclei_data['ap_reg_norm'] = (nuclei_data['ap_registered'] - minap) / (maxap - minap)
    
    miny = nuclei_data['yPos'].min()
    maxy = nuclei_data['yPos'].max()
    nuclei_data['yPos_norm'] = (nuclei_data['yPos'] - miny) / (maxy - miny)
    
    # Create DV bins (following 01_preprocess_eve_data.py logic)
    dv_edges = np.linspace(nuclei_data['yPos_norm'].min(), nuclei_data['yPos_norm'].max(), n_dv_bins + 1)
    for i in range(n_dv_bins):
        if i == 0:
            mask = nuclei_data['yPos_norm'].between(dv_edges[i], dv_edges[i+1], 'both')
        else:
            mask = nuclei_data['yPos_norm'].between(dv_edges[i], dv_edges[i+1], 'right')
        nuclei_data.loc[mask, 'yBin'] = i + 1
    
    # Create AP bins
    ap_edges = np.linspace(0, 1, n_ap_bins + 1)
    for i in range(n_ap_bins):
        if i == 0:
            mask = nuclei_data['ap_reg_norm'].between(ap_edges[i], ap_edges[i+1], 'both')
        else:
            mask = nuclei_data['ap_reg_norm'].between(ap_edges[i], ap_edges[i+1], 'right')
        nuclei_data.loc[mask, 'apBin'] = i + 1
    
    # Convert to 0-indexed for consistency with other scripts
    nuclei_data['apBin'] = (nuclei_data['apBin'] - 1).astype(int)
    nuclei_data['yBin'] = (nuclei_data['yBin'] - 1).astype(int)
    
    # Count nuclei per bin
    bin_counts = nuclei_data.groupby(['apBin', 'yBin']).size().reset_index(name='nuclei_count')
    
    # Ensure all bins are represented
    all_bins = pd.DataFrame([
        (ap, y) for ap in range(n_ap_bins) for y in range(n_dv_bins)
    ], columns=['apBin', 'yBin'])
    
    bin_counts = all_bins.merge(bin_counts, on=['apBin', 'yBin'], how='left')
    bin_counts['nuclei_count'] = bin_counts['nuclei_count'].fillna(0).astype(int)
    
    # Sort by apBin, then yBin for consistent ordering
    bin_counts = bin_counts.sort_values(['apBin', 'yBin']).reset_index(drop=True)
    
    return bin_counts['nuclei_count'].tolist()


def compute_nn_metrics_for_embryo(stripe_data: pd.DataFrame, embryo_id: int, k: int = 4, n_ap_bins: int = 5, n_dv_bins: int = 5) -> tuple:
    """
    Compute k-nearest neighbor metrics and bin nuclei counts for a single embryo.
    
    Follows methodology from exploreTranscriptionData.ipynb cell 9.
    
    Args:
        stripe_data: DataFrame with ap_registered, yPos columns for all embryos in stripe
        embryo_id: Embryo ID to process
        k: Number of nearest neighbors
        n_ap_bins: Number of bins along AP axis (default: 5)
        n_dv_bins: Number of bins along DV axis (default: 5)
    
    Returns:
        Tuple of (metrics_dict, normalized_data_df, bin_nuclei_counts)
    """
    embryo_data = stripe_data[stripe_data['embryo'] == embryo_id]

    if len(embryo_data) < 10: #less than 10 nuclei, skip
        return None, None, None
    
    # Normalize coordinates to [0, 1]
    x_min = stripe_data['ap_registered'].min()
    x_max = stripe_data['ap_registered'].max()
    y_min = stripe_data['yPos'].min()
    y_max = stripe_data['yPos'].max()
    
    embryo_data = embryo_data.copy()
    embryo_data['ap_norm'] = (embryo_data['ap_registered'] - x_min) / (x_max - x_min)
    embryo_data['yPos_norm'] = (embryo_data['yPos'] - y_min) / (y_max - y_min)
    
    # Compute k-nearest neighbors
    coords = embryo_data[['ap_norm', 'yPos_norm']].values
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(coords)
    distances, _ = nbrs.kneighbors(coords)
    
    # Remove first column (distance to self = 0)
    distances = distances[:, 1:]
    
    metrics = {
        'n_nuclei': len(embryo_data),
        'avg_nn_distance': np.mean(distances),
        'std_nn_distance': np.std(distances),
        'median_nn_distance': np.median(distances)
    }
    
    # Compute bin nuclei counts for this embryo
    bin_counts = compute_bin_nuclei_counts(embryo_data, n_ap_bins, n_dv_bins)
    
    return metrics, embryo_data[['ap_norm', 'yPos_norm']], bin_counts


def compute_validation_thresholds_for_stripe(
    pos_sum_df: pd.DataFrame, 
    stripe_name: str,
    ap_min: float, 
    ap_max: float,
    k: int = 4,
    n_ap_bins: int = 5,
    n_dv_bins: int = 5
) -> tuple:
    """
    Compute validation thresholds for a specific stripe.
    
    Args:
        pos_sum_df: DataFrame with all nuclei positions
        stripe_name: Name of stripe (e.g., 'stripe2')
        ap_min: Minimum AP coordinate for this stripe
        ap_max: Maximum AP coordinate for this stripe
        k: Number of nearest neighbors
        n_ap_bins: Number of bins along AP axis (default: 5)
        n_dv_bins: Number of bins along DV axis (default: 5)
    
    Returns:
        Tuple of (thresholds_dict, embryo_normalized_data_dict)
    """
    logger.info(f"Computing validation thresholds for {stripe_name} (AP: {ap_min:.3f}-{ap_max:.3f})")
    
    # Get unique embryo IDs
    embryo_ids = pos_sum_df['embryo'].unique()
    logger.info(f"  Found {len(embryo_ids)} embryos")
    
    results = []
    embryo_normalized_data = {}
    all_bin_counts = []
    
    for embryo_id in embryo_ids:
        # Filter for this stripe region
        stripe_data = pos_sum_df[
            (pos_sum_df['ap_registered'] > ap_min) & 
            (pos_sum_df['ap_registered'] < ap_max)
        ].copy()
        
        metrics, norm_data, bin_counts = compute_nn_metrics_for_embryo(stripe_data, embryo_id, k=k, n_ap_bins=n_ap_bins, n_dv_bins=n_dv_bins)
        if metrics:
            results.append({
                'embryo': embryo_id,
                **metrics
            })
            embryo_normalized_data[embryo_id] = {
                'norm_data': norm_data, 
                'metrics': metrics,
                'bin_counts': bin_counts
            }
            all_bin_counts.append(bin_counts)
            logger.info(f"    Embryo {embryo_id}: {metrics['n_nuclei']} nuclei, "
                       f"median NN={metrics['median_nn_distance']:.4f}, "
                       f"bin counts range=[{min(bin_counts)}, {max(bin_counts)}]")
    
    if not results:
        logger.warning(f"  No valid embryo data found for {stripe_name}")
        return None, None
    
    results_df = pd.DataFrame(results)
    
    # Compute aggregate k-NN statistics
    thresholds = {
        'expected_nuclei_count_mean': float(results_df['n_nuclei'].mean()),
        'expected_nuclei_count_std': float(results_df['n_nuclei'].std()),
        'expected_nn_distance_mean': float(results_df['avg_nn_distance'].mean()),
        'expected_nn_distance_median': float(results_df['median_nn_distance'].mean()),
        'min_nn_distance': float(results_df['median_nn_distance'].min()),
        'max_nn_distance': float(results_df['median_nn_distance'].max()),
        'k': k,
        'n_embryos_used': len(results)
    }
    
    # Compute bin count distribution statistics
    # Store per-embryo counts and combined distribution
    bin_count_distribution = {}
    
    # Per-embryo counts
    for embryo_id, data in embryo_normalized_data.items():
        bin_count_distribution[str(embryo_id)] = data['bin_counts']
    
    # All combined (flatten all embryo bin counts)
    all_counts_flat = [count for counts in all_bin_counts for count in counts]
    bin_count_distribution['all_combined'] = all_counts_flat
    
    # Statistics across bins (mean, std, q25, q75 per bin position)
    bin_counts_array = np.array(all_bin_counts)  # Shape: (n_embryos, 25)
    bin_count_distribution['stats'] = {
        'mean': bin_counts_array.mean(axis=0).tolist(),  # Mean count per bin position across embryos
        'std': bin_counts_array.std(axis=0).tolist(),    # Std per bin position
        'q25': np.percentile(bin_counts_array, 25, axis=0).tolist(),
        'q75': np.percentile(bin_counts_array, 75, axis=0).tolist()
    }
    
    thresholds['bin_count_distribution'] = bin_count_distribution
    
    logger.info(f"  {stripe_name} thresholds:")
    logger.info(f"    Expected nuclei: {thresholds['expected_nuclei_count_mean']:.1f} ± {thresholds['expected_nuclei_count_std']:.1f}")
    logger.info(f"    Expected median NN: {thresholds['expected_nn_distance_median']:.4f}")
    logger.info(f"    Valid range: [{thresholds['min_nn_distance']:.4f}, {thresholds['max_nn_distance']:.4f}]")
    logger.info(f"    Bin count stats: mean={np.mean(bin_count_distribution['stats']['mean']):.1f}, "
               f"overall range=[{np.min(bin_counts_array)}, {np.max(bin_counts_array)}]")
    
    return thresholds, embryo_normalized_data


def plot_stripe_nuclei_placement(stripe_name: str, embryo_data: dict, output_path: str):
    """
    Create multipanel figure showing normalized nuclei placement for all embryos in a stripe.
    
    Args:
        stripe_name: Name of the stripe
        embryo_data: Dict of {embryo_id: normalized_data_df}
        output_path: Path to save the figure
    """
    logger.info(f"Creating {stripe_name} nuclei placement figure at {output_path}")
    
    n_embryos = len(embryo_data)
    if n_embryos == 0:
        logger.warning(f"  No embryo data for {stripe_name}")
        return
    
    # Determine grid layout (prefer square-ish layout)
    ncols = int(np.ceil(np.sqrt(n_embryos)))
    nrows = int(np.ceil(n_embryos / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
    axes = axes.flatten() if n_embryos > 1 else [axes]
    
    plot_idx = 0
    for embryo_id, data in sorted(embryo_data.items()):
        ax = axes[plot_idx]
        
        norm_data = data['norm_data']
        metrics = data['metrics']
        median_nn = metrics['median_nn_distance']
        
        # Scatter plot of normalized nuclei positions
        ax.scatter(norm_data['ap_norm'], norm_data['yPos_norm'], 
                  s=10, alpha=0.6, color='steelblue')
        
        ax.set_xlabel('Normalized AP position', fontsize=10)
        ax.set_ylabel('Normalized DV position', fontsize=10)
        ax.set_title(f'{embryo_id}\n(n={len(norm_data)}, median NN={median_nn:.4f})', fontsize=11)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        plot_idx += 1
    
    # Hide unused subplots
    for idx in range(plot_idx, len(axes)):
        axes[idx].axis('off')
    
    # Add overall title for the stripe
    fig.suptitle(f'{stripe_name} - Nuclei Placement (Normalized Coordinates)', 
                fontsize=14, fontweight='bold', y=1.0)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  Saved {stripe_name} figure with {plot_idx} embryo panels")


def update_config_with_thresholds(config_path: str, validation_thresholds: dict):
    """
    Update config.yaml with computed validation thresholds.
    
    Args:
        config_path: Path to config.yaml
        validation_thresholds: Dictionary of thresholds per stripe
    """
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.default_flow_style = False
    
    # Load existing config
    config_path = Path(config_path)
    with open(config_path, 'r') as f:
        config = yaml.load(f)
    
    # Add validation thresholds section
    config['validation_thresholds'] = validation_thresholds
    
    # Add metadata
    config['validation_thresholds']['metadata'] = {
        'last_updated': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'source': 'Berrocal_2020 data',
        'script': '06_compute_validation_thresholds.py'
    }
    # Convert all top-level list values inside bin_count_distribution
    # to inline (flow-style) `CommentedSeq` so YAML emits them on one line.
    try:
        for stripe, stripe_thresh in config['validation_thresholds'].items():
            if stripe == 'metadata':
                continue
            if not isinstance(stripe_thresh, dict):
                continue
            bcd = stripe_thresh.get('bin_count_distribution')
            if not isinstance(bcd, dict):
                continue
            # Iterate over a static list of items to avoid runtime mutation issues
            for k, v in list(bcd.items()):
                if isinstance(v, list):
                    cs = CommentedSeq(v)
                    try:
                        cs.fa.set_flow_style()
                    except Exception:
                        try:
                            cs.yaml_set_flow_style()
                        except Exception:
                            pass
                    bcd[k] = cs
    except Exception as e:
        logger.debug(f"Could not convert lists to inline flow style: {e}")
    
    # Write updated config
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    
    logger.info(f"Updated {config_path} with validation thresholds")


def main():
    parser = argparse.ArgumentParser(
        description='Compute nuclei density validation thresholds from Berrocal_2020 data'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Path to Berrocal_2020 data CSV file'
    )
    parser.add_argument(
        '--config',
        required=True,
        help='Path to config.yaml'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to save multipanel nuclei placement figures (one per stripe)'
    )
    parser.add_argument(
        '--k',
        type=int,
        default=4,
        help='Number of nearest neighbors (default: 4)'
    )
    parser.add_argument(
        '--max-time',
        type=int,
        default=1200,
        help='Maximum time in seconds for fluorescence sum (default: 1200)'
    )
    parser.add_argument(
        '--n-ap-bins',
        type=int,
        default=5,
        help='Number of bins along AP axis (default: 5)'
    )
    parser.add_argument(
        '--n-dv-bins',
        type=int,
        default=5,
        help='Number of bins along DV axis (default: 5)'
    )
    
    args = parser.parse_args()
    
    # Load data
    df = load_berrocal_data(args.input)
    pos_sum_df = prepare_position_data(df, max_time_seconds=args.max_time)
    
    # Load config to get stripe ranges
    yaml = YAML()
    with open(args.config, 'r') as f:
        config = yaml.load(f)
    
    stripe_ranges = config.get('stripe_ranges', {})
    if not stripe_ranges:
        logger.error("No stripe_ranges found in config.yaml. Run 00_identify_stripe_ranges.py first.")
        return 1
    
    # Compute thresholds for each stripe
    validation_thresholds = {}
    all_stripe_data = {}
    for stripe_name, ranges in stripe_ranges.items():
        if isinstance(ranges, dict) and 'min' in ranges and 'max' in ranges:
            thresholds, embryo_data = compute_validation_thresholds_for_stripe(
                pos_sum_df,
                stripe_name,
                ranges['min'],
                ranges['max'],
                k=args.k,
                n_ap_bins=args.n_ap_bins,
                n_dv_bins=args.n_dv_bins
            )
            if thresholds:
                validation_thresholds[stripe_name] = thresholds
                all_stripe_data[stripe_name] = embryo_data
    
    # Update config
    if validation_thresholds:
        logger.info(f"\nComputed thresholds for {len(validation_thresholds)} stripes")
        update_config_with_thresholds(args.config, validation_thresholds)
        
        # Create one multipanel figure per stripe
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for stripe_name, embryo_data in all_stripe_data.items():
            output_path = output_dir / f"{stripe_name}_nuclei_placement.png"
            plot_stripe_nuclei_placement(stripe_name, embryo_data, str(output_path))
        
        return 0
    else:
        logger.error("No validation thresholds computed")
        return 1


if __name__ == '__main__':
    exit(main())
