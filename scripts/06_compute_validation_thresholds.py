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
from sklearn.neighbors import NearestNeighbors
from ruamel.yaml import YAML

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


def compute_nn_metrics_for_embryo(embryo_data: pd.DataFrame, k: int = 4) -> dict:
    """
    Compute k-nearest neighbor metrics for a single embryo.
    
    Follows methodology from exploreTranscriptionData.ipynb cell 9.
    
    Args:
        embryo_data: DataFrame with ap_registered, yPos columns
        k: Number of nearest neighbors
    
    Returns:
        Dictionary with nuclei count and NN distance metrics
    """
    if len(embryo_data) == 0:
        return None
    
    # Normalize coordinates to [0, 1]
    x_min = embryo_data['ap_registered'].min()
    x_max = embryo_data['ap_registered'].max()
    y_min = embryo_data['yPos'].min()
    y_max = embryo_data['yPos'].max()
    
    embryo_data = embryo_data.copy()
    embryo_data['ap_norm'] = (embryo_data['ap_registered'] - x_min) / (x_max - x_min)
    embryo_data['yPos_norm'] = (embryo_data['yPos'] - y_min) / (y_max - y_min)
    
    # Compute k-nearest neighbors
    coords = embryo_data[['ap_norm', 'yPos_norm']].values
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(coords)
    distances, _ = nbrs.kneighbors(coords)
    
    # Remove first column (distance to self = 0)
    distances = distances[:, 1:]
    
    return {
        'n_nuclei': len(embryo_data),
        'avg_nn_distance': np.mean(distances),
        'std_nn_distance': np.std(distances),
        'median_nn_distance': np.median(distances)
    }


def compute_validation_thresholds_for_stripe(
    pos_sum_df: pd.DataFrame, 
    stripe_name: str,
    ap_min: float, 
    ap_max: float,
    k: int = 4
) -> dict:
    """
    Compute validation thresholds for a specific stripe.
    
    Args:
        pos_sum_df: DataFrame with all nuclei positions
        stripe_name: Name of stripe (e.g., 'stripe2')
        ap_min: Minimum AP coordinate for this stripe
        ap_max: Maximum AP coordinate for this stripe
        k: Number of nearest neighbors
    
    Returns:
        Dictionary with validation thresholds
    """
    logger.info(f"Computing validation thresholds for {stripe_name} (AP: {ap_min:.3f}-{ap_max:.3f})")
    
    # Get unique embryo IDs
    embryo_ids = pos_sum_df['embryo'].unique()
    logger.info(f"  Found {len(embryo_ids)} embryos")
    
    results = []
    for embryo_id in embryo_ids:
        # Filter for this embryo and stripe region
        embryo_data = pos_sum_df[
            (pos_sum_df['embryo'] == embryo_id) & 
            (pos_sum_df['ap_registered'] > ap_min) & 
            (pos_sum_df['ap_registered'] < ap_max)
        ].copy()
        
        metrics = compute_nn_metrics_for_embryo(embryo_data, k=k)
        if metrics:
            results.append({
                'embryo': embryo_id,
                **metrics
            })
            logger.info(f"    Embryo {embryo_id}: {metrics['n_nuclei']} nuclei, "
                       f"median NN={metrics['median_nn_distance']:.4f}")
    
    if not results:
        logger.warning(f"  No valid embryo data found for {stripe_name}")
        return None
    
    results_df = pd.DataFrame(results)
    
    # Compute aggregate statistics
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
    
    logger.info(f"  {stripe_name} thresholds:")
    logger.info(f"    Expected nuclei: {thresholds['expected_nuclei_count_mean']:.1f} ± {thresholds['expected_nuclei_count_std']:.1f}")
    logger.info(f"    Expected median NN: {thresholds['expected_nn_distance_median']:.4f}")
    logger.info(f"    Valid range: [{thresholds['min_nn_distance']:.4f}, {thresholds['max_nn_distance']:.4f}]")
    
    return thresholds


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
    for stripe_name, ranges in stripe_ranges.items():
        if isinstance(ranges, dict) and 'min' in ranges and 'max' in ranges:
            thresholds = compute_validation_thresholds_for_stripe(
                pos_sum_df,
                stripe_name,
                ranges['min'],
                ranges['max'],
                k=args.k
            )
            if thresholds:
                validation_thresholds[stripe_name] = thresholds
    
    # Update config
    if validation_thresholds:
        update_config_with_thresholds(args.config, validation_thresholds)
        logger.info(f"\nComputed thresholds for {len(validation_thresholds)} stripes")
        return 0
    else:
        logger.error("No validation thresholds computed")
        return 1


if __name__ == '__main__':
    exit(main())
