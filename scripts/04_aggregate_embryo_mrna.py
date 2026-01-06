#!/usr/bin/env python3
"""
Aggregate mRNA counts from SASS position_data-intense.txt into spatial bins.

This script bins mRNA spot data into a 5×5 spatial grid matching the transcription
data binning strategy, computing average mRNA counts per nucleus in each bin.

Usage:
    python 04_aggregate_embryo_mrna.py <position_data_file> <stripe> <output_file>

The script uses identical spatial binning logic to 01_preprocess_eve_data.py to
ensure mRNA and transcription data are spatially aligned.
"""

import sys
import pandas as pd
import numpy as np
import yaml
from pathlib import Path


def load_stripe_ranges(config_path: str) -> dict:
    """Load stripe AP coordinate ranges from config.yaml."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config['stripe_ranges']


def load_position_data(filepath: str) -> pd.DataFrame:
    """
    Load SASS position_data-intense.txt file.
    
    Expected format: tab-separated with columns including:
    - nucx, nucy, nucz: nuclear positions
    - spotx, spoty, spotz: spot positions
    - nuc: nuclear ID
    - num_spots: number of spots per nucleus
    """
    df = pd.read_csv(filepath, sep='\t')
    return df


def bin_mrna_data(df: pd.DataFrame, ap_min: float, ap_max: float, 
                  n_ap_bins: int = 5, n_dv_bins: int = 5) -> pd.DataFrame:
    """
    Bin mRNA data into spatial grid matching transcription data binning.
    
    Strategy (matching 01_preprocess_eve_data.py):
    1. Filter nuclei to stripe AP range
    2. Normalize AP and DV coordinates to [0, 1]
    3. Create n_ap_bins × n_dv_bins spatial grid
    4. Compute average mRNA count per nucleus in each bin
    
    Args:
        df: DataFrame with nuclear positions and spot counts
        ap_min: Minimum AP coordinate for stripe
        ap_max: Maximum AP coordinate for stripe
        n_ap_bins: Number of bins along AP axis (default: 5)
        n_dv_bins: Number of bins along DV axis (default: 5)
    
    Returns:
        DataFrame with columns: apBin, yBin, avg_mrna_count
    """
    # Aggregate to one row per nucleus with total spot count
    nuc_data = df.groupby('nuc').agg({
        'nucx': 'first',  # AP coordinate (assuming nucx is AP)
        'nucy': 'first',  # DV coordinate (assuming nucy is DV)
        'num_spots': 'max'  # Total spots per nucleus
    }).reset_index()
    
    # Filter to stripe AP range
    nuc_data_filtered = nuc_data[
        (nuc_data['nucx'] >= ap_min) & 
        (nuc_data['nucx'] <= ap_max)
    ].copy()
    
    if len(nuc_data_filtered) == 0:
        raise ValueError(f"No nuclei found in AP range [{ap_min}, {ap_max}]")
    
    # Normalize coordinates to [0, 1]
    ap_range = ap_max - ap_min
    dv_min = nuc_data_filtered['nucy'].min()
    dv_max = nuc_data_filtered['nucy'].max()
    dv_range = dv_max - dv_min
    
    nuc_data_filtered['ap_norm'] = (nuc_data_filtered['nucx'] - ap_min) / ap_range
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
    
    return binned_data


def main():
    if len(sys.argv) != 4:
        print("Usage: python 04_aggregate_embryo_mrna.py <position_data_file> <stripe> <output_file>")
        sys.exit(1)
    
    position_data_file = sys.argv[1]
    stripe = sys.argv[2]
    output_file = sys.argv[3]
    
    # Load config to get stripe ranges
    config_path = Path(__file__).parent.parent / 'config.yaml'
    stripe_ranges = load_stripe_ranges(str(config_path))
    
    if stripe not in stripe_ranges:
        print(f"Error: Stripe '{stripe}' not found in config.yaml")
        sys.exit(1)
    
    ap_min = stripe_ranges[stripe]['min']
    ap_max = stripe_ranges[stripe]['max']
    
    print(f"Processing {position_data_file}")
    print(f"  Stripe: {stripe} (AP range: {ap_min} - {ap_max})")
    
    # Load and process data
    df = load_position_data(position_data_file)
    print(f"  Loaded {len(df)} spot records")
    
    binned_data = bin_mrna_data(df, ap_min, ap_max)
    print(f"  Binned into {len(binned_data)} spatial bins")
    
    # Write output in format expected by inference script (25 values, no header)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    binned_data['avg_mrna_count'].to_csv(
        output_file,
        index=False,
        header=False
    )
    
    print(f"  Wrote output to {output_file}")
    print(f"  Mean mRNA count: {binned_data['avg_mrna_count'].mean():.2f}")
    print(f"  Range: {binned_data['avg_mrna_count'].min():.2f} - {binned_data['avg_mrna_count'].max():.2f}")


if __name__ == '__main__':
    main()
