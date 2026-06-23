#!/usr/bin/env python
"""
Preprocess eve transcription data for spatial mRNA decay inference.

Filters nuclei by stripe membership, bins them spatially along the AP axis
using overlapping bins (50% overlap by default), and computes per-bin median
fluorescence traces. Outputs a CSV suitable for use as the transcription input
F(t) to the degradation inference scripts (02_infer_*).

Pipeline stage: Step 01 — runs after stripe range identification (00_*),
before inference (02_*).

Inputs (via CLI args, see argparse section at bottom):
    --input   : Path to filtered eve CSV (Berrocal_2020 / longform format)
    --config  : Path to config.yaml with stripe range definitions
    --stripe  : Stripe number to process (e.g. 2 or 3)
    --output  : Output CSV path (n_timepoints x n_bins, no header)

"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cmcrameri.cm as cmc
import argparse
from ruamel.yaml import YAML


def load_data(filtered_data_path):
    """Load filtered eve transcription data."""
    print(f"Loading data from: {filtered_data_path}")
    data_filtered = pd.read_csv(filtered_data_path, header=0)
    print(f"Loaded {len(data_filtered)} rows")
    return data_filtered


def extract_fluorescence_traces(data_filtered, max_time_seconds=1200):
    """
    Extract fluorescence traces for each nucleus.
    
    Args:
        data_filtered: Filtered eve data
        max_time_seconds: Maximum time in seconds (default 1200 = 20 minutes)
    
    Returns:
        DataFrame with nucleus_id and fluorescence over time
    """
    print(f"Extracting fluorescence traces up to {max_time_seconds} seconds...")
    
    # Pivot to get fluorescence traces per nucleus
    fluo_traces = data_filtered.pivot(index='nucleus_id', columns='time', values='fluo')
    fluo_traces = fluo_traces.fillna(0)
    fluo_traces.reset_index(inplace=True)
    
    # Select time points up to max_time_seconds (each timepoint is ~20s)
    max_timepoint_idx = int(max_time_seconds / 20) + 2  # +2 for index and nucleus_id columns
    fluo_traces = fluo_traces.iloc[:, 0:max_timepoint_idx]
    
    print(f"Extracted traces for {len(fluo_traces)} nuclei")
    return fluo_traces


def extract_nuclear_positions(data_filtered):
    """
    Extract median nuclear positions over time window.
    
    Returns:
        DataFrame with nucleus_id and median ap_registered, yPos
    """
    print("Extracting nuclear positions...")
    
    pos_data = data_filtered.drop(columns=[
        'particle_id', 'set_id', 'ap_raw', 'burst_stripe_id', 
        'time', 'fluo', 'inference_stripe_id', 'inference_fluo_id', 
        'v_state', 'v_fluo'
    ])
    
    # Assign median position over time window
    pos_data = pos_data.groupby('nucleus_id').median()
    pos_data.reset_index(inplace=True)
    
    return pos_data


def create_overlapping_bin_ranges(n_bins, min_val, max_val, overlap_fraction=0.5):
    """
    Create overlapping bin ranges.
    
    Args:
        n_bins: Number of bins
        min_val: Minimum value of the axis
        max_val: Maximum value of the axis
        overlap_fraction: Fraction of bin width to extend on each side (0.5 = 50%)
    
    Returns:
        List of (bin_min, bin_max) tuples for each bin (1-indexed)
    """
    bin_width = (max_val - min_val) / n_bins
    extension = bin_width * overlap_fraction
    
    bins = []
    for i in range(n_bins):
        bin_center = min_val + (i + 0.5) * bin_width
        
        # Extend on both sides, but clamp to min/max and only extend for non-edge bins
        if i == 0:  # Left edge bin - only extend right
            bin_min = min_val
            bin_max = min(max_val, bin_center + bin_width/2 + extension)
        elif i == n_bins - 1:  # Right edge bin - only extend left
            bin_min = max(min_val, bin_center - bin_width/2 - extension)
            bin_max = max_val
        else:  # Middle bins - extend both directions
            bin_min = bin_center - bin_width/2 - extension
            bin_max = bin_center + bin_width/2 + extension
        
        bins.append((bin_min, bin_max))
    
    return bins


def filter_stripe2_nuclei(pos_data, ap_min, ap_max, stripe):
    """
    Filter nuclei in stripe and adjacent interstripes.
    
    Args:
        pos_data: Position data for all nuclei
        ap_min: Minimum ap_registered value (default 0.33)
        ap_max: Maximum ap_registered value (default 0.45)
    
    Returns:
        Filtered position data for stripe
    """
    print(f"Filtering {stripe} nuclei (AP range: {ap_min} - {ap_max})...")
    
    pos_data_str2 = pos_data[
        (pos_data['ap_registered'] > ap_min) & 
        (pos_data['ap_registered'] < ap_max)
    ]
    
    print(f"Found {len(pos_data_str2)} nuclei in {stripe}")
    return pos_data_str2


def bin_nuclei_spatially(nuclei_in_str2, n_ap_bins, n_dv_bins, 
                         dv_crop_min=-1.0, dv_crop_max=-1.0):
    """
    Bin nuclei into spatial grid.
    
    Args:
        nuclei_in_str2: Nuclei in stripe 2
        n_ap_bins: Number of bins along AP axis (default 5)
        n_dv_bins: Number of bins along DV axis (default 5)
        dv_crop_min: Optional DV axis crop minimum (normalized, default -1.0 means no crop)
        dv_crop_max: Optional DV axis crop maximum (normalized, default -1.0 means no crop)
    
    Returns:
        DataFrame with apBin and yBin assignments
    """
    print(f"Binning nuclei into {n_ap_bins}x{n_dv_bins} spatial grid...")
    
    # Normalize AP and DV positions
    minap = nuclei_in_str2['ap_registered'].min()
    maxap = nuclei_in_str2['ap_registered'].max()
    nuclei_in_str2 = nuclei_in_str2.copy()
    nuclei_in_str2['ap_reg_norm'] = (
        (nuclei_in_str2['ap_registered'] - minap) / (maxap - minap)
    )
    
    miny = nuclei_in_str2['yPos'].min()
    maxy = nuclei_in_str2['yPos'].max()
    nuclei_in_str2['yPos_norm'] = (
        (nuclei_in_str2['yPos'] - miny) / (maxy - miny)
    )
    
    # Optional DV cropping (Jenny's note: may not be needed, eve pattern is DV-invariant)
    if dv_crop_min != -1.0 and dv_crop_max != -1.0:
        print(f"Cropping DV axis: {dv_crop_min} - {dv_crop_max}")
        nuclei_in_str2 = nuclei_in_str2[
            (nuclei_in_str2['yPos_norm'] > dv_crop_min) & 
            (nuclei_in_str2['yPos_norm'] < dv_crop_max)
        ]
    
    # Create bins: overlapping on AP (50%), non-overlapping on DV
    ap_overlap = 0.5
    dv_overlap = 0.0
    print(f"Using overlapping bins: {ap_overlap*100:.0f}% AP, {dv_overlap*100:.0f}% DV")

    dv_min = nuclei_in_str2['yPos_norm'].min()
    dv_max = nuclei_in_str2['yPos_norm'].max()
    dv_bins = create_overlapping_bin_ranges(n_dv_bins, dv_min, dv_max, dv_overlap)

    ap_bins = create_overlapping_bin_ranges(n_ap_bins, 0, 1, ap_overlap)
    
    # Assign nuclei to all overlapping bins they fall within
    result_rows = []
    for _, nucleus in nuclei_in_str2.iterrows():
        ap_norm = nucleus['ap_reg_norm']
        dv_norm = nucleus['yPos_norm']
        nucleus_id = nucleus['nucleus_id']
        embryo = nucleus['embryo']
        
        # Find all bins this nucleus belongs to
        for ap_idx, (ap_min, ap_max) in enumerate(ap_bins):
            if ap_min <= ap_norm <= ap_max:
                for dv_idx, (dv_min_bin, dv_max_bin) in enumerate(dv_bins):
                    if dv_min_bin <= dv_norm <= dv_max_bin:
                        result_rows.append({
                            'nucleus_id': nucleus_id,
                            'embryo': embryo,
                            'apBin': ap_idx + 1,  # 1-indexed
                            'yBin': dv_idx + 1     # 1-indexed
                        })
    
    result_df = pd.DataFrame(result_rows)
    total_assignments = len(result_df)
    unique_nuclei = len(nuclei_in_str2)
    print(f"Created {total_assignments} bin assignments for {unique_nuclei} nuclei "
          f"(avg {total_assignments/unique_nuclei:.1f} bins per nucleus)")
    
    return result_df


def compute_local_averages(binned_nuclear_ids, fluo_traces):
    """
    Compute mean fluorescence traces for each spatial bin.
    
    Returns:
        DataFrame with apBin, yBin, and mean fluorescence traces
    """
    print("Computing locally averaged fluorescence traces...")
    
    # Merge fluorescence traces with bin assignments
    fluo_traces_with_bins = binned_nuclear_ids.merge(fluo_traces, on='nucleus_id')
    
    # Compute mean per spatial bin
    locally_averaged = fluo_traces_with_bins.drop(columns=['nucleus_id', 'embryo'])
    locally_averaged = locally_averaged.groupby(['apBin', 'yBin']).mean() #TODO median or mean?####################################################################################################
    locally_averaged = locally_averaged.reset_index()
    
    print(f"Computed {len(locally_averaged)} locally averaged traces")
    return locally_averaged


def save_outputs(locally_averaged, output_path, output_path_no_ids):
    """Save locally averaged traces with and without bin IDs."""
    print(f"Saving locally averaged traces to: {output_path}")
    locally_averaged.to_csv(output_path, index=False)
    
    # Also save version without apBin/yBin for Julia script
    print(f"Saving traces without IDs to: {output_path_no_ids}")
    locally_averaged.drop(columns=['apBin', 'yBin']).to_csv(
        output_path_no_ids, index=False, header=False
    )


def plot_heatmap(locally_averaged, output_plot_path):
    """Generate and save summed fluorescence heatmap."""
    print(f"Generating heatmap: {output_plot_path}")
    
    # Compute sum fluorescence
    time_cols = [col for col in locally_averaged.columns 
                 if col not in ['apBin', 'yBin']]
    sum_fluo = locally_averaged[time_cols].sum(axis=1)
    locally_averaged_with_sum = locally_averaged.copy()
    locally_averaged_with_sum['sum_fluo'] = sum_fluo
    
    # Pivot and plot
    heatmap_data = locally_averaged_with_sum.pivot(
        index='yBin', columns='apBin', values='sum_fluo'
    )
    
    plt.figure(figsize=(7, 6))
    sns.heatmap(heatmap_data, annot=True, fmt='.0f', cmap=cmc.davos, 
                cbar_kws={'label': 'Sum Fluorescence'})#, yticklabels=heatmap_data.index[::-1])
    ax = plt.gca()
    ax.invert_yaxis()
    plt.title('Locally Averaged Transcription Activity')
    plt.xlabel('AP Bin')
    plt.ylabel('DV Bin')
    plt.tight_layout()
    plt.savefig(output_plot_path, dpi=300)
    plt.close()
    heatmap_data.to_csv(output_plot_path.replace('.pdf', '_data.csv'))


def main():
    parser = argparse.ArgumentParser( description='Preprocess eve transcription data for spatial mRNA decay analysis')
    parser.add_argument( '--input', required=True, help='Path to filtered eve data CSV')
    parser.add_argument( '--output', required=True, help='Output path for locally averaged traces')
    parser.add_argument( '--output-no-ids', required=True, help='Output path for traces without bin IDs (for Julia)')
    parser.add_argument( '--plot', required=True, help='Output path for heatmap plot')
    parser.add_argument( '--config', type=str, help='Path to config YAML (if provided, reads stripe AP ranges from here)')
    parser.add_argument( '--stripe', type=str, help='Stripe identifier (stripe2, stripe3, etc.)')
    parser.add_argument( '--ap-min', type=float, help='Minimum AP coordinate for stripe (required if --config not provided)')
    parser.add_argument( '--ap-max', type=float, help='Maximum AP coordinate for stripe (required if --config not provided)')
    parser.add_argument( '--dv-min', type=float, default=-1.0, help='Minimum DV coordinate for stripe (-1.0 means no crop, default: -1.0)')
    parser.add_argument( '--dv-max', type=float, default=-1.0, help='Maximum DV coordinate for stripe (-1.0 means no crop, default: -1.0)')
    parser.add_argument( '--n-ap-bins', type=int, default=5, help='Number of AP bins (default: 5)')
    parser.add_argument( '--n-dv-bins', type=int, default=5, help='Number of DV bins (default: 5)')
    parser.add_argument( '--max-time', type=int, default=1200, help='Maximum time in seconds (default: 1200 = 20 min)')
    
    args = parser.parse_args()
    
    # Determine AP range: either from config file or from command-line arguments
    if args.config and args.stripe:
        # Read stripe ranges from config file using ruamel.yaml for round-trip safety
        yaml = YAML()
        with open(args.config, 'r') as f:
            config = yaml.load(f)

        if 'stripe_ranges' not in config or args.stripe not in config['stripe_ranges']:
            raise ValueError(f"Stripe {args.stripe} not found in config file {args.config}")

        ap_min = config['stripe_ranges'][args.stripe]['min']
        ap_max = config['stripe_ranges'][args.stripe]['max']
        print(f"Read AP range from config: [{ap_min:.4f}, {ap_max:.4f}]")
    elif args.ap_min is not None and args.ap_max is not None:
        # Use command-line arguments
        ap_min = args.ap_min
        ap_max = args.ap_max
    else:
        raise ValueError("Must provide either (--config and --stripe) or (--ap-min and --ap-max)")
    
    # Print full command to reproduce the run
    print("\n=== Command to reproduce this script run ===")
    print(" ".join(sys.argv))
    
    # Log stripe information
    stripe_name = args.stripe if args.stripe else f"AP[{ap_min:.2f}, {ap_max:.2f}]"
    print(f"\n=== Processing {stripe_name} ===")
    print(f"AP range: [{ap_min:.4f}, {ap_max:.4f}]")
    
    # Run pipeline
    data_filtered = load_data(args.input)
    fluo_traces = extract_fluorescence_traces(data_filtered, args.max_time)
    pos_data = extract_nuclear_positions(data_filtered)
    pos_data_str2 = filter_stripe2_nuclei(pos_data, ap_min, ap_max, args.stripe)
    
    # Merge position and fluorescence data
    pos_sum_df = fluo_traces.merge(pos_data_str2, on='nucleus_id')
    nuclei_in_str2 = pos_sum_df[['nucleus_id', 'embryo', 'ap_registered', 'yPos']]
    
    # Bin nuclei spatially
    binned_nuclear_ids = bin_nuclei_spatially(
        nuclei_in_str2, 
        args.n_ap_bins, 
        args.n_dv_bins,
        dv_crop_min=args.dv_min,
        dv_crop_max=args.dv_max
    )
    
    # Compute local averages
    locally_averaged = compute_local_averages(binned_nuclear_ids, fluo_traces)
    
    # Save outputs
    save_outputs(locally_averaged, args.output, args.output_no_ids)
    plot_heatmap(locally_averaged, args.plot)
    
    print("\n=== Preprocessing complete ===")
    print(f"Stripe: {stripe_name}")
    print(f"Generated {len(locally_averaged)} locally averaged traces")
    print(f"Total bins: {args.n_ap_bins} AP × {args.n_dv_bins} DV = {args.n_ap_bins * args.n_dv_bins}")


if __name__ == '__main__':
    main()
