#!/usr/bin/env python3
"""
Identify stripe ranges from eve expression data.

This script analyzes the spatial distribution of eve expression along the
anterior-posterior (AP) axis to identify the 7 characteristic eve stripes.
It automatically detects peaks in fluorescence intensity and determines the
AP coordinate ranges for each stripe.

The detected stripe ranges are written to config.yaml for use by the analysis pipeline.

Usage:
    python scripts/identify_stripe_ranges.py --input data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv \
                                             --config config.yaml \
                                             --output results/figures/intermediate/transcription/stripe_identification.png
"""

import argparse
import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.signal as signal
import matplotlib.pyplot as plt
from ruamel.yaml import YAML

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_and_prepare_data(input_path, max_time_seconds=1200, time_step=20):
    """
    Load eve expression data and compute summed fluorescence per nucleus.
    
    Parameters
    ----------
    input_path : str or Path
        Path to input CSV file with eve expression data
    max_time_seconds : int
        Maximum time to include in fluorescence sum (default: 1200s = 20min)
    time_step : int
        Time resolution in seconds (default: 20s)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with nucleus_id, ap_registered, and sum_fluo_early columns
    """
    logger.info(f"Loading data from {input_path}")
    data_filtered = pd.read_csv(input_path, header=0)
    
    # Pivot to get fluorescence traces per nucleus
    fluo_traces = data_filtered.pivot(index='nucleus_id', columns='time', values='fluo')
    fluo_traces = fluo_traces.fillna(0)
    fluo_traces.reset_index(inplace=True)
    
    # Sum fluorescence before max_time (early expression period)
    max_time_col = int(max_time_seconds / time_step) + 2  # +2 for nucleus_id column + 1-indexing
    fluo_traces['sum_fluo_early'] = fluo_traces.iloc[:, 1:max_time_col].sum(axis=1)
    
    # Get nuclear positions (median over time)
    pos_data = data_filtered[['nucleus_id', 'ap_registered']].groupby('nucleus_id').median()
    pos_data.reset_index(inplace=True)
    
    # Merge position and fluorescence data
    result_df = fluo_traces[['nucleus_id', 'sum_fluo_early']].merge(pos_data, on='nucleus_id')
    
    logger.info(f"Loaded data for {len(result_df)} nuclei")
    return result_df


def bin_data_by_ap(data_df, bin_size=0.01):
    """
    Bin nuclei by AP position and compute median fluorescence per bin.
    
    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame with ap_registered and sum_fluo_early columns
    bin_size : float
        Size of AP bins (default: 0.01)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with ap_bin_center and median_fluo columns
    """
    logger.info(f"Binning data with bin size {bin_size}")
    
    min_ap = data_df['ap_registered'].min()
    max_ap = data_df['ap_registered'].max()
    bins = np.arange(min_ap, max_ap + bin_size, bin_size)
    
    data_df['ap_bin'] = pd.cut(data_df['ap_registered'], bins)
    binned_data = data_df.groupby('ap_bin')['sum_fluo_early'].median().reset_index()
    binned_data['ap_bin_center'] = binned_data['ap_bin'].apply(lambda x: x.mid)
    
    return binned_data[['ap_bin_center', 'sum_fluo_early']].rename(
        columns={'sum_fluo_early': 'median_fluo'}
    )


def smooth_fluorescence_data(fluo_data, method='savgol', **kwargs):
    """
    Apply smoothing to fluorescence data for peak detection.
    
    Parameters
    ----------
    fluo_data : np.ndarray
        1D array of fluorescence values
    method : str
        Smoothing method: 'savgol', 'moving_average', or 'none'
    **kwargs
        Additional parameters for smoothing method
    
    Returns
    -------
    np.ndarray
        Smoothed fluorescence data
    str
        Description of smoothing method used
    """
    if method == 'none':
        return fluo_data, "No smoothing"
    
    if method == 'savgol':
        window = kwargs.get('window_length', 11)
        polyorder = kwargs.get('polyorder', 3)
        try:
            smoothed = signal.savgol_filter(fluo_data, window_length=window, polyorder=polyorder)
            return smoothed, f"Savitzky-Golay (window={window}, poly={polyorder})"
        except np.linalg.LinAlgError:
            logger.warning(f"Savitzky-Golay failed with window={window}, trying moving average")
            method = 'moving_average'
    
    if method == 'moving_average':
        window = kwargs.get('window', 5)
        smoothed = np.convolve(fluo_data, np.ones(window)/window, mode='same')
        return smoothed, f"Moving average (window={window})"
    
    logger.warning(f"Unknown smoothing method '{method}', using original data")
    return fluo_data, "Original data (no smoothing)"


def detect_stripe_peaks_and_ranges(binned_data, prominence=300):
    """
    Detect stripe peaks and compute AP coordinate ranges.
    
    Parameters
    ----------
    binned_data : pd.DataFrame
        DataFrame with ap_bin_center and median_fluo columns
    prominence : float
        Minimum prominence for peak detection (default: 300)
    
    Returns
    -------
    list of dict
        List of stripe information dicts with keys: stripe_num, min, max, center
    np.ndarray
        Smoothed fluorescence data used for detection
    np.ndarray
        Indices of detected peaks
    np.ndarray
        Indices of detected troughs
    """
    logger.info(f"Detecting peaks with prominence threshold {prominence}")
    
    # Smooth the data
    smoothed_fluo, smooth_method = smooth_fluorescence_data(
        binned_data['median_fluo'].values,
        method='savgol',
        window_length=11,
        polyorder=3
    )
    logger.info(f"Smoothing method: {smooth_method}")
    
    # Find peaks (stripe centers) and troughs (inter-stripes)
    peaks, _ = signal.find_peaks(smoothed_fluo, prominence=prominence)
    troughs, _ = signal.find_peaks(-smoothed_fluo, prominence=prominence)
    
    logger.info(f"Detected {len(peaks)} peaks and {len(troughs)} troughs")
    
    # Extract AP coordinates
    ap_centers = binned_data['ap_bin_center'].values
    
    # For each peak, find boundaries based on nearest troughs
    stripe_ranges = []
    for i, peak_idx in enumerate(peaks):
        peak_ap = ap_centers[peak_idx]
        
        # Find nearest left trough
        left_troughs = troughs[troughs < peak_idx]
        if len(left_troughs) > 0:
            left_trough_idx = left_troughs[-1]
            left_ap = ap_centers[left_trough_idx]
        else:
            left_ap = ap_centers[0]  # Start of data
        
        # Find nearest right trough
        right_troughs = troughs[troughs > peak_idx]
        if len(right_troughs) > 0:
            right_trough_idx = right_troughs[0]
            right_ap = ap_centers[right_trough_idx]
        else:
            right_ap = ap_centers[-1]  # End of data
        
        # Create symmetric range around peak
        half_width = min(peak_ap - left_ap, right_ap - peak_ap)
        left_cutoff = peak_ap - half_width
        right_cutoff = peak_ap + half_width
        
        stripe_ranges.append({
            'stripe_num': i + 1,
            'min': float(left_cutoff),
            'max': float(right_cutoff),
            'center': float(peak_ap)
        })
    
    return stripe_ranges, smoothed_fluo, peaks, troughs


def plot_stripe_identification(binned_data, smoothed_fluo, peaks, troughs, 
                                stripe_ranges, output_path):
    """
    Create diagnostic plot showing detected stripes.
    
    Parameters
    ----------
    binned_data : pd.DataFrame
        DataFrame with ap_bin_center and median_fluo columns
    smoothed_fluo : np.ndarray
        Smoothed fluorescence data
    peaks : np.ndarray
        Indices of detected peaks
    troughs : np.ndarray
        Indices of detected troughs
    stripe_ranges : list of dict
        Detected stripe ranges
    output_path : str or Path
        Path to save the diagnostic plot
    """
    logger.info(f"Creating diagnostic plot at {output_path}")
    
    ap_centers = binned_data['ap_bin_center'].values
    raw_fluo = binned_data['median_fluo'].values
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot raw and smoothed data
    ax.plot(ap_centers, raw_fluo, 'o', alpha=0.3, markersize=3, label='Raw data')
    ax.plot(ap_centers, smoothed_fluo, 'b-', linewidth=2, label='Smoothed')
    
    # Mark peaks (stripe centers)
    ax.plot(ap_centers[peaks], smoothed_fluo[peaks], 'r^', markersize=10, 
            label='Stripe centers', zorder=5)
    
    # Mark troughs (inter-stripes)
    ax.plot(ap_centers[troughs], smoothed_fluo[troughs], 'gv', markersize=8, 
            label='Inter-stripes', zorder=5)
    
    # Shade stripe ranges
    for stripe in stripe_ranges:
        ax.axvspan(stripe['min'], stripe['max'], alpha=0.2, color='red')
        # Label stripe number
        ax.text(stripe['center'], ax.get_ylim()[1] * 0.95, 
                f"S{stripe['stripe_num']}", 
                ha='center', va='top', fontsize=10, fontweight='bold')
    
    ax.set_xlabel('AP Position (registered)', fontsize=12)
    ax.set_ylabel('Median Fluorescence (early expression)', fontsize=12)
    ax.set_title('Eve Stripe Identification', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    logger.info(f"Diagnostic plot saved to {output_path}")


def write_stripe_ranges_to_config(stripe_ranges, config_path, prominence, bin_size):
    """
    Update config.yaml with detected stripe ranges.
    
    Parameters
    ----------
    stripe_ranges : list of dict
        Detected stripe ranges
    config_path : str or Path
        Path to config.yaml file
    prominence : float
        Prominence threshold used for detection
    bin_size : float
        Bin size used for AP binning
    """
    logger.info(f"Updating configuration file: {config_path}")
    
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.default_flow_style = False
    yaml.width = 4096  # Prevent line wrapping
    
    # Load existing config
    config_path = Path(config_path)
    if config_path.exists():
        with open(config_path, 'r') as f:
            config = yaml.load(f)
    else:
        config = {}
    
    # Create stripe_ranges section
    stripe_ranges_dict = {}
    for stripe in stripe_ranges:
        stripe_name = f"stripe{stripe['stripe_num']}"
        stripe_ranges_dict[stripe_name] = {
            'center': round(stripe['center'], 4),
            'min': round(stripe['min'], 4),
            'max': round(stripe['max'], 4)
        }
    
    config['stripe_ranges'] = stripe_ranges_dict
    
    # Add metadata
    config['stripe_ranges_metadata'] = {
        'last_updated': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'detection_method': 'peak detection with Savitzky-Golay smoothing',
        'prominence_threshold': prominence,
        'bin_size': bin_size
    }
    
    # Write updated config
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    
    logger.info(f"Wrote {len(stripe_ranges)} stripe ranges to {config_path}")
    
    # Print summary
    logger.info("Detected stripe ranges:")
    for stripe in stripe_ranges:
        logger.info(f"  Stripe {stripe['stripe_num']}: "
                   f"AP [{stripe['min']:.4f}, {stripe['max']:.4f}], "
                   f"center = {stripe['center']:.4f}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Identify eve stripe ranges from expression data'
    )
    parser.add_argument(
        '--input',
        type=str,
        default='data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv',
        help='Input CSV file with eve expression data'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config.yaml',
        help='Config YAML file to update with stripe ranges'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='results/figures/intermediate/transcription/stripe_identification.png',
        help='Output path for diagnostic plot'
    )
    parser.add_argument(
        '--bin-size',
        type=float,
        default=0.01,
        help='Bin size for AP axis (default: 0.01)'
    )
    parser.add_argument(
        '--prominence',
        type=float,
        default=300,
        help='Minimum prominence for peak detection (default: 300)'
    )
    parser.add_argument(
        '--max-time',
        type=int,
        default=1200,
        help='Maximum time in seconds for fluorescence sum (default: 1200)'
    )
    
    args = parser.parse_args()
    
    # Load and prepare data
    data_df = load_and_prepare_data(args.input, max_time_seconds=args.max_time)
    
    # Bin data by AP position
    binned_data = bin_data_by_ap(data_df, bin_size=args.bin_size)
    
    # Detect stripe peaks and ranges
    stripe_ranges, smoothed_fluo, peaks, troughs = detect_stripe_peaks_and_ranges(
        binned_data, prominence=args.prominence
    )
    
    # Create diagnostic plot
    plot_stripe_identification(
        binned_data, smoothed_fluo, peaks, troughs, stripe_ranges, args.output
    )
    
    # Write to config file
    write_stripe_ranges_to_config(stripe_ranges, args.config, args.prominence, args.bin_size)
    
    logger.info("Stripe identification complete!")


if __name__ == '__main__':
    main()