#!/usr/bin/env python3
"""
Validate nuclei density in processed mRNA data against Berrocal_2020 benchmarks.

This script computes k=4 nearest neighbor distances for nuclei in normalized
coordinates and validates against stripe-specific thresholds computed from
Berrocal_2020 data.

Validation thresholds are loaded from config.yaml (computed by 06_compute_validation_thresholds.py).

Usage:
    python 05_validate_nuclei_density.py <position_data_file> <stripe> <config_yaml> <output_report>
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.neighbors import NearestNeighbors
from ruamel.yaml import YAML


def load_position_data(filepath: str) -> pd.DataFrame:
    """Load SASS position_data-intense.txt file."""
    df = pd.read_csv(filepath, sep='\t')
    return df


def compute_nuclei_density_metrics(df: pd.DataFrame, k: int = 4) -> dict:
    """
    Compute nuclei density metrics using k-nearest neighbors.
    
    Args:
        df: DataFrame with nuclear positions (nucx, nucy columns)
        k: Number of nearest neighbors to compute
    
    Returns:
        Dictionary with validation metrics
    """
    # Aggregate to one row per nucleus
    nuc_data = df.groupby('nuc').agg({
        'nucx': 'first',
        'nucy': 'first'
    }).reset_index()
    
    n_nuclei = len(nuc_data)
    
    if n_nuclei == 0:
        raise ValueError("No nuclei found in position data")
    
    # Normalize coordinates to [0, 1]
    x_min = nuc_data['nucx'].min()
    x_max = nuc_data['nucx'].max()
    y_min = nuc_data['nucy'].min()
    y_max = nuc_data['nucy'].max()
    
    nuc_data['nucx_norm'] = (nuc_data['nucx'] - x_min) / (x_max - x_min)
    nuc_data['nucy_norm'] = (nuc_data['nucy'] - y_min) / (y_max - y_min)
    
    # Compute k-nearest neighbors in normalized coordinates
    coords = nuc_data[['nucx_norm', 'nucy_norm']].values
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


def load_validation_thresholds(config_path: str, stripe: str) -> dict:
    """
    Load validation thresholds for a specific stripe from config.yaml.
    
    Args:
        config_path: Path to config.yaml
        stripe: Stripe name (e.g., 'stripe2', 'stripe3')
    
    Returns:
        Dictionary with validation thresholds
    """
    yaml = YAML()
    with open(config_path, 'r') as f:
        config = yaml.load(f)
    
    validation_thresholds = config.get('validation_thresholds', {})
    if stripe not in validation_thresholds:
        raise ValueError(
            f"No validation thresholds found for {stripe} in config.yaml. "
            f"Run 06_compute_validation_thresholds.py first."
        )
    
    return validation_thresholds[stripe]


def validate_metrics(metrics: dict, thresholds: dict) -> tuple[bool, list[str]]:
    """
    Validate metrics against stripe-specific Berrocal_2020 benchmarks.
    
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


def main():
    if len(sys.argv) != 5:
        print("Usage: python 05_validate_nuclei_density.py <position_data_file> <stripe> <config_yaml> <output_report>")
        sys.exit(1)
    
    position_data_file = sys.argv[1]
    stripe = sys.argv[2]
    config_yaml = sys.argv[3]
    output_report = sys.argv[4]
    
    print(f"Validating nuclei density for {position_data_file}")
    print(f"  Stripe: {stripe}")
    
    # Load validation thresholds
    try:
        thresholds = load_validation_thresholds(config_yaml, stripe)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # Load and compute metrics
    df = load_position_data(position_data_file)
    print(f"  Loaded {len(df)} spot records")
    
    k = thresholds.get('k', 4)
    metrics = compute_nuclei_density_metrics(df, k=k)
    
    # Validate
    passed, warnings = validate_metrics(metrics, thresholds)
    
    # Print results
    print("\nNuclei Density Metrics:")
    print(f"  Number of nuclei: {metrics['n_nuclei']}")
    print(f"  Average k={k} NN distance: {metrics['avg_nn_distance']:.4f}")
    print(f"  Median k={k} NN distance: {metrics['median_nn_distance']:.4f}")
    print(f"  Std k={k} NN distance: {metrics['std_nn_distance']:.4f}")
    print(f"  Range: {metrics['min_nn_distance']:.4f} - {metrics['max_nn_distance']:.4f}")
    
    print(f"\nValidation against Berrocal_2020 benchmarks ({stripe}):")
    print(f"  Expected nuclei count: ~{thresholds['expected_nuclei_count_mean']:.1f}")
    print(f"  Expected median NN distance: ~{thresholds['expected_nn_distance_median']:.4f}")
    print(f"  Valid range: {thresholds['min_nn_distance']:.4f} - {thresholds['max_nn_distance']:.4f}")
    
    if warnings:
        print("\nValidation messages:")
        for warning in warnings:
            print(f"  {warning}")
    
    if passed:
        print("\n✓ VALIDATION PASSED")
        status = "PASS"
    else:
        print("\n✗ VALIDATION FAILED - Nuclei density outside expected range")
        status = "FAIL"
    
    # Write report
    output_path = Path(output_report)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_report, 'w') as f:
        f.write(f"Nuclei Density Validation Report\n")
        f.write(f"=" * 80 + "\n\n")
        f.write(f"Input file: {position_data_file}\n")
        f.write(f"Stripe: {stripe}\n")
        f.write(f"Status: {status}\n\n")
        
        f.write(f"Metrics:\n")
        f.write(f"  Number of nuclei: {metrics['n_nuclei']}\n")
        f.write(f"  Average k={k} NN distance: {metrics['avg_nn_distance']:.4f}\n")
        f.write(f"  Median k={k} NN distance: {metrics['median_nn_distance']:.4f}\n")
        f.write(f"  Std k={k} NN distance: {metrics['std_nn_distance']:.4f}\n")
        f.write(f"  Range: {metrics['min_nn_distance']:.4f} - {metrics['max_nn_distance']:.4f}\n\n")
        
        f.write(f"Validation Criteria ({stripe}):\n")
        f.write(f"  Expected nuclei count: ~{thresholds['expected_nuclei_count_mean']:.1f}\n")
        f.write(f"  Expected median NN distance: ~{thresholds['expected_nn_distance_median']:.4f}\n")
        f.write(f"  Valid range: {thresholds['min_nn_distance']:.4f} - {thresholds['max_nn_distance']:.4f}\n")
        f.write(f"  Based on {thresholds['n_embryos_used']} Berrocal_2020 embryos\n\n")
        
        if warnings:
            f.write(f"Messages:\n")
            for warning in warnings:
                f.write(f"  {warning}\n")
    
    print(f"\nReport written to {output_report}")
    
    # Exit with error code if validation failed
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
