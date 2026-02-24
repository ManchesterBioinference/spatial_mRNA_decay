#!/usr/bin/env python3
"""
Create a spot-level expression density plot from SASS position_data-intense.txt.

This script only generates the expression density figure and does not run any
other aggregation or validation logic from 04_aggregate_and_validate_mrna.py.

Outputs both:
- Original (raw) spot counts per X-bin
- Smoothed spot counts per X-bin

Usage (simple):
    python scripts/09_center_mrna_data.py <position_data_file> <output_png>

Usage (legacy-compatible with script 04 argument shape):
    python scripts/09_center_mrna_data.py \
        <position_data_file> <stripe> <embryo_id> <config_yaml> \
        <output_file> <validation_report> <figure_file_or_dir> <no_groupByNuclei>

In legacy mode, output is written as:
    <figure_parent>/<embryo_id>_sass_formodel_expression_density.png
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def compute_spot_density(df: pd.DataFrame, bin_size: float = 0.75) -> pd.DataFrame:
    """Bin spot positions along X and count spots per bin."""
    if 'spotx' not in df.columns:
        raise ValueError("Input data must contain 'spotx' column.")

    x = df['spotx'].dropna()
    if x.empty:
        raise ValueError("No valid values found in 'spotx' column.")

    min_x = x.min()
    max_x = x.max()

    if np.isclose(min_x, max_x):
        bins = np.array([min_x - bin_size / 2, min_x + bin_size / 2])
    else:
        bins = np.arange(min_x, max_x + bin_size, bin_size)
        if len(bins) < 2:
            bins = np.array([min_x, max_x + bin_size])

    binned = pd.cut(x, bins=bins, include_lowest=True)
    density = binned.value_counts(sort=False).reset_index()
    density.columns = ['x_bin', 'spot_count']
    density['x_bin_center'] = density['x_bin'].apply(lambda interval: interval.mid)

    return density[['x_bin_center', 'spot_count']]


def smooth_series(values: np.ndarray, window: int = 5) -> np.ndarray:
    """Smooth 1D values with a centered moving average."""
    n = len(values)
    if n == 0:
        return values

    window = max(1, int(window))
    if window % 2 == 0:
        window += 1
    if window > n:
        window = n if n % 2 == 1 else max(1, n - 1)

    if window <= 1:
        return values.copy()

    return pd.Series(values).rolling(window=window, center=True, min_periods=1).mean().to_numpy()


def plot_expression_density(df: pd.DataFrame, output_path: Path, bin_size: float = 0.75, smooth_window: int = 5) -> None:
    """Plot original and smoothed expression density curves."""
    density = compute_spot_density(df, bin_size=bin_size)
    density['spot_count_smooth'] = smooth_series(density['spot_count'].to_numpy(), window=smooth_window)

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(
        density['x_bin_center'].to_numpy(),
        density['spot_count'].to_numpy(),
        'o-',
        linewidth=1.5,
        markersize=4,
        color='0.6',
        alpha=0.85,
        label='Original spot count per bin'
    )

    ax.plot(
        density['x_bin_center'].to_numpy(),
        density['spot_count_smooth'].to_numpy(),
        '-',
        linewidth=2.5,
        color='darkorange',
        label=f'Smoothed (moving average, window={smooth_window})'
    )

    ax.set_xlabel('X Position', fontsize=12)
    ax.set_ylabel('Number of Spots', fontsize=12)
    ax.set_title('Expression Density - Spot Level', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()



def parse_args(argv: list[str]) -> tuple[Path, Path]:
    """
    Parse CLI arguments.

    Returns:
        (position_data_file, output_plot_file)
    """
    if len(argv) == 3:
        # Simple mode: input + output plot path
        return Path(argv[1]), Path(argv[2])

    if len(argv) == 9:
        # Legacy-compatible mode using script-04 style args
        position_data_file = Path(argv[1])
        embryo_id = argv[3]
        fig_file_or_dir = Path(argv[7])
        fig_dir = fig_file_or_dir.parent
        output_plot = fig_dir / f"{embryo_id}_sass_formodel_expression_density.png"
        return position_data_file, output_plot

    raise ValueError(
        "Usage:\n"
        "  python scripts/09_center_mrna_data.py <position_data_file> <output_png>\n"
        "or legacy-compatible:\n"
        "  python scripts/09_center_mrna_data.py "
        "<position_data_file> <stripe> <embryo_id> <config_yaml> "
        "<output_file> <validation_report> <figure_file_or_dir> <no_groupByNuclei>"
    )



def main() -> None:
    try:
        position_data_file, output_plot = parse_args(sys.argv)
    except ValueError as e:
        print(e)
        sys.exit(1)

    print("\n=== Command to reproduce this script run ===")
    print(" ".join(sys.argv))

    if not position_data_file.exists():
        print(f"ERROR: input file not found: {position_data_file}")
        sys.exit(1)

    df = pd.read_csv(position_data_file, sep='\t')
    print(f"Loaded {len(df)} spot records from {position_data_file}")

    plot_expression_density(df, output_plot, smooth_window=35)
    print(f"Wrote expression density plot to {output_plot}")


if __name__ == '__main__':
    main()
