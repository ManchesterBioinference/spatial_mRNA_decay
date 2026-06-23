#!/usr/bin/env python3
"""
Apply a lateral shift to stripe AP coordinate ranges in a config YAML.

For each stripe in stripe_ranges, the shift amount is computed as:

    shift = (max - min) / (2 * n_ap_bins)

This moves the entire window (min, center, max) by ±half a bin width,
representing a sensitivity analysis to check whether data alignment between
Berrocal (transcription) and smFISH (mRNA) imaging data affects inferred
degradation rates.

The full window is translated rigidly — min, max, and center all shift by
the same amount, preserving window width and shape.

Usage
-----
    python scripts/apply_stripe_shift.py \\
        --input-config  results_1200/config.yaml \\
        --output-config shiftRight/results_1200/config.yaml \\
        --direction     right \\
        --n-ap-bins     5

Arguments
---------
--input-config   Path to source config.yaml (must have stripe_ranges populated
                 by identify_stripe_ranges before running this script).
--output-config  Destination path for the shifted config.yaml.
--direction      'left'  → shift toward lower AP coordinate (anterior)
                 'right' → shift toward higher AP coordinate (posterior)
--n-ap-bins      Number of AP bins used in the analysis (default: 5).
                 Determines bin width: bin_width = (max - min) / n_ap_bins.
                 Shift = bin_width / 2.

Outputs
-------
A YAML file identical to the input except:
- stripe_ranges: each stripe's min, max, center shifted by ±(bin_width/2)
- stripe_ranges_metadata: extended with shift_direction, shift_amounts,
  shift_applied_at fields for full reproducibility tracing.

Example (stripe2 with n_ap_bins=5)
------------------------------------
  Original : min=0.3525  center=0.3975  max=0.4425  → width=0.09, shift=0.009
  shiftRight: min=0.3615  center=0.4065  max=0.4515
  shiftLeft : min=0.3435  center=0.3885  max=0.4335
"""

import argparse
import copy
import os
import sys
from datetime import datetime

import yaml


def apply_shift(config: dict, direction: str, n_ap_bins: int) -> dict:
    """Return a deep copy of config with stripe AP ranges shifted by half a bin.

    Parameters
    ----------
    config : dict
        Parsed YAML config dict containing a populated ``stripe_ranges`` section.
    direction : str
        ``'left'`` moves the window toward lower AP values (anterior direction).
        ``'right'`` moves the window toward higher AP values (posterior direction).
    n_ap_bins : int
        Number of AP bins; used to compute bin_width = (max - min) / n_ap_bins.

    Returns
    -------
    dict
        Modified config where each stripe's min / max / center are shifted by
        ±(bin_width / 2).  A ``shift_metadata`` block is added inside
        ``stripe_ranges_metadata`` recording direction, per-stripe amounts,
        and timestamp.

    Raises
    ------
    ValueError
        If ``stripe_ranges`` is missing or empty in the input config.
    """
    cfg = copy.deepcopy(config)
    stripe_ranges = cfg.get("stripe_ranges", {})
    if not stripe_ranges:
        raise ValueError(
            "Config does not contain 'stripe_ranges'.  "
            "Run identify_stripe_ranges from the main Snakefile first, "
            "then pass the populated results_{max_time}/config.yaml here."
        )

    sign = 1.0 if direction == "right" else -1.0
    shift_amounts: dict[str, float] = {}

    for stripe, bounds in stripe_ranges.items():
        # Skip any non-dict entries that may appear under stripe_ranges
        if not isinstance(bounds, dict):
            continue

        ap_min = float(bounds["min"])
        ap_max = float(bounds["max"])
        bin_width = (ap_max - ap_min) / n_ap_bins
        shift = sign * bin_width / 2.0

        cfg["stripe_ranges"][stripe]["min"] = round(ap_min + shift, 6)
        cfg["stripe_ranges"][stripe]["max"] = round(ap_max + shift, 6)
        cfg["stripe_ranges"][stripe]["center"] = round(float(bounds["center"]) + shift, 6)
        shift_amounts[stripe] = round(shift, 6)

    # Extend stripe_ranges_metadata with shift provenance
    metadata = cfg.get("stripe_ranges_metadata", {})
    metadata["shift_direction"] = direction
    metadata["shift_amounts"] = shift_amounts
    metadata["shift_applied_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    cfg["stripe_ranges_metadata"] = metadata

    return cfg


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Shift stripe AP coordinate ranges by half a bin width for sensitivity analysis."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input-config",
        required=True,
        metavar="PATH",
        help="Source config.yaml with stripe_ranges populated (e.g. results_1200/config.yaml).",
    )
    parser.add_argument(
        "--output-config",
        required=True,
        metavar="PATH",
        help="Destination path for the shifted config.yaml.",
    )
    parser.add_argument(
        "--direction",
        required=True,
        choices=["left", "right"],
        help=(
            "'left'  → shift toward lower AP / anterior direction. "
            "'right' → shift toward higher AP / posterior direction."
        ),
    )
    parser.add_argument(
        "--n-ap-bins",
        type=int,
        default=5,
        metavar="N",
        help="Number of AP bins used in the analysis (default: 5).",
    )
    args = parser.parse_args()

    # ── Load source config ──────────────────────────────────────────────────
    if not os.path.isfile(args.input_config):
        print(f"ERROR: Input config not found: {args.input_config}", file=sys.stderr)
        sys.exit(1)

    with open(args.input_config) as fh:
        config = yaml.safe_load(fh)

    # ── Apply shift ─────────────────────────────────────────────────────────
    shifted = apply_shift(config, args.direction, args.n_ap_bins)

    # ── Write output ────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output_config)), exist_ok=True)
    with open(args.output_config, "w") as fh:
        yaml.dump(shifted, fh, default_flow_style=False, sort_keys=False)

    # ── Report ──────────────────────────────────────────────────────────────
    print(f"Shifted config ({args.direction}) written to: {args.output_config}")
    metadata = shifted.get("stripe_ranges_metadata", {})
    for stripe, amt in metadata.get("shift_amounts", {}).items():
        orig = config["stripe_ranges"][stripe]
        new = shifted["stripe_ranges"][stripe]
        print(
            f"  {stripe:12s}: center {orig['center']:.4f} → {new['center']:.4f}  "
            f"(shift {amt:+.6f})"
        )


if __name__ == "__main__":
    main()
