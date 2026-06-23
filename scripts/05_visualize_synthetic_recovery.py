#!/usr/bin/env python
"""
Visualise MCMC parameter recovery from synthetic (in silico) data.

Produces two panels that demonstrate how well the spatial inference model
recovers the manually-set degradation rates and transcription scaling factor
from noise-corrupted synthetic data.

Panel E — Posterior distributions with true-value markers
    One subplot per parameter (D[0]..D[4], gamma, sigma).
    - Grey shaded histogram of the posterior samples.
    - Yellow vertical line at the true parameter value.
    - Purple vertical line at the posterior mode (histogram peak).

Panel F — Recovery summary scatter plot
    Compares true vs. inferred (mode) values for all estimated parameters.
    Error bars show the 94% highest-density interval (HDI).
    A black dashed diagonal line marks perfect recovery (true == inferred).
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
import seaborn as sns


# ---------------------------------------------------------------------------
# Posterior utilities
# ---------------------------------------------------------------------------

def compute_mode(samples, n_bins=200):
    """
    Estimate the mode of a continuous distribution via a Gaussian KDE.

    Args:
        samples : 1-D array of MCMC samples
        n_bins  : number of evaluation points for the KDE

    Returns:
        float : estimated mode
    """
    samples = np.asarray(samples, dtype=np.float64)
    samples = samples[np.isfinite(samples)]
    if len(samples) < 2:
        return float(np.nanmedian(samples))
    try:
        kde = gaussian_kde(samples, bw_method="scott")
        x_grid = np.linspace(samples.min(), samples.max(), n_bins)
        return float(x_grid[np.argmax(kde(x_grid))])
    except Exception:
        return float(np.nanmedian(samples))


def compute_hdi(samples, credible_mass=0.94):
    """
    Compute the Highest Density Interval (HDI) for a 1-D sample array.

    Args:
        samples       : 1-D array of MCMC samples
        credible_mass : probability mass to include (default 0.94)

    Returns:
        (lower, upper) tuple of HDI bounds
    """
    samples = np.sort(np.asarray(samples, dtype=np.float64))
    n = len(samples)
    interval_idx_inc = int(np.floor(credible_mass * n))
    n_intervals = n - interval_idx_inc
    if n_intervals < 1:
        return float(samples[0]), float(samples[-1])
    interval_width = samples[interval_idx_inc:] - samples[:n_intervals]
    min_idx = int(np.argmin(interval_width))
    return float(samples[min_idx]), float(samples[min_idx + interval_idx_inc])


# ---------------------------------------------------------------------------
# Load chain
# ---------------------------------------------------------------------------

def load_chain(chain_path, n_ap_bins=5):
    """
    Load MCMC chain CSV produced by 02_infer_degradation_rates_spatial.py.

    Expected columns: D[0]..D[n_ap_bins-1], gamma, sigma

    Returns:
        dict mapping parameter name -> 1-D numpy array of samples
    """
    df = pd.read_csv(chain_path)
    samples = {}

    for i in range(n_ap_bins):
        col = f"D[{i}]"
        if col in df.columns:
            samples[col] = df[col].values

    for col in ("gamma", "sigma"):
        if col in df.columns:
            samples[col] = df[col].values

    return samples


# ---------------------------------------------------------------------------
# Panel E — posterior distributions
# ---------------------------------------------------------------------------

def plot_panel_e(samples, true_params, output_path):
    """
    Panel E: posterior histograms with yellow (true) and purple (mode) lines.

    Args:
        samples     : dict of parameter name -> array of MCMC samples
        true_params : dict of parameter name -> true scalar value
        output_path : path to save the PDF
    """
    param_names = list(samples.keys())[:-1]
    n_params = len(param_names)
    n_cols = 3
    n_rows = 6 #int(np.ceil(n_params / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, 1.5 * n_rows))#,sharex=True)

    # find the global min and max across all parameters for consistent x-axis limits. this must exclude the sigma parameter, which can be orders of magnitude larger than the D and gamma parameters and would make the x-axis limits too wide for the others.
    all_samples = np.concatenate([s for pname, s in samples.items() if pname != "sigma"])
    global_max = np.max(all_samples)
    colors = plt.cm.viridis(np.linspace(0, 1, 4))

    for row, pname in enumerate(param_names):
        s = samples[pname]
        #split s into 4 equal parts for plotting the traces in different colors, to help visually separate them. this is a bit hacky but it helps make the trace plots more readable when there are many samples. the colors will cycle through the 4 parts, so if there are more than 4 parts it will repeat the colors.
        s1, s2, s3, s4 = np.array_split(s, 4)
        if pname[0] == "D":
            name = f"D{int(pname[2])+1}"
        else:
            name = pname

        # Left column: Iteration/caterpillar plot
        ax_left = axes[row, 0]
        ax_left.plot(s1, alpha=0.8, color = colors[3])
        ax_left.plot(s2, alpha=0.8, color = colors[1])
        ax_left.plot(s3, alpha=0.8, color = colors[2])
        ax_left.plot(s4, alpha=0.8, color = colors[0])
        ax_left.set_title(name)
        if row == 2:
            ax_left.set_ylabel('Sample value', fontsize=13)
        if row == n_rows - 1:
            ax_left.set_xlabel('Iteration', fontsize=13)
        ax_left.grid(True, alpha=0.3)
        
        # middle column: Density plot
        ax_right = axes[row, 1]
        sns.kdeplot(s1, ax=ax_right, alpha=0.7, color=colors[3])
        sns.kdeplot(s2, ax=ax_right, alpha=0.7, color=colors[1])
        sns.kdeplot(s3, ax=ax_right, alpha=0.7, color=colors[2])
        sns.kdeplot(s4, ax=ax_right, alpha=0.7, color=colors[0])
        ax_right.set_title(name)
        if row == 2:
            ax_right.set_ylabel('Density', fontsize=13)
        else:
            ax_right.set_ylabel('')
        if row == n_rows - 1:
            ax_right.set_xlabel('Parameter value', fontsize=13)
        ax_right.grid(True, alpha=0.3)
        #ax_right.set_xlim(0, global_max+0.1)

        ax = axes[row,2]
        # Histogram
        ax.hist(s, bins=60, edgecolor='#cd96cd', color='#d8bfd8', alpha=0.7, 
                density=True, label=f"{name}")

        # True value
        if pname in true_params:
            ax.axvline(true_params[pname], color="gold", linewidth=2.0,
                       label=f"True: {true_params[pname]:.4g}")

        # Median
        median_val = np.nanmedian(s)
        ax.axvline(median_val, color='#332288', linewidth=2.0,# linestyle="--",
                   label=f"Median: {median_val:.4g}")

        ax.set_title(name)
        if name != "D5":
            ax.legend(loc="upper right")
        else:
            ax.legend(loc="upper left")
        ax.set_xlim(0, global_max+0.1)
        if row == 2:
            ax.set_ylabel('Frequency', fontsize=13)
        if row == n_rows - 1:
            ax.set_xlabel('Parameter value', fontsize=13)

    # Hide unused subplots
    for ax in axes[n_params:]:
        ax.set_visible(False)

    fig.suptitle("MCMC Chains                  Parameter posterior distributions", y=1.01, fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Panel E saved to: {output_path}")


# ---------------------------------------------------------------------------
# Panel F — recovery summary
# ---------------------------------------------------------------------------

def plot_panel_f(samples, true_params, output_path, credible_mass=0.94):
    """
    Panel F: true vs. inferred (mode) scatter with 94% HDI error bars.

    Includes only the D[0]..D[4] and gamma parameters (sigma is noise, not
    a target of the parameter recovery test).

    Args:
        samples       : dict of parameter name -> array of MCMC samples
        true_params   : dict of parameter name -> true scalar value
        output_path   : path to save the PDF
        credible_mass : HDI mass for error bars (default 0.94)
    """
    # Collect parameters to display (skip sigma)
    target_params = [k for k in samples.keys() if k != "sigma"]

    true_vals, mode_vals, lower_vals, upper_vals, labels = [], [], [], [], []

    for pname in target_params:
        if pname not in true_params:
            continue
        s = samples[pname]
        mode_val = compute_mode(s)
        lo, hi = compute_hdi(s, credible_mass)
        true_vals.append(true_params[pname])
        mode_vals.append(mode_val)
        lower_vals.append(mode_val - lo)
        upper_vals.append(hi - mode_val)
        labels.append(pname)

    true_vals = np.array(true_vals)
    mode_vals = np.array(mode_vals)
    lower_vals = np.array(lower_vals)
    upper_vals = np.array(upper_vals)

    # Use log scale when values span multiple orders of magnitude
    span = true_vals.max() / true_vals.min() if true_vals.min() > 0 else 1
    use_log = span > 20

    fig, ax = plt.subplots(figsize=(5.5, 5))

    # Error bars
    ax.errorbar(
        true_vals, mode_vals,
        yerr=[lower_vals, upper_vals],
        fmt="none",
        ecolor="purple", elinewidth=1.2, capsize=4, alpha=0.8,
        label=f"{int(credible_mass * 100)}% HDI"
    )

    # Points
    ax.scatter(true_vals, mode_vals, color="purple", s=60, zorder=5)

    # Labels next to points
    for x, y, lbl in zip(true_vals, mode_vals, labels):
        ax.annotate(lbl, (x, y), textcoords="offset points",
                    xytext=(6, 3), fontsize=7)

    # Perfect-recovery diagonal
    all_vals = np.concatenate([true_vals, mode_vals])
    lo_diag = all_vals.min() * 0.8
    hi_diag = all_vals.max() * 1.3
    ax.plot([lo_diag, hi_diag], [lo_diag, hi_diag],
            "k--", linewidth=1.2, label="Perfect recovery")

    if use_log:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel("True value", fontsize=10)
    ax.set_ylabel("Inferred value (posterior mode)", fontsize=10)
    ax.set_title("Panel F — Parameter recovery summary", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Panel F saved to: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Visualise MCMC parameter recovery from synthetic mRNA data"
    )
    parser.add_argument("--chain", required=True,
                        help="MCMC chain CSV from spatial inference")
    parser.add_argument("--D-values", type=float, nargs=5, required=True,
                        metavar=("D0", "D1", "D2", "D3", "D4"),
                        help="True degradation rates for AP bins 0-4 (min^-1)")
    parser.add_argument("--gamma", type=float, required=True,
                        help="True gamma (transcription scaling factor)")
    parser.add_argument("--noise-level", type=float, default=None,
                        help="True sigma / noise level (shown on sigma posterior, optional)")
    parser.add_argument("--n-ap-bins", type=int, default=5,
                        help="Number of AP bins (default: 5)")
    parser.add_argument("--output-posteriors", required=True,
                        help="Output path for Panel E (posterior distributions PDF)")
    parser.add_argument("--output-summary", required=True,
                        help="Output path for Panel F (recovery summary PDF)")

    args = parser.parse_args()

    # Load chain
    samples = load_chain(args.chain, args.n_ap_bins)
    if not samples:
        raise RuntimeError(f"No recognised parameter columns found in: {args.chain}")
    print(f"Loaded {len(next(iter(samples.values())))} samples for "
          f"{len(samples)} parameters from: {args.chain}")

    # Build true-parameter dictionary
    true_params = {}
    for i, d in enumerate(args.D_values):
        true_params[f"D[{i}]"] = d
    true_params["gamma"] = args.gamma
    if args.noise_level is not None:
        true_params["sigma"] = args.noise_level

    print(f"\nTrue parameters:")
    for k, v in true_params.items():
        print(f"  {k} = {v}")

    print(f"\nPosterior modes:")
    for pname, s in samples.items():
        mode = compute_mode(s)
        lo, hi = compute_hdi(s)
        true_v = true_params.get(pname, None)
        tag = f"  (true: {true_v:.4g})" if true_v is not None else ""
        print(f"  {pname}: mode={mode:.4g}  94%HDI=[{lo:.4g}, {hi:.4g}]{tag}")

    import os
    os.makedirs(os.path.dirname(args.output_posteriors), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_summary), exist_ok=True)

    # Panel E
    plot_panel_e(samples, true_params, args.output_posteriors)

    # Panel F
    plot_panel_f(samples, true_params, args.output_summary)


if __name__ == "__main__":
    main()
