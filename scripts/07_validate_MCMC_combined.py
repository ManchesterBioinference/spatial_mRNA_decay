#!/usr/bin/env python3
"""
Combined posterior predictive check for all four MCMC degradation rate models.

Produces a single scatter plot with all four models (spatial, null, age, biphasic)
overlaid, colored by model, with per-model lines of best fit and an x=y reference
line showing what perfect predictions would look like.

Usage:
    python 07_validate_MCMC_combined.py \\
        --spatial-chain  results/stripe2/e6/chains/degradation_chain.csv \\
        --null-chain     results/stripe2/e6_null/chains/degradation_chain.csv \\
        --age-chain      results/stripe2/e6_age/chains/degradation_chain.csv \\
        --biphasic-chain results/stripe2/e6_biphasic/chains/degradation_chain.csv \\
        --transcription  data/processed_transcription_data/transcription_traces_no_ids_stripe2_1200.csv \\
        --mrna           results/data/processed_mRNA_data_stripe2/e6_sass_formodel.csv \\
        --n-ap-bins 5 --n-dv-bins 5 \\
        --output results/comparison/stripe2/e6/posterior_predictive_check_combined.pdf
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid


# ---------------------------------------------------------------------------
# Solver functions (must match inference scripts exactly)
# ---------------------------------------------------------------------------

def solve_spatial(D_scalar, gamma, F_values, t_array):
    """
    Analytical solution for spatial (constant-D per AP bin) model.

    m(T) = gamma * integral( F(s) * exp(D*(s - T)) ds )
    """
    D_scalar = float(np.clip(D_scalar, 1e-6, 10.0))
    gamma = float(gamma)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)

    t_final = t_array[-1]
    exp_term = np.exp(D_scalar * (t_array - t_final))
    return gamma * trapezoid(F_values * exp_term, t_array)


def solve_analytical(D_age, gamma, F_values, t_array):
    """
    Age-dependent degradation model.

    m(T) = gamma * integral( F(t) * S(T-t) dt )
    where S(age) = exp(-cumsum(D_age) * dt)

    MUST match the numerical integration used in the inference scripts:
    - Rectangular rule (cumsum) for cumulative hazard
    - Trapezoidal rule for the convolution integral
    """
    D_age = np.asarray(D_age, dtype=np.float64)
    gamma = float(gamma)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)

    dt = t_array[1] - t_array[0]
    cumulative_hazard = np.cumsum(D_age) * dt
    survival_prob = np.exp(-cumulative_hazard)
    survival_profile_reversed = survival_prob[::-1]

    return gamma * trapezoid(F_values * survival_profile_reversed, t_array)


# ---------------------------------------------------------------------------
# Per-model prediction
# ---------------------------------------------------------------------------

def predict(chain_path, F_data, t_array, n_ap_bins, n_dv_bins):
    """
    Compute posterior-mean mRNA predictions for a single model.

    Returns
    -------
    m_pred : ndarray, shape (n_ap_bins * n_dv_bins,)
    """
    df = pd.read_csv(chain_path)
    n_timepoints = F_data.shape[1]

    D_cols = sorted(
        [c for c in df.columns if c.startswith("D[") and c.endswith("]")],
        key=lambda x: int(x.split("[")[1].split("]")[0]),
    )
    n_ages = len(D_cols)
    D_means_raw = np.array([df[c].mean() for c in D_cols])
    gamma_mean = df["gamma"].mean()

    is_spatial = n_ages == n_ap_bins
    is_null = n_ages == 1

    if is_null:
        D_age = np.full(n_timepoints, D_means_raw[0])
    elif not is_spatial:
        D_age = D_means_raw

    m_pred = []
    for i in range(n_ap_bins):
        for j in range(n_dv_bins):
            idx = i * n_dv_bins + j
            F_trace = F_data[idx, :]
            if is_spatial:
                val = solve_spatial(D_means_raw[i], gamma_mean, F_trace, t_array)
            else:
                val = solve_analytical(D_age, gamma_mean, F_trace, t_array)
            m_pred.append(val)

    return np.array(m_pred)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Combined PPC scatter plot for all four degradation models"
    )
    parser.add_argument("--spatial-chain",  required=True)
    parser.add_argument("--null-chain",     required=True)
    parser.add_argument("--age-chain",      required=True)
    parser.add_argument("--biphasic-chain", required=True)
    parser.add_argument("--transcription",  required=True)
    parser.add_argument("--mrna",           required=True)
    parser.add_argument("--n-ap-bins",  type=int, required=True)
    parser.add_argument("--n-dv-bins",  type=int, required=True)
    parser.add_argument("--output", default="posterior_predictive_check_combined.pdf")
    args = parser.parse_args()

    n_ap_bins = args.n_ap_bins
    n_dv_bins = args.n_dv_bins

    # Load shared data
    F_data = pd.read_csv(args.transcription, header=None).values
    m_obs  = pd.read_csv(args.mrna, header=None).values.flatten()
    n_timepoints = F_data.shape[1]
    t_array = np.arange(0, n_timepoints * 20, 20) / 60.0

    # Compute predictions for all four models
    models = {
        "Constant": args.null_chain,
        "Spatial":  args.spatial_chain,
        "Delayed":   args.age_chain,
        "Biphasic": args.biphasic_chain,
    }
    colors = {
        "Constant": "#ff7f0e",  # mpl orange
        "Spatial":  "#1f77b4",  # mpl blue
        "Delayed":   "#2ca02c",  # mpl green
        "Biphasic": "#9467bd",  # mpl purple
    }

    results = {}
    for name, chain_path in models.items():
        m_pred = predict(chain_path, F_data, t_array, n_ap_bins, n_dv_bins)
        residuals = m_obs - m_pred
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((m_obs - m_obs.mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        results[name] = {"m_pred": m_pred, "r2": r2}
        print(f"{name}: R² = {r2:.3f}")

    # Build combined figure
    MM = 1/25.4
    fig, ax = plt.subplots(figsize=(70*MM, 70*MM))

    all_vals = np.concatenate(
        [m_obs] + [r["m_pred"] for r in results.values()]
    )
    lo = max(0.0, np.nanmin(all_vals) * 0.95)
    hi = np.nanmax(all_vals) * 1.05

    # x=y reference line
    ax.plot([lo, hi], [lo, hi], color="red", linestyle="--", linewidth=1.5,
            alpha=0.7, zorder=1, label=None)

    for name, res in results.items():
        m_pred = res["m_pred"]
        r2     = res["r2"]
        color  = colors[name]

        # Scatter
        ax.scatter(m_obs, m_pred, color=color, edgecolors="k", linewidths=0.4,
                   s=15, alpha=0.75, zorder=2, label=f"{name}")# (R² = {r2:.3f})")

        # Line of best fit
        coeffs = np.polyfit(m_obs, m_pred, 1)
        x_fit = np.linspace(lo, hi, 200)
        ax.plot(x_fit, np.polyval(coeffs, x_fit), color=color,
                linewidth=1.5, zorder=3)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Observed mRNA (smFISH)", fontsize=12)
    ax.set_ylabel("Predicted mRNA", fontsize=12)
    ax.set_title("Posterior Predictive Check", fontsize=13)
    ax.legend(fontsize=7, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"\nCombined PPC plot saved to {args.output}")


if __name__ == "__main__":
    main()
