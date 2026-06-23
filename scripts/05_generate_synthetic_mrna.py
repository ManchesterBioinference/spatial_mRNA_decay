#!/usr/bin/env python
"""
Generate synthetic (in silico) mRNA data for parameter recovery validation.

Uses the real stripe2 MS2 transcription traces as the forward model input and
manually-specified degradation rates to produce synthetic mRNA counts. Gaussian
noise is added to mimic smFISH measurement variability.

Forward model (identical to 02_infer_degradation_rates_spatial.py):
    m(T) = gamma * integral( F(t) * exp(D*(t - T)) dt )   [trapezoidal rule]

Spatial structure:
    25 spatial bins = 5 AP positions × 5 DV traces per AP position.
    One D value is assigned per AP position; all 5 DV traces in that AP position
    use the same D (matching the inference model structure).

Outputs:
    --output-mrna        : 25 noisy values, no header (inference-ready format)
    --output-data-table  : 3-column CSV (simulated_data, noise, noisy_mRNA) for reference
    --output-figure      : Panel C — row 1: true D[0..4] as 1×5 heatmap;
                                      row 2: noisy synthetic mRNA as 5×5 heatmap
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmcrameri.cm as cmc
from scipy.integrate import trapezoid


# ---------------------------------------------------------------------------
# Forward model
# ---------------------------------------------------------------------------

def solve_ode_analytical(D, gamma, F_values, t_array):
    """
    Analytical solution to dm/dt = gamma*F(t) - D*m with m(0) = 0.

    m(T) = gamma * integral_0^T F(s) * exp(D*(s - T)) ds

    Numerically stable formulation using exp(D*(t - T)) which stays <= 1
    for t <= T and D >= 0.

    Matches solve_ode_analytical() in 02_infer_degradation_rates_spatial.py.

    Args:
        D      : degradation rate (min^-1), scalar
        gamma  : transcription scaling factor, scalar
        F_values: transcription values at each time point, shape (n_timepoints,)
        t_array : time points in minutes, shape (n_timepoints,)

    Returns:
        float : expected mRNA concentration at final time point
    """
    D = float(np.clip(D, 1e-6, 10.0))
    gamma = float(gamma)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)

    t_final = t_array[-1]
    exp_term = np.exp(D * (t_array - t_final))      # exp(D*(t-T)) ∈ [0, 1]
    integral = trapezoid(F_values * exp_term, t_array)
    m_final = gamma * integral

    return float(m_final) if np.isfinite(m_final) else 0.0


def generate_synthetic_mrna(F_data, D_values, gamma, noise_level, noise_seed,
                             n_ap_bins=5, n_dv_bins=5):
    """
    Generate synthetic mRNA for all 25 spatial bins.

    Structure mirrors the inference model:
        - Row ordering in F_data: AP bin 0 DV 0-4, AP bin 1 DV 0-4, ..., AP bin 4 DV 0-4
        - D_values[i] is applied to all 5 DV traces in AP bin i
        - Noise: iid N(0, noise_level); negative values clipped to 0

    Args:
        F_data     : transcription data, shape (25, n_timepoints)
        D_values   : degradation rates per AP bin, shape (5,)
        gamma      : transcription scaling factor
        noise_level: std of Gaussian noise (absolute mRNA units)
        noise_seed : random seed for reproducibility
        n_ap_bins  : number of AP bins (default 5)
        n_dv_bins  : number of DV bins per AP position (default 5)

    Returns:
        simulated  : noise-free expected mRNA, shape (25,)
        noise      : sampled noise, shape (25,)
        noisy_mrna : simulated + noise, clipped >= 0, shape (25,)
        t_array    : time points in minutes, shape (n_timepoints,)
    """
    n_timepoints = F_data.shape[1]

    # Time axis: 0 to max_time (seconds) sampled every 20 s, converted to minutes
    # Interval = 20 s = 1/3 min; 61 points span 0 to 20 min
    dt_seconds = 20.0
    t_seconds = np.arange(n_timepoints) * dt_seconds
    t_array = t_seconds / 60.0   # convert to minutes (same convention as inference)

    simulated = np.zeros(n_ap_bins * n_dv_bins)

    for ap_idx in range(n_ap_bins):
        D = D_values[ap_idx]
        for dv_idx in range(n_dv_bins):
            row = ap_idx * n_dv_bins + dv_idx
            F_trace = F_data[row, :]
            simulated[row] = solve_ode_analytical(D, gamma, F_trace, t_array)

    rng = np.random.default_rng(noise_seed)
    noise = rng.normal(0, noise_level, size=simulated.shape)
    noisy_mrna = np.maximum(simulated + noise, 0.0)

    return simulated, noise, noisy_mrna, t_array


# ---------------------------------------------------------------------------
# Panel C visualisation
# ---------------------------------------------------------------------------

def plot_panel_c(D_values, noisy_mrna, output_path, n_ap_bins=5, n_dv_bins=5):
    """
    Panel C: two-row figure.

    Row 1 — True D values displayed as a 1×5 heatmap (one cell per AP bin).
    Row 2 — Synthetic (noisy) mRNA displayed as a 5×5 heatmap (AP × DV).

    Args:
        D_values  : true degradation rates, shape (5,)
        noisy_mrna: noisy mRNA values to display, shape (25,)
        output_path: path to save PDF
        n_ap_bins : 5
        n_dv_bins : 5
    """
    mrna_grid = noisy_mrna.reshape(n_ap_bins, n_dv_bins)

    fig = plt.figure(figsize=(3, 4))
    gs = gridspec.GridSpec(
        2, 1,
        height_ratios=[1, 5],
        hspace=0.5,
        left=0.12, right=0.92, top=0.93, bottom=0.08
    )

    # ---- Row 1: True D heatmap (1 × 5) ----
    ax_d = fig.add_subplot(gs[0])
    D_grid = np.array(D_values).reshape(1, -1)
    im_d = ax_d.imshow(D_grid, aspect="auto", cmap=cmc.nuuk_r)
    ax_d.set_yticks([])
    ax_d.set_xticks(range(n_ap_bins))
    ax_d.set_xticklabels([f"D{i+1}" for i in range(n_ap_bins)], fontsize=11)
    ax_d.set_title(r"$\mathit{in\ silico}$ D values", fontsize=12, pad=4)
    for j, d in enumerate(D_values):
        if j == 0:
            ax_d.text(j, 0, f"{d:.3f}", ha="center", va="center", fontsize=8, color="black", fontweight="bold")
        elif j in [1,2]:
            ax_d.text(j, 0, f"{d:.2f}", ha="center", va="center", fontsize=8, color="black", fontweight="bold")
        else:
            ax_d.text(j, 0, f"{d:.1f}", ha="center", va="center", fontsize=8, color="white" if d > D_grid.max() * 0.5 else "black", fontweight="bold")
    cb_d = fig.colorbar(im_d, ax=ax_d, orientation="vertical",
                        pad=0.02, fraction=0.04)
    cb_d.set_label("D (min⁻¹)", fontsize=10)

    # ---- Row 2: Synthetic mRNA heatmap (5 AP × 5 DV) ----
    ax_m = fig.add_subplot(gs[1])
    im_m = ax_m.imshow(mrna_grid.T, aspect="auto", cmap=cmc.lipari, origin="lower")
    ax_m.set_xlabel("AP bin", fontsize=11)
    ax_m.set_ylabel("DV bin", fontsize=11)
    #ax_m.set_yticks(range(n_dv_bins))
    #ax_m.set_xticks(range(n_ap_bins))
    ax_m.set_yticklabels([f"" for j in range(n_dv_bins)], fontsize=7)
    ax_m.set_xticklabels([f"" for i in range(n_ap_bins)], fontsize=7)
    #ax_m.set_yticklabels([f""DV {j}" for j in range(n_dv_bins)], fontsize=7)
    #ax_m.set_xticklabels([f""AP {i}" for i in range(n_ap_bins)], fontsize=7)
    ax_m.set_title(r"$\mathit{in\ silico}$ mRNA data", fontsize=12, pad=4)
    # for i in range(n_ap_bins):
    #     for j in range(n_dv_bins, -1): #reverse order for correct orientation
    #         val = mrna_grid[i, j]
    #         ax_m.text(j, i, f"{val:.2f}", ha="center", va="center",
    #                  fontsize=6, color="white" if val > mrna_grid.max() * 0.5 else "black")
    cb_m = fig.colorbar(im_m, ax=ax_m, orientation="vertical",
                        pad=0.02, fraction=0.04)
    cb_m.set_label("mRNA / nucleus", fontsize=10)

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Panel C saved to: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic mRNA data for parameter recovery validation"
    )
    parser.add_argument("--transcription", required=True,
                        help="Transcription traces CSV (25 rows × n_timepoints, no header)")
    parser.add_argument("--output-mrna", required=True,
                        help="Output path for inference-ready mRNA CSV (25 values, no header)")
    parser.add_argument("--output-data-table", required=True,
                        help="Output path for 3-column reference CSV "
                             "(simulated_data, noise, noisy_mRNA)")
    parser.add_argument("--output-figure", required=True,
                        help="Output path for Panel C figure (PDF)")
    parser.add_argument("--gamma", type=float, required=True,
                        help="True transcription scaling factor gamma")
    parser.add_argument("--D-values", type=float, nargs=5, required=True,
                        metavar=("D0", "D1", "D2", "D3", "D4"),
                        help="True degradation rates for AP bins 0-4 (min^-1)")
    parser.add_argument("--noise-level", type=float, default=1.0,
                        help="Gaussian noise std (absolute mRNA units, default: 1.0)")
    parser.add_argument("--noise-seed", type=int, default=42,
                        help="Random seed for noise generation (default: 42)")
    parser.add_argument("--n-ap-bins", type=int, default=5,
                        help="Number of AP bins (default: 5)")
    parser.add_argument("--n-dv-bins", type=int, default=5,
                        help="Number of DV bins per AP position (default: 5)")

    args = parser.parse_args()

    # Load transcription data
    print(f"Loading transcription data: {args.transcription}")
    F_data = pd.read_csv(args.transcription, header=None).values
    print(f"  Transcription data shape: {F_data.shape}")

    # Generate synthetic mRNA
    print(f"\nGenerating synthetic mRNA:")
    print(f"  gamma     = {args.gamma}")
    print(f"  D_values  = {args.D_values}")
    print(f"  noise_level = {args.noise_level}")
    print(f"  noise_seed  = {args.noise_seed}")

    simulated, noise, noisy_mrna, t_array = generate_synthetic_mrna(
        F_data, args.D_values, args.gamma,
        args.noise_level, args.noise_seed,
        args.n_ap_bins, args.n_dv_bins
    )

    print(f"\nSynthetic mRNA summary:")
    print(f"  Clean values range: [{simulated.min():.4f}, {simulated.max():.4f}]")
    print(f"  Noise range:        [{noise.min():.4f}, {noise.max():.4f}]")
    print(f"  Noisy values range: [{noisy_mrna.min():.4f}, {noisy_mrna.max():.4f}]")

    # Warn if noise is large relative to signal
    sig_mean = simulated.mean()
    if sig_mean > 0 and args.noise_level / sig_mean > 0.3:
        print(f"\nWARNING: noise_level ({args.noise_level:.3f}) is >30% of mean signal "
              f"({sig_mean:.4f}). Consider reducing noise_level in config.yaml.")

    import os

    # Save inference-ready mRNA CSV (25 values, no header, one per row)
    os.makedirs(os.path.dirname(args.output_mrna), exist_ok=True)
    pd.DataFrame(noisy_mrna).to_csv(args.output_mrna, index=False, header=False)
    print(f"\nSaved inference-ready mRNA to: {args.output_mrna}")

    # Save 3-column reference table
    os.makedirs(os.path.dirname(args.output_data_table), exist_ok=True)
    df_table = pd.DataFrame({
        "simulated_data": simulated,
        "noise": noise,
        "noisy_mRNA": noisy_mrna
    })
    df_table.to_csv(args.output_data_table, index=False)
    print(f"Saved data table to: {args.output_data_table}")

    # Panel C figure
    os.makedirs(os.path.dirname(args.output_figure), exist_ok=True)
    plot_panel_c(args.D_values, noisy_mrna, args.output_figure,
                 args.n_ap_bins, args.n_dv_bins)


if __name__ == "__main__":
    main()
