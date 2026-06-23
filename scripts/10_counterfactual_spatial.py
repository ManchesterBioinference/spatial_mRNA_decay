#!/usr/bin/env python3
"""
Counterfactual analysis for the spatial mRNA decay model.

Tests whether spatially-varying degradation rates are necessary to explain the
observed mRNA stripe pattern by running two counterfactual scenarios:

  1. Fast-everywhere: all AP bins use the *edge* degradation rate
     (arithmetic mean of the two outermost bins' posterior means).
     Shows that fast degradation everywhere prevents the high mRNA levels
     seen at the stripe centre.

  2. Slow-everywhere: all AP bins use the *centre* degradation rate
     (posterior mean of the minimum-D bin, identified dynamically so the
     script works even when the slowest-decay region is off-centre).
     Shows that slow degradation everywhere would produce an overly wide stripe.

The centre bin is identified as argmin(D_means), not assumed to be the
geometric middle bin.  The edge rate is the average of the first and last AP
bin posterior means, so the result is robust to asymmetric stripe positions.

Requires a fitted *spatial* model chain (D[0] … D[n_ap_bins-1] columns).

Usage
-----
    python 10_counterfactual_spatial.py \\
        --chain results_1200/stripe2/e7um/chains/degradation_chain.csv \\
        --transcription data/processed_transcription_data/transcription_traces_no_ids_stripe2_1200.csv \\
        --mrna results_1200/data/processed_mRNA_data_stripe2/e7um_sass_formodel.csv \\
        --n-ap-bins 3 --n-dv-bins 5 \\
        --output results_1200/stripe2/e7um/figures/counterfactual_spatial.pdf
"""

import argparse
import sys
from pathlib import Path

import cmcrameri.cm as cmc
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid


# ── Core solver ───────────────────────────────────────────────────────────────

def solve_spatial(D_scalar: float, gamma: float, F_values: np.ndarray, t_array: np.ndarray) -> float:
    """
    Analytical solution for the spatial (constant-D) model: dm/dt = γ·F(t) − D·m.

    Identical to solve_ode_analytical in 02_infer_degradation_rates_spatial.py:
        m(T) = γ · ∫ F(s) · exp(D·(s − T)) ds

    Parameters
    ----------
    D_scalar : constant degradation rate for this AP bin (scalar, min⁻¹)
    gamma    : transcription scaling factor (scalar)
    F_values : transcription trace at each time point (array, length T)
    t_array  : time points in minutes (array, length T)

    Returns
    -------
    Predicted mRNA at the final time point (scalar).
    """
    D_scalar = float(np.clip(D_scalar, 1e-6, 10.0))
    gamma = float(gamma)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)

    t_final = t_array[-1]
    exp_term = np.exp(D_scalar * (t_array - t_final))
    return gamma * trapezoid(F_values * exp_term, t_array)


def predict_mrna(
    D_per_bin: np.ndarray,
    gamma: float,
    F_data: np.ndarray,
    t_array: np.ndarray,
    n_ap_bins: int,
    n_dv_bins: int,
) -> np.ndarray:
    """
    Predict mRNA for all AP × DV traces given one D value per AP bin.

    Parameters
    ----------
    D_per_bin : degradation rate assigned to each AP bin (shape: n_ap_bins)
    gamma     : shared transcription scaling factor
    F_data    : transcription matrix (shape: n_ap_bins*n_dv_bins × n_timepoints)
    t_array   : time points in minutes
    n_ap_bins, n_dv_bins : spatial bin counts

    Returns
    -------
    Flat predicted mRNA array (length n_ap_bins * n_dv_bins), ordered AP-major.
    """
    m_pred = np.empty(n_ap_bins * n_dv_bins, dtype=np.float64)
    for i in range(n_ap_bins):
        for j in range(n_dv_bins):
            idx = i * n_dv_bins + j
            m_pred[idx] = solve_spatial(D_per_bin[i], gamma, F_data[idx, :], t_array)
    return m_pred


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Counterfactual spatial mRNA decay: "
            "edge-speed-everywhere vs. centre-speed-everywhere"
        )
    )
    parser.add_argument("--chain",         required=True, help="Spatial model MCMC chain CSV")
    parser.add_argument("--transcription", required=True, help="Transcription traces CSV (no IDs)")
    parser.add_argument("--mrna",          required=True, help="Observed mRNA CSV (flat, for-model file)")
    parser.add_argument("--n-ap-bins",     type=int, required=True, help="Number of AP bins")
    parser.add_argument("--n-dv-bins",     type=int, required=True, help="Number of DV bins")
    parser.add_argument("--output", default="counterfactual_spatial.pdf",
                        help="Output PNG path")
    args = parser.parse_args()

    n_ap = args.n_ap_bins
    n_dv = args.n_dv_bins

    # ── 1. Load data ──────────────────────────────────────────────────────────
    df_samples = pd.read_csv(args.chain)
    F_data     = pd.read_csv(args.transcription, header=None).values          # (n_ap*n_dv, T)
    m_obs_flat = pd.read_csv(args.mrna, header=None).values.flatten()         # (n_ap*n_dv,)

    n_timepoints = F_data.shape[1]
    t_array = np.arange(0, n_timepoints * 20, 20) / 60.0                      # seconds → minutes

    # ── 2. Verify spatial model chain ─────────────────────────────────────────
    D_cols = sorted(
        [c for c in df_samples.columns if c.startswith("D[") and c.endswith("]")],
        key=lambda x: int(x[2:-1]),
    )
    n_d_cols = len(D_cols)
    if n_d_cols != n_ap:
        print(
            f"ERROR: Chain has {n_d_cols} D column(s) but --n-ap-bins={n_ap}.\n"
            f"This script requires a spatial model chain (one D per AP bin).\n"
            f"Detected columns: {D_cols}",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── 3. Posterior medians ──────────────────────────────────────────────
    D_medians    = np.array([float(np.median(df_samples[c].values)) for c in D_cols])
    gamma_median = float(np.median(df_samples["gamma"]))

    # Centre bin = AP bin with slowest decay (minimum D — not necessarily the middle)
    centre_bin = int(np.argmin(D_medians))
    D_slow     = float(D_medians[centre_bin])
    D_edge     = float((D_medians[0] + D_medians[-1]) / 2.0)

    print("Posterior median D per AP bin:")
    for i, d in enumerate(D_medians):
        tags = []
        if i == centre_bin:
            tags.append("← CENTRE (min D)")
        if i == 0 or i == n_ap - 1:
            tags.append("[edge]")
        tag_str = "  " + "  ".join(tags) if tags else ""
        print(f"  AP bin {i + 1}: D = {d:.4f} min⁻¹  (t½ = {np.log(2) / d:.1f} min){tag_str}")

    print(
        f"\nEdge D (avg of AP bin 1 and AP bin {n_ap}): "
        f"{D_edge:.4f} min⁻¹  (t½ = {np.log(2) / D_edge:.1f} min)"
    )
    print(
        f"Centre D (AP bin {centre_bin + 1}):              "
        f"{D_slow:.4f} min⁻¹  (t½ = {np.log(2) / D_slow:.1f} min)"
    )
    print(f"Gamma:                               {gamma_median:.4f}")

    # ── 4. Three prediction scenarios ─────────────────────────────────────────
    # A – Fitted: each AP bin uses its own inferred D
    m_fitted = predict_mrna(D_medians,              gamma_median, F_data, t_array, n_ap, n_dv)
    # B – Fast everywhere: all bins use the edge (faster) rate
    m_fast   = predict_mrna(np.full(n_ap, D_edge), gamma_median, F_data, t_array, n_ap, n_dv)
    # C – Slow everywhere: all bins use the centre (slowest) rate
    m_slow   = predict_mrna(np.full(n_ap, D_slow), gamma_median, F_data, t_array, n_ap, n_dv)

    # ── 5. Fit statistics (fitted scenario only) ───────────────────────────────
    residuals = m_obs_flat - m_fitted
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((m_obs_flat - m_obs_flat.mean()) ** 2))
    r2     = 1.0 - ss_res / ss_tot

    print(f"\nFitted model — R² = {r2:.3f}, RMSE = {rmse:.2f}")

    # ── 6. Reshape to (n_ap, n_dv) grids ────────────────────────────────────
    obs_grid    = m_obs_flat.reshape(n_ap, n_dv)
    fitted_grid = m_fitted.reshape(n_ap, n_dv)
    fast_grid   = m_fast.reshape(n_ap, n_dv)
    slow_grid   = m_slow.reshape(n_ap, n_dv)

    # ── 7. Figure ─────────────────────────────────────────────────────────────
    # Layout (2 rows × 4 cols, heatmaps in cols 0-1 of both rows):
    #   Row 1: Observed heatmap | Fitted heatmap   | D-profile bar chart (cols 2-3)
    #   Row 2: Fast heatmap     | Slow heatmap     | Mean-mRNA line plot (cols 2-3)
    fig = plt.figure(figsize=(22, 13))
    gs  = gridspec.GridSpec(
        2, 3,
        height_ratios=[1.0, 1.0],
        hspace=0.50,
        wspace=0.35,
    )

    # Shared colour scale across all four heatmaps
    all_values = np.concatenate([obs_grid, fitted_grid, fast_grid, slow_grid])
    vmax = float(all_values.max()) * 1.05
    cmap = cmc.lipari
    imshow_kw = dict(
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=vmax,
    )

    ap_labels = [f"AP {i + 1}" for i in range(n_ap)]
    dv_labels = [f"DV {j + 1}" for j in range(n_dv)]

    # (grid_row, grid_col, data, title, colour)
    heatmap_specs = [
        (0, 0, obs_grid,    "Observed mRNA",                                              "black"),
        (0, 1, fitted_grid, f"Fitted  (spatial model)\nR² = {r2:.2f}, RMSE = {rmse:.1f}", "#1a7abf"),
        (1, 0, fast_grid,   f"Counterfactual: fast everywhere\n"
                             f"D_edge = {D_edge:.3f} min⁻¹  (t½ = {np.log(2)/D_edge:.0f} min)",
                                                                                           "#c0392b"),
        (1, 1, slow_grid,   f"Counterfactual: slow everywhere\n"
                             f"D_centre = {D_slow:.3f} min⁻¹  (t½ = {np.log(2)/D_slow:.0f} min)",
                                                                                           "#27ae60"),
    ]

    heatmap_axes = []
    img_ref = None
    for grid_row, grid_col, grid, title, _ in heatmap_specs:
        ax = fig.add_subplot(gs[grid_row, grid_col])
        heatmap_axes.append(ax)
        img = ax.imshow(grid.T, **imshow_kw)
        if img_ref is None:
            img_ref = img

        ax.set_xticks(range(n_ap))
        ax.set_xticklabels(ap_labels, fontsize=9)
        ax.set_yticks(range(n_dv))
        ax.set_yticklabels(dv_labels, fontsize=9)
        ax.set_xlabel("AP bin", fontsize=10)
        if grid_col == 0:
            ax.set_ylabel("DV bin", fontsize=10)
        ax.set_title(title, fontsize=10, pad=7)

#        # Mark the identified centre AP bin with a white dashed line
#        ax.axvline(x=centre_bin, color="white", linewidth=1.8, linestyle="--", alpha=0.85,
#                   label=f"Centre AP bin ({centre_bin + 1})")
        if grid_row == 0 and grid_col == 0:
            ax.legend(fontsize=7, loc="upper right", framealpha=0.6)

    # Shared colourbar to the right of the 2×2 heatmap block (cols 0-1, both rows)
    cbar = fig.colorbar(img_ref, ax=heatmap_axes, shrink=0.7, pad=0.04, aspect=35)
    cbar.set_label("mRNA count", fontsize=10)

    # ── Row 2, right: mean mRNA per AP bin (line plot) ────────────────────────
    ax_line = fig.add_subplot(gs[1, 2])

    ap_x = np.arange(n_ap)
    line_panels = [
        (obs_grid,    "Observed",              "black",   "-",  "o"),
        (fitted_grid, "Fitted (spatial model)", "#1a7abf", "-",  "s"),
        (fast_grid,   "Fast everywhere",        "#c0392b", "--", "^"),
        (slow_grid,   "Slow everywhere",        "#27ae60", ":",  "D"),
    ]
    for grid, label, colour, ls, marker in line_panels:
        ap_means = grid.mean(axis=1)   # average over DV bins → shape (n_ap,)
        ax_line.plot(
            ap_x, ap_means,
            color=colour, linestyle=ls, marker=marker,
            linewidth=2, markersize=7, label=label,
        )

    ax_line.axvline(x=centre_bin, color="grey", linewidth=1.2, linestyle="--",
                    alpha=0.6, label=f"Centre AP bin ({centre_bin + 1})")
    ax_line.set_xticks(ap_x)
    ax_line.set_xticklabels(ap_labels, fontsize=10)
    ax_line.set_xlabel("AP bin", fontsize=11)
    ax_line.set_ylabel("Mean mRNA count (avg over DV bins)", fontsize=11)
    ax_line.set_title("Mean mRNA per AP bin across scenarios", fontsize=10)
    ax_line.legend(fontsize=9, loc="upper right")
    ax_line.grid(True, alpha=0.3)
    # Expand top so the legend floats above the data
    y_max_line = max(
        g.mean(axis=1).max()
        for g in [obs_grid, fitted_grid, fast_grid, slow_grid]
    )
    ax_line.set_ylim(bottom=0, top=y_max_line * 1.45)

    # ── Row 1, right: D-profile bar chart ─────────────────────────────────────
    ax_d = fig.add_subplot(gs[0, 2])

    bar_x  = np.arange(n_ap)
    bar_colours = ["#4c72b0"] * n_ap

    # Colour edge bins red, centre bin green
    bar_colours[0]         = "#c0392b"
    bar_colours[-1]        = "#c0392b"
    bar_colours[centre_bin] = "#27ae60"

    ax_d.bar(
        bar_x, D_medians,
        color=bar_colours,
        edgecolor="black",
        linewidth=0.8,
        alpha=0.85,
        zorder=3,
    )

    # Reference lines
    ax_d.axhline(
        D_edge, color="#c0392b", linestyle="--", linewidth=1.8, zorder=4,
        label=f"Edge D (avg bins 1 & {n_ap}) = {D_edge:.4f} min⁻¹",
    )
    ax_d.axhline(
        D_slow, color="#27ae60", linestyle=":",  linewidth=1.8, zorder=4,
        label=f"Centre D (bin {centre_bin + 1}) = {D_slow:.4f} min⁻¹",
    )

    # Half-life annotations above each bar
    y_top = D_medians.max()
    for i, d in enumerate(D_medians):
        hl = np.log(2) / d
        ax_d.text(
            i, d + y_top * 0.03,
            f"t½={hl:.0f} min",
            ha="center", va="bottom", fontsize=8, color="black",
        )

    ax_d.set_xticks(bar_x)
    ax_d.set_xticklabels(ap_labels, fontsize=10)
    ax_d.set_xlabel("AP bin", fontsize=11)
    ax_d.set_ylabel("Degradation rate D (min⁻¹)", fontsize=11)
    ax_d.set_title(
        "Inferred D per AP bin\n"
        "(dashed = edge avg, dotted = centre min)",
        fontsize=10,
    )
    ax_d.legend(fontsize=9, loc="upper right")
    ax_d.grid(True, axis="y", alpha=0.3, zorder=0)
    # Expand top to clear the half-life text annotations AND the legend
    ax_d.set_ylim(bottom=0, top=y_top * 1.25)

    # ── Save figure ───────────────────────────────────────────────────────────
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"\nCounterfactual figure saved → {args.output}")

    # ── Save figure data to CSVs ──────────────────────────────────────────────
    out_path   = Path(args.output)
    # Derive figure_data/ as a sibling of figures/
    fig_data_dir = out_path.parent.parent / "figure_data"
    fig_data_dir.mkdir(parents=True, exist_ok=True)

    # Build a short tag from the stripe and embryo parts of the path, e.g. stripe4_e8_9
    path_parts  = out_path.parts          # (..., 'stripe4', 'e8_9um', 'figures', 'file.pdf')
    stripe_tag  = path_parts[-4]          # e.g. stripe4
    embryo_tag  = path_parts[-3].removesuffix("um")  # e.g. e8_9
    suffix      = f"{stripe_tag}_{embryo_tag}"

    # 1. Heatmap grids (n_ap × n_dv)
    heatmap_exports = [
        (obs_grid,    "observed"),
        (fitted_grid, "fitted"),
        (fast_grid,   "fast_everywhere"),
        (slow_grid,   "slow_everywhere"),
    ]
    for grid, tag in heatmap_exports:
        df_hm = pd.DataFrame(
            grid,
            index   = ap_labels,
            columns = dv_labels,
        )
        df_hm.index.name = "ap_bin"
        csv_path = fig_data_dir / f"heatmap_{tag}_{suffix}.csv"
        df_hm.to_csv(csv_path)
        print(f"  Data saved → {csv_path}")

    # 2. D-profile (bar chart data)
    df_d = pd.DataFrame({
        "ap_bin":        ap_labels,
        "D_median":      D_medians,
        "half_life_min": np.log(2) / D_medians,
    })
    csv_path = fig_data_dir / f"dprofile_{suffix}.csv"
    df_d.to_csv(csv_path, index=False)
    print(f"  Data saved → {csv_path}")

    # 3. AP-mean line plot data
    df_ap = pd.DataFrame({
        "ap_bin":        ap_labels,
        "observed_mean":       obs_grid.mean(axis=1),
        "fitted_mean":         fitted_grid.mean(axis=1),
        "fast_everywhere_mean": fast_grid.mean(axis=1),
        "slow_everywhere_mean": slow_grid.mean(axis=1),
    })
    csv_path = fig_data_dir / f"ap_means_{suffix}.csv"
    df_ap.to_csv(csv_path, index=False)
    print(f"  Data saved → {csv_path}")


if __name__ == "__main__":
    main()
