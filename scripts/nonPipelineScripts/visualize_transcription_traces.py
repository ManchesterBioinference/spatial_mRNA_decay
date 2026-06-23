import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--transcription', required=True, help='Path to transcription CSV (e.g., data/processed_transcription_data/transcription_traces_stripe2.csv)')
    parser.add_argument('--out_dir', default='results/figures/intermediate/transcription', help='Output directory for saving plots')
    args = parser.parse_args()

    # Load data - assuming headers exist based on your description
    # If the file has no headers, we skip the first two columns (apBin, yBin)
    df = pd.read_csv(args.transcription)
    
    # Extract time column names (numeric headers from 0 to 1200) and sort them numerically
    time_cols = sorted([c for c in df.columns if c not in ['apBin', 'yBin']], key=lambda x: float(x))
    time_points = np.array([float(c) for c in time_cols])
    
    # Get unique bins dynamically
    unique_ap = sorted(df['apBin'].unique())
    unique_dv = sorted(df['yBin'].unique())
    num_ap = len(unique_ap)
    num_dv = len(unique_dv)
    
    # Colors for AP groups: blue, orange, green, red, purple (cycle if more)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    # Line styles for DV bins: solid, dashed, dashdot, dotted, custom (cycle if more)
    line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
    
    # Create subplots for each (AP, DV) bin in a grid (rows=DV, columns=AP for spatial intuition)
    fig, axes = plt.subplots(num_dv, num_ap, figsize=(num_ap * 2, num_dv * 2), sharex=True, sharey=True)
    if num_ap == 1 and num_dv == 1:
        axes = np.array([[axes]])  # Handle single subplot case
    elif num_ap == 1:
        axes = axes.reshape(-1, 1)
    elif num_dv == 1:
        axes = axes.reshape(1, -1)
    
    for j, dv in enumerate(unique_dv):
        for i, ap in enumerate(unique_ap):
            ax = axes[j, i]
            # Filter for the specific bin
            bin_data = df[(df['apBin'] == ap) & (df['yBin'] == dv)]
            
            if not bin_data.empty:
                row = bin_data.iloc[0]  # Assume one row per bin
                ax.plot(time_points, row[time_cols].values.astype(float), 
                        color=colors[i % len(colors)], linestyle=line_styles[j % len(line_styles)],
                        alpha=0.8)
            
            ax.set_title(f'AP {ap}, DV {dv}', fontsize=10)
            ax.grid(True, alpha=0.3)
    
    # Set axis labels
    for j in range(num_dv):
        axes[j, 0].set_ylabel('Intensity (F)', fontsize=8)
    for i in range(num_ap):
        axes[-1, i].set_xlabel('Time (seconds)', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, 'transcription_traces_check.pdf'), dpi=300)
    print(f"Trace visualization saved to {os.path.join(args.out_dir, 'transcription_traces_check.pdf')}")
    
    # -------------------------
    # SINGLE PLOT: All traces in one panel
    # -------------------------
    fig_single, ax_single = plt.subplots(1, 1, figsize=(12, 8))
    
    for j, dv in enumerate(unique_dv):
        for i, ap in enumerate(unique_ap):
            # Filter for the specific bin
            bin_data = df[(df['apBin'] == ap) & (df['yBin'] == dv)]
            
            if not bin_data.empty:
                row = bin_data.iloc[0]  # Assume one row per bin
                ax_single.plot(time_points, row[time_cols].values.astype(float), 
                               color=colors[i % len(colors)], linestyle=line_styles[j % len(line_styles)],
                               alpha=0.8, label=f'AP {ap}, DV {dv}')
    
    ax_single.set_title('All Transcription Traces')
    ax_single.set_xlabel('Time (seconds)')
    ax_single.set_ylabel('Intensity (F)')
    ax_single.grid(True, alpha=0.3)
    ax_single.legend(loc='upper left', fontsize='small', ncol=3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, 'transcription_traces_all.pdf'), dpi=300)
    print(f"All traces in one plot saved to {os.path.join(args.out_dir, 'transcription_traces_all.pdf')}")
    
    # -------------------------
    # GROUPED BY AP: Subplots for each AP bin, with all DV lines
    # -------------------------
    fig_grouped, axes_grouped = plt.subplots(num_ap, 1, figsize=(10, num_ap * 3), sharex=True, sharey=True)
    if num_ap == 1:
        axes_grouped = [axes_grouped]
    
    for i, ap in enumerate(unique_ap):
        ax = axes_grouped[i]
        for j, dv in enumerate(unique_dv):
            # Filter for the specific bin
            bin_data = df[(df['apBin'] == ap) & (df['yBin'] == dv)]
            
            if not bin_data.empty:
                row = bin_data.iloc[0]  # Assume one row per bin
                ax.plot(time_points, row[time_cols].values.astype(float), 
                        color=colors[i % len(colors)], linestyle=line_styles[j % len(line_styles)],
                        alpha=0.8, label=f'DV {dv}')
        
        ax.set_title(f'Transcription Traces: AP {ap}')
        ax.set_ylabel('Intensity (F)')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize='small', ncol=2)
    
    axes_grouped[-1].set_xlabel('Time (seconds)')
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, 'transcription_traces_grouped_ap.pdf'), dpi=300)
    print(f"Traces grouped by AP saved to {os.path.join(args.out_dir, 'transcription_traces_grouped_ap.pdf')}")
    
    # -------------------------
    # CUMULATIVE LINE PLOTS (additive across time)
    # For each bin, show cumulative sum from t0 to current time
    # -------------------------
    fig_cum, axes_cum = plt.subplots(num_dv, num_ap, figsize=(num_ap * 2, num_dv * 2), sharex=True, sharey=True)
    if num_ap == 1 and num_dv == 1:
        axes_cum = np.array([[axes_cum]])
    elif num_ap == 1:
        axes_cum = axes_cum.reshape(-1, 1)
    elif num_dv == 1:
        axes_cum = axes_cum.reshape(1, -1)
    
    for j, dv in enumerate(unique_dv):
        for i, ap in enumerate(unique_ap):
            ax = axes_cum[j, i]
            # Filter for the specific bin
            bin_data = df[(df['apBin'] == ap) & (df['yBin'] == dv)]
            
            if not bin_data.empty:
                row = bin_data.iloc[0]
                # Compute cumulative sum across time for this trace
                intensity_values = row[time_cols].values.astype(float)
                cumulative_values = np.cumsum(intensity_values)
                ax.plot(time_points, cumulative_values, 
                        color=colors[i % len(colors)], linestyle=line_styles[j % len(line_styles)],
                        alpha=0.8)
            
            ax.set_title(f'Cumulative: AP {ap}, DV {dv}', fontsize=10)
            ax.grid(True, alpha=0.3)
    
    # Set axis labels
    for j in range(num_dv):
        axes_cum[j, 0].set_ylabel('Cumulative Intensity (a.u.)', fontsize=8)
    for i in range(num_ap):
        axes_cum[-1, i].set_xlabel('Time (seconds)', fontsize=8)
    
    plt.tight_layout()
    out_cum = os.path.join(args.out_dir, 'transcription_traces_cumulative.pdf')
    plt.savefig(out_cum, dpi=300)
    print(f"Cumulative line traces saved to {out_cum}")
    
    # -------------------------
    # SINGLE PLOT: All cumulative traces in one panel
    # -------------------------
    fig_cum_single, ax_cum_single = plt.subplots(1, 1, figsize=(12, 8))
    
    for j, dv in enumerate(unique_dv):
        for i, ap in enumerate(unique_ap):
            # Filter for the specific bin
            bin_data = df[(df['apBin'] == ap) & (df['yBin'] == dv)]
            
            if not bin_data.empty:
                row = bin_data.iloc[0]
                # Compute cumulative sum across time for this trace
                intensity_values = row[time_cols].values.astype(float)
                cumulative_values = np.cumsum(intensity_values)
                ax_cum_single.plot(time_points, cumulative_values, 
                                   color=colors[i % len(colors)], linestyle=line_styles[j % len(line_styles)],
                                   alpha=0.8, label=f'AP {ap}, DV {dv}')
    
    ax_cum_single.set_title('All Cumulative Transcription Traces')
    ax_cum_single.set_xlabel('Time (seconds)')
    ax_cum_single.set_ylabel('Cumulative Intensity (a.u.)')
    ax_cum_single.grid(True, alpha=0.3)
    ax_cum_single.legend(loc='upper left', fontsize='small', ncol=3)
    
    plt.tight_layout()
    out_cum_all = os.path.join(args.out_dir, 'transcription_traces_cumulative_all.pdf')
    plt.savefig(out_cum_all, dpi=300)
    print(f"All cumulative traces in one plot saved to {out_cum_all}")
    
    # -------------------------
    # GROUPED BY AP: Cumulative subplots for each AP bin, with all DV lines
    # -------------------------
    fig_cum_grouped, axes_cum_grouped = plt.subplots(num_ap, 1, figsize=(10, num_ap * 3), sharex=True, sharey=True)
    if num_ap == 1:
        axes_cum_grouped = [axes_cum_grouped]
    
    for i, ap in enumerate(unique_ap):
        ax = axes_cum_grouped[i]
        for j, dv in enumerate(unique_dv):
            # Filter for the specific bin
            bin_data = df[(df['apBin'] == ap) & (df['yBin'] == dv)]
            
            if not bin_data.empty:
                row = bin_data.iloc[0]
                # Compute cumulative sum across time for this trace
                intensity_values = row[time_cols].values.astype(float)
                cumulative_values = np.cumsum(intensity_values)
                ax.plot(time_points, cumulative_values, 
                        color=colors[i % len(colors)], linestyle=line_styles[j % len(line_styles)],
                        alpha=0.8, label=f'DV {dv}')
        
        ax.set_title(f'Cumulative Transcription: AP {ap}')
        ax.set_ylabel('Cumulative Intensity (a.u.)')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize='small', ncol=2)
    
    axes_cum_grouped[-1].set_xlabel('Time (seconds)')
    plt.tight_layout()
    out_cum_grouped = os.path.join(args.out_dir, 'transcription_traces_cumulative_grouped_ap.pdf')
    plt.savefig(out_cum_grouped, dpi=300)
    print(f"Cumulative traces grouped by AP saved to {out_cum_grouped}")
    
#    # -------------------------
#    # ADDITIVE IMAGE (cumulative across time)
#    # For each AP bin, build a matrix (rows = yBin sorted, cols = time) then cumulative-sum across time
#    # and display as an image so intensity is additive as time progresses.
#    # -------------------------
#    fig2, axes2 = plt.subplots(5, 1, figsize=(10, 15), sharex=True)
#    vmin = 0.0
#    # Compute a robust vmax across all bins to use a consistent color scale
#    global_max = 0.0
#    for i in range(1, 6):
#        bin_data = df[df['apBin'] == float(i)]
#        if bin_data.empty:
#            continue
#        matrix = bin_data[time_cols].astype(float).values
#        cum_matrix = np.cumsum(matrix, axis=1)
#        global_max = max(global_max, np.nanmax(cum_matrix))
#    if global_max == 0:
#        global_max = 1.0
#
#    for i in range(1, 6):
#        ax = axes2[i-1]
#        bin_data = df[df['apBin'] == float(i)]
#        if bin_data.empty:
#            ax.set_visible(False)
#            continue
#        # sort by yBin so rows map to DV position
#        bin_data = bin_data.sort_values('yBin')
#        y_bins = bin_data['yBin'].values
#        matrix = bin_data[time_cols].astype(float).values
#        cum_matrix = np.cumsum(matrix, axis=1)
#        
#        # Display cumulative (additive) image: y axis = DV bins, x axis = time
#        im = ax.imshow(cum_matrix, aspect='auto', cmap='magma',
#                       vmin=vmin, vmax=global_max,
#                       extent=[time_points[0], time_points[-1], y_bins.min()-0.5, y_bins.max()+0.5],
#                       origin='lower', interpolation='nearest')
#        ax.set_ylabel(f'AP {i} DV Bin')
#        ax.set_title(f'Additive (cumulative) Intensity — AP Bin {i}')
#        ax.grid(False)
#
#    axes2[-1].set_xlabel('Time (seconds)')
#    # shared colorbar on the right
#    cbar = fig2.colorbar(im, ax=axes2, orientation='vertical', fraction=0.02, pad=0.02)
#    cbar.set_label('Cumulative Intensity (a.u.)')
#    plt.tight_layout()
#    out_add = os.path.join(args.out_dir, 'transcription_traces_additive.pdf')
#    plt.savefig(out_add, dpi=300)
#    print(f"Additive cumulative image saved to {out_add}")

if __name__ == "__main__":
    main()