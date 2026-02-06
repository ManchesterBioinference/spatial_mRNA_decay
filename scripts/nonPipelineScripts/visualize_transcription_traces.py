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
    
    # Create subplots for each AP Bin (line traces)
    fig, axes = plt.subplots(5, 1, figsize=(10, 15), sharex=True, sharey=True)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i in range(1, 6):
        ax = axes[i-1]
        # Filter for the specific AP bin
        bin_data = df[df['apBin'] == float(i)]
        
        for idx, row in bin_data.iterrows():
            ax.plot(time_points, row[time_cols].values.astype(float), 
                    alpha=0.7, label=f'DV Bin {row["yBin"]}')
        
        ax.set_title(f'Transcription Traces: AP Bin {i}')
        ax.set_ylabel('Intensity (F)')
        ax.grid(True, alpha=0.3)
        if i == 1:
            ax.legend(loc='upper right', fontsize='small', ncol=2)

    axes[-1].set_xlabel('Time (seconds)')
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, 'transcription_traces_check.png'), dpi=300)
    print(f"Trace visualization saved to {os.path.join(args.out_dir, 'transcription_traces_check.png')}")
    
    # -------------------------
    # CUMULATIVE LINE PLOTS (additive across time)
    # For each trace, show cumulative sum from t0 to current time
    # -------------------------
    fig_cum, axes_cum = plt.subplots(5, 1, figsize=(10, 15), sharex=True, sharey=True)
    
    for i in range(1, 6):
        ax = axes_cum[i-1]
        # Filter for the specific AP bin
        bin_data = df[df['apBin'] == float(i)]
        
        for idx, row in bin_data.iterrows():
            # Compute cumulative sum across time for this trace
            intensity_values = row[time_cols].values.astype(float)
            cumulative_values = np.cumsum(intensity_values)
            ax.plot(time_points, cumulative_values, 
                    alpha=0.7, label=f'DV Bin {row["yBin"]}')
        
        ax.set_title(f'Cumulative Transcription: AP Bin {i}')
        ax.set_ylabel('Cumulative Intensity (a.u.)')
        ax.grid(True, alpha=0.3)
        if i == 1:
            ax.legend(loc='upper left', fontsize='small', ncol=2)

    axes_cum[-1].set_xlabel('Time (seconds)')
    plt.tight_layout()
    out_cum = os.path.join(args.out_dir, 'transcription_traces_cumulative.png')
    plt.savefig(out_cum, dpi=300)
    print(f"Cumulative line traces saved to {out_cum}")
    
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
#    out_add = os.path.join(args.out_dir, 'transcription_traces_additive.png')
#    plt.savefig(out_add, dpi=300)
#    print(f"Additive cumulative image saved to {out_add}")

if __name__ == "__main__":
    main()