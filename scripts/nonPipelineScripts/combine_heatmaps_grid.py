import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import glob
import re
import sys

def main():
    stripe = sys.argv[1] if len(sys.argv) > 1 else 'stripe2'
    # Define embryos (columns)
    if stripe == 'stripe2':
        embryos = ['e7um','e8_9um','e9_10um']
    else: #if stripe == 'stripe3':
        embryos = ['e7um_4', 'e8_9um', 'e9_10um']
    
    # Identify results directories (rows)
    # Using the list from the user's workspace but filtering for results_XXX
    results_dirs = []
    for d in os.listdir('.'):
        if os.path.isdir(d) and re.match(r'^results_\d+$', d):
            results_dirs.append(d)
    
    # Sort numerically based on the XXX part
    results_dirs.sort(key=lambda x: int(re.search(r'\d+', x).group()))
    
    if not results_dirs:
        print("No results_XXX directories found.")
        return

    n_rows = len(results_dirs)
    n_cols = len(embryos)
    
    print(f"Creating grid of {n_rows} rows and {n_cols} columns...")
    
    # Set up the figure
    # We'll adjust the figure size based on the number of rows/cols
    # Assuming each heatmap is roughly rectangular
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), 
                             gridspec_kw={'wspace': 0.05, 'hspace': 0.05})
    
    # Ensure axes is 2D even if 1 row/col
    if n_rows == 1:
        axes = [axes]
    if n_cols == 1:
        axes = [[ax] for ax in axes]
        
    for i, res_dir in enumerate(results_dirs):
        for j, embryo in enumerate(embryos):
            ax = axes[i][j]
            img_path = os.path.join(res_dir, stripe, embryo, 'figures', 'halflife_heatmap.png')
            
            if os.path.exists(img_path):
                img = mpimg.imread(img_path)
                ax.imshow(img)
            else:
                ax.text(0.5, 0.5, 'Missing', ha='center', va='center')
            
            ax.set_xticks([])
            ax.set_yticks([])
            
            # Label columns on the first row
            if i == 0:
                ax.set_title(embryo, fontsize=12)
            
            # Label rows on the first column
            if j == 0:
                ax.set_ylabel(res_dir, rotation=0, labelpad=40, verticalalignment='center', fontsize=10)

    plt.tight_layout()
    output_path = f'results/combined_halflife_heatmaps_grid_{stripe}.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    print(f"Combined image saved to {output_path}")

if __name__ == "__main__":
    main()
