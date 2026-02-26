import os
import csv
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import re

STRIPES = {
    'stripe2': ['e7um', 'e8_9um', 'e9_10um'],
    'stripe3': ['e7um_4', 'e8_9um', 'e9_10um'],
    'stripe4': ['e7um_4', 'e8_9um', 'e9_10um'],
}

FIGURES = [
    'loo_compare.png',
    'pareto_k_by_ap_bin_summary.png',
]


def has_loo_warning(img_dir):
    """Return True if loo_comparison.csv in img_dir has any True in the 'warning' column."""
    csv_path = os.path.join(img_dir, 'loo_comparison.csv')
    if not os.path.exists(csv_path):
        return False
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('warning', '').strip().lower() == 'true':
                return True
    return False


def find_results_dirs():
    results_dirs = [
        d for d in os.listdir('.')
        if os.path.isdir(d) and re.match(r'^results_\d+$', d)
    ]
    results_dirs.sort(key=lambda x: int(re.search(r'\d+', x).group()))
    return results_dirs


def make_grid(results_dirs, stripe, embryos, figure_filename):
    n_rows = len(results_dirs)
    n_cols = len(embryos)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 4, n_rows * 3),
        gridspec_kw={'wspace': 0.05, 'hspace': 0.05},
    )

    # Normalise axes to always be a 2-D list
    if n_rows == 1:
        axes = [axes]
    if n_cols == 1:
        axes = [[ax] for ax in axes]

    for i, res_dir in enumerate(results_dirs):
        for j, embryo in enumerate(embryos):
            ax = axes[i][j]
            img_path = os.path.join(
                res_dir, 'comparison', stripe, embryo, figure_filename
            )

            if os.path.exists(img_path):
                img = mpimg.imread(img_path)
                ax.imshow(img)
            else:
                ax.text(0.5, 0.5, 'Missing', ha='center', va='center',
                        transform=ax.transAxes, fontsize=10)

            ax.set_xticks([])
            ax.set_yticks([])

            # Red border if this is the pareto_k figure and loo_comparison.csv
            # has any warning == True in the same directory.
            img_dir = os.path.join(res_dir, 'comparison', stripe, embryo)
            if (figure_filename == 'pareto_k_by_ap_bin_summary.png'
                    and has_loo_warning(img_dir)):
                for spine in ax.spines.values():
                    spine.set_edgecolor('red')
                    spine.set_linewidth(3)
                    spine.set_visible(True)

            if i == 0:
                ax.set_title(embryo, fontsize=12)

            if j == 0:
                ax.set_ylabel(res_dir, rotation=0, labelpad=50,
                              verticalalignment='center', fontsize=10)

    figure_stem = os.path.splitext(figure_filename)[0]
    output_path = os.path.join('results', f'combined_{figure_stem}_grid_{stripe}.png')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    results_dirs = find_results_dirs()
    if not results_dirs:
        print("No results_XXX directories found.")
        return

    print(f"Found {len(results_dirs)} results directories: {results_dirs[0]} … {results_dirs[-1]}")

    for stripe, embryos in STRIPES.items():
        for figure_filename in FIGURES:
            print(f"Building grid: {stripe} / {figure_filename}")
            make_grid(results_dirs, stripe, embryos, figure_filename)


if __name__ == "__main__":
    main()
