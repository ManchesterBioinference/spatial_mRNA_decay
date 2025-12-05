---
name: scientific-visualization
description: Create publication-quality scientific figures with matplotlib, seaborn, and plotly. Includes multi-panel layouts, error bars, significance markers, colorblind-safe palettes, and journal-specific export (PDF/EPS/TIFF). Use when creating figures for manuscripts, presentations, or any research visualization.
---

# Scientific Visualization

> Create publication-ready figures that are clear, accurate, and accessible.

## When to Use

- Creating plots or visualizations for manuscripts
- Preparing figures for journal submission
- Ensuring figures are colorblind-friendly and accessible
- Making multi-panel figures with consistent styling
- Exporting figures at correct resolution and format
- Creating figures in `manuscript/figures/figN/`
- During the ANALYSIS or WRITING phases

## Quick Start: Publication Figure

```python
import matplotlib.pyplot as plt
import numpy as np

# Publication settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# Create figure (single column = 3.5 inches / 89mm)
fig, ax = plt.subplots(figsize=(3.5, 2.5))

# Plot with colorblind-safe colors
x = np.linspace(0, 10, 100)
ax.plot(x, np.sin(x), color='#0072B2', label='sin(x)')
ax.plot(x, np.cos(x), color='#D55E00', label='cos(x)')

# Proper labeling with units
ax.set_xlabel('Time (s)')
ax.set_ylabel('Amplitude (mV)')
ax.legend(frameon=False)

# Clean style
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Save
fig.savefig('manuscript/figures/fig1/figure1.pdf')
fig.savefig('manuscript/figures/fig1/figure1.png', dpi=300)
plt.close()
```

---

## Colorblind-Safe Palettes

### Okabe-Ito Palette (Recommended)

```python
# Best for categorical data - distinguishable by all color vision types
OKABE_ITO = [
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermillion
    '#CC79A7',  # Reddish Purple
    '#000000',  # Black
]

plt.rcParams['axes.prop_cycle'] = plt.cycler(color=OKABE_ITO)
```

### For Continuous Data

```python
# Use perceptually uniform colormaps
cmap = 'viridis'     # Default, good for most cases
cmap = 'plasma'      # High contrast
cmap = 'cividis'     # Colorblind-safe diverging

# AVOID: 'jet', 'rainbow', 'hsv' - these are NOT colorblind-safe
```

### For Diverging Data

```python
# Colorblind-safe diverging colormaps
cmap = 'RdBu_r'     # Red-Blue (reversed so blue = low)
cmap = 'PuOr'       # Purple-Orange
cmap = 'BrBG'       # Brown-Blue-Green

# Always center diverging maps at meaningful zero
sns.heatmap(data, cmap='RdBu_r', center=0)
```

---

## Journal Figure Sizes

### Common Widths (in inches)

| Journal | Single Column | Double Column |
|---------|---------------|---------------|
| Nature | 3.5" (89mm) | 7.2" (183mm) |
| Science | 2.2" (55mm) | 6.9" (175mm) |
| Cell | 3.3" (85mm) | 7.0" (178mm) |
| PLOS | 5.2" | 7.5" |
| Default | 3.5" | 7.0" |

### Resolution Requirements

| Content Type | Minimum DPI |
|--------------|-------------|
| Line art (graphs) | 600-1200 (or vector) |
| Photographs | 300-600 |
| Combination | 600 |

### File Formats

| Format | Use For | Notes |
|--------|---------|-------|
| PDF | Primary | Vector, preserves quality |
| EPS | Journals | Vector, legacy support |
| TIFF | Images | Lossless raster |
| PNG | Web/screen | Lossless, smaller than TIFF |
| JPEG | NEVER | Lossy compression artifacts |

---

## Multi-Panel Figures

### Using GridSpec

```python
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from string import ascii_uppercase

# Create figure with custom layout
fig = plt.figure(figsize=(7, 5))
gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.4)

# Create panels
ax_a = fig.add_subplot(gs[0, 0])      # Top left
ax_b = fig.add_subplot(gs[0, 1:])     # Top right (spans 2 columns)
ax_c = fig.add_subplot(gs[1, 0])      # Bottom left
ax_d = fig.add_subplot(gs[1, 1])      # Bottom middle
ax_e = fig.add_subplot(gs[1, 2])      # Bottom right

# Add panel labels
axes = [ax_a, ax_b, ax_c, ax_d, ax_e]
for i, ax in enumerate(axes):
    ax.text(-0.15, 1.05, ascii_uppercase[i], transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top')
    
    # Example: clean up each axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.savefig('manuscript/figures/fig1/figure1.pdf')
```

### Using Subplots

```python
fig, axes = plt.subplots(2, 2, figsize=(7, 6))
axes = axes.flatten()

for i, ax in enumerate(axes):
    ax.text(-0.1, 1.05, ascii_uppercase[i], transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top')

plt.tight_layout()
```

---

## Common Plot Types

### Bar Plot with Error Bars

```python
import matplotlib.pyplot as plt
import numpy as np

groups = ['Control', 'Treatment A', 'Treatment B']
means = [10.2, 15.8, 18.3]
sems = [1.2, 1.5, 1.8]  # Standard error of mean
n = [20, 20, 20]

fig, ax = plt.subplots(figsize=(3.5, 3))

bars = ax.bar(groups, means, yerr=sems, capsize=4, 
              color=['#0072B2', '#E69F00', '#D55E00'],
              edgecolor='black', linewidth=0.5)

ax.set_ylabel('Response (AU)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add individual data points (transparency)
for i, (group, n_i) in enumerate(zip(groups, n)):
    x_jitter = np.random.normal(i, 0.05, n_i)
    y_data = np.random.normal(means[i], sems[i]*np.sqrt(n_i), n_i)
    ax.scatter(x_jitter, y_data, alpha=0.4, s=15, c='black')

plt.tight_layout()
```

### Line Plot with Confidence Bands

```python
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 50)
y = np.sin(x) + np.random.normal(0, 0.2, 50)
y_smooth = np.sin(x)
ci = 0.3  # Confidence interval width

fig, ax = plt.subplots(figsize=(4, 3))

# Confidence band
ax.fill_between(x, y_smooth - ci, y_smooth + ci, alpha=0.3, color='#0072B2')
# Line
ax.plot(x, y_smooth, color='#0072B2', linewidth=2)
# Data points
ax.scatter(x, y, s=20, alpha=0.5, color='#0072B2')

ax.set_xlabel('Time (s)')
ax.set_ylabel('Signal (AU)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
```

### Box Plot with Significance Markers

```python
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(4, 4))

# Box plot
sns.boxplot(data=df, x='group', y='value', palette='colorblind', ax=ax)

# Add individual points
sns.stripplot(data=df, x='group', y='value', color='black', 
              alpha=0.3, size=4, ax=ax)

# Add significance bar
y_max = df['value'].max()
ax.plot([0, 1], [y_max*1.1, y_max*1.1], 'k-', linewidth=1)
ax.text(0.5, y_max*1.12, '***', ha='center', fontsize=10)

ax.set_ylabel('Response (AU)')
sns.despine()
```

### Heatmap with Proper Colorbar

```python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Correlation matrix example
fig, ax = plt.subplots(figsize=(5, 4))

# Mask upper triangle
mask = np.triu(np.ones_like(corr_matrix, dtype=bool))

sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f',
            cmap='RdBu_r', center=0, square=True,
            linewidths=0.5, cbar_kws={'shrink': 0.8},
            vmin=-1, vmax=1, ax=ax)

ax.set_title('Correlation Matrix')
plt.tight_layout()
```

---

## Seaborn for Statistical Plots

```python
import seaborn as sns

# Apply publication style
sns.set_theme(style='ticks', context='paper', font_scale=1.1)
sns.set_palette('colorblind')

# Violin plot with split
fig, ax = plt.subplots(figsize=(4, 3))
sns.violinplot(data=df, x='timepoint', y='expression',
               hue='treatment', split=True, inner='quartile', ax=ax)
ax.set_ylabel('Gene Expression (AU)')
sns.despine()

# Regression plot with CI
fig, ax = plt.subplots(figsize=(4, 3))
sns.regplot(data=df, x='predictor', y='outcome', 
            ci=95, scatter_kws={'alpha': 0.5}, ax=ax)
sns.despine()
```

---

## Figure Export Workflow

### For RA Projects

```python
def save_figure(fig, fig_num, name):
    """Save figure in multiple formats following RA conventions."""
    import os
    
    # Create figure directory
    fig_dir = f'manuscript/figures/fig{fig_num}'
    os.makedirs(fig_dir, exist_ok=True)
    
    # Save in multiple formats
    fig.savefig(f'{fig_dir}/{name}.pdf', format='pdf', bbox_inches='tight')
    fig.savefig(f'{fig_dir}/{name}.png', format='png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{fig_dir}/{name}.tiff', format='tiff', dpi=300, bbox_inches='tight')
    
    print(f"Saved: {fig_dir}/{name}.{{pdf,png,tiff}}")
    
# Usage
save_figure(fig, 1, 'treatment_comparison')
```

### Caption Template

Create `manuscript/figures/fig1/caption.md`:

```markdown
**Figure 1. Treatment comparison across experimental groups.**

(A) Response levels in control (n=20), Treatment A (n=20), and Treatment B (n=20) 
groups. Data shown as mean ± SEM with individual data points. ***p < 0.001, 
one-way ANOVA with Tukey's post-hoc test.

(B) Time course of response following treatment. Shaded regions indicate 95% 
confidence intervals. Treatment A (orange) shows significantly faster response 
than control (blue), p < 0.01.
```

---

## Accessibility Checklist

- [ ] **Colorblind-safe palette** (Okabe-Ito or viridis)
- [ ] **Works in grayscale** (add patterns/markers if needed)
- [ ] **Font size ≥ 6pt** at final print size
- [ ] **All axes labeled** with units
- [ ] **Error bars defined** in caption (SEM, SD, or CI)
- [ ] **No 3D effects** or chartjunk
- [ ] **Legend readable** and positioned appropriately
- [ ] **Resolution ≥ 300 DPI** for raster, or vector format

---

## Integration with RA Workflow

1. Create figure directory: `manuscript/figures/figN/`
2. Save figure files: `figure.pdf`, `figure.png`
3. Create caption: `caption.md`
4. Reference in `manuscript/results.md`
5. Log activity in `.research/logs/activity.md`
