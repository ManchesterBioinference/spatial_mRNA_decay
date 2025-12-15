# Scripts

This directory contains all scripts for the research project, including both pipeline analysis scripts and utility scripts.

> **Note**: Research Assistant tools (like transcription) are in the `tools/` directory, not here.
> See `tools/README.md` for transcription and other RA utilities.

## Purpose

Use this directory for:

- **Pipeline scripts**: Core analysis scripts tracked by DVC (numbered for order)
- **Data acquisition**: Scripts to download or fetch raw data
- **Environment setup**: Scripts to configure the development environment
- **Utilities**: Helper scripts for common tasks
- **Automation**: CI/CD scripts, batch processing, etc.

## Organization

```
scripts/
├── 01_preprocess.py       # Pipeline: Data cleaning
├── 02_features.py         # Pipeline: Feature engineering
├── 03_train.py            # Pipeline: Model training
├── 04_evaluate.py         # Pipeline: Evaluation
├── 05_figures.py          # Pipeline: Figure generation
├── data/                  # Data acquisition scripts
│   └── download_data.py
├── utils/                 # Utility functions
│   └── helpers.py
└── setup/                 # Environment setup
    └── install_deps.sh
```

## Pipeline Scripts (DVC Tracked)

Numbered scripts (e.g., `01_preprocess.py`) are part of the reproducible DVC pipeline:

- Must be referenced in `dvc.yaml`
- Should use parameters from `params.yaml`
- Outputs go to `results/` directory

## Documentation Standards

Every script should have:

1. **Module docstring**: Purpose, inputs, outputs
2. **Function docstrings**: Args, returns, examples
3. **Inline comments**: Explain non-obvious logic
4. **No hardcoded paths**: Use params.yaml or CLI arguments

### Example Pipeline Script

```python
"""
Preprocessing module for [project name].

This script handles data cleaning and normalization.

Inputs:
    - data/raw/input.csv
    
Outputs:
    - data/processed/cleaned.csv

Usage:
    python scripts/01_preprocess.py
    
Author: [Name]
Date: [Date]
"""

import argparse
import pandas as pd
# ...
```

## Stripe Range Identification

The `identify_stripe_ranges.py` script automatically detects the AP coordinate ranges for all 7 eve stripes from expression data. This must be run before the analysis pipeline to populate `config.yaml` with stripe boundaries.

### Running Stripe Identification

```bash
# Generate stripe ranges and diagnostic plot
python scripts/identify_stripe_ranges.py \
    --input data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv \
    --config config.yaml \
    --output results/figures/stripe_identification.png
```

**What it does:**
1. Bins nuclei by AP position (0.01 bin size)
2. Computes median fluorescence per bin (early expression period)
3. Applies Savitzky-Golay smoothing
4. Detects peaks (stripe centers) and troughs (inter-stripes)
5. Calculates symmetric ranges around each peak
6. Writes results to `config.yaml` in `stripe_ranges` section
7. Creates diagnostic plot showing detected stripes

**Output in config.yaml:**
```yaml
stripe_ranges:
  stripe1:
    center: 0.1234
    min: 0.0987
    max: 0.1481
  stripe2:
    center: 0.3900
    min: 0.3300
    max: 0.4500
  # ... through stripe7
```

### How Other Scripts Use Stripe Ranges

**In Snakefile:**
```python
# Stripe ranges loaded from config
STRIPE_RANGES = config.get("stripe_ranges", {})
AP_MIN_LOOKUP = {stripe: STRIPE_RANGES[stripe]["min"] for stripe in STRIPES}
AP_MAX_LOOKUP = {stripe: STRIPE_RANGES[stripe]["max"] for stripe in STRIPES}

# Passed to preprocessing script
params:
    ap_min=lambda wildcards: AP_MIN_LOOKUP[wildcards.stripe],
    ap_max=lambda wildcards: AP_MAX_LOOKUP[wildcards.stripe]
```

**In preprocessing script:**
```bash
python scripts/01_preprocess_eve_data.py \
    --stripe stripe2 \
    --ap-min 0.33 \
    --ap-max 0.45 \
    # ... other arguments
```

### When to Regenerate

Re-run `identify_stripe_ranges.py` if:
- Using a different dataset
- Adjusting peak detection parameters (prominence threshold)
- Fine-tuning stripe boundaries for specific analysis

### Relationship to Pipeline Output Files

Before this update, files were named with encoded AP ranges:
- `transcription_traces_033045.csv` (meant AP [0.33, 0.45])

Now files use stripe names:
- `transcription_traces_stripe2.csv`
- `results/stripe2/chains/degradation_chain.csv`

This makes the analysis more intuitive and extensible to all 7 stripes.

## Notes

- Pipeline scripts use numbered prefixes to indicate order
- Utility scripts in subdirectories (data/, utils/, setup/) are not tracked in main pipeline
- Keep exploratory notebooks in a separate `notebooks/` directory if needed
- Always run `identify_stripe_ranges.py` before starting the main pipeline
````
