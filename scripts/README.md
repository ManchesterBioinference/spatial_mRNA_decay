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
├── 00_identify_stripe_ranges.py   # Detect AP ranges for eve stripes
├── 01_preprocess_eve_data.py      # Pipeline: Transcription data preprocessing
├── 02_infer_degradation_rates.py  # Pipeline: Bayesian inference
├── 03_visualize_results.py        # Pipeline: Results visualization
├── 04_aggregate_embryo_mrna.py    # Pipeline: mRNA spatial binning
├── 05_validate_nuclei_density.py  # Pipeline: QC validation
├── 06_compute_validation_thresholds.py  # Compute QC thresholds from data
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
    --output results/figures/intermediate/transcription/stripe_identification.png
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

## Bayesian Inference of Degradation Rates (Script 02)

The `02_infer_degradation_rates.py` script uses PyMC to fit a Bayesian ODE model inferring spatially-varying mRNA degradation rates across anterior-posterior bins.

### Model Overview

**ODE Model:**
```
dm/dt = γ*F(t) - D*m
```
Where:
- `m`: mRNA concentration
- `F(t)`: Time-varying transcription input
- `D`: Degradation rate (varies by spatial bin)
- `γ`: Transcription scaling factor

**Bayesian Priors (for reproducibility):**
- `D` (degradation rates): `LogNormal(mu=-2, sigma=1)`
  - Targets 1-10 min half-lives (ln(2)/D), prevents numerical overflow
- `γ` (transcription scaling): `HalfNormal(sigma=100)`
  - Conditioned for normalized data scales
- `σ` (observation noise): `InverseGamma(alpha=2, beta=3)`
  - Stable for likelihood with scaled data

### Running Inference

```bash
python scripts/02_infer_degradation_rates.py \
    --transcription data/processed_transcription_data/transcription_traces_no_ids_stripe3.csv \
    --mrna data/processed_mRNA_data_stripe3/e1_sass_formodel.csv \
    --output-chain results/stripe3/e1/chains/degradation_chain.csv \
    --output-trace results/stripe3/e1/figures/mcmc_trace.png \
    --n-samples 10000 \
    --n-chains 4 \
    --n-ap-bins 5 \
    --n-dv-bins 5
```

**Key Features:**
- Analytical ODE solution for numerical stability
- NUTS MCMC sampler with convergence diagnostics
- Data normalization to prevent overflow
- Outputs posterior samples, trace plots, and summary statistics

**Outputs:**
- MCMC chain CSV with posterior samples for D, γ, σ
- Trace plot for convergence checking
- Summary statistics including half-lives per bin

### Prior Rationale

Priors are chosen to:
1. **Match Biology**: Half-lives 1-10 min for dynamic developmental genes
2. **Ensure Stability**: Bounds prevent extreme values causing exp() overflow
3. **Center Appropriately**: mu=-2; LogNormal centered at -2.0 gives mode ~0.13 min^-1 (t1/2 ~ 5 min)
4. **Allow Flexibility**: Wide enough for data-driven inference

See script docstring and `.research/logs/activity.md` for detailed decision log.

## mRNA Processing Pipeline (Scripts 04-06)

### 06_compute_validation_thresholds.py

Computes stripe-specific nuclei density validation thresholds from Berrocal_2020 data.

```bash
python scripts/06_compute_validation_thresholds.py \
    --input data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv \
    --config config.yaml \
    --k 4 \
    --max-time 1200
```

**What it does:**

1. Analyzes Berrocal_2020 transcription data per stripe
2. Computes k=4 nearest neighbor distances for each embryo
3. Calculates expected nuclei count, median NN distance
4. Determines valid range (min/max across embryos)
5. Writes thresholds to config.yaml under `validation_thresholds`

**Output format in config.yaml:**

```yaml
validation_thresholds:
  stripe2:
    expected_nuclei_count_mean: 186.5
    expected_nn_distance_median: 0.163
    min_nn_distance: 0.136
    max_nn_distance: 0.201
    k: 4
    n_embryos_used: 4
  stripe3:
    # ... computed for stripe3
```

This must run after `00_identify_stripe_ranges.py` but before validation. The pipeline handles this automatically.

### 04_aggregate_embryo_mrna.py

Bins mRNA spot data from SASS position_data-intense.txt into 5×5 spatial grid matching transcription data binning.

```bash
python scripts/04_aggregate_embryo_mrna.py \
    data/Ali_embryos/stripe3/e1/position_data-intense.txt \
    stripe3 \
    data/processed_mRNA_data_stripe3/e1_sass_formodel.csv
```

**What it does:**

1. Loads SASS-processed spot-to-nucleus assignments
2. Filters nuclei to stripe AP range (from config.yaml)
3. Normalizes AP and DV coordinates to [0, 1]
4. Creates 5×5 spatial bins (matching transcription binning)
5. Computes average mRNA count per nucleus in each bin
6. Outputs 25 values (no header) for model input

### 05_validate_nuclei_density.py

Validates nuclei density against stripe-specific Berrocal_2020 benchmarks using k-nearest neighbors analysis.

```bash
python scripts/05_validate_nuclei_density.py \
    data/Ali_embryos/stripe3/e1/position_data-intense.txt \
    stripe3 \
    config.yaml \
    results/stripe3/validation/e1_nuclei_density_validation.txt
```

**What it does:**

1. Loads stripe-specific validation thresholds from config.yaml
2. Computes k=4 nearest neighbor distances in normalized coordinates
3. Validates median NN distance is within observed range (hard gate)
4. Checks nuclei count and other metrics (soft warnings)
5. Writes validation report
6. **Exits with error if validation fails** (blocks pipeline)

**Validation criteria (stripe-specific):**

- **Hard gate**: Median NN distance must be within [min, max] from Berrocal data
- **Soft warning**: Nuclei count should be within expected range
- **Soft warning**: Median NN distance should be close to expected value

Thresholds are dynamically computed from Berrocal_2020 data per stripe, not hardcoded.

## Notes

- Pipeline scripts use numbered prefixes to indicate order
- Utility scripts in subdirectories (data/, utils/, setup/) are not tracked in main pipeline
- Keep exploratory notebooks in a separate `notebooks/` directory if needed
- Always run `identify_stripe_ranges.py` before starting the main pipeline
- mRNA processing (scripts 04-05) runs automatically through Snakemake for stripe3+ data
````
