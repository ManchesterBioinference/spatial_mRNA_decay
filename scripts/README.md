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

## Notes

- Pipeline scripts use numbered prefixes to indicate order
- Utility scripts in subdirectories (data/, utils/, setup/) are not DVC tracked
- Keep exploratory notebooks in a separate `notebooks/` directory if needed
