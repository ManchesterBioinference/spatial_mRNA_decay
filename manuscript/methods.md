# Methods

<!-- 
Document your methodology here.
This section should enable complete reproduction of your work.

Run /write_methods after completing scripts to auto-generate from code.
Run /review_script first to ensure code is properly documented.
-->

## Overview

<!-- Brief summary of your analytical approach (1 paragraph) -->


## Data

### Data Sources

<!-- Where does your data come from? Access information? -->


### Data Preprocessing

<!-- Cleaning, filtering, transformation steps -->
<!-- Reference: pipeline/scripts/[preprocessing_script] -->


## Analysis Pipeline

### [Stage 1 Name]

<!-- Description of first analysis stage -->
<!-- Script: pipeline/scripts/[script_name] -->
<!-- Key parameters: -->


### [Stage 2 Name]

<!-- Continue for each pipeline stage -->


## Statistical Analysis

<!-- Tests used, significance thresholds, multiple testing correction -->


## Software and Reproducibility

<!-- 
List tools and versions used.
Include command for reproduction.
-->

All analyses were performed using:
- Python X.X with [packages]
- DVC X.X for pipeline management

The complete pipeline can be reproduced with:

```bash
# Clone repository
git clone [repository-url]
cd [project-name]

# Set up environment
uv sync  # or: conda env create -f environment.yml

# Run pipeline
dvc repro
```

---

## Code-to-Methods Reference

| Script | Section | Purpose |
|--------|---------|---------|
| | | |

<!-- This table is auto-populated by /write_methods -->
