# Data Directory

This directory contains all project data. Follow these guidelines to maintain reproducibility and data integrity.

## Directory Structure

```
data/
├── raw/           # Original, immutable source data
├── processed/     # Cleaned, transformed data (DVC-tracked)
└── .sensitive/    # Protected data (excluded from AI via .copilotignore)
```

## Data Management Rules

### 1. Raw Data is Immutable

**NEVER modify files in `data/raw/`**

- Raw data should be exactly as received/downloaded
- All transformations go to `data/processed/`
- Document data provenance (source, date acquired, version)

```python
# ❌ WRONG
raw_data = pd.read_csv("data/raw/samples.csv")
raw_data.dropna(inplace=True)
raw_data.to_csv("data/raw/samples.csv")  # Overwrites original!

# ✅ CORRECT
raw_data = pd.read_csv("data/raw/samples.csv")
processed_data = raw_data.dropna()
processed_data.to_csv("data/processed/samples_cleaned.csv")
```

### 2. Document Data Provenance

For each dataset in `data/raw/`, create a companion file:

```markdown
<!-- data/raw/samples_provenance.md -->
# Dataset: samples.csv

**Source**: [Database/Lab/Collaborator name]
**URL**: [Download link if applicable]
**Date acquired**: YYYY-MM-DD
**Version**: [If versioned dataset]
**Description**: [Brief description]
**Columns**: [List of columns and their meaning]
**Rows**: [Number of samples/observations]
**Notes**: [Any relevant context]
```

### 3. Tidy Data Principles

When processing data, aim for tidy format:

- Each variable in its own column
- Each observation in its own row
- Each type of observational unit in its own table

```python
# ❌ Wide format (harder to analyze)
# subject, day1_score, day2_score, day3_score

# ✅ Tidy/long format (easier to analyze)
# subject, day, score
```

### 4. Verify Data Integrity

Use checksums to verify raw data hasn't been corrupted:

```bash
# Generate checksum when first acquiring data
md5 data/raw/samples.csv > data/raw/samples.csv.md5

# Verify later
md5 -c data/raw/samples.csv.md5
```

Or in Python:
```python
import hashlib

def file_checksum(filepath):
    return hashlib.md5(open(filepath, 'rb').read()).hexdigest()

# Store this value in provenance documentation
print(file_checksum("data/raw/samples.csv"))
```

### 5. Use DVC for Large/Processed Data

Files in `data/processed/` should be tracked with DVC:

```bash
dvc add data/processed/large_dataset.csv
git add data/processed/large_dataset.csv.dvc
git commit -m "data: Add processed dataset"
```

### 6. Naming Conventions

```
YYYYMMDD_description_version.ext

Examples:
20241202_samples_raw_v01.csv
20241202_samples_normalized_v02.csv
```

- Use ISO dates for sortability
- Include version numbers
- Use underscores, no spaces
- Be descriptive but concise

## Sensitive Data

The `.sensitive/` subdirectory is:
- Excluded from git (via `.gitignore`)
- Excluded from AI assistant (via `.copilotignore`)
- For data that cannot be shared (PHI, proprietary, etc.)

See `data/.sensitive/README.md` for handling instructions.

## Data Lifecycle

```
1. ACQUIRE    → Place in data/raw/, create provenance doc
2. VALIDATE   → Generate checksum, verify integrity
3. PROCESS    → Transform to data/processed/ via pipeline
4. TRACK      → DVC add processed files
5. ANALYZE    → Reference processed data in scripts
6. ARCHIVE    → Include data statement in manuscript
```

## Data Availability Statement Template

For your manuscript:

```markdown
The data that support the findings of this study are available 
[in the repository / from the corresponding author upon reasonable request / 
from [source] at [URL]]. [Any restrictions on data access].
```
