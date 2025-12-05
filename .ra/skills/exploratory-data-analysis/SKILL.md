---
name: exploratory-data-analysis
description: Perform comprehensive exploratory data analysis on research data. Automatically analyze data structure, quality, distributions, and generate insights. Use when the user provides a dataset, asks to "explore data", "analyze this file", or needs to understand their data before formal analysis.
---

# Exploratory Data Analysis (EDA)

> Systematically explore and understand datasets before formal analysis.

## When to Use

- User provides a data file for analysis (CSV, Excel, HDF5, etc.)
- User asks to "explore", "analyze", or "summarize" data
- Starting the ANALYSIS phase of a research project
- Before running formal statistical tests
- Assessing data quality and completeness
- Understanding data distributions and relationships
- Identifying outliers and anomalies

## EDA Workflow

```
1. LOAD DATA     → Read file, check structure
2. SUMMARIZE     → Basic statistics, data types
3. QUALITY       → Missing values, outliers, duplicates
4. DISTRIBUTIONS → Visualize variable distributions
5. RELATIONSHIPS → Correlations, group comparisons
6. DOCUMENT      → Generate EDA report
```

---

## Step 1: Load and Inspect Data

```python
import pandas as pd
import numpy as np

# Load data
df = pd.read_csv('data.csv')  # Adjust for your file type

# Basic inspection
print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
print(f"\nColumn types:\n{df.dtypes}")
print(f"\nFirst few rows:\n{df.head()}")
print(f"\nMemory usage: {df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")
```

---

## Step 2: Summary Statistics

```python
# Numerical columns
print("Numerical Summary:")
print(df.describe().T[['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']])

# Categorical columns
for col in df.select_dtypes(include=['object', 'category']).columns:
    print(f"\n{col}:")
    print(df[col].value_counts().head(10))
```

### Key Statistics to Report

| Statistic | Purpose |
|-----------|---------|
| n (count) | Sample size, completeness |
| Mean | Central tendency |
| SD | Spread/variability |
| Min/Max | Range, potential outliers |
| Quartiles | Distribution shape |
| Unique values | Cardinality for categoricals |

---

## Step 3: Data Quality Assessment

### Missing Values

```python
# Missing value summary
missing = df.isnull().sum()
missing_pct = (missing / len(df) * 100).round(2)
missing_df = pd.DataFrame({
    'Missing': missing,
    'Percent': missing_pct
}).query('Missing > 0').sort_values('Percent', ascending=False)

print("Missing Values:")
print(missing_df)

# Visualize missing pattern
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 6))
sns.heatmap(df.isnull(), cbar=True, yticklabels=False, cmap='viridis')
plt.title('Missing Data Pattern')
plt.tight_layout()
plt.savefig('results/eda_missing_pattern.png', dpi=150)
```

### Outlier Detection

```python
def detect_outliers_iqr(data, column):
    """Detect outliers using IQR method."""
    Q1 = data[column].quantile(0.25)
    Q3 = data[column].quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    outliers = data[(data[column] < lower) | (data[column] > upper)]
    return outliers, lower, upper

# Check all numerical columns
for col in df.select_dtypes(include=[np.number]).columns:
    outliers, lower, upper = detect_outliers_iqr(df, col)
    if len(outliers) > 0:
        print(f"{col}: {len(outliers)} outliers ({len(outliers)/len(df)*100:.1f}%)")
        print(f"  Range: [{lower:.2f}, {upper:.2f}]")
```

### Duplicates

```python
# Check for duplicate rows
duplicates = df.duplicated().sum()
print(f"Duplicate rows: {duplicates} ({duplicates/len(df)*100:.1f}%)")

# Check for duplicate IDs (if applicable)
if 'id' in df.columns:
    dup_ids = df['id'].duplicated().sum()
    print(f"Duplicate IDs: {dup_ids}")
```

---

## Step 4: Distribution Analysis

### Numerical Variables

```python
import matplotlib.pyplot as plt
import seaborn as sns

numerical_cols = df.select_dtypes(include=[np.number]).columns

fig, axes = plt.subplots(len(numerical_cols), 2, figsize=(12, 4*len(numerical_cols)))

for i, col in enumerate(numerical_cols):
    # Histogram
    axes[i, 0].hist(df[col].dropna(), bins=30, edgecolor='black', alpha=0.7)
    axes[i, 0].set_title(f'{col} - Distribution')
    axes[i, 0].set_xlabel(col)
    axes[i, 0].set_ylabel('Frequency')
    
    # Box plot
    axes[i, 1].boxplot(df[col].dropna())
    axes[i, 1].set_title(f'{col} - Box Plot')
    axes[i, 1].set_ylabel(col)

plt.tight_layout()
plt.savefig('results/eda_distributions.png', dpi=150)
```

### Categorical Variables

```python
categorical_cols = df.select_dtypes(include=['object', 'category']).columns

for col in categorical_cols:
    plt.figure(figsize=(10, 4))
    df[col].value_counts().head(15).plot(kind='bar', edgecolor='black')
    plt.title(f'{col} - Value Counts')
    plt.xlabel(col)
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'results/eda_{col}_counts.png', dpi=150)
```

---

## Step 5: Relationship Analysis

### Correlation Matrix

```python
# Numerical correlations
corr_matrix = df.select_dtypes(include=[np.number]).corr()

plt.figure(figsize=(10, 8))
mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f', 
            cmap='RdBu_r', center=0, square=True)
plt.title('Correlation Matrix')
plt.tight_layout()
plt.savefig('results/eda_correlations.png', dpi=150)

# Identify strong correlations
strong_corr = []
for i in range(len(corr_matrix.columns)):
    for j in range(i+1, len(corr_matrix.columns)):
        if abs(corr_matrix.iloc[i, j]) > 0.7:
            strong_corr.append({
                'var1': corr_matrix.columns[i],
                'var2': corr_matrix.columns[j],
                'correlation': corr_matrix.iloc[i, j]
            })
if strong_corr:
    print("Strong correlations (|r| > 0.7):")
    for c in strong_corr:
        print(f"  {c['var1']} ↔ {c['var2']}: r = {c['correlation']:.3f}")
```

### Pairwise Scatter Plots

```python
# For key variables only (to avoid overwhelming output)
key_vars = ['var1', 'var2', 'var3']  # Adjust to your variables
sns.pairplot(df[key_vars], diag_kind='hist')
plt.savefig('results/eda_pairplot.png', dpi=150)
```

### Group Comparisons

```python
# If you have a grouping variable
if 'group' in df.columns:
    for col in df.select_dtypes(include=[np.number]).columns:
        plt.figure(figsize=(8, 5))
        df.boxplot(column=col, by='group')
        plt.title(f'{col} by Group')
        plt.suptitle('')  # Remove automatic title
        plt.tight_layout()
        plt.savefig(f'results/eda_{col}_by_group.png', dpi=150)
```

---

## Step 6: Generate EDA Report

### Report Template

```markdown
# Exploratory Data Analysis Report

**Dataset**: [filename]
**Date**: [date]
**Analyst**: [name]

## 1. Data Overview

- **Rows**: X
- **Columns**: Y
- **File size**: Z MB

## 2. Variable Summary

| Variable | Type | Non-Null | Unique | Mean | SD |
|----------|------|----------|--------|------|-----|
| var1 | float64 | 100 | 50 | 25.3 | 5.2 |
| ... | ... | ... | ... | ... | ... |

## 3. Data Quality

### Missing Values
- [List variables with missing data and percentages]

### Outliers
- [List variables with outliers detected]

### Duplicates
- [Number of duplicate rows]

## 4. Key Findings

1. **Finding 1**: Description
2. **Finding 2**: Description
3. **Finding 3**: Description

## 5. Recommendations

- [ ] Handle missing values in [variable] using [method]
- [ ] Consider transformation for [variable] (skewed distribution)
- [ ] Investigate outliers in [variable]
- [ ] Check data collection for [issue noted]

## 6. Next Steps

Based on this EDA, the following analyses are recommended:
1. [Recommended analysis 1]
2. [Recommended analysis 2]
```

---

## Integration with RA Workflow

### ANALYSIS Phase Connection

After completing EDA:
1. Document findings in `.research/logs/activity.md`
2. Update `tasks.md` with identified issues to address
3. Proceed to formal statistical analysis with `/statistical_analysis`
4. Save figures to `results/intermediate/` or `manuscript/figures/`

### Files to Create

| File | Location | Purpose |
|------|----------|---------|
| EDA report | `results/eda_report.md` | Document findings |
| Distribution plots | `results/intermediate/` | Quality check |
| Correlation matrix | `results/intermediate/` | Relationship overview |
| Missing data pattern | `results/intermediate/` | Data quality |

---

## Quick EDA Checklist

- [ ] Loaded data and verified structure
- [ ] Checked data types are correct
- [ ] Calculated summary statistics
- [ ] Identified and documented missing values
- [ ] Detected outliers
- [ ] Checked for duplicates
- [ ] Visualized distributions
- [ ] Examined correlations/relationships
- [ ] Documented key findings
- [ ] Listed recommended next steps
