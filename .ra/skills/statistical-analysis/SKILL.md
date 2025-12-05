---
name: statistical-analysis
description: Comprehensive statistical analysis toolkit for research. Conduct hypothesis tests (t-test, ANOVA, chi-square), regression, correlation, Bayesian stats, power analysis, assumption checks, and APA reporting. Use when the user asks about statistics, needs help analyzing data, or when writing methods sections that include statistical approaches.
---

# Statistical Analysis

> Comprehensive statistical testing, power analysis, and experimental design for reproducible research.

## When to Use

- Conducting statistical hypothesis tests (t-tests, ANOVA, chi-square)
- Performing regression or correlation analyses
- Running Bayesian statistical analyses
- Checking statistical assumptions and diagnostics
- Calculating effect sizes and conducting power analyses
- Reporting statistical results in APA format
- Planning experiments with proper power calculations
- Helping with the ANALYSIS phase of a research project

## Workflow Decision Tree

```
START
│
├─ Need to SELECT a statistical test?
│  └─ See "Test Selection Guide"
│
├─ Ready to check ASSUMPTIONS?
│  └─ See "Assumption Checking"
│
├─ Ready to run ANALYSIS?
│  └─ See "Running Statistical Tests"
│
└─ Need to REPORT results?
   └─ See "Reporting Results (APA)"
```

---

## Test Selection Guide

### Quick Reference: Choosing the Right Test

**Comparing Two Groups:**
| Data Type | Distribution | Design | Test |
|-----------|--------------|--------|------|
| Continuous | Normal | Independent | Independent t-test |
| Continuous | Non-normal | Independent | Mann-Whitney U |
| Continuous | Normal | Paired | Paired t-test |
| Continuous | Non-normal | Paired | Wilcoxon signed-rank |
| Binary | - | - | Chi-square / Fisher's exact |

**Comparing 3+ Groups:**
| Data Type | Distribution | Design | Test |
|-----------|--------------|--------|------|
| Continuous | Normal | Independent | One-way ANOVA |
| Continuous | Non-normal | Independent | Kruskal-Wallis |
| Continuous | Normal | Paired | Repeated measures ANOVA |
| Continuous | Non-normal | Paired | Friedman test |

**Relationships:**
| Analysis | Use Case | Test |
|----------|----------|------|
| Two continuous vars | Normal | Pearson correlation |
| Two continuous vars | Non-normal | Spearman correlation |
| Continuous outcome + predictor(s) | Prediction | Linear regression |
| Binary outcome + predictor(s) | Classification | Logistic regression |

---

## Assumption Checking

**ALWAYS check assumptions before interpreting test results.**

### Key Assumptions to Check

```python
import scipy.stats as stats
import numpy as np

# 1. Normality Test (Shapiro-Wilk)
stat, p = stats.shapiro(data)
print(f"Shapiro-Wilk: W={stat:.3f}, p={p:.3f}")
if p < 0.05:
    print("⚠️ Normality assumption violated - consider non-parametric test")

# 2. Homogeneity of Variance (Levene's test)
stat, p = stats.levene(group1, group2)
print(f"Levene's: F={stat:.3f}, p={p:.3f}")
if p < 0.05:
    print("⚠️ Variance assumption violated - use Welch's t-test")

# 3. Outlier Detection (IQR method)
Q1, Q3 = np.percentile(data, [25, 75])
IQR = Q3 - Q1
outliers = data[(data < Q1 - 1.5*IQR) | (data > Q3 + 1.5*IQR)]
print(f"Outliers detected: {len(outliers)}")
```

### What to Do When Assumptions Are Violated

| Assumption | Violation | Solution |
|------------|-----------|----------|
| Normality (mild, n>30) | Proceed | Parametric tests are robust |
| Normality (severe) | Transform | Use log/sqrt or non-parametric |
| Homogeneity of variance | t-test | Use Welch's t-test |
| Homogeneity of variance | ANOVA | Use Welch's ANOVA |
| Linearity (regression) | Violated | Add polynomial terms or use GAM |

---

## Running Statistical Tests

### Python Libraries

```python
import scipy.stats as stats       # Core statistical tests
import statsmodels.api as sm      # Regression, diagnostics
import pingouin as pg             # User-friendly testing
import numpy as np
import pandas as pd
```

### Common Analyses

#### T-Test with Complete Reporting

```python
import pingouin as pg

# Independent t-test with effect size
result = pg.ttest(group_a, group_b, correction='auto')
print(f"t({result['dof'].values[0]:.0f}) = {result['T'].values[0]:.2f}, "
      f"p = {result['p-val'].values[0]:.3f}, "
      f"d = {result['cohen-d'].values[0]:.2f}")
```

#### One-Way ANOVA with Post-Hoc

```python
import pingouin as pg

# ANOVA
aov = pg.anova(dv='score', between='group', data=df, detailed=True)
print(f"F = {aov['F'].values[0]:.2f}, p = {aov['p-unc'].values[0]:.3f}, "
      f"η²_p = {aov['np2'].values[0]:.3f}")

# Post-hoc if significant
if aov['p-unc'].values[0] < 0.05:
    posthoc = pg.pairwise_tukey(dv='score', between='group', data=df)
    print(posthoc[['A', 'B', 'diff', 'p-tukey']])
```

#### Linear Regression with Diagnostics

```python
import statsmodels.api as sm

# Fit model
X = sm.add_constant(predictors)
model = sm.OLS(outcome, X).fit()
print(model.summary())

# Key outputs
print(f"R² = {model.rsquared:.3f}, Adjusted R² = {model.rsquared_adj:.3f}")
print(f"F({model.df_model:.0f}, {model.df_resid:.0f}) = {model.fvalue:.2f}, p = {model.f_pvalue:.4f}")
```

#### Correlation with Confidence Intervals

```python
import pingouin as pg

# Pearson correlation with CI
result = pg.corr(x, y, method='pearson')
print(f"r = {result['r'].values[0]:.3f}, "
      f"p = {result['p-val'].values[0]:.3f}, "
      f"95% CI [{result['CI95%'].values[0][0]:.3f}, {result['CI95%'].values[0][1]:.3f}]")
```

---

## Effect Sizes

**Always report effect sizes alongside p-values.**

### Quick Reference: Effect Size Benchmarks

| Test | Effect Size | Small | Medium | Large |
|------|-------------|-------|--------|-------|
| T-test | Cohen's d | 0.20 | 0.50 | 0.80 |
| ANOVA | η²_p (partial eta²) | 0.01 | 0.06 | 0.14 |
| Correlation | r | 0.10 | 0.30 | 0.50 |
| Regression | R² | 0.02 | 0.13 | 0.26 |
| Chi-square | Cramér's V | 0.07 | 0.21 | 0.35 |

**Important**: These are guidelines only. Practical significance depends on context.

---

## Power Analysis

### A Priori Power Analysis (Before Study)

```python
from statsmodels.stats.power import tt_ind_solve_power, FTestAnovaPower

# T-test: Required n for d=0.5, power=0.80, alpha=0.05
n = tt_ind_solve_power(effect_size=0.5, alpha=0.05, power=0.80, ratio=1.0)
print(f"Required n per group: {n:.0f}")

# ANOVA: Required n for f=0.25, 3 groups
power_anova = FTestAnovaPower()
n = power_anova.solve_power(effect_size=0.25, ngroups=3, alpha=0.05, power=0.80)
print(f"Required n per group: {n:.0f}")
```

### Sensitivity Analysis (After Study)

```python
# What effect could we detect with n=50 per group?
detectable_d = tt_ind_solve_power(effect_size=None, nobs1=50, alpha=0.05, 
                                   power=0.80, ratio=1.0)
print(f"Minimum detectable effect: d = {detectable_d:.2f}")
```

---

## Reporting Results (APA Format)

### Templates for Common Tests

**Independent T-Test:**
```
Group A (n = 48, M = 75.2, SD = 8.5) scored significantly higher than
Group B (n = 52, M = 68.3, SD = 9.2), t(98) = 3.82, p < .001, d = 0.77,
95% CI [0.36, 1.18].
```

**One-Way ANOVA:**
```
A one-way ANOVA revealed a significant main effect of treatment on test
scores, F(2, 147) = 8.45, p < .001, η²_p = .10. Post hoc comparisons
using Tukey's HSD indicated that Condition A (M = 78.2, SD = 7.3) 
differed significantly from Condition B (M = 71.5, SD = 8.1, p = .002).
```

**Pearson Correlation:**
```
There was a significant positive correlation between study hours and
exam scores, r(98) = .45, p < .001, 95% CI [.28, .59].
```

**Multiple Regression:**
```
Multiple regression was conducted with exam scores as the outcome.
The model was significant, F(3, 146) = 45.2, p < .001, R² = .48.
Study hours (β = .35, p < .001) and prior GPA (β = .28, p < .001)
were significant predictors.
```

---

## Integration with RA Workflow

### During PLANNING Phase
- Help determine appropriate sample sizes with power analysis
- Suggest statistical approaches for research design

### During ANALYSIS Phase
- Run assumption checks on collected data
- Perform planned statistical analyses
- Generate effect sizes and confidence intervals

### During WRITING Phase
- Format results for methods and results sections
- Generate APA-formatted statistical reports
- Connect to `/write_methods` and `/write_results` skills

---

## Essential Reporting Elements

Always include:
1. **Descriptive statistics**: M, SD, n for all groups
2. **Test statistics**: Name, statistic value, df, exact p-value
3. **Effect sizes**: With confidence intervals when possible
4. **Assumption checks**: What was tested, results, any corrections
5. **All planned analyses**: Including non-significant findings
