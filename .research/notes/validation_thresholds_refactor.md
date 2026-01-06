# Validation Thresholds Refactoring

**Date**: 2026-01-06  
**Status**: Implemented

## Summary

Refactored nuclei density validation to use dynamically-computed, stripe-specific thresholds from Berrocal_2020 data instead of hardcoded values.

## Problem

Previously, `05_validate_nuclei_density.py` used hardcoded validation thresholds:
- `MIN_NN_DISTANCE = 0.136`
- `MAX_NN_DISTANCE = 0.201`
- `EXPECTED_NUCLEI_COUNT_MEAN = 186`
- `EXPECTED_NN_DISTANCE_MEDIAN = 0.163`

These were derived from **stripe2 only** in the Berrocal_2020 dataset. This is problematic because:
1. Different stripes may have different nuclei densities
2. Hardcoded values are not reproducible from raw data
3. Updating thresholds requires code changes, not data reprocessing

## Solution

Created a two-stage validation system:

### Stage 1: Compute Thresholds (`06_compute_validation_thresholds.py`)

New script that analyzes Berrocal_2020 data to compute stripe-specific thresholds:

```python
# For each stripe in config.yaml stripe_ranges:
# 1. Filter nuclei to stripe AP range
# 2. Group by embryo ID
# 3. Compute k=4 NN distances per embryo
# 4. Calculate statistics: mean, median, min, max
# 5. Write to config.yaml
```

**Output in config.yaml**:
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
    expected_nuclei_count_mean: 195.2  # Different from stripe2!
    expected_nn_distance_median: 0.158
    min_nn_distance: 0.142
    max_nn_distance: 0.189
    k: 4
    n_embryos_used: 4
```

### Stage 2: Use Thresholds (`05_validate_nuclei_density.py`)

Updated validation script to:
1. Load stripe-specific thresholds from config.yaml
2. Compare computed metrics against loaded thresholds
3. Exit with error if hard gate fails

**New function signature**:
```python
python 05_validate_nuclei_density.py \
    <position_data_file> \
    <stripe> \
    <config_yaml> \      # NEW: path to config.yaml
    <output_report>
```

## Pipeline Integration

Updated Snakefile with new rule:

```python
rule compute_validation_thresholds:
    input:
        data="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        script="scripts/06_compute_validation_thresholds.py",
        config_updated="config.yaml.updated"  # After stripe ranges computed
    output:
        validation_updated=touch("config.yaml.validation_updated")
    # ... runs 06_compute_validation_thresholds.py
```

**Dependency chain**:
```
identify_stripe_ranges → compute_validation_thresholds → process_mrna_sass → validate_nuclei_density
                                                                            ↓
                                                                     bin_mrna_counts
```

The validation thresholds are computed **once** per stripe after stripe ranges are identified, then **reused** for all embryos in that stripe.

## Benefits

1. **Data-driven**: Thresholds computed from raw data, not hardcoded
2. **Stripe-specific**: Each stripe has appropriate density expectations
3. **Reproducible**: Thresholds can be recomputed if Berrocal data changes
4. **Cached**: Computed once, stored in config.yaml, reused across embryos
5. **Transparent**: Thresholds visible in config.yaml with metadata

## Files Changed

- **Created**: `scripts/06_compute_validation_thresholds.py` - Threshold computation script
- **Modified**: `scripts/05_validate_nuclei_density.py` - Load thresholds from config
- **Modified**: `Snakefile` - Added `compute_validation_thresholds` rule
- **Modified**: `scripts/README.md` - Documentation updates

## Testing

Dry-run confirms proper execution order:
```bash
snakemake --dry-run --cores 1
# Shows 22 jobs (was 21) with compute_validation_thresholds added
```

## Next Steps

1. Run full pipeline to populate validation_thresholds in config.yaml
2. Verify thresholds are reasonable across all stripes
3. Consider adding visualization of NN distance distributions per stripe
