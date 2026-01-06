# Data Source Notes

## Stripe2 Data (from Jenny)

### Nuclear Cycle Stage Reference

The stripe2 mRNA data filenames use **nuclear cycle notation**, not embryo IDs:

| Filename Prefix | Nuclear Cycle | Development Stage | Notes |
|----------------|---------------|-------------------|-------|
| `e7_...` | Nuclear cycle 7 | ~2h 10min AEL | Earlier timepoint |
| `e8_9_...` | Cycle 8-9 transition | ~2h 20-30min AEL | **Primary analysis stage** |
| `e9_10_...` | Cycle 9-10 transition | ~2h 40-50min AEL | Later timepoint |

**AEL**: After Egg Laying

### Data Characteristics

- **Source**: Received pre-processed from Jenny
- **Processing**: Already run through SASS pipeline
- **Edge spots**: Already removed
- **Stripe centering**: Already centered
- **Original raw images**: Not currently available (attempting to locate)

The e8_9 stage (nuclear cycle 8-9 transition) is used for primary analysis because it represents mid-gastrulation when eve stripes are well-established and mRNA patterns are stable.

## Stripe3+ Data (from Ali)

### Embryo Numbering

Filenames use **embryo numbering** (e1, e2, e3, e4), not nuclear cycle stages.

All embryos were imaged at the same developmental stage (approximately nuclear cycle 8-9 transition, matching the e8_9 timepoint used for stripe2).

### Data Characteristics

- **Source**: Raw imaging data from Ali
- **Processing**: Raw Imaris exports in `data/Ali_embryos/{stripe}/{embryo}/`
- **Edge spots**: Manually removed before Imaris export
- **Stripe centering**: Manually centered before imaging
- **Pipeline**: Full SASS processing in Snakemake pipeline

## Naming Convention Summary

| Data Type | Stripe | Naming Convention | Example | Meaning |
|-----------|--------|-------------------|---------|---------|
| mRNA (Jenny) | stripe2 | Nuclear cycle | `e8_9_...` | Cycle 8-9 transition |
| mRNA (Ali) | stripe3+ | Embryo ID | `e1_...` | Embryo 1 |
| Transcription | All | AP range code | `033045` | AP 0.33-0.45 |

**Key Point**: The letter "e" followed by numbers means different things in different contexts. Always check the data source and documentation.
