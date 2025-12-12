# Project Handoff Notes

**From:** Jenny  
**Date:** December 2024  
**Project:** Eve transcription spatial mRNA decay analysis

## Data Processing Overview

### Data Filtering & Cleaning

- **Filtered data used** to remove:
  - Nuclei over-traveling during the time window
  - Nuclei tracked for only a few frames (avoids double counting)

### Data Registration & Alignment

- Used **registered AP data** to align data across different embryos
- Stripe positioned **perpendicular to AP axis**
- smiFISH data aligned the same way (AP axis perpendicular) on microscope for registration

### Stripe 2 Selection Criteria

**Initial approach:**
- Looked at min/max `ap_registered` values for nuclei in stripe 2
- Extended boundaries to include interstripes
- **Final filters: 0.33 and 0.45**

**Border optimization tested:**
- Tried posterior borders: 0.44, 0.45, 0.46
- Results consistent across these values
- **Final choice: 0.45** - best matched observed nuclei count from smFISH data (~10 nuclei)

**Validation:**
- Spot density consistent across domain
- No missing nuclei with 0 fluorescence in filtered interstripes

### Binning Strategy

**Spatial binning:**
- **5 bins across AP axis**
- **5 bins across DV axis**
- Total: **25 bins**

**Fluorescence calculation:**
- Mean fluorescent trace computed per bin

**Alignment validation:**
- Summed heatmaps: stripe center aligns well with smFISH data bins ✓

### smFISH Data Processing

**Spot detection & assignment:**
- Spots detected in Imaris
- Assigned to nuclei using sass script (Tom Minchington)
- Binned in X and Y (no Y-cropping needed)
- Average mRNAs per nucleus plotted per bin

## Analysis Pipeline

1. **Preprocessing:** Filtered transcription data → binned data (25 bins)
2. **Modeling:** Julia script processes transcription + averaged mRNA per cell data
3. **Visualization:** Python plotting script (Julia plotting "clunky", reproduced in Python)

## Important Decision: DV Cropping

**Initial concern:** Central third of DV axis cropping needed?

**Resolution (Jenny's follow-up):**
- ❌ **NOT needed** - single embryo plots show central section already correct for AP bins
- Data already covers whole DV axis appropriately
- Eve pattern **invariant over DV axis** for each AP bin
- May check for edge effects, but cropping shouldn't affect results
- ✅ Confirmed: current DV region is appropriate

## Files Added to Project

**Scripts:**
- `eve_transcription_data_processing_feb25.ipynb`
- `gp_fits_pbody_coloc_learning_noise.ipynb`
- `infer_D_across_stripe2.jl`
- `plotting_modelling_results.ipynb`
- `plotting_mRNA_heatmaps_oct24.ipynb`

**Data:**
- `data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520.csv`
- `data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv`
- `data/processed_transcription_data_feb25/*.csv`

## Next Steps / To Verify

- [ ] Confirm DV region captures full pattern (edge effect check)
- [ ] Document which specific embryo plots validated the central section decision
- [ ] Verify 0.45 AP boundary choice with final results
- [ ] Link binning parameters in `params.yaml` to this documentation

---

**Note:** Jenny away until 17 December 2024, monitoring email for urgent items.
