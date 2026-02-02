# Tasks

> Quick todos for this project. Updated by the Research Assistant and manually.
>
> Items here should take < 2 hours. Larger items → GitHub Issues.

## High Priority

- [X] Prototype overlapping-bin modeling and compare to current 5-bin outputs
- [X] Recompute validation thresholds with updated AP/DV bin counts
- [X] Clean up scripts/README.md - remove outdated mu=0.07 prior explanations (now using Julia-matched priors)
- [X] Read Jenny's paper on Eve stripe analysis (Before next meeting)
- [X] Contact Ali about SNFISH data processing methods (This week)
- [X] Contact Ali about mRNA counting methodology (This week)

## Normal Priority

- [x] Verify DV region captures full eve pattern (check for edge effects)
- [x] Document binning parameters (5×5 grid, 25 bins) in params.yaml
- [x] Link AP boundary choice (0.45) to params.yaml with justification
- [x] Process new Ali data for e8_9 and e9_10
- [x] Reassess time window alignment between public time-series and embryo staging
- [x] Update methods.md with spatial binning strategy from handoff notes
- [x] Compare 5x5 vs 5x7 results and save summary
- [x] Review combined half-life heatmaps grid at results/combined_halflife_heatmaps_grid.png
- [x] Get bin ratios between mRNA and transcription data for 1200
- [ ] Implement size adjustment normalization for failing embryos (Magnus meeting 2026-01-30)
- [ ] Process stripe 4 data when Ali provides it (Magnus meeting 2026-01-30)
- [ ] Reprocess stripe 3 with updated edge removal when Ali provides it (Magnus meeting 2026-01-30)

## Low Priority / Someday

## Completed

- [x] 2026.01.06 - Refactored validation thresholds to be data-driven and stripe-specific

- [x] Created `06_compute_validation_thresholds.py` to compute thresholds from Berrocal_2020 data
- [x] Updated `05_validate_nuclei_density.py` to load thresholds from config.yaml
- [x] Added `compute_validation_thresholds` rule to Snakefile
- [x] Thresholds now computed per stripe and cached in config.yaml

---

## From Meetings

### 2026.01.06 Ali & Logan Chat

**HIGH PRIORITY:**

- [X] Build automated stripe-alignment pipeline (currently manual; should align centers and handle multi-stripe analysis)

**MEDIUM PRIORITY:**

- [X] Remove Y-axis cropping filter from analysis (confirmed public data uses full height; simplify pipeline)
- [X] Process remaining embryo data on CSF (can handle up to 100+ GB memory per job; will be much faster than laptop)

---

*Use `/summarize_meeting` to automatically extract tasks from meeting notes.*
*Use `/plan_week` to organize tasks into a weekly plan.*
