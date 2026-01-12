# Tasks

> Quick todos for this project. Updated by the Research Assistant and manually.
> Items here should take < 2 hours. Larger items → GitHub Issues.

## High Priority
<!-- Do these first -->

- [ ] Clean up scripts/README.md - remove outdated mu=0.07 prior explanations (now using Julia-matched priors)
- [ ] Read Jenny's paper on Eve stripe analysis (Before next meeting)
- [ ] Contact Ali about SNFISH data processing methods (This week)
- [ ] Contact Ali about mRNA counting methodology (This week)

## Normal Priority
<!-- Do when high priority is clear -->

- [ ] Verify DV region captures full eve pattern (check for edge effects)
- [ ] Document binning parameters (5×5 grid, 25 bins) in params.yaml
- [ ] Link AP boundary choice (0.45) to params.yaml with justification
- [ ] Update methods.md with spatial binning strategy from handoff notes

## Low Priority / Someday
<!-- Nice to have, not urgent -->


## Completed
<!-- Move items here when done, with date -->

- [x] 2026.01.06 - Refactored validation thresholds to be data-driven and stripe-specific
  - Created `06_compute_validation_thresholds.py` to compute thresholds from Berrocal_2020 data
  - Updated `05_validate_nuclei_density.py` to load thresholds from config.yaml
  - Added `compute_validation_thresholds` rule to Snakefile
  - Thresholds now computed per stripe and cached in config.yaml

---

## From Meetings
<!-- Tasks extracted from meeting summaries appear here -->

### 2026.01.06 Ali & Logan Chat

**HIGH PRIORITY:**

- [ ] Build automated stripe-alignment pipeline (currently manual; should align centers and handle multi-stripe analysis)

**MEDIUM PRIORITY:**

- [ ] Remove Y-axis cropping filter from analysis (confirmed public data uses full height; simplify pipeline)
- [ ] Process remaining embryo data on CSF (can handle up to 100+ GB memory per job; will be much faster than laptop)

---

*Use `/summarize_meeting` to automatically extract tasks from meeting notes.*
*Use `/plan_week` to organize tasks into a weekly plan.*
