# Activity Log

> Running log of project activity. Updated after significant work sessions.
> The Research Assistant uses this to understand recent context.

---

## Quick Session Entry Template

```markdown
## [YYYY-MM-DD]

**Session focus**: [What you worked on]

**Accomplished**:
- [What was completed]

**Decisions made**:
- [Any choices or directions taken]

**Next steps**:
- [What to do next]

**Notes**:
- [Any other context]
```

---

## Lab Notebook Entry Template (For Experimental Work)

Use this format when documenting experiments, analyses, or pipeline runs:

```markdown
## [YYYY-MM-DD] - [Experiment/Analysis Name]

### Objective
[What you're trying to accomplish in 1-2 sentences]

### Protocol/Method
[Key steps or deviations from standard protocol]
- Used parameters: [list key settings from params.yaml]
- Modified from standard: [any changes]

### Observations
[What happened, including unexpected events]
- Runtime: [how long it took]
- Warnings/errors: [any issues encountered]

### Data/Outputs
[Where outputs are stored]
- Output files: `results/[filename]`
- Figures generated: `manuscript/figures/[figN]/`

### Preliminary Conclusions
[Initial interpretation - can be updated later]

### Next Steps
[What to do based on these results]

### Cross-References
- Related scripts: `pipeline/scripts/[name].py`
- Related params: `params.yaml` section `[section]`
- Related issue: #[issue number]
```

---

## [2026-01-06]

**Session focus**: Meeting with Ali on data alignment, pipeline refactoring for stripe-specific validation thresholds

**Accomplished**:
- **Meeting with Ali** (~40 min chat about data alignment and processing)
  - Confirmed nuclei density alignment: Ali's Embryo 2 perfectly matches Berrocal public data (44-45 nuclei, 0.16 spacing)
  - Identified Y-axis cropping as unnecessary (public data uses full stripe height)
  - Learned about edge spot removal methodology and SpotMe V2 workflow
  - Discussed Stripe 3 & 4 status (imaging complete, needs processing)
- **Refactored validation thresholds to be data-driven and stripe-specific**
  - Created `06_compute_validation_thresholds.py` to dynamically compute thresholds from Berrocal_2020 data
  - Updated `05_validate_nuclei_density.py` to load stripe-specific thresholds from config.yaml
  - Modified Snakefile to add `compute_validation_thresholds` rule in dependency chain
  - Documented entire refactor in `.research/notes/validation_thresholds_refactor.md`
- **Documentation updates**:
  - Enhanced `scripts/README.md` with detailed script descriptions and pipeline architecture
  - Updated READMEs for processed data directories
  - Documented meeting in `.research/meetings/transcripts/2026.01.06_Ali_Logan_Chat.md`
- **Pipeline improvements**:
  - Added SASS submodule integration to Snakefile
  - Created `sass.yml` environment for SASS processing
  - Updated `tasks.md` with new action items from Ali meeting

**Decisions made**:
- Validation thresholds should be stripe-specific and data-driven (not hardcoded)
- Can remove Y-axis cropping filter from pipeline (confirmed unnecessary)
- Will process remaining embryo data on CSF supercomputer (much faster than laptop)
- Need automated stripe-alignment pipeline (currently manual centering process)

**Next steps**:
- **HIGH PRIORITY**: Build automated stripe-alignment pipeline for multi-stripe analysis
- **HIGH PRIORITY**: Receive and process Stripe 3 data from Ali (expected tomorrow)
- **MEDIUM PRIORITY**: Remove Y-axis cropping filter from pipeline to simplify
- **MEDIUM PRIORITY**: Process remaining embryo data on CSF (100+ GB memory available)
- Review Jenny's paper on Eve stripe analysis (still pending)
- Contact Ali about SNFISH data processing methods (still pending)

**Notes**:
- Large uncommitted changes: 8 files modified (+1370, -19 lines)
- Significant refactoring but no commits made today
- Meeting transcript is 3200+ lines with detailed technical discussion
- Next meeting with team Thursday, Jan 9

---

## [2025-12-16]

**Session focus**: Received updated SNFISH data from Ali with corrected analysis methodology

**Accomplished**:

- Received new SNFISH data from Ali (`data/Ali_embryos/embryo1/`, `embryo2/`)
  - Data now processed according to Jenny's original methodology
  - Ali corrected his analysis based on more detailed instructions from Jenny
  - New file format: `position_data-intense.txt` with spot-level intensity measurements
- Explored transcription data in notebook (`notebooks/afterHandoff/exploreTranscriptionData.ipynb`)
  - Modified cells working with spatial peak detection and stripe identification
  - Currently viewing embryo-level filtering approach (AP range: 0.355-0.415)
- Started new EDA notebook for Ali's position data (`notebooks/eda_position_data_intense.ipynb`)

**Decisions made**:

- Need to review Ali's corrected data first thing tomorrow before integrating into pipeline
- Keep previous Berrocal data separate until confirming alignment between Ali's and Jenny's approaches

**Next steps**:

- **PRIORITY**: Review Ali's updated results tomorrow morning
  - Compare new `position_data-intense.txt` format with Jenny's previous processing
  - Verify that Ali's corrections align with Jenny's methodology
  - Check if this resolves reviewer concerns about data alignment
- Document differences between Ali's first pass and corrected analysis
- Integrate Ali's corrected data into preprocessing pipeline if validated

**Notes**:

- Ali got more detailed instructions from Jenny on proper SNFISH processing methodology
- This corrected analysis should better match Jenny's original work
- New data structure includes time series spot tracking with nuclear assignments
- Two embryos provided: embryo1 and embryo2 in `data/Ali_embryos/`

---

## [2025-12-15]

**Session focus**: Pipeline reorganization for multi-stripe analysis and mRNA processing integration

**Accomplished**:

- Reorganized pipeline to run on multiple stripes in structured way:
  - Stripe 2 (original focus) now runs through complete pipeline
  - Results stored in organized `results/stripe2/` directory structure
  - Generated new transcription heatmaps and stripe identification visualizations
- Received SASS script and processed mRNA data from Ali:
  - Edge spot removal methodology (need clarification from Jenny on rationale)
  - Raw data for mRNA processing now available
  - Plan to incorporate full SASS processing into pipeline for complete analysis chain
- Successfully ran end-to-end Snakemake pipeline with `snakemake --cores 4 --use-conda`
- Cleaned up legacy data structure (removed outdated `processed_mRNA_data/` files)
- Enhanced documentation:
  - Updated `scripts/README.md` with preprocessing step details
  - Improved `scripts/01_preprocess_eve_data.py` for stripe 2 compatibility

**Decisions made**:

- Keep stripe 2 results separately organized to enable comparison across stripes
- Prioritize incorporating SASS mRNA processing into pipeline once Jenny clarifies edge spot removal

**Next steps**:

- Commit changes with message describing pipeline reorganization and stripe 2 results
- **CRITICAL**: Hear back from Jenny on edge spot removal methodology and rationale
- Contact Ali about SNFISH data processing methods and mRNA counting (high priority)
- Read Jenny's paper on Eve stripe analysis (before next meeting)
- Plan integration of SASS processing into Snakemake pipeline once methodology is documented

**Notes**:

- Pipeline architecture is scalable across stripes
- Full analysis chain (SASS mRNA processing → spatial binning → model fitting) can now be contained within single pipeline
- Edge spot removal is critical preprocessing step; understanding methodology is prerequisite for reproducibility
- Verify DV region captures full eve pattern

**Notes**:

- Large refactor with 10 files changed (+253, -559 lines)
- No commits made yet - all work staged/modified
- Exploring transcription data notebook currently open
- Generated results include degradation chains and summary statistics for dataset 033045
- I now have a programmatic way to get the peak range. Hopefully we can work out a similar technique for the data that Ali has. This way we can be consistent in our alignment for all three stripes.

---

## [2024-12-12]

**Session focus**: Major project restructure - transition from DVC to Snakemake workflow management

**Accomplished**:

- Migrated workflow from DVC to Snakemake (removed dvc.yaml, params.yaml, created Snakefile)
- Consolidated configuration into config.yaml with improved organization
- Refactored all 3 analysis scripts (01_preprocess, 02_infer, 03_visualize) to:
  - Use config.yaml instead of params.yaml
  - Remove DVC-specific dependencies
  - Simplify command-line interfaces
  - Improve code organization
- Added specialized analysis environment (envs/analysis.yml) with PyMC dependencies
- Staged new research infrastructure files:
  - Handoff documentation
  - Meeting transcripts and audio
  - Project telos updates
  - Activity logging structure
- Organized data files including transcription traces and Berrocal 2020 dataset
- Updated tasks.md with specific action items from handoff

**Decisions made**:

- Chose Snakemake over DVC for better flexibility and explicit workflow definition
- Separated environment dependencies (base environment.yml vs analysis envs/analysis.yml)
- Maintained original script structure but streamlined I/O handling
- Kept spatial binning parameters in config (5×5 grid, 25 bins)

**Next steps**:

- Test Snakemake workflow end-to-end
- Commit current changes with descriptive message
- Read Jenny's paper on Eve stripe analysis
- Contact Ali about SNFISH data and mRNA counting methodology

---

## File Naming Convention

Use consistent naming for all project files:

```
YYYYMMDD_experiment_condition_version.ext

Examples:
20241202_qpcr_treatment_v01.csv
20241202_preprocessing_output_v02.csv
20241202_figure1_draft.png
```

**Rules:**

- Use ISO date format (YYYYMMDD) for sortability
- Use underscores, never spaces
- Include version numbers for iterations
- Keep names descriptive but concise

---

## Digital Documentation Standards

### Core Principles

1. **Date everything** - ISO format (YYYY-MM-DD) for sorting
2. **Never overwrite** - Create new versions instead of modifying
3. **Cross-reference** - Link notebook entries to scripts and data files
4. **Back up regularly** - Follow 3-2-1 rule: 3 copies, 2 media types, 1 offsite

### Version Control Integration

- Commit at least daily when actively working
- Reference notebook entries in commit messages
- Tag versions before major analyses: `git tag -a v1.0 -m "Pre-analysis checkpoint"`

---

<!-- 
Add new entries at the top.
Keep entries brief but informative.
The RA will prompt you to update this after significant sessions.
-->
