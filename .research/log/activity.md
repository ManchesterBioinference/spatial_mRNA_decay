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