# Project Telos

> This file is the source of truth for your project's purpose, aims, and current state.
> The Research Assistant (RA) reads this file to understand where you are and what you're trying to accomplish.

## Project Mission
<!-- One sentence: What question are you answering? What problem are you solving? 
-->

How does mRNA degradation vary across space in developing Drosophila embryos, and does this spatial variation contribute to pattern formation in gap genes like *eve*?

## Context
<!-- Is this part of a larger grant? Related to other projects? 
-->

**Parent Grant/Project**: Part of larger manuscript on *eve* stripe formation and gene regulation in Drosophila
**Specific Aim Addressed**: Understanding spatial/temporal aspects of mRNA dynamics; paper being revised based on reviewer comments
**Target Output**: Journal paper (currently in revision; withdrawn to address substantial reviewer feedback)
**Target Venue**: Related to original submission; resubmitting after revision
**Collaborators**: 
- **Magnus Rattray** (PI/supervisor) - Conceptual guidance, statistical methods
- **Jenny Pantin** (previous postdoc) - Developed core analysis pipeline, STAN implementation
- **Hilary** - Experimental work, MS2 live cell imaging, SNFISH experiments
- **Ali** - Data processing, SNFISH methodology, alignment expertise

## Specific Aims

<!-- 
Define 1-3 specific aims. Each aim should be testable and lead to concrete outputs.
Use the format below for each aim.
-->

### Aim 1: Quantify spatially-varying mRNA degradation rates in *eve* stripe (COMPLETED by Jenny)

- **Hypothesis**: mRNA degradation is higher between stripes and lower within stripes, enabling sharp stripe boundaries in *eve* expression patterns
- **Approach**: 
  1. Combine MS2 live cell imaging data (temporal dynamics of pre-mRNA) from Mike Eisen's lab with SNFISH data (mature mRNA snapshot) from our lab
  2. Align the two datasets spatially across the anterior-posterior axis
  3. Fit a mathematical model where degradation varies across 5 spatial regions (stripes + inter-stripe gaps)
  4. Use Bayesian inference (STAN) to estimate posterior distributions of degradation parameters
- **Success Criteria**: 
  - Successfully quantify degradation rates for one *eve* stripe
  - Demonstrate clear spatial pattern: higher degradation between stripes, lower inside
  - Validate that combining MS2 and SNFISH data is feasible
- **Status**: ✓ **COMPLETED** (by Jenny Pantin)
- **Key Scripts**: Jenny's original analysis code (to be documented and refactored in Aim 2)
- **Key Figures**: [Spatial degradation maps for stripe 1; degradation parameter estimates]
- **Notes**: 
  - Core model: dm/dt = s*p(t) - δ*m, where p(t) is pre-mRNA from MS2, m is mature mRNA from SNFISH
  - Eve has fast degradation (~5-10 min half-life)
  - Potential weakness: SNFISH and MS2 from different experiments; alignment quality critical (raised by reviewers)
  - New approach (vs. previous Gaussian process work): MS2 temporal resolution is fine-grained, so can sum bins directly instead of regressing between sparse time points
  - This work forms the foundation for paper currently being revised


### Aim 2: Reproduce Jenny's analysis in a reproducible DVC/Snakemake pipeline

- **Hypothesis**: A well-documented DVC/Snakemake pipeline will enable exact reproduction of Jenny's results and facilitate future extensions to other stripes and genes
- **Approach**: 
  1. Understand Jenny's existing analysis code and methodology
  2. Document all analysis steps with clear docstrings
  3. Structure as DVC or Snakemake pipeline stages: data alignment → spatial binning → STAN model fitting → posterior inference
  4. Parameterize all key choices (e.g., bin size, STAN parameters) in params.yaml
  5. Validate that pipeline reproduces Jenny's original results on stripe 1 data
- **Success Criteria**: 
  - All scripts have clear docstrings and are version-controlled
  - dvc.yaml or Snakefile captures full pipeline with clear dependency structure
  - params.yaml documents all tunable parameters
  - Running `dvc repro` or `snakemake` produces results matching Jenny's original analysis
  - Methods section accurately reflects current implementation
- **Status**: **IN PROGRESS** (current focus)
- **Key Scripts**: [Alignment, binning, STAN fitting, posterior visualization - to be refactored]
- **Key Figures**: [Reproducibility validation: comparison of pipeline output vs. Jenny's original results]
- **Notes**: 
  - Key challenge: understanding and documenting spatial alignment between SNFISH and MS2 (needs Ali's expertise)
  - This work is also training ground for Pontus' master's project and testing new research workflow tools
  - Must decide between DVC vs. Snakemake based on project needs


### Aim 3: Extend pipeline to additional *eve* stripes (3 and 4)

- **Hypothesis**: The reproducible pipeline from Aim 2 will enable robust estimation of degradation rates across multiple stripe replicates, demonstrating that spatial degradation patterns are consistent features of *eve* regulation
- **Approach**: 
  1. Use SNFISH data from Hilary for stripes 3 and 4
  2. Apply validated pipeline from Aim 2 to new stripe data
  3. Compare degradation parameter estimates across all three stripes (1, 3, 4)
  4. Assess robustness and reproducibility of spatial degradation patterns
- **Success Criteria**: 
  - Successfully process and analyze stripes 3 and 4 using pipeline from Aim 2
  - Demonstrate that degradation estimates are robust across stripe replicates
  - Show consistent spatial pattern across all three stripes: higher degradation between stripes, lower inside
  - Generate comparative figures showing degradation patterns across multiple stripes
- **Status**: **PLANNED** (waiting on completion of Aim 2)
- **Key Scripts**: Same pipeline as Aim 2, applied to new data
- **Key Figures**: [Multi-stripe comparison: degradation maps across stripes 1, 3, 4; statistical comparison of parameters]
- **Notes**: 
  - Timeline depends on Hilary's experimental schedule
  - Originally planned stripes 1, 5, 6 but changed to 1, 3, 4 based on Hilary's timeline
  - This will provide the replication needed to strengthen paper for resubmission

## Current Phase

<!-- Check the current phase. The RA uses this to guide suggestions. 
-->

- [x] **SETUP** - Environment, structure, git configured
- [x] **PLANNING** - Aims defined, literature reviewed, background drafted
- [ ] **DEVELOPMENT** - Pipeline being built, scripts documented
- [ ] **ANALYSIS** - Experiments running, results being generated
- [ ] **WRITING** - Drafting manuscript sections
- [ ] **REVIEW** - Polishing, reproducibility verification, submission prep

## Goals (Prioritized)

<!-- 
List current goals in priority order. G1 is highest priority.
Update these regularly as you progress.
-->

- G1: Reproduce Jenny's analysis on original eve stripe data to understand the pipeline
- G2: Obtain and integrate Hillary's new SNFISH data (stripes 3 and 4) and extend analysis
- G3: Work with Ali to understand and document the spatial alignment methodology
- G4: Build a complete reproducible DVC pipeline for multi-stripe analysis
- G5: Mentor Pontus on technical implementation while testing research workflow tools

## Risks / Blockers

<!-- 
What's blocking progress? What are you worried about?
The RA can help address these.
-->

- R1: Understanding existing analysis - Jenny's code and STAN implementation may have a learning curve; need to read paper and understand methods thoroughly first
- R2: Data access and preprocessing - Need to coordinate with Ali on SNFISH counting methodology and with Hilary on timeline for new stripe experiments (stripes 3 & 4)
- R3: Cross-experiment alignment quality - SNFISH and MS2 from different experiments; alignment is critical but noted as weak point by reviewers; need to validate alignment approaches
- R4: Scaling complexity - Current discrete spatial binning may not scale; future PDE-based approaches (Aim 2) may require different mathematical framework

## Key Decisions Log

<!-- 
Record important methodological or direction decisions here.
These inform the methods section and help maintain project memory.
-->

| Date | Decision | Rationale |
|------|----------|-----------|
| 2025-12-05 | Use STAN for Bayesian inference instead of custom MCMC | Cleaner inference without tedious tuning; allows posterior distributions of degradation parameters |
| 2025-12-05 | Combine MS2 (temporal) with SNFISH (snapshot) data from different experiments | Allows linking pre-mRNA dynamics to mature mRNA counts; Drosophila reproducibility justifies cross-experiment alignment |
| 2025-12-05 | Model degradation as piecewise constant across spatial regions (~5 bins) | Simpler first approach; future extension to continuous fields (PDEs) possible |
| 2025-12-05 | Withdraw paper and resubmit as new submission instead of revised submission | Allows addressing substantial reviewer comments without back-and-forth; likely same reviewers but cleaner process |
| 2025-12-05 | Focus on 3 stripe replicates (1, 3, 4 instead of 1, 5, 6) | Hillary's experimental timeline; provides replication for robustness |
| 2026-02-02 | No edge removal in data processing; use spot counts without aggregating by nuclei first for bin mean calculation | Both edge removal approaches (all edges vs only left/right) introduced biases; spots are smoother without artificial removal |

## Activity Log

<!-- 
Running log of significant activities. The RA may prompt you to update this.
Newest entries at top.
-->

### [YYYY-MM-DD]
- [What you worked on]
- [What was accomplished]
- [What's next]

---

*Last updated: [Date]*
*Current phase: [PHASE]*
