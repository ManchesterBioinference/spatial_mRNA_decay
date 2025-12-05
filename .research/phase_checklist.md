# Phase Checklist

> This file tracks your progress through research phases.
> The RA uses this to ensure you don't skip critical steps.

## Current Phase: SETUP

---

## SETUP Phase

**Goal**: Establish a solid, reproducible foundation for your project.

### Required (Must complete before PLANNING)

- [ ] Git repository initialized
- [ ] Remote repository connected (GitHub/GitLab)
- [ ] Environment manager chosen and configured
  - [ ] pyproject.toml / environment.yml / requirements.txt exists
- [ ] Directory structure created
- [ ] .gitignore configured
- [ ] .copilotignore configured (sensitive data excluded)
- [ ] DVC initialized (if using large data files)

### Recommended

- [ ] README.md has project overview
- [ ] Initial commit made
- [ ] Collaborator access configured (if applicable)

### Exit Criteria
All required items checked → Ready for PLANNING

---

## PLANNING Phase

**Goal**: Define clear aims and understand the field before building.

### Required (Must complete before DEVELOPMENT)

- [ ] Project mission defined in project_telos.md
- [ ] At least 1 specific aim defined with hypothesis
- [ ] Literature search completed on core concepts
  - [ ] .research/literature/ has at least one .md and .bib file
- [ ] Background section drafted (rough is fine)
  - [ ] manuscript/background.md has content
- [ ] Key methodological approach identified

### Recommended

- [ ] 2-3 aims defined
- [ ] Target venue/journal identified
- [ ] Collaborator roles clarified
- [ ] Timeline estimated

### Exit Criteria
Project aims are clear, you understand the field, background draft exists → Ready for DEVELOPMENT

---

## DEVELOPMENT Phase

**Goal**: Build your computational pipeline with documentation.

### Required (Must complete before ANALYSIS)

- [ ] At least one script in pipeline/scripts/
- [ ] All scripts have docstrings
- [ ] dvc.yaml has at least one stage defined
- [ ] params.yaml has configurable parameters
- [ ] Data acquisition method documented
- [ ] Methods section reflects current scripts
  - [ ] manuscript/methods.md updated

### Recommended

- [ ] Pipeline runs end-to-end (dvc repro works)
- [ ] Unit tests for critical functions
- [ ] README in pipeline/scripts/ explains workflow
- [ ] Regular commits with descriptive messages

### Exit Criteria
Pipeline is functional and documented → Ready for ANALYSIS

---

## ANALYSIS Phase

**Goal**: Run experiments and generate results.

### Required (Must complete before WRITING)

- [ ] Full pipeline executed successfully
- [ ] At least one figure generated
- [ ] Figure captions drafted
  - [ ] manuscript/figures/fig1/caption.md exists for each figure
- [ ] Key results identified and noted
- [ ] Methods section is complete and accurate

### Recommended

- [ ] Multiple figures for different aims
- [ ] Statistical tests run and results recorded
- [ ] Sensitivity analyses completed (if applicable)
- [ ] Results replicate on clean run (dvc repro from scratch)

### Exit Criteria
Figures exist with captions, methods are complete → Ready for WRITING

---

## WRITING Phase

**Goal**: Draft all manuscript sections.

### Required (Must complete before REVIEW)

- [ ] Background section complete (not just draft)
- [ ] Methods section complete and matches code
- [ ] Results section drafted from figures
- [ ] Discussion section drafted
- [ ] All figures have finalized captions
- [ ] Abstract drafted

### Recommended

- [ ] Co-authors have reviewed sections
- [ ] References formatted correctly
- [ ] Supplementary materials organized
- [ ] Target journal guidelines reviewed

### Exit Criteria
All sections drafted → Ready for REVIEW

---

## REVIEW Phase

**Goal**: Polish, verify reproducibility, prepare for submission.

### Checklist

- [ ] Full reproducibility test (clone repo, run pipeline, verify outputs)
- [ ] All figures regenerate correctly
- [ ] Methods match code exactly
- [ ] References verified (no hallucinated citations)
- [ ] Co-author approval obtained
- [ ] README has full reproduction instructions
- [ ] Code cleaned and commented
- [ ] Data availability statement prepared
- [ ] Submission materials gathered

### Journal-Specific

- [ ] Word/character limits checked
- [ ] Figure format requirements met
- [ ] Cover letter drafted
- [ ] Suggested reviewers identified (if required)

### Exit Criteria
Everything verified and polished → SUBMIT!

---

## Progress Summary

| Phase | Status | Started | Completed |
|-------|--------|---------|-----------|
| SETUP | ⬜ Not started | | |
| PLANNING | ⬜ Not started | | |
| DEVELOPMENT | ⬜ Not started | | |
| ANALYSIS | ⬜ Not started | | |
| WRITING | ⬜ Not started | | |
| REVIEW | ⬜ Not started | | |

Legend: ⬜ Not started | 🟨 In progress | ✅ Complete

---

*Last updated: [Date]*
