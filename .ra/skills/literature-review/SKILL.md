---
name: literature-review
description: Conduct comprehensive, systematic literature reviews using multiple databases (PubMed, bioRxiv, Semantic Scholar, OpenAlex). Creates documented searches, synthesizes findings thematically, verifies citations, and generates professional markdown reports with multiple citation styles (APA, Nature, Vancouver). Use when the user needs thorough literature research or types /deep_research.
---

# Systematic Literature Review

> Conduct rigorous, reproducible literature searches with verified citations.

## When to Use

- Conducting systematic literature review
- Synthesizing knowledge on a topic
- Writing the PLANNING phase literature review
- Populating `.research/literature/` directory
- Updating `manuscript/background.md`
- Meta-analysis or scoping reviews
- Finding prior work on techniques

## Core Workflow

```
1. PLAN         → Define question, scope, databases
2. SEARCH       → Systematic, documented searches
3. SCREEN       → Apply inclusion/exclusion criteria
4. EXTRACT      → Pull key data from sources
5. SYNTHESIZE   → Identify themes, patterns, gaps
6. VERIFY       → Validate all citations
7. DOCUMENT     → Create reproducible record
```

---

## Phase 1: Planning

### Define Research Question (PICO Framework)

For clinical/biomedical topics:
- **P**opulation: Who is studied?
- **I**ntervention: What treatment/exposure?
- **C**omparison: Versus what?
- **O**utcome: What is measured?

Example: "What is the efficacy of CRISPR-Cas9 (I) for treating sickle cell disease (P) compared to standard care (C)?"

### Develop Search Strategy

```markdown
## Search Strategy for: [Topic]

### Core Concepts
1. [Primary concept] - synonyms: [alternative terms]
2. [Secondary concept] - synonyms: [alternative terms]
3. [Methodological aspect]

### Search Queries
- "[concept1] AND [concept2]"
- "[method] in [application domain]"
- "review [topic]" (for overview papers)

### Databases to Search
- [ ] PubMed (biomedical)
- [ ] bioRxiv (preprints)
- [ ] arXiv (computational, physics)
- [ ] Semantic Scholar (cross-disciplinary)
- [ ] OpenAlex (open access)
- [ ] Google Scholar (broad coverage)

### Inclusion Criteria
- Date range: [start] to [end]
- Languages: [list]
- Study types: [list]

### Exclusion Criteria
- [criterion 1]
- [criterion 2]
```

---

## Phase 2: Systematic Search

### Documenting Searches

For each database, record:

```markdown
## Database: PubMed
- **Date searched**: YYYY-MM-DD
- **Date range**: [start] to [end]
- **Search string**: 
  ```
  ("CRISPR"[Title] OR "Cas9"[Title])
  AND ("sickle cell"[MeSH] OR "SCD"[Title/Abstract])
  ```
- **Results**: N articles
- **Notes**: [any filter applied]
```

### PubMed Search Tips

- Use MeSH terms: `"disease name"[MeSH]`
- Field tags: `[Title]`, `[Title/Abstract]`, `[Author]`
- Boolean: AND, OR, NOT
- Date filter: `2020:2024[Publication Date]`
- Article type: `"Review"[Publication Type]`

### bioRxiv/medRxiv

- Preprints - NOT peer-reviewed
- Check if subsequently published
- Note version and date

### Semantic Scholar / OpenAlex

- Good for citation networks
- Cross-disciplinary coverage
- Can find highly-cited papers

---

## Phase 3: Screening

### PRISMA Flow

```
Records identified (n = X)
├─ Duplicates removed (n = Y)
├─ Records screened (n = Z)
│  └─ Excluded by title/abstract (n = A)
├─ Full-text assessed (n = B)
│  └─ Excluded with reasons (n = C)
└─ Studies included (n = D)
```

### Screening Checklist

For each paper:
- [ ] Meets inclusion criteria
- [ ] Doesn't meet exclusion criteria
- [ ] Relevant to research question
- [ ] Sufficient quality (peer-reviewed, appropriate methods)

---

## Phase 4: Data Extraction

### Extract from Each Source

1. **Metadata**: Authors, year, journal, DOI
2. **Study design**: Methods, population, sample size
3. **Key findings**: Main results, effect sizes
4. **Limitations**: Noted by authors
5. **Relevance**: How it relates to your question
6. **Quality**: Assessment of rigor

### Literature Note Template

Create in `.research/literature/[topic].md`:

```markdown
# Literature Review: [Topic]

Date: YYYY-MM-DD
Search strategy: [link to search doc]

## Summary
[2-3 sentence overview of the literature landscape]

## Key Themes

### Theme 1: [Name]
[Synthesis across multiple papers - NOT paper-by-paper]

Key points:
- Point 1 (Author1, Year; Author2, Year)
- Point 2 (Author3, Year)

### Theme 2: [Name]
...

## Gaps and Controversies
- [Gap 1]: Limited research on...
- [Controversy]: Conflicting results regarding...

## Key Papers

### Author (Year) - Title
- **DOI**: 10.xxxx/xxxxx
- **Main finding**: [One sentence]
- **Relevance**: [How it relates to your work]
- **Limitations**: [Key caveats]

## References
[BibTeX or formatted references]
```

---

## Phase 5: Thematic Synthesis

### Do This ✓

Organize by themes, not by paper:
```
Multiple delivery approaches have been investigated. Viral vectors 
(AAV) showed high efficiency (65-85%)¹⁻¹⁵ but raised immunogenicity 
concerns³,⁷. Lipid nanoparticles demonstrated lower efficiency 
(40-60%) but improved safety¹⁶⁻²³.
```

### Don't Do This ✗

Paper-by-paper summary:
```
Smith (2020) found that viral vectors work well. Jones (2021) 
also used viral vectors. Wang (2022) tried nanoparticles...
```

### Synthesis Template

```markdown
## Thematic Synthesis

### Theme: [Name]

**Current understanding**: [Summary of what is known]

**Consensus areas**: [What most studies agree on]

**Controversies**: [Conflicting findings]

**Gaps**: [What's missing from the literature]

**Evidence quality**: [Overall assessment]
```

---

## Phase 6: Citation Verification

### CRITICAL: Verify All Citations

**NEVER fabricate citations.** Every claim must have a verifiable source.

### Verification Checklist

For each citation:
- [ ] DOI resolves correctly
- [ ] Author names match
- [ ] Year is correct
- [ ] Title matches
- [ ] Journal name correct

### Citation Styles

#### APA (7th Edition)
```
Smith, J. D., Johnson, M. L., & Williams, K. R. (2023). Title of article. 
Journal Name, 22(4), 301-318. https://doi.org/10.xxx/yyy
```

#### Nature
```
Smith, J. D., Johnson, M. L. & Williams, K. R. Title of article. 
Nat. Rev. Drug Discov. 22, 301-318 (2023).
```

#### Vancouver
```
Smith JD, Johnson ML, Williams KR. Title of article. 
Nat Rev Drug Discov. 2023;22(4):301-18.
```

---

## Phase 7: Documentation

### Output Files

| File | Location | Purpose |
|------|----------|---------|
| Search strategy | `.research/literature/search_[topic].md` | Reproducibility |
| Literature notes | `.research/literature/[topic].md` | Findings |
| BibTeX file | `.research/literature/[topic].bib` | Citations |
| PRISMA diagram | `.research/literature/prisma_[topic].md` | Flow |

### BibTeX Entry Template

```bibtex
@article{smith2023,
  author = {Smith, John D. and Johnson, Mary L. and Williams, Karen R.},
  title = {Title of the Article},
  journal = {Journal Name},
  year = {2023},
  volume = {22},
  number = {4},
  pages = {301--318},
  doi = {10.xxx/yyy},
  pmid = {12345678}
}
```

---

## Integration with RA Workflow

### PLANNING Phase Connection

Literature review is a gate for entering DEVELOPMENT:
- [ ] At least one literature search completed → ✓ Creates `.research/literature/` content
- [ ] background.md has at least a rough draft → Uses literature synthesis

### After Literature Review

1. Update `.research/literature/` with findings
2. Draft sections in `manuscript/background.md`
3. Add citations to `.research/literature/*.bib`
4. Log activity in `.research/logs/activity.md`
5. Update `tasks.md` with identified gaps

### Connection to Other Skills

- **→ `/write_background`**: Uses literature synthesis to draft background
- **→ `/hypothesis_generation`**: Literature informs hypothesis development
- **← `/deep_research`**: Alternative invocation for same functionality

---

## Quality Checklist

Before completing a literature review:

- [ ] Search strategy documented with all parameters
- [ ] Multiple databases searched (minimum 3)
- [ ] Inclusion/exclusion criteria clearly stated
- [ ] PRISMA flow documented
- [ ] Findings synthesized thematically (not paper-by-paper)
- [ ] All DOIs verified
- [ ] Citations formatted consistently
- [ ] Gaps and controversies identified
- [ ] Key papers annotated
- [ ] BibTeX file created

---

## Common Pitfalls

| Pitfall | Problem | Solution |
|---------|---------|----------|
| Single database | Misses relevant papers | Search 3+ databases |
| No documentation | Not reproducible | Record all search parameters |
| Paper-by-paper | No synthesis | Organize thematically |
| Unverified citations | Errors and fabrication | Check every DOI |
| Too broad | Thousands of results | Refine with specific terms |
| Too narrow | Misses relevant work | Include synonyms |
| Ignoring preprints | Misses latest findings | Include bioRxiv/medRxiv |
| No quality assessment | Treats all evidence equally | Assess study quality |
