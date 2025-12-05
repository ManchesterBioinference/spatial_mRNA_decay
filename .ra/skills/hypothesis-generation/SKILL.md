---
name: hypothesis-generation
description: Generate testable scientific hypotheses from observations, data, or research questions. Develop competing explanations, design experiments, and formulate predictions. Use during the PLANNING phase or when developing research aims.
---

# Scientific Hypothesis Generation

> Systematically develop testable explanations and experimental designs.

## When to Use

- Developing hypotheses from observations or preliminary data
- Filling out `.research/project_telos.md` aims section
- During the PLANNING phase
- Exploring competing explanations for phenomena
- Designing experiments to test research questions
- Generating predictions for research proposals

## Workflow

```
1. UNDERSTAND    → Clarify the phenomenon/question
2. RESEARCH      → Survey existing literature
3. SYNTHESIZE    → Integrate evidence and identify gaps
4. GENERATE      → Develop 3-5 competing hypotheses
5. EVALUATE      → Assess hypothesis quality
6. DESIGN        → Plan experimental tests
7. PREDICT       → Formulate testable predictions
```

---

## Step 1: Understand the Phenomenon

### Clarifying Questions

Ask these to define the research question:

1. **What is the core observation or pattern?**
   - What specifically needs to be explained?
   
2. **What is the scope?**
   - What are the boundaries of the phenomenon?
   - What is included/excluded?

3. **What is already known?**
   - Established facts vs. assumptions
   - Previous attempts at explanation

4. **What are the constraints?**
   - Methodological limitations
   - Time/resource constraints
   - Ethical considerations

### Example

```markdown
## Phenomenon Definition

**Observation**: Treatment X reduces tumor growth in mice, but only in 
animals with intact immune systems.

**Scope**: Focus on solid tumors in murine models. Excludes metastasis 
and blood cancers.

**Known**: Treatment X has no direct cytotoxic effect on cancer cells 
in vitro.

**Question**: How does Treatment X reduce tumor growth in an 
immune-dependent manner?
```

---

## Step 2: Literature Research

Before generating hypotheses, ground them in existing evidence.

### Search Strategy

1. **Identify key concepts** from the research question
2. **List synonyms** and related terms
3. **Search systematically** (use `/deep_research` if needed)
4. **Document findings** for later reference

### Integration with RA

Use `/deep_research` to gather literature:

```
/deep_research Treatment X mechanism of action cancer immunity
```

This creates documented summaries in `.research/literature/`

---

## Step 3: Synthesize Existing Evidence

### Evidence Summary Template

```markdown
## Literature Synthesis

### What is Established
1. [Established fact 1 with citation]
2. [Established fact 2 with citation]

### Current Theories
1. [Theory A]: [Brief description] (Supporting: [refs], Contradicting: [refs])
2. [Theory B]: [Brief description]

### Knowledge Gaps
1. [Gap 1]: No studies have examined...
2. [Gap 2]: Conflicting results regarding...

### Relevant Mechanisms from Related Systems
1. [Analogous system 1]: [What can be learned]
2. [Analogous system 2]: [Potential parallel]
```

---

## Step 4: Generate Competing Hypotheses

Develop 3-5 distinct hypotheses that could explain the phenomenon.

### Hypothesis Requirements

Each hypothesis must be:
- **Mechanistic**: Explains *how* and *why*, not just *what*
- **Testable**: Can be evaluated empirically
- **Distinguishable**: Different from other hypotheses in testable ways
- **Evidence-based**: Grounded in existing knowledge

### Strategies for Generation

| Strategy | Description | Example |
|----------|-------------|---------|
| **Analogical** | Apply mechanisms from similar systems | "Similar to how X works in Y system" |
| **Mechanistic decomposition** | Break down into component processes | "Step 1 leads to Step 2 which causes..." |
| **Level shifting** | Consider different scales | "At the molecular level..." vs "At the systems level..." |
| **Contradiction exploration** | What if the opposite were true? | "What if X inhibits rather than activates?" |
| **Integration** | Combine known mechanisms in new ways | "If A and B act together..." |

### Example Hypotheses

```markdown
## Competing Hypotheses

### H1: T-cell Activation Hypothesis
Treatment X enhances T-cell activation through direct binding to 
checkpoint receptors, leading to increased tumor infiltration and 
cytotoxicity.

*Mechanism*: X → checkpoint binding → T-cell activation → tumor killing
*Key prediction*: T-cell depletion would abolish the effect

### H2: Dendritic Cell Priming Hypothesis
Treatment X stimulates dendritic cell maturation and antigen presentation,
leading to enhanced adaptive immune response against tumor antigens.

*Mechanism*: X → DC maturation → antigen presentation → T-cell priming
*Key prediction*: Effect would require intact antigen presentation

### H3: Tumor Microenvironment Remodeling Hypothesis
Treatment X alters the immunosuppressive tumor microenvironment by
depleting regulatory T cells or MDSCs, allowing existing immune 
responses to act.

*Mechanism*: X → Treg/MDSC depletion → reduced immunosuppression
*Key prediction*: Effect correlates with reduction in immunosuppressive cells

### H4: Innate Immune Activation Hypothesis
Treatment X activates innate immune cells (NK cells, macrophages) that
directly kill tumor cells and enhance adaptive immunity.

*Mechanism*: X → innate activation → tumor killing + cytokine release
*Key prediction*: Early innate response precedes adaptive response
```

---

## Step 5: Evaluate Hypothesis Quality

### Quality Criteria

| Criterion | Question | Score (1-5) |
|-----------|----------|-------------|
| **Testability** | Can it be empirically tested? | |
| **Falsifiability** | What would disprove it? | |
| **Parsimony** | Is it the simplest explanation? | |
| **Explanatory power** | How much does it explain? | |
| **Scope** | What range of observations does it cover? | |
| **Consistency** | Does it fit established principles? | |
| **Novelty** | Does it offer new insights? | |

### Evaluation Template

```markdown
## Hypothesis Evaluation

| Hypothesis | Testability | Falsifiability | Parsimony | Explanatory Power | Priority |
|------------|-------------|----------------|-----------|-------------------|----------|
| H1: T-cell | 5 | 5 | 4 | 4 | **High** |
| H2: DC | 4 | 4 | 3 | 4 | Medium |
| H3: TME | 5 | 5 | 4 | 3 | Medium |
| H4: Innate | 4 | 4 | 4 | 3 | Low |

**Strongest hypothesis**: H1 (most testable, clear predictions)
**Alternative to test**: H3 (could explain H1 results)
```

---

## Step 6: Design Experimental Tests

### For Each Hypothesis, Define:

1. **Key experiment**: The critical test
2. **Controls**: What comparisons are needed
3. **Methods**: How would you measure outcomes
4. **Expected results**: If hypothesis is correct
5. **Alternative outcomes**: What other results could mean

### Experimental Design Template

```markdown
## Experimental Design: Testing H1

### Critical Experiment
Deplete CD8+ T cells using anti-CD8 antibodies before Treatment X 
administration.

### Experimental Groups
1. Treatment X + isotype control antibody (n=10)
2. Treatment X + anti-CD8 antibody (n=10)
3. Vehicle + isotype control (n=10)
4. Vehicle + anti-CD8 antibody (n=10)

### Primary Outcome
Tumor volume at day 14 post-treatment

### Expected Results (if H1 correct)
- Group 1 shows reduced tumor growth
- Group 2 shows NO reduction (similar to Group 3)
- This would demonstrate T-cell dependence

### Alternative Outcomes
- If Group 2 still shows reduction → Effect is T-cell independent
  → Support for H3 or H4
- If Group 4 shows increased growth → T cells contribute to 
  baseline control → Consider immunocompetent models
```

---

## Step 7: Formulate Testable Predictions

### Prediction Requirements

Good predictions are:
- **Specific**: Clear, measurable outcomes
- **Quantitative** (when possible): Expected magnitude or direction
- **Conditional**: Specify under what conditions
- **Distinguishing**: Differentiate between hypotheses

### Prediction Format

```markdown
## Predictions

### If H1 is correct:
1. CD8+ T cell infiltration will increase >2-fold after Treatment X
2. T cell depletion will abolish >80% of tumor reduction effect
3. PD-1/PD-L1 blockade will enhance Treatment X efficacy synergistically

### If H2 is correct:
1. Dendritic cell maturation markers (CD80, CD86) will increase
2. Antigen presentation blockade will eliminate the effect
3. The effect will require 7+ days (time for adaptive response)

### Discriminating Predictions:
- H1 predicts rapid effect (days); H2 predicts delayed effect (weeks)
- H1 predicts T cell depletion is sufficient; H2 predicts DC depletion 
  is also required
```

---

## Integration with RA Workflow

### Output to project_telos.md

The generated hypotheses should be added to `.research/project_telos.md`:

```markdown
### Aim 1: Determine the mechanism of Treatment X efficacy

- **Hypothesis**: Treatment X enhances T-cell activation through 
  checkpoint receptor binding, leading to increased tumor infiltration 
  and cytotoxicity. (H1 from hypothesis generation)
- **Alternative**: The effect may be mediated by tumor microenvironment 
  remodeling (H3).
- **Approach**: T-cell depletion experiments followed by immune profiling
- **Success Criteria**: Identify the critical immune cell population 
  required for Treatment X efficacy
- **Status**: Not started
```

### Phase Gate Contribution

This skill helps complete the PLANNING phase requirements:
- [ ] Project aims are defined ← Hypothesis generation contributes here
- [ ] At least one literature search completed
- [ ] background.md has at least a rough draft

---

## Hypothesis Generation Checklist

- [ ] Phenomenon clearly defined and bounded
- [ ] Literature searched and synthesized
- [ ] 3-5 competing hypotheses generated
- [ ] All hypotheses are mechanistic and testable
- [ ] Quality evaluation completed
- [ ] Experimental designs outlined
- [ ] Predictions formulated and distinguish between hypotheses
- [ ] Added to project_telos.md
