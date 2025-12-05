---
name: monthly-review
description: Conduct a monthly review of project alignment and progress toward aims. Use when the user types /monthly_review, at the start of a new month, after completing a major milestone, or when uncertain about project direction.
---

# Monthly Review

> Conduct a monthly review of project alignment and progress toward aims.

## When to Use
- First week of each month (RA will prompt)
- After completing a major milestone
- When feeling uncertain about direction
- Before reporting to PI/collaborators

## Execution Steps

### 1. Gather Context

Read these files:
- `.research/project_telos.md` - Aims and goals
- `.research/phase_checklist.md` - Phase progress
- `.research/logs/weekly/*.md` - All weekly reviews this month
- `.research/logs/monthly/` - Previous monthly reviews
- `manuscript/` - Manuscript progress
- Git log for the month

### 2. Generate Monthly Review

```markdown
# Monthly Review: [MONTH YEAR]

## Executive Summary

**Project**: [Project name from project_telos.md]
**Current Phase**: [PHASE]
**Overall Health**: 🟢 On Track / 🟡 Needs Attention / 🔴 At Risk

## Progress Against Aims

### Aim 1: [Title]
- **Status**: [Not started / In progress / Complete]
- **Progress this month**: 
  - [Specific accomplishment]
  - [Specific accomplishment]
- **Blockers**: [None / List]
- **On track for completion?**: [Yes/No - by when]

### Aim 2: [Title]
- **Status**: 
- **Progress this month**: 
- **Blockers**: 
- **On track?**: 

### Aim 3: [Title]
[Continue pattern]

## Phase Progress

| Phase | Last Month | This Month | Change |
|-------|------------|------------|--------|
| SETUP | ✅ | ✅ | - |
| PLANNING | 80% | ✅ | Completed |
| DEVELOPMENT | 20% | 60% | +40% |
| ANALYSIS | 0% | 10% | Started |
| WRITING | 0% | 0% | - |
| REVIEW | 0% | 0% | - |

## Manuscript Status

| Section | Status | Last Updated |
|---------|--------|--------------|
| Background | [Draft/Revision/Complete] | [Date] |
| Methods | [Draft/Revision/Complete] | [Date] |
| Results | [Draft/Revision/Complete] | [Date] |
| Discussion | [Not started/Draft/Complete] | [Date] |

## Key Metrics

### Productivity Dashboard
| Metric | Target | Achieved | YTD Total |
|--------|--------|----------|----------|
| Weekly reviews completed | 4 | | |
| Active coding days | 15 | | |
| Commits | 20 | | |
| Scripts documented | | | |
| Figures generated | | | |
| Manuscript sections drafted | | | |

### Pipeline
- DVC stages defined: [N]
- Pipeline runs successfully: [Yes/No]
- Figures generated: [N]

### Off-Track Indicators

Check for these warning signs:

- [ ] Behind schedule by >2 weeks on milestones
- [ ] Same blocker appearing across multiple weeks
- [ ] Scope expanding without adjustment
- [ ] Key decisions being deferred
- [ ] Documentation falling behind

## Next Month Goals

### Must Complete
1. [Critical goal]
2. [Critical goal]

### Should Complete
1. [Important goal]
2. [Important goal]

### Stretch Goals
1. [Ambitious goal]

## Action Items

- [ ] [Specific action from this review]
- [ ] [Specific action]

---

*Review completed: [TIMESTAMP]*
*Next monthly review: [DATE]*
```

### 3. Save Review

Save to `.research/logs/monthly/YYYY-MM.md`

### 4. Suggest Next Steps

```
Monthly review saved!

Based on your progress:

A) Adjust project timeline or scope
   [If off-track indicators detected]

B) Update aims or methodology
   [If direction has shifted]

C) Plan next month's priorities
   Run /plan_week for immediate focus

D) Continue with current work
   You're on track!

What would you like to focus on?
```

## Related Skills

- `weekly-review` - Tactical progress
- `quarterly-review` - Strategic alignment
- `next` - Get next suggestion
