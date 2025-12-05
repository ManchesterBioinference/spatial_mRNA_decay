---
name: weekly-review
description: Conduct a weekly review of project progress and plan for the upcoming week. Use when the user types /weekly_review, when it's Monday and no review exists for this week, or at the end of the work week. Aggregates daily activity entries and identifies patterns.
---

# Weekly Review

> Conduct a weekly review of project progress and plan for the upcoming week.

## When to Use
- Every Monday (RA will prompt via /next)
- End of work week
- When feeling stuck or overwhelmed
- Before planning the next week

## Execution Steps

### 1. Gather Context

Read these files:
- `.research/project_telos.md` - Current aims and phase
- `.research/phase_checklist.md` - Phase progress
- `.research/logs/activity.md` - This week's activity
- `.research/logs/weekly/` - Previous reviews
- `tasks.md` - Pending tasks
- Git log for the past week

### 2. Generate Weekly Review

```markdown
# Weekly Review: [Week of DATE]

## What Happened This Week

### Accomplishments ✅
<!-- Auto-generated from git commits and activity.md -->
- [Accomplishment 1]
- [Accomplishment 2]
- [Accomplishment 3]

### Challenges Encountered
<!-- From context or ask user -->
- [Challenge 1]
- [Challenge 2]

### Time Spent
<!-- Estimate based on activity -->
- [Task category]: ~X hours
- [Task category]: ~X hours

## Project Health Check

### Phase Progress
- **Current Phase**: [PHASE]
- **Phase completion**: [X]% (based on checklist items)
- **On track?**: [Yes/No/At risk]

### Aims Status
| Aim | Status | Progress This Week |
|-----|--------|-------------------|
| Aim 1 | In progress | [Brief update] |
| Aim 2 | Not started | No change |

### Risks/Blockers
- [Current blockers from project_telos.md]

### Key Metrics
| Metric | Target | Actual | Status |
|--------|--------|--------|--------|
| Experiments/analyses completed | | | ✅/❌ |
| Scripts documented | | | ✅/❌ |
| Figures generated | | | ✅/❌ |
| Writing progress (words/sections) | | | ✅/❌ |
| Commits made | | | ✅/❌ |

## Patterns & Observations

### What Worked Well
- [Pattern to continue]

### What Could Improve
- [Area for improvement]

### Energy/Focus
<!-- Ask user or infer -->
- Peak productivity days: [Days]
- Challenging days: [Days]

## Looking Ahead

### Priorities for Next Week
1. [Top priority]
2. [Second priority]
3. [Third priority]

### Potential Blockers
- [Anticipated challenge]

### Questions to Resolve
- [Open question that needs answering]

---

*Review completed: [TIMESTAMP]*
*Next weekly review due: [DATE]*
```

### 3. Interactive Components

Ask the user:
```
Let's complete your weekly review. A few questions:

1. What are you most proud of accomplishing this week?
2. What was your biggest challenge or frustration?
3. Is there anything blocking your progress right now?
4. On a scale of 1-5, how focused did you feel this week?
5. What's the ONE thing that would make next week a success?
```

### 4. Save Review

Save to `.research/logs/weekly/YYYY-MM-DD.md`

### 5. Update Project State

- Update activity.md with review summary
- Flag any blockers in project_telos.md
- Update phase checklist if progress made

### 6. Suggest Next Steps

```
Weekly review saved to .research/logs/weekly/[DATE].md

Based on your review, suggested next steps:

A) 🗓️ Plan next week's work
   Run /plan_week to create a focused plan

B) 🚧 Address blockers
   [Specific blocker] needs attention

C) 📝 Continue with [top priority]
   Ready to dive in?

What would you like to do?
```

## Patterns to Watch For

### Warning Signs
Flag if detected:
- No commits in >7 days
- Same blocker appearing multiple weeks
- Phase not progressing
- Declining focus scores
- Growing task backlog
- **Motivation decline**: Reduced engagement, avoidance behaviors
- **Scope creep**: Project expanding without planned adjustment

### Recovery Strategies (When Off-Track)

If warning signs are present:

1. **Stop and assess** - Don't just push harder; pause to understand
2. **Identify root cause** - Is it technical, motivational, unclear direction, or external?
3. **Adjust scope/timeline** - Be realistic, not optimistic
4. **Communicate** - Tell collaborators, PI, or stakeholders early
5. **Document the decision** - Write down what changed and why
6. **Set clear checkpoints** - Define how you'll know if recovery is working

### Positive Patterns
Celebrate:
- Consistent commits
- Phase progression
- Blockers resolved
- Completed aims

## Related Skills

- `plan-week` - Create detailed weekly plan
- `monthly-review` - Zoom out for monthly perspective
- `next` - Get next suggestion

## Notes

- Keep reviews brief (15-20 minutes)
- Honesty > perfection
- Patterns emerge over time - review past weeks
- Celebrate progress, no matter how small
