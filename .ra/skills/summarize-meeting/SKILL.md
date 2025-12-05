---
name: summarize-meeting
description: Extract action items, decisions, and key points from meeting transcripts. Automatically routes items to tasks.md or GitHub Issues based on complexity. Use when the user types /summarize_meeting or after running /transcribe.
---

# Meeting Summary Generator

> Extract action items, decisions, and key points from meeting transcripts.
> Automatically routes items to tasks.md or GitHub Issues.

## Usage
```
/summarize_meeting [transcript_file]
/summarize_meeting .research/meetings/transcripts/2024-12-02-lab-meeting.md
```

## When to Use
- After running /transcribe
- On any meeting transcript or notes
- To process handwritten meeting notes (type them first)

## Prerequisites
- Transcript exists in `.research/meetings/transcripts/` folder
- Transcript is in markdown format

## Execution Steps

### 1. Load Transcript

Read the meeting transcript and project context:
- `.research/meetings/transcripts/[filename].md` - The transcript
- `.research/project_telos.md` - Project aims (for context)
- `tasks.md` - Current tasks (avoid duplicates)

### 2. Extract Key Information

Analyze transcript for:
- **Decisions made**
- **Action items** (who, what, when)
- **Questions raised**
- **Key insights or ideas**
- **Follow-up needed**

### 3. Generate Meeting Summary

Append summary to transcript or create separate file:

```markdown
---

# Meeting Summary

## Key Decisions
<!-- Decisions that were made during the meeting -->
1. [Decision 1] - [Brief context]
2. [Decision 2] - [Brief context]

## Action Items

### Tasks (< 2 hours, single implementation)
<!-- These will be added to tasks.md -->

| Item | Owner | Due | Priority |
|------|-------|-----|----------|
| [Task description] | [Name/@you] | [Date/ASAP] | [High/Med/Low] |
| [Task description] | [Name/@you] | [Date] | [Priority] |

### Issues (> 2 hours, needs tracking)
<!-- These will become GitHub Issues -->

1. **[Issue title]**
   - Description: [What needs to be done]
   - Why: [Why this is needed]
   - Complexity: [Estimate]
   - Labels: [Suggested labels]

2. **[Issue title]**
   - Description: [Details]
   - Why: [Rationale]

## Key Insights
<!-- Important points or ideas worth remembering -->
- [Insight 1]
- [Insight 2]

## Open Questions
<!-- Questions that weren't resolved -->
- [ ] [Question 1] - Needs: [Who/what to resolve]
- [ ] [Question 2] - Needs: [Who/what to resolve]

## Follow-up Needed
<!-- Things to discuss or check on later -->
- [Follow-up item]

## Next Meeting
<!-- If discussed -->
- **Date**: [If scheduled]
- **Agenda items**: [If mentioned]

---

*Summary generated: [Timestamp]*
```

### 4. Task vs Issue Classification

Apply this heuristic:

| Criteria | → Task (tasks.md) | → Issue (GitHub) |
|----------|-------------------|------------------|
| Estimated time | < 2 hours | > 2 hours |
| Scope | Single action | Multiple steps |
| Branching | Not needed | Needs own branch |
| Comparison | No | Comparing alternatives |
| Documentation | Not needed | Should be tracked |
| Project direction | Doesn't change | May change direction |

**When uncertain, ask:**
```
I found this action item: "[Item description]"

This could be:
A) A quick task (< 2 hours, add to tasks.md)
B) A larger issue (needs GitHub Issue for tracking)

Which fits better? (Or provide more context)
```

### 5. Update tasks.md

Add new tasks with meeting reference:

```markdown
## From Meeting: 2024-12-02-lab-meeting

- [ ] [Task 1] (Due: [date])
- [ ] [Task 2] (Due: [date])
- [ ] [Task 3]
```

### 6. Create GitHub Issues

For items classified as issues, offer to create:

```
I identified 2 items that should be GitHub Issues:

1. "Compare SVM vs Random Forest for classification"
   - Would require testing both approaches
   - Results should be documented for paper
   
2. "Implement alternative normalization method"
   - Needs research into options
   - May change downstream pipeline

Create these as GitHub Issues? (Y/n)
```

If yes, create issues with:
- Clear title
- Description from meeting context
- Labels (if determinable)
- Reference to meeting transcript

### 7. Post-Summary Actions

```
Meeting summarized!

Summary added to: .research/meetings/transcripts/2024-12-02-lab-meeting.md
Tasks added: 3 new items in tasks.md
Issues to create: 2 (awaiting confirmation)

Next steps:
A) Create the GitHub Issues
B) Review and prioritize new tasks
C) Update project_telos.md with decisions made
D) Continue with other work

What would you like to do?
```

## Example Output

```markdown
# Meeting Summary

## Key Decisions
1. Use random forest as primary classifier (SVM as comparison)
2. Deadline for analysis: end of month
3. Weekly check-ins moving to Tuesdays

## Action Items

### Tasks
| Item | Owner | Due | Priority |
|------|-------|-----|----------|
| Fix axis labels on Figure 2 | @you | Dec 4 | Low |
| Send PI the draft methods section | @you | Dec 3 | High |
| Update README with new instructions | @you | Dec 5 | Med |

### Issues
1. **Compare SVM vs Random Forest performance**
   - Description: Run both classifiers with same CV scheme, compare metrics
   - Why: Reviewer may ask about method choice
   - Complexity: ~4-6 hours
   - Labels: analysis, methodology

## Key Insights
- PI suggested looking at recent paper by Smith et al. on normalization
- Feature importance might be more interesting than just accuracy

## Open Questions
- [ ] Which normalization method should we use? - Needs: literature review
- [ ] Include supplementary figures in main text? - Needs: check journal guidelines
```

## Related Skills

- `transcribe` - Generate transcript from audio
- `weekly-review` - See tasks in context of weekly work
- `plan-week` - Incorporate new tasks into weekly plan
- `next` - Get next suggestion

## Notes

- Review summaries for accuracy - AI may misinterpret discussion
- Action items should have clear owners
- Link issues back to meeting for context
- Consider who said what for attribution
