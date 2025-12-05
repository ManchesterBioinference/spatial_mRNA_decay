---
name: task
description: Rapid task capture without breaking flow. Use when the user types /task or wants to quickly add a todo item to tasks.md. Supports priority flags (!high, !low). Use this instead of /note for action items.
---

# Rapid Task Capture

> Rapid task capture without breaking flow.

## When to Use
- When a task pops into your head while working
- After identifying something that needs to be done
- When someone assigns you something in a meeting
- Any time you need to remember to do something later

## Syntax

```
/task [task description]
/task !high [task description]
/task !low [task description]
```

**Examples:**
```
/task Fix the data loader bug
/task !high Review PI's comments before meeting
/task !low Refactor the plotting utilities
/task Update documentation for preprocessing script
```

## Execution Steps

### 1. Parse Input

Extract priority and task text:

| Input | Priority | Task Text |
|-------|----------|-----------|
| `/task Fix bug` | Normal | Fix bug |
| `/task !high Fix bug` | High | Fix bug |
| `/task !low Fix bug` | Low | Fix bug |

Priority flags must be at the start of the task text.

### 2. Add to tasks.md

Append task to the appropriate section in `tasks.md`:

| Priority | Section |
|----------|---------|
| `!high` | `## High Priority` |
| (default) | `## Normal Priority` |
| `!low` | `## Low Priority / Someday` |

Format:
```markdown
- [ ] [Task description]
```

### 3. Confirm (Brief)

Single-line confirmation:

```
✅ Added to [priority] tasks.
```

Examples:
```
✅ Added to high priority tasks.
✅ Added to tasks.
✅ Added to low priority tasks.
```

## Behavior Notes

### Speed is Priority
- No confirmation dialogs
- No follow-up questions  
- Minimal output
- Get the user back to work immediately

### Task Format
Tasks are added as unchecked markdown checkboxes:
```markdown
- [ ] Task description here
```

### Append Position
New tasks are added at the **end** of the relevant priority section, before any blank lines or the next section.

### No Deduplication
If the user adds the same task twice, add it twice. They may have a reason, and checking for duplicates slows things down.

## Edge Cases

### Empty task
```
/task
```
Response:
```
What's the task? Usage: /task [description] or /task !high [description]
```

### Task looks like it should be an issue
If task description suggests >2 hours work (mentions "implement", "build", "create new", "research"):
```
✅ Added to tasks. (This sounds substantial - consider creating a GitHub issue for tracking)
```

Only suggest this for obviously large tasks.

### Task with context
User can add any context they want:
```
/task Fix preprocessing bug - see note from 10:30am
```
Just add it as-is.

## Example Session

```
User: /task Fix the edge case in validate_input.py

RA: ✅ Added to tasks.

User: /task !high Respond to reviewer comments

RA: ✅ Added to high priority tasks.

User: /task !low Someday refactor the whole pipeline

RA: ✅ Added to low priority tasks.
```

**Resulting tasks.md:**
```markdown
## High Priority
- [ ] Respond to reviewer comments

## Normal Priority
- [ ] Fix the edge case in validate_input.py

## Low Priority / Someday
- [ ] Someday refactor the whole pipeline
```

## Related Skills
- `note` - For observations/thoughts (not action items)
- `plan-week` - Organize tasks into weekly plan
- `summarize-meeting` - Extracts tasks from meeting notes automatically
