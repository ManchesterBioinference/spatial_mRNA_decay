# Research Assistant Template

> A cloneable project template that transforms VS Code + GitHub Copilot into a domain-specific Research Assistant (RA) for computational researchers.

## What Is This?

This is not a standalone tool or extension. It's a **project template** containing context files, slash commands, and configuration that guide VS Code's GitHub Copilot to behave as a specialized research assistant.

Clone this template for each new research project. The RA will:
- Guide you through research phases (planning → development → analysis → writing)
- Help write manuscript sections from your actual code and data
- Ensure reproducibility through proper documentation and DVC pipelines
- Track tasks and action items from meetings
- Conduct regular reviews to keep you on track

## Quick Start

```bash
# Clone the template
git clone https://github.com/[username]/research-assistant-template my-new-project
cd my-new-project

# Remove template's git history
rm -rf .git
git init

# Open in VS Code
code .
```

Then open GitHub Copilot Chat and type:
```
/next
```

The RA will guide you through setup.

## Core Philosophy

**"Feed the Beast"** - The RA is only as good as the context you give it.

This template provides:
- **Context files** that tell the RA about your project and preferences
- **Slash commands** that define specific workflows
- **Templates** for manuscript sections, figure captions, and reviews
- **Phase gates** that prevent skipping critical reproducibility steps

## Key Features

### 🎯 Phase-Based Guidance
Projects progress through: SETUP → PLANNING → DEVELOPMENT → ANALYSIS → WRITING → REVIEW

The RA tracks where you are and won't let you skip critical steps.

### 📝 Manuscript Integration
Write your methods section as you code. Generate results from figures. Draft background from literature reviews.

### 🔍 Literature Research
Deep research with verified citations. No hallucinated references.

### 📋 Task Management
Extract action items from meeting transcripts. Route to tasks.md or GitHub Issues based on complexity.

### 📊 Regular Reviews
Weekly project check-ins. Monthly alignment checks. Quarterly mission reviews.

## Commands

| Command | Purpose |
|---------|---------|
| `/next` | **Start here.** Get suggested next steps based on project state. |
| `/deep_research [topic]` | Literature search with verified citations |
| `/write_background` | Draft background section from literature |
| `/write_methods` | Generate methods from documented scripts |
| `/write_results` | Draft results from figures and captions |
| `/review_script [path]` | Check code documentation quality |
| `/weekly_review` | Weekly progress review |
| `/monthly_review` | Monthly alignment check |
| `/quarterly_review` | Quarterly mission review |
| `/plan_week` | Create focused weekly plan |
| `/transcribe [file]` | Transcribe meeting audio |
| `/summarize_meeting [file]` | Extract action items from transcript |

## Directory Structure

```
├── .ra/                               # RA tool framework
│   ├── copilot-instructions.md        # RA brain and personality
│   ├── commands/                      # Slash command definitions
│   └── tools/                         # RA utilities (transcribe.py, etc.)
├── .research/
│   ├── project_telos.md               # Your project's aims and state
│   ├── phase_checklist.md             # Progress tracking
│   ├── literature/                    # Research outputs with citations
│   ├── meetings/                      # Meeting recordings and transcripts
│   │   ├── audio/                     # Audio files (.m4a, .mp3, .wav)
│   │   └── transcripts/               # Transcript markdown files
│   └── logs/                          # Activity and review logs
├── data/
│   ├── raw/                           # Immutable source data
│   ├── processed/                     # Derived data
│   └── .sensitive/                    # Excluded from AI (IRB data, etc.)
├── scripts/                           # Your analysis code
├── manuscript/
│   ├── background.md
│   ├── methods.md
│   ├── results.md
│   ├── discussion.md
│   └── figures/                       # Figures with caption files
├── meetings/                          # Audio files and transcripts
├── dvc.yaml                           # Reproducible pipeline
├── params.yaml                        # Pipeline parameters
├── tasks.md                           # Quick todos
└── researcher_telos_template.md       # Your profile (moves to ~/.researchAssistant/)
```

## Two-Tier Context

**RA Framework** (`.ra/`):
- `copilot-instructions.md` - RA brain and personality
- `commands/` - Slash command definitions
- `tools/` - RA utilities (transcribe.py, etc.)
- **Committed to git** - gets cloned with every project
- **Readable by Copilot** - Copilot needs these to function as RA
- **Don't modify** - part of the framework you're cloning

**User-level** (`~/.researchAssistant/`):
- `researcher_telos.md` - Your preferences, productivity patterns, strengths
- Persists across all projects
- Set up once, used everywhere

**Project-level** (`.research/`):
- `project_telos.md` - This project's aims, phase, goals
- `phase_checklist.md` - Progress through research phases
- Specific to each project

## Privacy

Sensitive data can be placed in `data/.sensitive/`. This directory is:
- Listed in `.copilotignore` - AI cannot read these files
- Listed in `.gitignore` - Files are not committed

## Requirements

- VS Code with GitHub Copilot enabled
- Python 3.10+ (or R, if preferred)
- DVC for pipeline management
- Git for version control

## Contributing

This template is designed to be forked and customized. If you develop useful commands or improvements, please consider contributing back.

## License

MIT License - Use freely for your research.

---

*Inspired by the "Feed the Beast" philosophy and Telos-based context management.*