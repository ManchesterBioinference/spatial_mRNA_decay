# Project Title

> A reproducible research project powered by the Research Assistant (RA)

[![DVC](https://img.shields.io/badge/DVC-tracked-blue)](https://dvc.org)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Quick Start

1. **Clone this repository**
   ```bash
   git clone [repository-url]
   cd [project-name]
   ```

2. **Set up environment**
   ```bash
   # Create and activate conda environment
   conda env create -f environment.yml
   conda activate research-assistant
   
   # Copy environment variables template
   cp .env.example .env
   # Edit .env to add your HuggingFace token (optional, for speaker diarization)
   ```

3. **Run the pipeline**
   ```bash
   dvc repro
   ```

4. **Get started with the Research Assistant**
   Open in VS Code with GitHub Copilot enabled, then type:
   ```
   /next
   ```

## Environment

This project uses a **conda environment** to ensure all dependencies are available cross-platform.

### What's Included

| Package | Purpose |
|---------|---------|
| `ffmpeg` | Audio format conversion for transcription |
| `pytorch` | ML framework for speaker diarization |
| `dvc` | Data version control and pipelines |
| `snakemake` | Alternative workflow management |
| `faster-whisper` | Speech-to-text (99+ languages) |
| `pyannote.audio` | Speaker diarization (optional) |

### First-Time Model Download

On first use, Whisper models are downloaded automatically (~2GB for `small` model):

```bash
# Test transcription (downloads model on first run)
python tools/transcribe.py --help
```

## Project Overview

**Mission**: [Describe your research question/goal]

**Status**: [Current phase: SETUP / PLANNING / DEVELOPMENT / ANALYSIS / WRITING / REVIEW]

**Target Output**: [Journal paper / Thesis chapter / Tool]

## Repository Structure

```
.
├── .ra/
│   ├── copilot-instructions.md   # Research Assistant configuration
│   ├── commands/                 # Slash command definitions
│   └── tools/                    # RA utilities (transcribe.py, etc.)
├── .research/
│   ├── project_telos.md          # Project aims and state
│   ├── phase_checklist.md        # Progress tracking
│   ├── literature/               # Literature reviews and citations
│   ├── meetings/                 # Meeting recordings and transcripts
│   │   ├── audio/                # Audio files (.m4a, .mp3, .wav)
│   │   └── transcripts/          # Transcript markdown files
│   └── logs/                     # Activity and review logs
├── data/
│   ├── raw/                      # Original, immutable data
│   ├── processed/                # Cleaned/transformed data
│   └── .sensitive/               # Private data (not tracked)
├── scripts/                      # Analysis and utility scripts
├── results/                      # Pipeline outputs and metrics
│   ├── intermediate/             # Intermediate stage outputs
│   ├── final/                    # Final analysis results
│   ├── metrics/                  # DVC-tracked metrics
│   └── logs/                     # Execution logs
├── manuscript/
│   ├── background.md
│   ├── methods.md
│   ├── results.md
│   ├── discussion.md
│   └── figures/                  # Figures with captions
├── dvc.yaml                      # Pipeline definition
├── params.yaml                   # Pipeline parameters
├── tasks.md                      # Current tasks
└── README.md                     # This file
```

## Using the Research Assistant

This project template includes context files that guide VS Code's GitHub Copilot to act as a Research Assistant (RA). The RA helps you:

- **Stay on track**: Guides you through research phases
- **Write better**: Drafts manuscript sections from your work
- **Stay organized**: Tracks tasks and action items
- **Stay reproducible**: Ensures proper documentation

### Key Commands

| Command | Purpose |
|---------|---------|
| `/next` | Get suggested next steps (catch-all for beginners) |
| `/deep_research [topic]` | Literature search with citations |
| `/write_background` | Draft background from literature |
| `/write_methods` | Document scripts as methods |
| `/write_results` | Draft results from figures |
| `/review_script [path]` | Check code documentation |
| `/weekly_review` | Weekly progress review |
| `/plan_week` | Create weekly plan |
| `/transcribe [file]` | Transcribe meeting audio |
| `/summarize_meeting [file]` | Extract action items |

### Transcribing Meetings

Record meetings and transcribe them for documentation:

```bash
# Transcribe a recording
python tools/transcribe.py meetings/recording.m4a

# Or use VS Code chat command
/transcribe meetings/recording.m4a

# Use larger model for better accuracy
python tools/transcribe.py --model large-v3 meetings/recording.m4a
```

**Speaker Diarization**: To identify who is speaking, set up a HuggingFace token. See `tools/README.md` for instructions.

### Getting Started

1. Open VS Code with GitHub Copilot enabled
2. Type `/next` in the chat
3. Follow the RA's guidance to set up your project

## Reproducing This Work

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- Git
- (Optional) NVIDIA GPU with CUDA for faster transcription

### Full Reproduction

```bash
# Clone and enter repository
git clone [repository-url]
cd [project-name]

# Set up environment
conda env create -f environment.yml
conda activate research-assistant

# Copy environment config
cp .env.example .env

# Pull data (if using DVC remote)
dvc pull

# Run complete pipeline
dvc repro

# Figures will be generated in manuscript/figures/
```

### Verify Outputs

```bash
# Check pipeline status
dvc status

# Compare metrics
dvc metrics show
```

## Data

**Source**: [Describe where data comes from]

**Access**: [How to obtain the data]

**Sensitive data**: Files in `data/.sensitive/` are not tracked and must be obtained separately.

## Citation

If you use this work, please cite:

```bibtex
@article{author2024title,
  author = {Author, First and Author, Second},
  title = {Title of the Paper},
  journal = {Journal Name},
  year = {2024},
  doi = {10.xxxx/xxxxx}
}
```

## License

[Choose appropriate license]

## Contact

[Your name and contact information]

---

*This project uses the [Research Assistant Template](https://github.com/[link-to-template])*
