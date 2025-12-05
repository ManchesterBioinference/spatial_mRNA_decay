---
name: transcribe
description: Transcribe audio files from meetings into text documents using Whisper. Use when the user types /transcribe, has a new audio recording, or when RA detects new audio files in meetings/audio/. Supports speaker diarization with pyannote.
---

# Audio Transcription

> Transcribe audio files from meetings into text documents.

## Usage
```
/transcribe [filename]
/transcribe .research/meetings/audio/2024-12-02-lab-meeting.m4a
/transcribe .research/meetings/audio/  # Transcribe all untranscribed audio in directory
```

## When to Use
- After recording a meeting, seminar, or discussion
- When RA detects new audio files in meetings/audio/ folder
- Before running /summarize_meeting

## Supported Formats
- .m4a, .mp3, .wav, .webm, .mp4 (audio track)
- .ogg, .flac

## Execution

The command runs:

```bash
conda run -n research-assistant python .ra/skills/transcribe/scripts/transcribe.py [filename or .research/meetings/audio/]
```

**Behavior:**
- If `[filename]` provided: transcribe that specific audio file
- If no filename (or `.research/meetings/audio/` specified): automatically detect all audio files without transcripts and process them
- If transcript already exists for a file: skip it
- Output saves to `.research/meetings/transcripts/[same-name].md`

## Post-Transcription Options

```
Transcription complete! 

A) Run /summarize_meeting to extract action items and create tasks
B) Open transcript to review manually first
C) Continue with other work

What would you like to do?
```

## Quality Notes

### Improving Transcription Quality
- Use good microphone/recording quality
- Minimize background noise
- Speak clearly and at moderate pace
- Identify speakers at start if possible

### Limitations
- Speaker diarization may be imperfect
- Technical terms may need manual correction
- Timestamps are approximate

## Configuration

Environment variables (optional):
- `WHISPER_MODEL`: Model size (default: "small", options: tiny, base, small, medium, large-v3)
- `WHISPER_LANGUAGE`: Force language (default: auto-detect)
- `HF_TOKEN`: HuggingFace token for speaker diarization

## Related Skills

- `summarize-meeting` - Extract action items from transcript
- `next` - Get next suggestion

## Notes

- Raw transcripts may contain errors - review before citing
- Keep original audio files as source of truth
- Transcripts are for internal use, not publication
