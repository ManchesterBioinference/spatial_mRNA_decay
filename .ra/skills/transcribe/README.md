# Research Assistant Tools

This directory contains utilities that are part of the Research Assistant (RA) template, separate from your research project scripts.

## Transcription Tool

Transcribe meeting recordings with speaker diarization:

```bash
# Transcribe a single file
python tools/transcribe.py meetings/recording.m4a

# Transcribe all untranscribed audio in meetings/
python tools/transcribe.py meetings/

# Use a larger model for better accuracy (slower)
python tools/transcribe.py --model large-v3 meetings/recording.m4a

# Specify language (default: auto-detect)
python tools/transcribe.py --language ja meetings/recording.m4a

# Skip speaker diarization (faster)
python tools/transcribe.py --no-diarization meetings/recording.m4a
```

### Whisper Model Sizes

| Model | RAM | 45 min audio (CPU) | Notes |
|-------|-----|-------------------|-------|
| `tiny` | ~1 GB | ~15 min | Fastest, good for drafts |
| `base` | ~1 GB | ~25 min | Fast |
| `small` | ~2 GB | ~1 hour | **Default** - good balance |
| `medium` | ~5 GB | ~2-3 hours | Better accuracy |
| `large-v3` | ~10 GB | ~4-6 hours | Best accuracy |
| `turbo` | ~6 GB | ~30 min | English-optimized |

### Speaker Diarization Setup (Optional)

Speaker diarization identifies who is speaking in a recording. To enable:

1. **Create HuggingFace account**: https://huggingface.co/join

2. **Accept model terms** (required for each model):
   - https://hf.co/pyannote/speaker-diarization-3.1
   - https://hf.co/pyannote/segmentation-3.0

3. **Create access token**: https://hf.co/settings/tokens
   - Click "New token"
   - Name: "research-assistant" (or any name)
   - Type: Read
   - Copy the token

4. **Add to your `.env` file**:
   ```
   HF_TOKEN=hf_your_token_here
   ```

Without `HF_TOKEN`, transcription works normally but without speaker labels.

### Speaker Management

After transcribing with diarization, you can manage speaker profiles:

```bash
# List known speakers
python tools/transcribe.py --list-speakers

# Rename an unknown speaker
python tools/transcribe.py --rename-speaker UNKNOWN_20241202_SPEAKER_0 "Alice"

# Delete a speaker profile
python tools/transcribe.py --delete-speaker "Old Name"
```

---

## Other Tools

Additional RA tools will be added here as needed.
