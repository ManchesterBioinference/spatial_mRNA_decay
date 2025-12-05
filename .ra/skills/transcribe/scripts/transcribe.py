#!/usr/bin/env python3
"""
Transcribe audio files using Whisper with optional speaker diarization.

This script is designed to work with the Research Assistant /transcribe skill.
It uses faster-whisper for cross-platform, multi-language transcription and
optionally pyannote.audio for speaker diarization.

Usage:
    python .ra/skills/transcribe/scripts/transcribe.py .research/meetings/audio/recording.m4a
    python .ra/skills/transcribe/scripts/transcribe.py .research/meetings/audio/  # Process all untranscribed audio
    python .ra/skills/transcribe/scripts/transcribe.py --model large-v3 --language en .research/meetings/audio/recording.m4a

For speaker diarization, set HF_TOKEN environment variable. See .ra/skills/transcribe/README.md.
"""

import os
import sys
import glob
import subprocess
import argparse
import json
import tempfile
from datetime import datetime
from pathlib import Path

# Load environment variables from .env file
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass  # dotenv is optional

from tqdm import tqdm
import numpy as np

# Configuration from environment
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent.parent  # .ra/skills/transcribe/scripts -> project root

MEETINGS_AUDIO_DIR = PROJECT_ROOT / '.research' / 'meetings' / 'audio'
MEETINGS_TRANSCRIPTS_DIR = PROJECT_ROOT / '.research' / 'meetings' / 'transcripts'
SPEAKER_DB_FILE = PROJECT_ROOT / '.research' / 'speaker_profiles.json'

# Whisper configuration (can be overridden via CLI or environment)
DEFAULT_MODEL = os.environ.get('WHISPER_MODEL', 'small')
DEFAULT_LANGUAGE = os.environ.get('WHISPER_LANGUAGE', None)  # None = auto-detect
DEFAULT_COMPUTE_TYPE = os.environ.get('WHISPER_COMPUTE_TYPE', 'auto')

# HuggingFace token for pyannote models
HF_TOKEN = os.environ.get('HF_TOKEN', None)

# Similarity threshold for speaker recognition (0.0 to 1.0, higher = stricter)
SPEAKER_SIMILARITY_THRESHOLD = 0.75

# Supported audio formats
AUDIO_EXTENSIONS = {'.m4a', '.mp3', '.wav', '.webm', '.mp4', '.ogg', '.flac'}

# Lazy-loaded models
_whisper_model = None
_diarization_pipeline = None
_embedding_model = None

# Check for optional dependencies
DIARIZATION_AVAILABLE = False
try:
    import torch
    from pyannote.audio import Pipeline
    from pyannote.audio.pipelines.utils.hook import ProgressHook
    DIARIZATION_AVAILABLE = True
except ImportError:
    pass


def get_whisper_model(model_size: str = None, compute_type: str = None):
    """Get or initialize the Whisper model."""
    global _whisper_model
    
    model_size = model_size or DEFAULT_MODEL
    compute_type = compute_type or DEFAULT_COMPUTE_TYPE
    
    if _whisper_model is None:
        try:
            from faster_whisper import WhisperModel
            
            print(f"Loading Whisper model '{model_size}'...")
            
            # Determine device and compute type
            device = "cpu"
            if compute_type == "auto":
                # Auto-detect best settings
                if DIARIZATION_AVAILABLE and torch.cuda.is_available():
                    device = "cuda"
                    compute_type = "float16"
                elif DIARIZATION_AVAILABLE and torch.backends.mps.is_available():
                    # MPS (Apple Silicon) - use CPU with int8 for now
                    # faster-whisper doesn't fully support MPS yet
                    device = "cpu"
                    compute_type = "int8"
                else:
                    device = "cpu"
                    compute_type = "int8"
            
            print(f"  Device: {device}, Compute type: {compute_type}")
            
            _whisper_model = WhisperModel(
                model_size,
                device=device,
                compute_type=compute_type
            )
            
            print(f"  Model loaded successfully")
            
        except Exception as e:
            print(f"Error: Failed to load Whisper model: {e}")
            print("Make sure faster-whisper is installed: pip install faster-whisper")
            sys.exit(1)
    
    return _whisper_model


def get_diarization_pipeline():
    """Get or initialize the speaker diarization pipeline."""
    global _diarization_pipeline
    
    if not DIARIZATION_AVAILABLE:
        return None
    
    if _diarization_pipeline is None:
        if HF_TOKEN is None:
            print("\nWarning: HF_TOKEN environment variable not set.")
            print("To enable speaker diarization:")
            print("  1. Create account at: https://huggingface.co/join")
            print("  2. Accept model terms at: https://hf.co/pyannote/speaker-diarization-3.1")
            print("  3. Accept model terms at: https://hf.co/pyannote/segmentation-3.0")
            print("  4. Create token at: https://hf.co/settings/tokens")
            print("  5. Add to .env file: HF_TOKEN=hf_your_token_here")
            print("\nContinuing without speaker diarization...\n")
            return None
        
        try:
            print("Loading speaker diarization model...")
            
            # Fix for PyTorch 2.6+ weights_only default change
            from pyannote.audio.core.task import Specifications, Problem, Resolution
            torch.serialization.add_safe_globals([Specifications, Problem, Resolution])
            
            _diarization_pipeline = Pipeline.from_pretrained(
                "pyannote/speaker-diarization-3.1",
                use_auth_token=HF_TOKEN
            )
            
            # Send to GPU if available
            if torch.cuda.is_available():
                _diarization_pipeline.to(torch.device("cuda"))
                print("  Using CUDA GPU for diarization")
            elif torch.backends.mps.is_available():
                _diarization_pipeline.to(torch.device("mps"))
                print("  Using Apple Silicon GPU for diarization")
            else:
                print("  Using CPU for diarization (this may be slow)")
                
        except Exception as e:
            print(f"Warning: Failed to load diarization pipeline: {e}")
            print("Continuing without speaker diarization...\n")
            return None
    
    return _diarization_pipeline


def get_embedding_model():
    """Get or initialize the speaker embedding model for recognition."""
    global _embedding_model
    
    if not DIARIZATION_AVAILABLE:
        return None
    
    if _embedding_model is None:
        if HF_TOKEN is None:
            return None
        
        try:
            from pyannote.audio import Model, Inference
            
            print("Loading speaker embedding model...")
            
            # Temporarily allow unsafe loading for pyannote models
            original_load = torch.load
            torch.load = lambda *args, **kwargs: original_load(*args, **{**kwargs, 'weights_only': False})
            
            try:
                model = Model.from_pretrained(
                    "pyannote/embedding",
                    use_auth_token=HF_TOKEN
                )
            finally:
                torch.load = original_load
            
            _embedding_model = Inference(model, window="whole")
            
            # Send to GPU if available
            if torch.cuda.is_available():
                _embedding_model.to(torch.device("cuda"))
            elif torch.backends.mps.is_available():
                _embedding_model.to(torch.device("mps"))
                
        except Exception as e:
            print(f"Warning: Failed to load embedding model: {e}")
            return None
    
    return _embedding_model


def get_audio_duration(file_path: Path) -> float:
    """Get audio duration in seconds using ffprobe."""
    try:
        result = subprocess.run(
            ['ffprobe', '-v', 'error', '-show_entries', 'format=duration',
             '-of', 'default=noprint_wrappers=1:nokey=1', str(file_path)],
            capture_output=True,
            text=True
        )
        if result.returncode == 0 and result.stdout.strip():
            return float(result.stdout.strip())
    except Exception:
        pass
    return None


def format_timestamp(seconds: float) -> str:
    """Format seconds as HH:MM:SS."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def transcribe_audio(file_path: Path, model_size: str = None, 
                     language: str = None, compute_type: str = None) -> tuple:
    """
    Transcribe audio file using Whisper.
    
    Returns:
        tuple: (segments_list, detected_language, duration)
        Each segment is (start_time, end_time, text)
    """
    model = get_whisper_model(model_size, compute_type)
    duration = get_audio_duration(file_path)
    duration_str = f" ({duration/60:.1f} minutes)" if duration else ""
    
    print(f"Transcribing: {file_path.name}{duration_str}")
    
    try:
        # Transcribe with word-level timestamps for diarization alignment
        segments, info = model.transcribe(
            str(file_path),
            language=language,
            beam_size=5,
            word_timestamps=True,
            vad_filter=True,  # Filter out non-speech
            vad_parameters=dict(min_silence_duration_ms=500)
        )
        
        # Collect segments with timestamps
        transcript_segments = []
        for segment in segments:
            transcript_segments.append((
                segment.start,
                segment.end,
                segment.text.strip()
            ))
        
        detected_lang = info.language
        print(f"  Detected language: {detected_lang} (probability: {info.language_probability:.2f})")
        print(f"  Transcribed {len(transcript_segments)} segments")
        
        return transcript_segments, detected_lang, duration
        
    except Exception as e:
        print(f"Error: Transcription failed: {e}")
        return None, None, None


def perform_diarization(file_path: Path) -> list:
    """
    Perform speaker diarization on an audio file.
    
    Returns:
        list: [(start_time, end_time, speaker_label), ...]
    """
    pipeline = get_diarization_pipeline()
    if pipeline is None:
        return None
    
    print(f"Performing speaker diarization...")
    
    try:
        import torchaudio
        
        # Convert to WAV if needed (torchaudio handles most formats but WAV is safest)
        if file_path.suffix.lower() in {'.m4a', '.mp3', '.aac', '.ogg', '.flac', '.webm', '.mp4'}:
            print("  Converting audio to WAV format...")
            with tempfile.NamedTemporaryFile(suffix='.wav', delete=False) as tmp:
                tmp_path = tmp.name
            
            result = subprocess.run(
                ['ffmpeg', '-y', '-i', str(file_path), '-ar', '16000', '-ac', '1', tmp_path],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                print(f"  Warning: FFmpeg conversion failed: {result.stderr[:200]}")
                return None
            
            waveform, sample_rate = torchaudio.load(tmp_path)
            os.unlink(tmp_path)
        else:
            waveform, sample_rate = torchaudio.load(str(file_path))
        
        # Convert to mono if needed
        if waveform.shape[0] > 1:
            waveform = waveform.mean(dim=0, keepdim=True)
        
        # Create audio dict for pyannote
        audio_dict = {"waveform": waveform, "sample_rate": sample_rate}
        
        # Run diarization with progress
        with ProgressHook() as hook:
            diarization = pipeline(audio_dict, hook=hook)
        
        # Extract speaker segments
        segments = []
        for turn, _, speaker in diarization.itertracks(yield_label=True):
            segments.append((turn.start, turn.end, speaker))
        
        num_speakers = len(set(s[2] for s in segments))
        print(f"  Identified {num_speakers} speakers in {len(segments)} segments")
        
        return segments
        
    except Exception as e:
        print(f"  Warning: Diarization failed: {e}")
        return None


def combine_transcript_with_diarization(transcript_segments: list, 
                                         diarization_segments: list,
                                         speaker_mapping: dict = None) -> str:
    """
    Combine transcript with speaker labels from diarization.
    
    Returns formatted transcript with speaker annotations.
    """
    if not transcript_segments:
        return ""
    
    if not diarization_segments:
        # No diarization - just return plain transcript
        lines = []
        for start, end, text in transcript_segments:
            lines.append(f"[{format_timestamp(start)}] {text}")
        return '\n\n'.join(lines)
    
    # Create speaker name mapping
    if speaker_mapping is None:
        speaker_ids = sorted(set(seg[2] for seg in diarization_segments))
        speaker_mapping = {sid: f"SPEAKER_{i+1}" for i, sid in enumerate(speaker_ids)}
    
    result_lines = []
    current_speaker = None
    current_text = []
    current_start = None
    
    for t_start, t_end, text in transcript_segments:
        # Find the speaker with most overlap for this segment
        best_speaker = None
        best_overlap = 0
        
        for d_start, d_end, speaker in diarization_segments:
            overlap_start = max(t_start, d_start)
            overlap_end = min(t_end, d_end)
            overlap = max(0, overlap_end - overlap_start)
            
            if overlap > best_overlap:
                best_overlap = overlap
                best_speaker = speaker
        
        speaker_name = speaker_mapping.get(best_speaker, "UNKNOWN") if best_speaker else "UNKNOWN"
        
        # Group consecutive text from same speaker
        if speaker_name == current_speaker:
            current_text.append(text)
        else:
            # Output previous speaker's text
            if current_speaker and current_text:
                timestamp = format_timestamp(current_start)
                combined_text = ' '.join(current_text)
                result_lines.append(f"[{timestamp}] [{current_speaker}]:\n{combined_text}")
            
            current_speaker = speaker_name
            current_text = [text]
            current_start = t_start
    
    # Don't forget the last segment
    if current_speaker and current_text:
        timestamp = format_timestamp(current_start)
        combined_text = ' '.join(current_text)
        result_lines.append(f"[{timestamp}] [{current_speaker}]:\n{combined_text}")
    
    return '\n\n'.join(result_lines)


def load_speaker_database() -> dict:
    """Load saved speaker profiles from disk."""
    if SPEAKER_DB_FILE.exists():
        try:
            with open(SPEAKER_DB_FILE, 'r') as f:
                data = json.load(f)
                return {
                    name: {
                        'embedding': np.array(profile['embedding']),
                        'description': profile.get('description', ''),
                        'sample_count': profile.get('sample_count', 1)
                    }
                    for name, profile in data.items()
                }
        except Exception as e:
            print(f"Warning: Failed to load speaker database: {e}")
    return {}


def save_speaker_database(speaker_db: dict):
    """Save speaker profiles to disk."""
    SPEAKER_DB_FILE.parent.mkdir(parents=True, exist_ok=True)
    
    data = {
        name: {
            'embedding': profile['embedding'].tolist(),
            'description': profile.get('description', ''),
            'sample_count': profile.get('sample_count', 1)
        }
        for name, profile in speaker_db.items()
    }
    
    try:
        with open(SPEAKER_DB_FILE, 'w') as f:
            json.dump(data, f, indent=2)
    except Exception as e:
        print(f"Warning: Failed to save speaker database: {e}")


def cosine_similarity(a, b) -> float:
    """Calculate cosine similarity between two embeddings."""
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))


def generate_transcript_markdown(file_path: Path, transcript_text: str,
                                  detected_language: str, duration: float,
                                  model_used: str, has_diarization: bool) -> str:
    """Generate markdown-formatted transcript document."""
    date_str = datetime.now().strftime('%Y-%m-%d %H:%M')
    duration_str = f"{duration/60:.1f} minutes" if duration else "Unknown"
    
    # Extract date from filename if present (e.g., 2024-12-02-meeting.m4a)
    recording_date = "See filename"
    name_parts = file_path.stem.split('-')
    if len(name_parts) >= 3:
        try:
            # Check if first three parts look like a date
            int(name_parts[0])
            int(name_parts[1])
            int(name_parts[2])
            recording_date = f"{name_parts[0]}-{name_parts[1]}-{name_parts[2]}"
        except ValueError:
            pass
    
    md = f"""# Meeting Transcript: {file_path.stem}

**Date**: {recording_date}
**Duration**: {duration_str}
**Participants**: {"See transcript" if has_diarization else "Speaker diarization not available"}
**Audio file**: {file_path.name}

---

## Transcript

{transcript_text}

---

## Metadata

- **Transcribed**: {date_str}
- **Model**: Whisper {model_used}
- **Language**: {detected_language or "auto-detected"}
- **Speaker diarization**: {"Yes" if has_diarization else "No"}

## Next Steps

Run `/summarize_meeting {file_path.stem}.md` to:
- Extract action items
- Identify key decisions
- Create tasks and issues

---

*Raw transcript - not yet summarized*
"""
    return md


def find_untranscribed_audio(directory: Path) -> list:
    """Find audio files in directory that don't have matching .md transcripts."""
    audio_files = []
    
    for ext in AUDIO_EXTENSIONS:
        audio_files.extend(directory.glob(f'*{ext}'))
        audio_files.extend(directory.glob(f'*{ext.upper()}'))
    
    # Filter out files that already have transcripts
    untranscribed = []
    for audio_file in audio_files:
        # Check if transcript exists in the transcripts directory
        transcript_path = MEETINGS_TRANSCRIPTS_DIR / audio_file.with_suffix('.md').name
        if not transcript_path.exists():
            untranscribed.append(audio_file)
    
    return sorted(untranscribed)


def process_file(file_path: Path, model_size: str = None, language: str = None,
                 compute_type: str = None, enable_diarization: bool = True,
                 enable_recognition: bool = True) -> bool:
    """
    Process a single audio file: transcribe and optionally diarize.
    
    Returns True if successful, False otherwise.
    """
    print(f"\n{'='*60}")
    print(f"Processing: {file_path.name}")
    print('='*60)
    
    # Transcribe
    transcript_segments, detected_language, duration = transcribe_audio(
        file_path, model_size, language, compute_type
    )
    
    if transcript_segments is None:
        print(f"Error: Failed to transcribe {file_path.name}")
        return False
    
    # Optionally perform diarization
    diarization_segments = None
    speaker_mapping = None
    
    if enable_diarization and DIARIZATION_AVAILABLE:
        diarization_segments = perform_diarization(file_path)
        
        # TODO: Add speaker recognition if enabled
        # For now, just use generic speaker labels
    
    # Combine transcript with diarization
    transcript_text = combine_transcript_with_diarization(
        transcript_segments, diarization_segments, speaker_mapping
    )
    
    # Generate markdown document
    model_used = model_size or DEFAULT_MODEL
    markdown = generate_transcript_markdown(
        file_path, transcript_text, detected_language, duration,
        model_used, diarization_segments is not None
    )
    
    # Save transcript
    # Transcripts go to .research/meetings/transcripts/ with the same stem as the audio
    output_path = MEETINGS_TRANSCRIPTS_DIR / file_path.with_suffix('.md').name
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(markdown)
        print(f"\n✓ Transcript saved: {output_path}")
        return True
    except Exception as e:
        print(f"Error: Failed to save transcript: {e}")
        return False


def list_speakers():
    """List all speakers in the database."""
    speaker_db = load_speaker_database()
    
    if not speaker_db:
        print("No speakers in database.")
        print("Process audio files with diarization to build speaker profiles.")
        return
    
    print(f"\nSpeaker Database ({len(speaker_db)} profiles):")
    print("-" * 60)
    
    named = [(n, p) for n, p in speaker_db.items() if not n.startswith('UNKNOWN_')]
    unknown = [(n, p) for n, p in speaker_db.items() if n.startswith('UNKNOWN_')]
    
    for name, profile in sorted(named) + sorted(unknown):
        samples = profile.get('sample_count', 1)
        desc = profile.get('description', '')
        print(f"  {name}")
        print(f"    Samples: {samples}, Description: {desc[:50]}")
    
    if unknown:
        print(f"\nTip: Rename unknown speakers with: --rename-speaker OLD_NAME NEW_NAME")


def rename_speaker(old_name: str, new_name: str) -> bool:
    """Rename a speaker in the database."""
    speaker_db = load_speaker_database()
    
    if old_name not in speaker_db:
        print(f"Error: Speaker '{old_name}' not found.")
        print(f"Available speakers: {', '.join(speaker_db.keys())}")
        return False
    
    if new_name in speaker_db:
        print(f"Error: Speaker '{new_name}' already exists. Use --merge-speakers to combine.")
        return False
    
    speaker_db[new_name] = speaker_db.pop(old_name)
    speaker_db[new_name]['description'] = f'Renamed from {old_name}'
    save_speaker_database(speaker_db)
    print(f"Renamed '{old_name}' to '{new_name}'")
    return True


def delete_speaker(name: str) -> bool:
    """Delete a speaker from the database."""
    speaker_db = load_speaker_database()
    
    if name not in speaker_db:
        print(f"Error: Speaker '{name}' not found.")
        return False
    
    del speaker_db[name]
    save_speaker_database(speaker_db)
    print(f"Deleted speaker '{name}'")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Transcribe audio files with optional speaker diarization.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s .research/meetings/audio/recording.m4a              # Transcribe single file
  %(prog)s .research/meetings/audio/                            # Process all untranscribed audio
  %(prog)s --model large-v3 recording.m4a      # Use larger model for better accuracy
  %(prog)s --language ja recording.m4a         # Specify Japanese language
  %(prog)s --no-diarization recording.m4a      # Skip speaker identification

Model sizes (speed vs accuracy tradeoff):
  tiny   - Fastest, ~1GB RAM, good for quick drafts
  base   - Fast, ~1GB RAM
  small  - Balanced (default), ~2GB RAM, good multilingual
  medium - Slower, ~5GB RAM, better accuracy
  large-v3 - Slowest, ~10GB RAM, best accuracy
  turbo  - Fast like base, accuracy like large (English-optimized)

For 45 min audio on CPU: tiny ~15min, small ~1hr, large-v3 ~4-6hr
        """
    )
    
    # Input
    parser.add_argument('input', nargs='?', default=None,
                        help='Audio file or directory to process')
    
    # Model options
    parser.add_argument('--model', '-m', default=None,
                        choices=['tiny', 'base', 'small', 'medium', 'large-v3', 'turbo'],
                        help=f'Whisper model size (default: {DEFAULT_MODEL})')
    parser.add_argument('--language', '-l', default=None,
                        help='Language code (e.g., en, ja, de). Default: auto-detect')
    parser.add_argument('--compute-type', default=None,
                        choices=['auto', 'int8', 'float16', 'float32'],
                        help='Compute type for inference (default: auto)')
    
    # Diarization options
    parser.add_argument('--no-diarization', action='store_true',
                        help='Disable speaker diarization (faster)')
    parser.add_argument('--no-recognition', action='store_true',
                        help='Disable speaker recognition (still labels speakers)')
    
    # Speaker database management
    parser.add_argument('--list-speakers', action='store_true',
                        help='List all speakers in the database')
    parser.add_argument('--rename-speaker', nargs=2, metavar=('OLD', 'NEW'),
                        help='Rename a speaker in the database')
    parser.add_argument('--delete-speaker', metavar='NAME',
                        help='Delete a speaker from the database')
    
    args = parser.parse_args()
    
    # Handle speaker database commands
    if args.list_speakers:
        list_speakers()
        return
    
    if args.rename_speaker:
        rename_speaker(args.rename_speaker[0], args.rename_speaker[1])
        return
    
    if args.delete_speaker:
        delete_speaker(args.delete_speaker)
        return
    
    # Determine input files
    if args.input is None:
        # Default to .research/meetings/audio/ directory
        input_path = MEETINGS_AUDIO_DIR
    else:
        input_path = Path(args.input)
        if not input_path.is_absolute():
            input_path = Path.cwd() / input_path
    
    if not input_path.exists():
        print(f"Error: Path not found: {input_path}")
        sys.exit(1)
    
    # Collect files to process
    if input_path.is_file():
        if input_path.suffix.lower() not in AUDIO_EXTENSIONS:
            print(f"Error: Not a supported audio file: {input_path}")
            print(f"Supported formats: {', '.join(sorted(AUDIO_EXTENSIONS))}")
            sys.exit(1)
        files_to_process = [input_path]
    else:
        # Directory - find untranscribed audio
        files_to_process = find_untranscribed_audio(input_path)
        
        if not files_to_process:
            print(f"No untranscribed audio files found in: {input_path}")
            print("(Audio files with existing .md transcripts are skipped)")
            return
        
        print(f"Found {len(files_to_process)} untranscribed audio file(s):")
        for f in files_to_process:
            print(f"  - {f.name}")
        print()
    
    # Process files
    model_size = args.model or DEFAULT_MODEL
    language = args.language or DEFAULT_LANGUAGE
    compute_type = args.compute_type or DEFAULT_COMPUTE_TYPE
    enable_diarization = not args.no_diarization
    enable_recognition = not args.no_recognition
    
    # Print configuration
    print("Configuration:")
    print(f"  Model: {model_size}")
    print(f"  Language: {language or 'auto-detect'}")
    print(f"  Diarization: {'enabled' if enable_diarization else 'disabled'}")
    if enable_diarization and not DIARIZATION_AVAILABLE:
        print("  (Note: pyannote.audio not installed, diarization unavailable)")
    elif enable_diarization and not HF_TOKEN:
        print("  (Note: HF_TOKEN not set, diarization unavailable)")
    
    success_count = 0
    for file_path in files_to_process:
        if process_file(file_path, model_size, language, compute_type,
                       enable_diarization, enable_recognition):
            success_count += 1
    
    print(f"\n{'='*60}")
    print(f"Completed: {success_count}/{len(files_to_process)} files transcribed successfully")
    
    if success_count > 0:
        print("\nNext steps:")
        print("  Run /summarize_meeting to extract action items from transcripts")


if __name__ == "__main__":
    main()
