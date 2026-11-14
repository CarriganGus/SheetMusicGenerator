# Sheet Music Generator

EECS 351 (Digital Signal Processing) term project, University of Michigan.

Transcribes audio recordings into musical notes (pitch, onset, offset) using
DSP techniques (STFT, spectral flux onset detection, harmonic product spectrum),
and evaluates the results against MIDI ground truth.

## Team
- Gus Carrigan
- _TBD_
- _TBD_
- _TBD_

## Repository Layout
```
matlab/       MATLAB implementation
python/       Python implementation
data/audio/   WAV recordings (not tracked - see Data below)
data/midi/    MIDI ground-truth files (tracked)
results/      Figures and evaluation output
docs/         Report and website material
```

## Data
Audio files are too large for the repo and are gitignored. Place WAV files in
`data/audio/` with the same base name as their matching MIDI file in
`data/midi/` (e.g. `c_major_scale.wav` <-> `c_major_scale.mid`).

Data sources: _TBD_

## How to Run
### MATLAB
_TBD_

### Python
_TBD_

## Evaluation
Detected notes are compared to MIDI ground truth using note-level precision,
recall, and F1 (onset tolerance +/-50 ms, pitch tolerance +/-50 cents).
