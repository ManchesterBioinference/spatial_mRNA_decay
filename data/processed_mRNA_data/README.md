# Processed mRNA Data

## Files in This Directory

### mRNA Count Data (for model fitting)
- `e7_edgespotsremoved_sass_formodel.csv`
- `e8_9_edgespotsremoved_sass_formodel.csv` ← **Used in current pipeline**
- `e9_10_edgespotsremoved_sass_formodel.csv`

## File Naming Convention

**The "e" prefix numbers refer to nuclear cycle stages, NOT embryo IDs.**

| Filename | Nuclear Cycle Stage | Notes |
|----------|---------------------|-------|
| `e7_...` | Nuclear cycle 7 | Earlier stage |
| `e8_9_...` | Nuclear cycle 8-9 transition | **Used in primary analysis** |
| `e9_10_...` | Nuclear cycle 9-10 transition | Later stage |

These represent different developmental timepoints in the Drosophila embryo.

## Data Format

**Format**: CSV with **no header row**
- One value per line
- 25 values total (one per spatial bin)
- Values represent: **average mRNA count per nucleus** in each spatial bin

**Spatial organization**:
- Values correspond to 5×5 spatial grid (5 AP bins × 5 DV bins)
- Organized by spatial bin (matching transcription data binning)

## Processing Details

**Processing steps applied**:
1. **Edge spots removed**: mRNA spots near nuclear boundaries excluded to avoid overcounting
2. **SASS method**: Single-cell Automated Segmentation and Spot-calling
3. **Averaged per bin**: Mean mRNA count computed for all nuclei in each spatial bin

**Source**: smFISH imaging data (single-molecule fluorescence in situ hybridization)

## Usage in Pipeline

These data serve as the **observed mRNA counts** in the Bayesian inference model:

```python
# Model: dm/dt = γ*F(t) - D*m
# Where:
#   F(t) = transcription trace (from processed_transcription_data_feb25/)
#   m = observed mRNA (from this directory)
#   D = degradation rate (to be inferred)
```

The inference fits the ODE model to match these observed mRNA values.

## Why Multiple Nuclear Cycle Stages?

Different nuclear cycle stages were analyzed to:
1. Understand how degradation rates change during development
2. Test model robustness across developmental timepoints
3. Compare early vs. late nuclear cycle dynamics

**For primary analysis, we use `e8_9_edgespotsremoved_sass_formodel.csv`** (nuclear cycle 8-9 transition).

## Related Files

- **Transcription data**: `data/processed_transcription_data_feb25/`
- **Pipeline config**: `config.yaml`
- **Inference script**: `scripts/02_infer_degradation_rates.py`

## Contact

For questions about mRNA counting methodology, contact Ali or refer to the SASS methods paper.
