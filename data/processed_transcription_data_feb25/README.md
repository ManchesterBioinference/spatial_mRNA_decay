# Processed Transcription Data (February 2025)

## Important: File Naming Convention

**IMPORTANT: Different naming conventions for different data types**

### Transcription Data Naming (This Directory)

**The numeric suffixes in filenames encode AP coordinate ranges, NOT embryo IDs.**

| Filename Suffix | AP Range Encoding | Meaning |
|----------------|-------------------|---------|
| `033044` | min=0.33, max=0.44 | AP coordinates from 0.33 to 0.44 |
| `033045` | min=0.33, max=0.45 | AP coordinates from 0.33 to 0.45 (**PRIMARY ANALYSIS**) |
| `033046` | min=0.33, max=0.46 | AP coordinates from 0.33 to 0.46 |

**Example**: 
- `locally_averaged_transcription_traces_feb25_033045_jan25.csv` 
- This file contains data filtered to AP range [0.33, 0.45]
- The "033045" is **not an embryo ID**

### mRNA Data Naming (See data/processed_mRNA_data_stripe2/README.md)

mRNA data uses different conventions depending on source:
- **Stripe2**: Uses nuclear cycle notation (e7, e8_9, e9_10) - pre-processed by Jenny
- **Stripe3+**: Uses embryo numbering (e1, e2, e3, e4) - processed from Ali's raw images

## Files in This Directory

### Locally Averaged Transcription Traces (with bin IDs)
- `locally_averaged_transcription_traces_feb25_033044_jan25.csv`
- `locally_averaged_transcription_traces_feb25_033045_jan25.csv` ← **Primary analysis file**
- `locally_averaged_transcription_traces_feb25_033046_jan25.csv`

**Format**: 
- Columns: `apBin`, `yBin`, followed by timepoint columns (0, 20, 40, ..., 1200 seconds)
- Each row represents a spatial bin's average transcription trace
- 25 rows total (5 AP bins × 5 DV bins)

### Transcription Traces Without IDs
- `noIDs_locally_averaged_transcription_traces_feb25_033044_jan25.csv`
- `noIDs_locally_averaged_transcription_traces_feb25_033045_jan25.csv` ← **Primary analysis file**
- `noIDs_locally_averaged_transcription_traces_feb25_033046_jan25.csv`

**Format**:
- No header row
- No bin ID columns (just timepoint values)
- Used as input to inference scripts (Julia/Python)

### Heatmap Visualizations
- `033044.png` - Spatial heatmap for AP range [0.33, 0.44]
- `033045.png` - Spatial heatmap for AP range [0.33, 0.45] ← **Primary**
- `033046.png` - Spatial heatmap for AP range [0.33, 0.46]

## Processing Details

**Source Data**: `data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv`

**Processing Steps**:
1. Filter nuclei to specified AP coordinate range
2. Bin spatially into 5×5 grid (5 AP bins × 5 DV bins)
3. Compute mean fluorescence trace per spatial bin
4. Extract timepoints 0-1200 seconds (61 timepoints at 20s intervals)

**Spatial Grid**:
- **AP axis**: 5 bins spanning the specified AP range
- **DV axis**: 5 bins spanning dorsal-ventral extent
- **Total bins**: 25 (5 × 5)

## Why Multiple AP Ranges?

These different AP ranges were generated to:
1. Test sensitivity of inference to stripe 2 boundary definition
2. Compare results with slightly wider/narrower spatial windows
3. Assess whether including more interstripe region affects degradation rate estimates

**For primary analysis, use `033045` files** (AP range 0.33-0.45).

## Usage in Pipeline

The Snakemake pipeline (`Snakefile`) uses these files as input to the inference stage:

```yaml
rule infer_degradation_rates:
    input:
        transcription="data/processed_transcription_data_feb25/noIDs_locally_averaged_transcription_traces_feb25_{ap_range}_jan25.csv"
```

The `{ap_range}` wildcard is replaced with `033044`, `033045`, or `033046` depending on which analysis is being run.

## Related Configuration

See `config.yaml` for AP range definitions:

```yaml
ap_ranges:
  - id: "033044"
    min: 0.33
    max: 0.44
  - id: "033045"      # PRIMARY
    min: 0.33
    max: 0.45
  - id: "033046"
    min: 0.33
    max: 0.46
```

## Contact

If you have questions about these data files, contact the project maintainer or refer to:
- `scripts/01_preprocess_eve_data.py` - Processing script
- `.research/meetings/` - Meeting notes discussing AP boundary choices
