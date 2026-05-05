# Methods

## Overview

To quantify spatially varying *eve* mRNA degradation in *Drosophila* embryos, we integrated (i) MS2 live-cell imaging measurements of *eve* transcription dynamics (Berrocal et al., 2020; eLife doi:10.7554/eLife.61635) with (ii) single-molecule fluorescence *in situ* hybridisation (smFISH) measurements of mature *eve* mRNA abundance. The analysis pipeline was implemented in Python and executed via Snakemake to produce a fully reproducible workflow from raw inputs to inferred degradation parameters and diagnostic figures.

Unless stated otherwise, time was represented in minutes and transcription traces were sampled at 20 s resolution (i.e., $\Delta t = 20/60$ min). For sensitivity analyses of transcription-window choice, the pipeline was executed across multiple maximum transcription integration windows (`max_time_seconds`), producing per-window result directories (`results_{max_time}/`).

## Data

### MS2 transcription data (Berrocal 2020)

MS2 fluorescence traces were obtained from Berrocal et al. (2020) and downloaded from the eLife supplementary archive (`https://cdn.elifesciences.org/articles/61635/elife-61635-supp1-v2.zip`). The supplementary notebook and CSVs were programmatically filtered by executing the initial notebook cells used in the original analysis (script: `scripts/00_download_berrocal_data.py`), producing a filtered longform table (`eve_data_longform_w_nuclei_060520_FILTERED.csv`) containing per-nucleus fluorescence over time together with registered anterior–posterior (AP) positions.

### smFISH mature mRNA data

**Stripe 2.** Stripe 2 smFISH spot-count data were provided pre-processed by Jenny Pantin in a format consistent with the SASS output (see below), with embryos at multiple developmental stages (e.g., e7, e8_9, e9_10).

**Stripes 3 and 4.** For stripes 3 and 4, raw Imaris exports were provided (directories `data/Ali_embryos/{stripe}/{embryo}/nuclei_Statistics/` and `.../spots_Statistics/`). Prior to export, images had been manually centered on the target stripe and edge spots removed.

## Stripe identification and spatial binning

### Automated stripe window detection

Stripe AP windows were derived from the Berrocal MS2 dataset (script: `scripts/00_identify_stripe_ranges.py`). For each nucleus, fluorescence was summarized using a windowed sum over the final four frames ending at a configurable timepoint (`max_time_seconds`, in 20 s steps), which emphasized the peak stripe signal while reducing contributions from very early/late timepoints. Nuclei were binned along AP using bins of width 0.005 (Snakefile parameterization), and the median windowed fluorescence was computed per AP bin.

The resulting 1D AP fluorescence profile was smoothed using a Savitzky–Golay filter (window length 17, polynomial order 3) and stripe centers were detected using prominence-based peak finding (`scipy.signal.find_peaks`) with a relative prominence threshold of 0.1. For each peak, the nearest left and right “zero-crossing” positions were defined as the closest AP bins where the smoothed profile fell to $\le 0$. Stripe boundaries were set symmetrically around the peak center using the smaller of the left/right half-widths, with an optional buffer term (Snakefile parameter `widthBuffer`, set to 0.0 for the current pipeline). Detected stripe ranges were written into the per-window config copy (`results_{max_time}/config.yaml`) for downstream preprocessing.

### Transcription preprocessing and binning

Transcription traces were filtered to nuclei within the detected stripe AP window and then binned into a $5\times 5$ AP×DV spatial grid (script: `scripts/01_preprocess_eve_data.py`). Within each stripe, AP and DV coordinates were normalized to $[0,1]$ based on the observed coordinate range. AP bins were defined with 50% overlap (to reduce discretization artifacts along AP), and DV bins were defined without overlap.

Each nucleus was assigned to all AP×DV bins it fell within under this overlapping-binning scheme. For each spatial bin, the mean MS2 fluorescence trace across assigned nuclei was computed at each timepoint, producing 25 binned transcription traces $F_{a,d}(t)$ (AP bin $a=1..5$, DV bin $d=1..5$). These traces served as the transcription input for degradation inference.

### smFISH spot-to-nucleus assignment (SASS)

For stripes with raw Imaris exports (stripes 3 and 4), spots were assigned to nuclei using the SASS pipeline (`external/sass/spotMe_v2.py` invoked via the Snakemake rule `process_mrna_sass`). This produced a `position_data.txt` table containing spot and nucleus coordinates and identifiers, which was used for downstream aggregation and QC.

### smFISH centering, QC, and binning

smFISH data were aggregated into a $5\times 5$ AP×DV grid using `scripts/04_aggregate_and_validate_mrna.py` (Snakemake rule `bin_mrna_counts`). In the current pipeline, binning was performed in “spot mode” (`no_groupByNuclei=True` in the Snakefile): spatial bin values corresponded to total spot counts per bin rather than mean per-nucleus spot counts. This choice matched project decisions to avoid nucleus-level aggregation during bin-mean construction.

Before binning, the AP axis was automatically centered using the smoothed spot-density peak along the AP coordinate (`spotx`). Spot positions were binned at 0.75-unit resolution and smoothed with a centered moving average (window 25); if the density peak was offset from the field-of-view center, one AP edge was trimmed so the peak aligned with the center.

To ensure comparable nuclei density across datasets, we applied two QC checks derived from the Berrocal MS2 data (script: `scripts/06_compute_validation_thresholds.py`):

1. **Nuclei density check:** nuclei positions were evaluated using a $k$-nearest-neighbor metric ($k=4$), and the median neighbor distance was required to fall within a stripe-specific [min, max] range learned from Berrocal embryos.
2. **Bin nuclei count check:** nuclei were binned using the same overlapping AP/non-overlapping DV scheme as the transcription preprocessing, and the resulting 25-element bin-count distribution was compared to the stripe-specific Berrocal reference distribution.

If the mRNA field-of-view was **too dense** (median NN distance below the stripe minimum), the AP axis was iteratively trimmed symmetrically from both ends (5% per end per iteration, up to 10 iterations) until the NN metric passed. If the mRNA field-of-view was **too sparse** (median NN distance above the stripe maximum), the pipeline iteratively narrowed the stripe AP window in the per-window config file (5% per end per iteration, up to 4 iterations) and recomputed both transcription traces and QC thresholds for consistency. After centering and any QC-driven adjustments, spot counts were binned into the $5\times 5$ grid; missing bins were filled with zero.

## Bayesian inference of degradation

### Common observation structure

For each embryo, the transcription preprocessing yielded 25 binned transcription traces (5 AP × 5 DV), and the smFISH preprocessing yielded a corresponding 25-element vector of binned spot counts. All degradation models were fit separately for each embryo.

Across models, expected mature mRNA at the final timepoint $T$ was computed as a history integral of the transcription input, weighted by a survival kernel derived from the degradation process. Numerical integration used a trapezoidal rule with $\Delta t = 20/60$ min.

Observed binned mRNA values $m_{\mathrm{obs},k}$ were modeled with independent Gaussian noise:

$$m_{\mathrm{obs},k} \sim \mathcal{N}(\hat m_k,\,\sigma)$$

where $k$ indexes the 25 spatial bins.

### Spatial AP-binned constant-degradation model

The primary model allowed degradation to vary along AP while remaining constant across DV within each AP bin (script: `scripts/02_infer_degradation_rates_spatial.py`). For each AP bin $a$, a degradation rate $D_a$ was inferred and shared across the five DV bins within that AP position. Expected mRNA at the final timepoint $T$ for a given trace was computed from the analytical solution of

$$\frac{dm}{dt} = \gamma F(t) - D_a m, \quad m(0)=0$$

which yields

$$\hat m(T) = \gamma \int_0^T F(t)\,\exp\{-D_a\,(T-t)\}\,dt.$$

Priors (as implemented):

- $D_a \sim \mathrm{LogNormal}(\mu=-2.0,\,\sigma=1.0)$ independently for $a=1..5$
- $\gamma \sim \mathrm{HalfNormal}(\sigma=100.0)$
- $\sigma \sim \mathrm{InverseGamma}(\alpha=2,\,\beta=3)$

### Null constant-degradation model

As a baseline, a single constant degradation rate $D_0$ was shared across all spatial bins (script: `scripts/02_infer_degradation_rates_exponential_null.py`). Expected mRNA used the same kernel form with a single rate parameter. Priors (as implemented):

- $D_0 \sim \mathrm{HalfNormal}(\sigma=1.0)$
- $\gamma \sim \mathrm{InverseGamma}(\alpha=2,\,\beta=3)$
- $\sigma \sim \mathrm{InverseGamma}(\alpha=2,\,\beta=3)$

### Age-dependent degradation (Gaussian random walk)

To capture degradation that depends on molecular age rather than spatial position, we inferred an age-indexed degradation profile $D(\tau)$ shared across all spatial bins (script: `scripts/02_infer_degradation_rates_delayed.py`). Survival probability was defined as

$$S(\tau) = \exp\left\{-\int_0^{\tau} D(u)\,du\right\}$$

and expected mRNA at $T$ was computed by the convolution

$$\hat m(T) = \gamma \int_0^T F(t)\,S(T-t)\,dt.$$

The log degradation profile was modeled as a Gaussian random walk:

- $\log D(\tau_k) \sim \mathrm{GaussianRandomWalk}(\sigma=0.1)$ with initial distribution $\mathcal{N}(-3.5, 0.5)$ on the first step.
- $\gamma \sim \mathrm{HalfNormal}(\sigma=500.0)$
- $\sigma \sim \mathrm{InverseGamma}(\alpha=2,\,\beta=3)$

### Biphasic age-dependent degradation via poly-A protection

We also fit a mechanistic age-dependent model in which degradation is modulated by poly-A tail shortening and loss of protection (script: `scripts/02_infer_degradation_rates_biphasic.py`). Poly-A length as a function of age was modeled as $N_A(\tau) = \max\{0, N_{A,0} - r\tau\}$, and a protection factor was defined as $p(\tau) = \tanh(\beta N_A(\tau))$. The resulting degradation profile was

$$D(\tau) = D_{\mathrm{protected}} + D_{\mathrm{unprotected}}\,[1 - p(\tau)].$$

Expected mRNA was computed by the same convolution form as the age-dependent model above. Priors (as implemented):

- $N_{A,0} \sim \mathcal{N}(60, 10)$
- $r \sim \mathrm{HalfNormal}(\sigma=10.0)$
- $\beta \sim \mathcal{N}(0.096, 0.005)$
- $D_{\mathrm{protected}} \sim \mathrm{HalfNormal}(\sigma=0.1)$
- $D_{\mathrm{unprotected}} \sim \mathrm{HalfNormal}(\sigma=1.0)$
- $\gamma \sim \mathrm{HalfNormal}(\sigma=500.0)$
- $\sigma \sim \mathrm{InverseGamma}(\alpha=2,\,\beta=3)$

**Prior distributions:**

| Parameter          | Spatial                  | Null               | Delayed                                | Biphasic            |
| ------------------ | ------------------------ | ------------------ | -------------------------------------- | ------------------- |
| *D* / *D*₀         | LogNormal(−2, 1) per bin | LogNormal(−2, 1)   | GRW on log *D*, init Normal(−3.5, 0.5) | —                   |
| *γ*                | HalfNormal(100)          | HalfNormal(100)    | HalfNormal(500)                        | HalfNormal(500)     |
| *σ* (noise)        | InverseGamma(2, 3)       | InverseGamma(2, 3) | InverseGamma(2, 3)                     | InverseGamma(2, 3)  |
| *N*_A0             | —                        | —                  | —                                      | Normal(70, 10)      |
| Deadenylation rate | —                        | —                  | —                                      | HalfNormal(2.0)     |
| *β*                | —                        | —                  | —                                      | Normal(0.096, 0.02) |
| *D*_protected      | —                        | —                  | —                                      | HalfNormal(0.05)    |
| *D*_unprotected    | —                        | —                  | —                                      | HalfNormal(0.5)     |

### MCMC sampling

All Bayesian models were fit using the No-U-Turn Sampler (NUTS) in PyMC. For the spatial model, sampling used 4 chains with 10,000 draws per chain, target acceptance 0.9, and random seed 42 (configurable via the Snakemake config for sample and chain counts). Sampling was repeated independently for each embryo.

## Model checking and comparison

Posterior predictive checks were performed by recomputing $\hat m(T)$ using posterior mean parameter values and comparing predicted against observed binned mRNA values (script: `scripts/07_validate_MCMC_results.py`). Fit quality was summarized using $R^2$ and root mean squared error (RMSE).

Model comparison was performed using Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO) implemented in ArviZ (script: `scripts/08_compare_models_loo_ppc.py`). Per-observation log-likelihood values were reconstructed from posterior samples and used to compute expected log predictive density (ELPD); Pareto-$k$ diagnostics were used to identify influential observations.

## Software and reproducibility

The workflow was orchestrated with Snakemake (pipeline definition: `Snakefile`) and executed in a conda environment (`envs/analysis.yml`) using Python 3.11. Key packages included NumPy, pandas, SciPy, scikit-learn, PyMC, and ArviZ.

The complete analysis can be reproduced by cloning the repository (including the SASS submodule) and running Snakemake with conda support:

```bash
git clone https://github.com/ManchesterBioinference/spatial_mRNA_decay
cd spatial_mRNA_decay
git submodule update --init
conda env create -f envs/analysis.yml
conda activate spatial-mrna-decay
snakemake --cores 4 --use-conda
```
