# Spatial mRNA Decay

Infers spatially-varying mRNA degradation rates from *eve* stripe transcription and smFISH data using Bayesian (PyMC) inference. Four degradation models are compared: spatial ODE, null (constant), age-dependent (Gaussian random walk), and biphasic (poly-A).

## Setup

```bash
# Create and activate the Python environment
uv sync
```

## Pipelines

### Main analysis (`Snakefile`)

Processes transcription data, bins mRNA counts, runs MCMC inference for all four models, and produces posterior predictive checks and model comparisons.

```bash
# Dry-run (preview)
snakemake --cores 1 --dry-run

# Local run
snakemake --cores 4

# SLURM (multicore partition, 36 cores)
sbatch slurm.sh
```

### AP-shift sensitivity (`Snakefile_shift_sensitivity`)

Repeats the full analysis with the stripe AP window shifted ±half a bin width in each direction, to verify that spatial misalignment between the Berrocal transcription data and smFISH imaging does not bias the inferred degradation rates.

**Prerequisite:** the main Snakefile must have completed stripe identification first (`results_{max_time}/config.yaml.updated` must exist).

```bash
# Dry-run (preview)
snakemake --snakefile Snakefile_shift_sensitivity --cores 1 --dry-run

# Local run
snakemake --snakefile Snakefile_shift_sensitivity --cores 4
```

Outputs are written to `shiftLeft/` and `shiftRight/` directories, keeping them separate from the main results.