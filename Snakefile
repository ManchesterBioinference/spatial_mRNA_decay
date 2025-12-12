"""
Snakemake pipeline for spatial mRNA decay analysis.

Based on Jenny's handoff instructions for eve transcription data processing.

Pipeline stages:
1. Preprocess eve transcription data (filter, bin spatially, compute local averages)
2. Infer spatially-varying degradation rates using Bayesian inference
3. Visualize results (posterior distributions and spatial heatmaps)

Usage:
    conda run -n research-assistant snakemake --cores 4 --use-conda   # run the full pipeline
    conda run -n research-assistant snakemake --cores 4 --use-conda  --dry-run  # preview what will run
"""

import os

# Configuration
configfile: "config.yaml"

# AP threshold ranges to analyze
# These encode min/max AP coordinates: "033045" = AP range [0.33, 0.45]
AP_RANGES = config.get("ap_ranges", [
    {"id": "033045", "min": 0.33, "max": 0.45}
])

# Extract just the IDs for wildcard expansion
AP_RANGE_IDS = [ap_range["id"] for ap_range in AP_RANGES]

# Create lookup dictionaries for min/max values
AP_MIN_LOOKUP = {ap_range["id"]: ap_range["min"] for ap_range in AP_RANGES}
AP_MAX_LOOKUP = {ap_range["id"]: ap_range["max"] for ap_range in AP_RANGES}

# Spatial binning parameters
N_AP_BINS = config.get("n_ap_bins", 5)
N_DV_BINS = config.get("n_dv_bins", 5)
DV_MIN = config.get("dv_min", -1.0)
DV_MAX = config.get("dv_max", -1.0)
MAX_TIME = config.get("max_time_seconds", 1200)

# MCMC parameters
N_MCMC_SAMPLES = config.get("n_mcmc_samples", 10000)
N_MCMC_CHAINS = config.get("n_mcmc_chains", 4)

# All output files
rule all:
    input:
        # Preprocessing outputs
        expand("data/processed_transcription_data/transcription_traces_{ap_range}.csv", ap_range=AP_RANGE_IDS),
        expand("data/processed_transcription_data/transcription_traces_no_ids_{ap_range}.csv", ap_range=AP_RANGE_IDS),
        expand("results/figures/transcription_heatmap_{ap_range}.png", ap_range=AP_RANGE_IDS),
        
        # Inference outputs
        expand("results/chains/degradation_chain_{ap_range}.csv", ap_range=AP_RANGE_IDS),
        expand("results/figures/mcmc_trace_{ap_range}.png", ap_range=AP_RANGE_IDS),
        
        # Visualization outputs
        expand("results/figures/degradation_posteriors_{ap_range}.png", ap_range=AP_RANGE_IDS),
        expand("results/figures/halflife_heatmap_{ap_range}.png", ap_range=AP_RANGE_IDS),
        expand("results/summary_statistics_{ap_range}.csv", ap_range=AP_RANGE_IDS)


rule preprocess_eve_data:
    """
    Preprocess eve transcription data: filter nuclei in stripe 2, bin spatially (5×5 grid),
    and compute locally averaged fluorescence traces.
    
    The wildcard {ap_range} encodes the AP coordinate range:
    - "033044" = AP range [0.33, 0.44]
    - "033045" = AP range [0.33, 0.45] 
    - "033046" = AP range [0.33, 0.46]
    """
    input:
        data="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        script="scripts/01_preprocess_eve_data.py"
    output:
        traces="data/processed_transcription_data/transcription_traces_{ap_range}.csv",
        traces_no_ids="data/processed_transcription_data/transcription_traces_no_ids_{ap_range}.csv",
        heatmap="results/figures/transcription_heatmap_{ap_range}.png"
    params:
        ap_min=lambda wildcards: AP_MIN_LOOKUP[wildcards.ap_range],
        ap_max=lambda wildcards: AP_MAX_LOOKUP[wildcards.ap_range],
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS,
        dv_min=DV_MIN,
        dv_max=DV_MAX,
        max_time=MAX_TIME
    conda:
        "envs/analysis.yml"
    log:
        "results/logs/preprocess_{ap_range}.log"
    shell:
        """
        python scripts/01_preprocess_eve_data.py \
            --input {input.data} \
            --output {output.traces} \
            --output-no-ids {output.traces_no_ids} \
            --plot {output.heatmap} \
            --ap-min {params.ap_min} \
            --ap-max {params.ap_max} \
            --n-ap-bins {params.n_ap_bins} \
            --n-dv-bins {params.n_dv_bins} \
            --dv-min {params.dv_min} \
            --dv-max {params.dv_max} \
            --max-time {params.max_time} \
            2>&1 | tee {log}
        """


rule infer_degradation_rates:
    """
    Infer spatially-varying mRNA degradation rates using Bayesian inference.
    
    Fits ODE model: dm/dt = γ*F(t) - D*m
    where F(t) is transcription input and D is degradation rate (varies by AP position).
    
    Uses preprocessed data:
    - Transcription: 25 bins × 61 timepoints (5 AP × 5 DV spatial grid)
    - mRNA: 25 values (average mRNA/nucleus per spatial bin from smFISH)
    
    Model infers 5 degradation rates (one per AP position), each applied to all 5 DV bins.
    
    The wildcard {ap_range} specifies which AP coordinate range was used for preprocessing.
    """
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{ap_range}.csv",
        mrna="data/processed_mRNA_data/e8_9_edgespotsremoved_sass_formodel.csv",
        script="scripts/02_infer_degradation_rates.py"
    output:
        chain="results/chains/degradation_chain_{ap_range}.csv",
        trace_plot="results/figures/mcmc_trace_{ap_range}.png"
    params:
        n_samples=N_MCMC_SAMPLES,
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    conda:
        "envs/analysis.yml"
    threads: N_MCMC_CHAINS
    log:
        "results/logs/inference_{ap_range}.log"
    shell:
        """
        python scripts/02_infer_degradation_rates.py \
            --transcription {input.transcription} \
            --mrna {input.mrna} \
            --output-chain {output.chain} \
            --output-trace {output.trace_plot} \
            --n-samples {params.n_samples} \
            --n-chains {params.n_chains} \
            --n-ap-bins {params.n_ap_bins} \
            --n-dv-bins {params.n_dv_bins} \
            2>&1 | tee {log}
        """


rule visualize_results:
    """
    Visualize inference results: degradation rate posteriors and spatial heatmaps of mRNA half-life.
    
    Creates:
    - Posterior distributions for 5 AP-specific degradation rates
    - Spatial heatmap showing half-life variation along AP axis
    
    The wildcard {ap_range} specifies which AP coordinate range was analyzed.
    """
    input:
        chain="results/chains/degradation_chain_{ap_range}.csv",
        script="scripts/03_visualize_results.py"
    output:
        posteriors="results/figures/degradation_posteriors_{ap_range}.png",
        heatmap="results/figures/halflife_heatmap_{ap_range}.png",
        summary="results/summary_statistics_{ap_range}.csv"
    params:
        n_bins=N_AP_BINS  # Number of AP bins (5 degradation rates inferred)
    conda:
        "envs/analysis.yml"
    log:
        "results/logs/visualize_{ap_range}.log"
    shell:
        """
        python scripts/03_visualize_results.py \
            --chain {input.chain} \
            --posteriors-plot {output.posteriors} \
            --heatmap-plot {output.heatmap} \
            --summary {output.summary} \
            --n-bins {params.n_bins} \
            2>&1 | tee {log}
        """


rule clean:
    """Remove all generated files."""
    shell:
        """
        rm -rf data/processed_transcription_data/transcription_traces_*.csv
        rm -rf results/chains/*.csv
        rm -rf results/figures/*.png
        rm -rf results/summary_statistics_*.csv
        rm -rf results/logs/*.log
        """
