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

# Stripes to analyze
STRIPES = config.get("stripes", ["stripe2"])

# Stripe AP coordinate ranges (loaded from config)
STRIPE_RANGES = config.get("stripe_ranges", {})

# Create lookup dictionaries for min/max/center values by stripe
AP_MIN_LOOKUP = {stripe: STRIPE_RANGES[stripe]["min"] for stripe in STRIPES if stripe in STRIPE_RANGES}
AP_MAX_LOOKUP = {stripe: STRIPE_RANGES[stripe]["max"] for stripe in STRIPES if stripe in STRIPE_RANGES}
AP_CENTER_LOOKUP = {stripe: STRIPE_RANGES[stripe]["center"] for stripe in STRIPES if stripe in STRIPE_RANGES}

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
        # Stripe identification (prerequisite for all analyses)
        "results/figures/stripe_identification.png",
        
        # Preprocessing outputs (per stripe - each stripe has its own AP range)
        expand("data/processed_transcription_data/transcription_traces_{stripe}.csv", stripe=STRIPES),
        expand("data/processed_transcription_data/transcription_traces_no_ids_{stripe}.csv", stripe=STRIPES),
        expand("results/figures/transcription_heatmap_{stripe}.png", stripe=STRIPES),
        
        # Inference outputs (per stripe)
        expand("results/{stripe}/chains/degradation_chain.csv", stripe=STRIPES),
        expand("results/{stripe}/figures/mcmc_trace.png", stripe=STRIPES),
        
        # Visualization outputs (per stripe)
        expand("results/{stripe}/figures/degradation_posteriors.png", stripe=STRIPES),
        expand("results/{stripe}/figures/halflife_heatmap.png", stripe=STRIPES),
        expand("results/{stripe}/summary_statistics.csv", stripe=STRIPES)


rule identify_stripe_ranges:
    """
    Identify AP coordinate ranges for all 7 eve stripes from expression data.
    
    This rule must run before preprocessing to populate config.yaml with stripe_ranges.
    It detects peaks in fluorescence intensity along the AP axis and calculates
    symmetric boundaries around each peak.
    
    Outputs:
    - Updated config.yaml with stripe_ranges section
    - Diagnostic plot showing detected stripes and boundaries
    """
    input:
        data="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        script="scripts/00_identify_stripe_ranges.py"
    output:
        plot="results/figures/stripe_identification.png",
        config_updated=touch("config.yaml.updated")  # Timestamp file to track config updates
    params:
        config_file="config.yaml",
        bin_size=0.01,
        prominence=200,
        max_time=MAX_TIME
    conda:
        "envs/analysis.yml"
    log:
        "results/logs/identify_stripe_ranges.log"
    shell:
        """
        python scripts/00_identify_stripe_ranges.py \
            --input {input.data} \
            --config {params.config_file} \
            --output {output.plot} \
            --bin-size {params.bin_size} \
            --prominence {params.prominence} \
            --max-time {params.max_time} \
            2>&1 | tee {log}
        """


rule preprocess_eve_data:
    """
    Preprocess eve transcription data: filter nuclei to specified stripe, bin spatially (5×5 grid),
    and compute locally averaged fluorescence traces.
    
    The wildcard {stripe} determines which stripe to analyze (stripe2, stripe3, stripe4, etc.).
    AP coordinate ranges are loaded from config.yaml stripe_ranges section.
    """
    input:
        data="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        script="scripts/01_preprocess_eve_data.py",
        config_updated="config.yaml.updated"  # Ensures stripe ranges are identified first
    output:
        traces="data/processed_transcription_data/transcription_traces_{stripe}.csv",
        traces_no_ids="data/processed_transcription_data/transcription_traces_no_ids_{stripe}.csv",
        heatmap="results/figures/transcription_heatmap_{stripe}.png"
    params:
        stripe=lambda wildcards: wildcards.stripe,
        ap_min=lambda wildcards: AP_MIN_LOOKUP[wildcards.stripe],
        ap_max=lambda wildcards: AP_MAX_LOOKUP[wildcards.stripe],
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS,
        dv_min=DV_MIN,
        dv_max=DV_MAX,
        max_time=MAX_TIME
    conda:
        "envs/analysis.yml"
    log:
        "results/logs/preprocess_{stripe}.log"
    shell:
        """
        python scripts/01_preprocess_eve_data.py \
            --input {input.data} \
            --output {output.traces} \
            --output-no-ids {output.traces_no_ids} \
            --plot {output.heatmap} \
            --stripe {params.stripe} \
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
    - Transcription: 25 bins × 61 timepoints (5 AP × 5 DV spatial grid) - STRIPE-SPECIFIC
    - mRNA: 25 values (average mRNA/nucleus per spatial bin from smFISH) - STRIPE-SPECIFIC
    
    Model infers 5 degradation rates (one per AP position), each applied to all 5 DV bins.
    
    Wildcards:
    - {stripe}: Which eve stripe to analyze (stripe2, stripe3, stripe4)
    """
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}.csv",
        mrna="data/processed_mRNA_data_{stripe}/e8_9_edgespotsremoved_sass_formodel.csv",
        script="scripts/02_infer_degradation_rates.py"
    output:
        chain="results/{stripe}/chains/degradation_chain.csv",
        trace_plot="results/{stripe}/figures/mcmc_trace.png"
    params:
        n_samples=N_MCMC_SAMPLES,
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    conda:
        "envs/analysis.yml"
    threads: N_MCMC_CHAINS
    log:
        "results/{stripe}/logs/inference.log"
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
    
    Wildcards:
    - {stripe}: Which eve stripe was analyzed (stripe2, stripe3, stripe4)
    """
    input:
        chain="results/{stripe}/chains/degradation_chain.csv",
        script="scripts/03_visualize_results.py"
    output:
        posteriors="results/{stripe}/figures/degradation_posteriors.png",
        heatmap="results/{stripe}/figures/halflife_heatmap.png",
        summary="results/{stripe}/summary_statistics.csv"
    params:
        n_bins=N_AP_BINS  # Number of AP bins (5 degradation rates inferred)
    conda:
        "envs/analysis.yml"
    log:
        "results/{stripe}/logs/visualize.log"
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
        rm -rf results/stripe*/chains/*.csv
        rm -rf results/stripe*/figures/*.png
        rm -rf results/stripe*/summary_statistics*.csv
        rm -rf results/stripe*/logs/*.log
        rm -rf results/figures/transcription_heatmap_*.png
        rm -rf results/figures/stripe_identification.png
        """
