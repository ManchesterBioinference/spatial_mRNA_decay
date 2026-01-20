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

# Embryos per stripe (for mRNA processing)
EMBRYOS = config.get("embryos", {
    "stripe3": ["e1", "e2", "e3", "e4"]
})

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

# Helper function to get all stripe/embryo combinations
def get_embryo_outputs(pattern,ignoreStripes=[]):
    """Generate output paths for all stripe/embryo combinations."""
    outputs = []
    for stripe in STRIPES:
        if stripe in ignoreStripes:
            continue
        if stripe in EMBRYOS:
            for embryo in EMBRYOS[stripe]:
                outputs.append(pattern.format(stripe=stripe, embryo=embryo))
    return outputs

# All output files
rule all:
    input:
        # Stripe identification (prerequisite for all analyses)
        "results/figures/intermediate/transcription/stripe_identification.png",
        
        # Validation thresholds and QC figures (directory containing all detected stripes)
        "results/figures/intermediate/transcription/nucleiDistributions",
        
        # Preprocessing outputs (per stripe - each stripe has its own AP range)
        expand("data/processed_transcription_data/transcription_traces_{stripe}.csv", stripe=STRIPES),
        expand("data/processed_transcription_data/transcription_traces_no_ids_{stripe}.csv", stripe=STRIPES),
        expand("results/figures/intermediate/transcription/transcription_heatmap_{stripe}.png", stripe=STRIPES),
        
        # mRNA processing and validation outputs (per stripe and embryo)
        get_embryo_outputs("data/Ali_embryos/{stripe}/{embryo}/time_data/position_data-intense.txt",['stripe2']),
        get_embryo_outputs("data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",['stripe2']),
        get_embryo_outputs("results/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_heatmap.png",['stripe2']),
        get_embryo_outputs("results/figures/intermediate/mrna/{stripe}/{embryo}_bin_count_ridge.png",['stripe2']),
        get_embryo_outputs("results/{stripe}/validation/{embryo}_nuclei_density_validation.txt",['stripe2']),
        
        # Inference outputs (per embryo)
        get_embryo_outputs("results/{stripe}/{embryo}/chains/degradation_chain.csv"),
        get_embryo_outputs("results/{stripe}/{embryo}/figures/mcmc_trace.png"),
        
        # Visualization outputs (per embryo)
        get_embryo_outputs("results/{stripe}/{embryo}/figures/degradation_posteriors.png"),
        get_embryo_outputs("results/{stripe}/{embryo}/figures/halflife_heatmap.png"),
        get_embryo_outputs("results/{stripe}/{embryo}/summary_statistics.csv")


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
        plot="results/figures/intermediate/transcription/stripe_identification.png",
        config_updated=touch("config.yaml.updated")  # Timestamp file to track config updates
    params:
        config_file="config.yaml",
        bin_size=0.01,
        prominence=200,
        max_time=MAX_TIME,
        widthBuffer=0.0
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
            --widthBuffer {params.widthBuffer} \
            2>&1 | tee {log}
        """


rule compute_validation_thresholds:
    """
    Compute nuclei density validation thresholds from Berrocal_2020 data.
    
    This rule analyzes the transcription data to compute expected nuclei density
    metrics for each stripe. These thresholds are used as QC gates when processing
    mRNA data from imaging.
    
    Thresholds are computed per stripe and stored in config.yaml, including:
    - k-nearest neighbor distance statistics
    - Expected nuclei count per stripe
    - Bin nuclei count distributions (per embryo and combined)
    
    Also generates multipanel figures (one per stripe) showing normalized nuclei 
    placement for each embryo in that stripe. Figures are created for ALL stripes
    detected in the Berrocal data, not just those in STRIPES config.
    """
    input:
        data="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        script="scripts/06_compute_validation_thresholds.py",
        config_updated="config.yaml.updated"  # Ensure stripe ranges are computed first
    output:
        validation_updated=touch("config.yaml.validation_updated"),  # Track validation threshold updates
        figures_dir=directory("results/figures/intermediate/transcription/nucleiDistributions")
    params:
        config_file="config.yaml",
        output_dir="results/figures/intermediate/transcription/nucleiDistributions",
        k=4,
        max_time=MAX_TIME,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    conda:
        "envs/analysis.yml"
    log:
        "results/logs/compute_validation_thresholds.log"
    shell:
        """
        python scripts/06_compute_validation_thresholds.py \
            --input {input.data} \
            --config {params.config_file} \
            --output-dir {params.output_dir} \
            --k {params.k} \
            --max-time {params.max_time} \
            --n-ap-bins {params.n_ap_bins} \
            --n-dv-bins {params.n_dv_bins} \
            2>&1 | tee {log}
        """


rule process_mrna_sass:
    """
    Run SASS spotMe_v2.py to process raw Imaris imaging data.
    
    This rule assigns mRNA spots to nuclei using the SASS (Single-cell Automated 
    Segmentation and Spot-calling) method.
    
    IMPORTANT: Raw imaging data has already undergone manual preprocessing:
    - Images manually centered on the target stripe
    - Edge spots manually removed before export from Imaris
    
    Input directories must contain:
    - nuclei_Statistics/: CSV files with nuclear measurements from Imaris
    - spots_Statistics/: CSV files with mRNA spot measurements from Imaris
    
    Wildcards:
    - {stripe}: Eve stripe (stripe3, stripe4, etc.)
    - {embryo}: Embryo ID (e1, e2, e3, e4, etc.)
    
    """
    input:
        nuclei_dir="data/Ali_embryos/{stripe}/{embryo}/nuclei_Statistics",
        spots_dir="data/Ali_embryos/{stripe}/{embryo}/spots_Statistics",
        #config_updated="config.yaml.updated",
        #validation_updated="config.yaml.validation_updated"  # Ensure thresholds computed
    output:
        position_data="data/Ali_embryos/{stripe}/{embryo}/time_data/position_data-intense.txt"
    conda:
        "envs/analysis.yml"
    log:
        "results/{stripe}/logs/sass_{embryo}.log"
    shell:
        """
        # spotMe_v2.py expects the parent directory containing nuclei_Statistics and spots_Statistics
        embryo_dir=$(dirname {input.nuclei_dir})
        
        python external/sass/spotMe_v2.py "$embryo_dir" \
            2>&1 | tee {log}
        """


rule bin_mrna_counts:
    """
    Aggregate mRNA counts and validate against Berrocal_2020 benchmarks.
    
    This rule:
    1. Bins mRNA spot data into a 5×5 spatial grid matching transcription data
    2. Computes k=4 nearest neighbor metrics for nuclei density validation
    3. Computes bin nuclei counts and validates against transcription distribution
    4. Generates QC visualizations (heatmap, scatter, density, ridge plot)
    5. Writes validation report
    
    Output format: 25 values (no header), one per spatial bin, representing
    average mRNA count per nucleus in each bin.
    
    HARD GATE: Pipeline will fail if validation fails (k-NN distance outside
    expected range OR bin counts significantly deviate from transcription data).
    
    Wildcards:
    - {stripe}: Eve stripe
    - {embryo}: Embryo ID
    """
    input:
        position_data="data/Ali_embryos/{stripe}/{embryo}/time_data/position_data-intense.txt",
        script="scripts/04_aggregate_and_validate_mrna.py",
        validation_updated="config.yaml.validation_updated"  # Ensure thresholds exist
    output:
        binned_mrna="data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        validation_report="results/{stripe}/validation/{embryo}_nuclei_density_validation.txt",
        heatmap="results/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_heatmap.png",
        spatial_xy="results/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_spatial_xy.png",
        density="results/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_expression_density.png",
        nuc_density="results/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_nuclei_expression_density.png",
        ridge_plot="results/figures/intermediate/mrna/{stripe}/{embryo}_bin_count_ridge.png"
    params:
        config_file="config.yaml",
        no_groupByNuclei=True
    conda:
        "envs/analysis.yml"
    log:
        "results/{stripe}/logs/bin_mrna_{embryo}.log"
    shell:
        """
        python {input.script} \
            {input.position_data} \
            {wildcards.stripe} \
            {wildcards.embryo} \
            {params.config_file} \
            {output.binned_mrna} \
            {output.validation_report} \
            {output.heatmap} \
            {params.no_groupByNuclei} \
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
        heatmap="results/figures/intermediate/transcription/transcription_heatmap_{stripe}.png"
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
    - mRNA: 25 values (average mRNA/nucleus per spatial bin from smFISH) - EMBRYO-SPECIFIC
    
    Model infers 5 degradation rates (one per AP position), each applied to all 5 DV bins.
    
    Each embryo is analyzed independently to capture biological variability.
    
    Wildcards:
    - {stripe}: Which eve stripe to analyze (stripe2, stripe3, stripe4)
    - {embryo}: Which embryo to analyze (e1, e2, e3, e4)
    """
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}.csv",
        mrna="data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/02_infer_degradation_rates.py"
    output:
        chain="results/{stripe}/{embryo}/chains/degradation_chain.csv",
        trace_plot="results/{stripe}/{embryo}/figures/mcmc_trace.png"
    params:
        n_samples=N_MCMC_SAMPLES,
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    conda:
        "envs/analysis.yml"
    threads: N_MCMC_CHAINS
    log:
        "results/{stripe}/{embryo}/logs/inference.log"
    shell:
        """
        python scripts/02_infer_degradation_rates_new.py \
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


rule test_julia_inference:
    """
    TEMPORARY COMPARISON: Run original Julia inference script to compare with Python version.
    
    This rule runs the original infer_D_across_stripe2.jl script with the same inputs
    as the Python version. Useful for validating the Python implementation.
    
    TO REMOVE THIS COMPARISON:
    - Delete this rule from Snakefile
    - Run: snakemake clean_julia_comparison
    
    REQUIREMENTS:
    - Load Julia BEFORE running snakemake: module load apps/binapps/julia
    - Julia packages will be installed automatically on first run (~10-15 min)
    
    USAGE:
    - module load apps/binapps/julia
    - conda run -n research-assistant snakemake -c4 results/stripe3/e1/julia_comparison/degradation_chain.csv
    """
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}.csv",
        mrna="data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/fromJenny/infer_D_across_stripe2.jl"
    output:
        chain="results/{stripe}/{embryo}/julia_comparison/degradation_chain.csv",
        script_copy="results/{stripe}/{embryo}/julia_comparison/infer_D_modified.jl"
    log:
        "results/{stripe}/{embryo}/julia_comparison/inference_julia.log"
    shell:
        """
        # Check if Julia is available
        if ! command -v julia &> /dev/null; then
            echo "ERROR: Julia not found. Please load Julia module first:"
            echo "  module load apps/binapps/julia"
            echo "Then rerun snakemake"
            exit 1
        fi
        
        # Install Julia packages if not already installed (suppress verbose output)
        # julia envs/setup_julia_packages.jl 2>&1 | head -n 50

        # Create output directory
        mkdir -p $(dirname {output.chain})
        
        # Copy script and inject file paths using sed
        cp {input.script} {output.script_copy}
        
        # Replace empty CSV read paths with actual data paths (preserve closing parens)
        sed -i 's|CSV.File(""; header=false))|CSV.File("{input.transcription}"; header=false))|' {output.script_copy}
        sed -i 's|CSV.File(""; header=false))|CSV.File("{input.mrna}"; header=false))|' {output.script_copy}
        
        # Replace empty CSV write path with output path
        sed -i 's|CSV.write("",df)|CSV.write("{output.chain}",df)|' {output.script_copy}
        
        # Run Julia script from its output directory (plots will save there)
        cd $(dirname {output.script_copy})
        julia infer_D_modified.jl 2>&1 | tee $(basename {log})
        
        # Move log to correct location
        mv $(basename {log}) {log}
        """


rule visualize_results:
    """
    Visualize inference results: degradation rate posteriors and spatial heatmaps of mRNA half-life.
    
    Creates:
    - Posterior distributions for 5 AP-specific degradation rates
    - Spatial heatmap showing half-life variation along AP axis
    
    Each embryo is visualized independently to capture biological variability.
    
    Wildcards:
    - {stripe}: Which eve stripe was analyzed (stripe2, stripe3, stripe4)
    - {embryo}: Which embryo was analyzed (e1, e2, e3, e4)
    """
    input:
        chain="results/{stripe}/{embryo}/chains/degradation_chain.csv",
        script="scripts/03_visualize_results.py"
    output:
        posteriors="results/{stripe}/{embryo}/figures/degradation_posteriors.png",
        heatmap="results/{stripe}/{embryo}/figures/halflife_heatmap.png",
        summary="results/{stripe}/{embryo}/summary_statistics.csv"
    params:
        n_bins=N_AP_BINS  # Number of AP bins (5 degradation rates inferred)
    conda:
        "envs/analysis.yml"
    log:
        "results/{stripe}/{embryo}/logs/visualize.log"
    shell:
        """
        python scripts/03_visualize_results.py \
            --chain {input.chain} \
            --posteriors-plot {output.posteriors} \
            --heatmap-plot {output.heatmap} \
            --summary {output.summary} \
            --n-ap-bins {params.n_bins} \
            2>&1 | tee {log}
        """


rule clean:
    """Remove all generated files."""
    shell:
        """
        rm -rf data/processed_transcription_data/transcription_traces_*.csv
        rm -rf results/stripe*/*/chains/*.csv
        rm -rf results/stripe*/*/figures/*.png
        rm -rf results/stripe*/*/summary_statistics*.csv
        rm -rf results/stripe*/*/logs/*.log
        rm -rf results/stripe*/validation/*.txt
        rm -rf results/figures/intermediate/transcription/transcription_heatmap_*.png
        rm -rf results/figures/intermediate/transcription/stripe_identification.png
        rm -rf results/figures/intermediate/mrna/*/*.png
        """

rule clean_julia_comparison:
    """Remove Julia comparison outputs only."""
    shell:
        """
        rm -rf results/*/*/julia_comparison/
        """
