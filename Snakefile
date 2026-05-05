"""
Snakemake pipeline for spatial mRNA decay analysis.

Based on Jenny's handoff instructions for eve transcription data processing.

Pipeline stages:
1. Preprocess eve transcription data (filter, bin spatially, compute local averages)
2. Infer spatially-varying degradation rates using Bayesian inference
3. Visualize results (posterior distributions and spatial heatmaps)

Usage:
    source .venv/bin/activate && snakemake --cores 4         # run the full pipeline
    source .venv/bin/activate && snakemake --cores 4 --dry-run  # preview what will run
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
#MAX_TIME = config.get("max_time_seconds", 1200)
# Loop over max_time values from 1200 to 3000 by 100
# These will be used to create results_{max_time}/... output trees
MAX_TIME_VALUES = [1200] #list(range(1100, 1301, 20))# + list(range(1200, 1801, 100))

# MCMC parameters
N_MCMC_SAMPLES = config.get("n_mcmc_samples", 10000)
N_MCMC_CHAINS = config.get("n_mcmc_chains", 4)

# Helper function to get all stripe/embryo combinations
def get_embryo_outputs(pattern, ignoreStripes=[], max_times=None):
    """Generate output paths for all stripe/embryo combinations.

    If `max_times` is provided, pattern may include a `{max_time}` placeholder
    which will be filled for each value in `max_times`.
    """
    outputs = []
    for stripe in STRIPES:
        if stripe in ignoreStripes:
            continue
        if stripe in EMBRYOS:
            for embryo in EMBRYOS[stripe]:
                if max_times is None:
                    outputs.append(pattern.format(stripe=stripe, embryo=embryo))
                else:
                    for mt in max_times:
                        outputs.append(pattern.format(max_time=mt, stripe=stripe, embryo=embryo))
    return outputs

# All output files
rule all:
    input:
        # Per-max_time config files (prerequisite for everything)
        expand("results_{max_time}/config.yaml", max_time=MAX_TIME_VALUES),
        "data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        
        # Stripe identification (prerequisite for all analyses)
        expand("results_{max_time}/figures/intermediate/transcription/stripe_identification.pdf", max_time=MAX_TIME_VALUES),

        # Validation thresholds and QC figures (directory containing all detected stripes)
        expand("results_{max_time}/figures/intermediate/transcription/nucleiDistributions", max_time=MAX_TIME_VALUES),

        # Preprocessing outputs (per stripe - each stripe has its own AP range)
        expand("data/processed_transcription_data/transcription_traces_{stripe}_{max_time}.csv", stripe=STRIPES, max_time=MAX_TIME_VALUES),
        expand("data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv", stripe=STRIPES, max_time=MAX_TIME_VALUES),

        # mRNA processing and validation outputs (per stripe and embryo)
        get_embryo_outputs("data/Ali_embryos/{stripe}/{embryo}/time_data/position_data.txt", []),
        get_embryo_outputs("results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv", [], max_times=MAX_TIME_VALUES),

        # results under different max_time values
        get_embryo_outputs("results_{max_time}/figures/intermediate/mrna/{stripe}/{embryo}_bin_count_ridge.pdf", [], max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/validation/{embryo}_nuclei_density_validation.txt", [], max_times=MAX_TIME_VALUES),

        # Spatial model inference outputs
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/figures/mcmc_trace.pdf", max_times=MAX_TIME_VALUES),

        # Null (exponential) model inference outputs
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_null/chains/degradation_chain.csv", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_null/figures/mcmc_trace.pdf", max_times=MAX_TIME_VALUES),

        # Delayed (GRW) model inference outputs
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_age/chains/degradation_chain.csv", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_age/figures/mcmc_trace.pdf", max_times=MAX_TIME_VALUES),

        # Biphasic (poly-A) model inference outputs
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_biphasic/chains/degradation_chain.csv", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_biphasic/figures/mcmc_trace.pdf", max_times=MAX_TIME_VALUES),

        # Visualization outputs (spatial model)
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/figures/degradation_posteriors.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/figures/halflife_heatmap.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/figures/spatial_overview.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/summary_statistics.csv", max_times=MAX_TIME_VALUES),

        # Visualization outputs (null, age, biphasic models)
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_null/figures/degradation_vs_age.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_null/summary_statistics.csv", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_age/figures/degradation_vs_age.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_age/summary_statistics.csv", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_biphasic/figures/degradation_vs_age.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_biphasic/summary_statistics.csv", max_times=MAX_TIME_VALUES),
        # Validation outputs (all four models)
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/figures/posterior_predictive_check.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_null/figures/posterior_predictive_check.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_age/figures/posterior_predictive_check.pdf", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}_biphasic/figures/posterior_predictive_check.pdf", max_times=MAX_TIME_VALUES),

        # Combined PPC plot (all four models in one figure)
        get_embryo_outputs("results_{max_time}/comparison/{stripe}/{embryo}/posterior_predictive_check_combined.pdf", max_times=MAX_TIME_VALUES),

        # Model comparison and LOO diagnostics
        get_embryo_outputs("results_{max_time}/comparison/{stripe}/{embryo}/comparison_summary.txt", max_times=MAX_TIME_VALUES),
        get_embryo_outputs("results_{max_time}/comparison/{stripe}/{embryo}/loo_diagnostics_summary.txt", max_times=MAX_TIME_VALUES),

        # Simulation
        get_embryo_outputs("results_{max_time}/{stripe}/{embryo}/figures/counterfactual_spatial.pdf", max_times=[1200])

        #expand("results_{max_time}/figures/intermediate/transcription/transcription_heatmap_{stripe}.pdf", stripe=STRIPES, max_time=MAX_TIME_VALUES),
        #get_embryo_outputs("results_{max_time}/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_heatmap.pdf", [], max_times=MAX_TIME_VALUES),



rule copy_config:
    """
    Copy base config.yaml to results_{max_time}/ directory.
    
    Each max_time value gets its own config copy to avoid conflicts when
    rules modify the config file (e.g., adding stripe_ranges or validation_thresholds).
    
    This rule must run before any other rules that read or modify config.yaml.
    """
    input:
        "config.yaml"
    output:
        "results_{max_time}/config.yaml"
    log:
        "results_{max_time}/logs/copy_config.log"
    shell:
        """
        mkdir -p $(dirname {output})
        cp {input} {output} 2>&1 | tee {log}
        """


rule download_berrocal_data:
    """
    Download and process Berrocal_2020 supplementary data using Berrocal filtering notebook.

    Downloads the eLife 61635 supplementary zip from:
        https://cdn.elifesciences.org/articles/61635/elife-61635-supp1-v2.zip

    Extracts only Data/ CSVs and Berrocal_2020.ipynb (Figures/ and Movies/
    are discarded).  Then executes the first 17 notebook cells plus a new
    save cell to produce the FILTERED CSV used by all downstream rules.
    """
    output:
        filtered_csv="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        raw_csv="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520.csv",
        notebook="data/Berrocal_2020/Berrocal_2020.ipynb"
    log:
        "logs/download_berrocal_data.log"
    shell:
        """
        mkdir -p $(dirname {log})
        MPLBACKEND=Agg uv run scripts/00_download_berrocal_data.py \
            --output-dir data/Berrocal_2020 \
            2>&1 | tee {log}
        """


rule identify_stripe_ranges:
    """
    Identify AP coordinate ranges for all 7 eve stripes from expression data.
    
    This rule must run before preprocessing to populate config.yaml with stripe_ranges.
    It detects peaks in fluorescence intensity along the AP axis and calculates
    symmetric boundaries around each peak.
    
    Outputs:
    - Updated per-max_time config.yaml with stripe_ranges section
    - Diagnostic plot showing detected stripes and boundaries
    """
    input:
        data="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        config=ancient("results_{max_time}/config.yaml"),
        script="scripts/00_identify_stripe_ranges.py"
    output:
        plot="results_{max_time}/figures/intermediate/transcription/stripe_identification.pdf",
        config_updated=touch("results_{max_time}/config.yaml.updated")  # Timestamp file to track config updates
    params:
        config_file=lambda wildcards: f"results_{wildcards.max_time}/config.yaml",
        bin_size=0.005,
        relativeProminence=0.1,
        max_time=lambda wildcards: int(wildcards.max_time),
        widthBuffer=0.0,
        window_length=17,
        smooth_method='savgol'#moving_average' # '
    log:
        "results_{max_time}/logs/identify_stripe_ranges.log"
    shell:
        """
        MPLBACKEND=Agg python scripts/00_identify_stripe_ranges.py \
            --input {input.data} \
            --config {params.config_file} \
            --output {output.plot} \
            --bin-size {params.bin_size} \
            --relativeProminence {params.relativeProminence} \
            --max-time {params.max_time} \
            --widthBuffer {params.widthBuffer} \
            --window-length {params.window_length} \
            --smooth-method {params.smooth_method} \
            2>&1 | tee {log}
        """


rule compute_validation_thresholds:
    """
    Compute nuclei density validation thresholds from Berrocal_2020 data.
    
    This rule analyzes the transcription data to compute expected nuclei density
    metrics for each stripe. These thresholds are used as QC gates when processing
    mRNA data from imaging.
    
    Thresholds are computed per stripe and stored in the per-max_time config.yaml, including:
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
        config_updated="results_{max_time}/config.yaml.updated"  # Ensure stripe ranges are computed first
    output:
        validation_updated=touch("results_{max_time}/config.yaml.validation_updated"),  # Track validation threshold updates
        figures_dir=directory("results_{max_time}/figures/intermediate/transcription/nucleiDistributions")
    params:
        config_file=lambda wildcards: f"results_{wildcards.max_time}/config.yaml",
        output_dir="results_{max_time}/figures/intermediate/transcription/nucleiDistributions",
        k=4,
        max_time=lambda wildcards: int(wildcards.max_time),
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    log:
        "results_{max_time}/logs/compute_validation_thresholds.log"
    shell:
        """
        MPLBACKEND=Agg python scripts/06_compute_validation_thresholds.py \
            --input {input.data} \
            --config {params.config_file} \
            --output-dir {params.output_dir} \
            --k {params.k} \
            --max-time {params.max_time} \
            --n-ap-bins {params.n_ap_bins} \
            --n-dv-bins {params.n_dv_bins} \
            2>&1 | tee {log}
        """

rule initialize_submodule:
    """
    Initialize the git submodule to ensure external/sass/spotMe_v2.py is available.
    """
    output:
        "external/sass/spotMe_v2.py"
    shell:
        "git submodule update --init"


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
        script="external/sass/spotMe_v2.py"
        #config_updated="config.yaml.updated",
        #validation_updated="config.yaml.validation_updated"  # Ensure thresholds computed
    output:
        position_data="data/Ali_embryos/{stripe}/{embryo}/time_data/position_data.txt"
    log:
        "results/{stripe}/logs/sass_{embryo}.log"
    shell:
        """
        # spotMe_v2.py expects the parent directory containing nuclei_Statistics and spots_Statistics
        embryo_dir=$(dirname {input.nuclei_dir})
        
        python external/sass/spotMe_v2.py "$embryo_dir" \
            2>&1 | tee {log}
        """
        # if {wildcard.stripe} == 'stripe2':
        #     touch {output.position_data}
        # fi


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
        position_data="data/Ali_embryos/{stripe}/{embryo}/time_data/position_data.txt",
        script="scripts/04_aggregate_and_validate_mrna.py",
        validation_updated="results_{max_time}/config.yaml.validation_updated"  # Ensure thresholds exist
    output:
        binned_mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        validation_report="results_{max_time}/{stripe}/validation/{embryo}_nuclei_density_validation.txt",
        heatmap="results_{max_time}/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_heatmap.pdf",
        spatial_xy="results_{max_time}/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_spatial_xy.pdf",
        density="results_{max_time}/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_expression_density.pdf",
        nuc_density="results_{max_time}/figures/intermediate/mrna/{stripe}/{embryo}_sass_formodel_nuclei_expression_density.pdf",
        ridge_plot="results_{max_time}/figures/intermediate/mrna/{stripe}/{embryo}_bin_count_ridge.pdf"
    params:
        config_file=lambda wildcards: f"results_{wildcards.max_time}/config.yaml",
        no_groupByNuclei=True
    log:
        "results_{max_time}/{stripe}/logs/bin_mrna_{embryo}.log"
    shell:
        """
        MPLBACKEND=Agg python {input.script} \
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
    AP coordinate ranges are loaded from per-max_time config.yaml stripe_ranges section.
    """
    input:
        data="data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv",
        config=ancient("results_{max_time}/config.yaml"),
        script="scripts/01_preprocess_eve_data.py",
        config_updated="results_{max_time}/config.yaml.updated"  # Ensures stripe ranges are identified first
    output:
        traces="data/processed_transcription_data/transcription_traces_{stripe}_{max_time}.csv",
        traces_no_ids="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        heatmap="results_{max_time}/figures/intermediate/transcription/transcription_heatmap_{stripe}.pdf"
    params:
        config_file=lambda wildcards: f"results_{wildcards.max_time}/config.yaml",
        stripe=lambda wildcards: wildcards.stripe,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS,
        dv_min=DV_MIN,
        dv_max=DV_MAX,
        max_time=lambda wildcards: int(wildcards.max_time)
    log:
        "results_{max_time}/logs/preprocess_{stripe}.log"
    shell:
        """
        MPLBACKEND=Agg python scripts/01_preprocess_eve_data.py \
            --input {input.data} \
            --config {params.config_file} \
            --output {output.traces} \
            --output-no-ids {output.traces_no_ids} \
            --plot {output.heatmap} \
            --stripe {params.stripe} \
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
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/02_infer_degradation_rates_spatial.py"
    output:
        chain="results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv",
        trace_plot="results_{max_time}/{stripe}/{embryo}/figures/mcmc_trace.pdf"
    params:
        n_samples=N_MCMC_SAMPLES,
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    threads: N_MCMC_CHAINS
    log:
        "results_{max_time}/{stripe}/{embryo}/logs/inference.log"
    shell:
        """
        python scripts/02_infer_degradation_rates_spatial.py \
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


rule infer_null_constant_degradation:
    """
    Infer constant (null model) mRNA degradation rate using Bayesian inference.
    
    This rule fits a NULL MODEL where degradation rate is CONSTANT across all
    molecular ages and spatial positions. Serves as baseline for model comparison.
    
    Model: D(age) = D₀ (constant for all ages and positions)
    
    Uses same preprocessed data as age-dependent model:
    - Transcription: 25 bins × timepoints (5 AP × 5 DV spatial grid)
    - mRNA: 25 values (average mRNA/nucleus per spatial bin)
    
    Outputs single constant degradation rate D₀ shared by all observations.
    Results are saved to {stripe}/{embryo}_null/ directory to distinguish
    from age-dependent model outputs.
    
    Wildcards:
    - {stripe}: Which eve stripe to analyze (stripe2, stripe3, stripe4)
    - {embryo}: Which embryo to analyze (e1, e2, e3, e4)
    """
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/02_infer_degradation_rates_exponential_null.py"
    output:
        chain="results_{max_time}/{stripe}/{embryo}_null/chains/degradation_chain.csv",
        trace_plot="results_{max_time}/{stripe}/{embryo}_null/figures/mcmc_trace.pdf"
    params:
        n_samples=N_MCMC_SAMPLES,
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    threads: N_MCMC_CHAINS
    log:
        "results_{max_time}/{stripe}/{embryo}_null/logs/inference.log"
    shell:
        """
        python {input.script} \
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


rule infer_delayed_age_degradation:
    """
    Infer age-dependent mRNA degradation rates using Bayesian inference (GRW model).

    Fits a model where D varies as a Gaussian Random Walk over molecular age bins.
    Results are saved to {stripe}/{embryo}_age/ directory.

    Wildcards:
    - {stripe}: Which eve stripe to analyze
    - {embryo}: Which embryo to analyze (base name, without suffix)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/02_infer_degradation_rates_delayed.py"
    output:
        chain="results_{max_time}/{stripe}/{embryo}_age/chains/degradation_chain.csv",
        trace_plot="results_{max_time}/{stripe}/{embryo}_age/figures/mcmc_trace.pdf"
    params:
        n_samples=N_MCMC_SAMPLES,
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    threads: N_MCMC_CHAINS
    log:
        "results_{max_time}/{stripe}/{embryo}_age/logs/inference.log"
    shell:
        """
        python {input.script} \\
            --transcription {input.transcription} \\
            --mrna {input.mrna} \\
            --output-chain {output.chain} \\
            --output-trace {output.trace_plot} \\
            --n-samples {params.n_samples} \\
            --n-chains {params.n_chains} \\
            --n-ap-bins {params.n_ap_bins} \\
            --n-dv-bins {params.n_dv_bins} \\
            2>&1 | tee {log}
        """


rule infer_biphasic_degradation:
    """
    Infer biphasic mRNA degradation rates using Bayesian inference (mechanistic poly-A model).

    Fits a mechanistic model capturing slow decay while poly-A protected, then fast decay
    after decapping. Results are saved to {stripe}/{embryo}_biphasic/ directory.

    Wildcards:
    - {stripe}: Which eve stripe to analyze
    - {embryo}: Which embryo to analyze (base name, without suffix)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/02_infer_degradation_rates_biphasic.py"
    output:
        chain="results_{max_time}/{stripe}/{embryo}_biphasic/chains/degradation_chain.csv",
        trace_plot="results_{max_time}/{stripe}/{embryo}_biphasic/figures/mcmc_trace.pdf"
    params:
        n_samples=N_MCMC_SAMPLES,
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    threads: N_MCMC_CHAINS
    log:
        "results_{max_time}/{stripe}/{embryo}_biphasic/logs/inference.log"
    shell:
        """
        python {input.script} \\
            --transcription {input.transcription} \\
            --mrna {input.mrna} \\
            --output-chain {output.chain} \\
            --output-trace {output.trace_plot} \\
            --n-samples {params.n_samples} \\
            --n-chains {params.n_chains} \\
            --n-ap-bins {params.n_ap_bins} \\
            --n-dv-bins {params.n_dv_bins} \\
            2>&1 | tee {log}
        """


rule visualize_results:
    """
    Visualize spatial (binned) degradation inference results.

    Creates:
    - Per-bin posterior histograms of degradation rate D
    - AP×DV heatmap of mRNA half-lives
    - Combined spatial overview (histograms + heatmap)

    Each embryo is visualized independently to capture biological variability.

    Wildcards:
    - {stripe}: Which eve stripe was analyzed (stripe2, stripe3, stripe4)
    - {embryo}: Which embryo was analyzed (e8_9um, e7um, e9_10um, etc.)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        chain="results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv",
        script="scripts/03_visualize_results.py"
    output:
        degradation="results_{max_time}/{stripe}/{embryo}/figures/degradation_posteriors.pdf",
        halflife="results_{max_time}/{stripe}/{embryo}/figures/halflife_heatmap.pdf",
        overview="results_{max_time}/{stripe}/{embryo}/figures/spatial_overview.pdf",
        summary="results_{max_time}/{stripe}/{embryo}/summary_statistics.csv"
    log:
        "results_{max_time}/{stripe}/{embryo}/logs/visualize.log"
    shell:
        """
        python scripts/03_visualize_results.py \
            --chain {input.chain} \
            --degradation-plot {output.degradation} \
            --halflife-plot {output.halflife} \
            --overview-plot {output.overview} \
            --summary {output.summary} \
            --model-type spatial \
            2>&1 | tee {log}
        """


rule visualize_null_results:
    """Visualize null constant degradation inference results."""
    input:
        chain="results_{max_time}/{stripe}/{embryo}_null/chains/degradation_chain.csv",
        script="scripts/03_visualize_results.py"
    output:
        degradation="results_{max_time}/{stripe}/{embryo}_null/figures/degradation_vs_age.pdf",
        halflife="results_{max_time}/{stripe}/{embryo}_null/figures/halflife_vs_age.pdf",
        overview="results_{max_time}/{stripe}/{embryo}_null/figures/overview.pdf",
        summary="results_{max_time}/{stripe}/{embryo}_null/summary_statistics.csv"
    log:
        "results_{max_time}/{stripe}/{embryo}_null/logs/visualize.log"
    shell:
        """
        python scripts/03_visualize_results.py \\
            --chain {input.chain} \\
            --degradation-plot {output.degradation} \\
            --halflife-plot {output.halflife} \\
            --overview-plot {output.overview} \\
            --summary {output.summary} \\
            2>&1 | tee {log}
        """


rule visualize_age_results:
    """Visualize delayed-age (GRW) degradation inference results."""
    input:
        chain="results_{max_time}/{stripe}/{embryo}_age/chains/degradation_chain.csv",
        script="scripts/03_visualize_results.py"
    output:
        degradation="results_{max_time}/{stripe}/{embryo}_age/figures/degradation_vs_age.pdf",
        halflife="results_{max_time}/{stripe}/{embryo}_age/figures/halflife_vs_age.pdf",
        overview="results_{max_time}/{stripe}/{embryo}_age/figures/overview.pdf",
        summary="results_{max_time}/{stripe}/{embryo}_age/summary_statistics.csv"
    log:
        "results_{max_time}/{stripe}/{embryo}_age/logs/visualize.log"
    shell:
        """
        python scripts/03_visualize_results.py \\
            --chain {input.chain} \\
            --degradation-plot {output.degradation} \\
            --halflife-plot {output.halflife} \\
            --overview-plot {output.overview} \\
            --summary {output.summary} \\
            2>&1 | tee {log}
        """


rule visualize_biphasic_results:
    """Visualize biphasic degradation inference results."""
    input:
        chain="results_{max_time}/{stripe}/{embryo}_biphasic/chains/degradation_chain.csv",
        script="scripts/03_visualize_results.py"
    output:
        degradation="results_{max_time}/{stripe}/{embryo}_biphasic/figures/degradation_vs_age.pdf",
        halflife="results_{max_time}/{stripe}/{embryo}_biphasic/figures/halflife_vs_age.pdf",
        overview="results_{max_time}/{stripe}/{embryo}_biphasic/figures/overview.pdf",
        summary="results_{max_time}/{stripe}/{embryo}_biphasic/summary_statistics.csv"
    log:
        "results_{max_time}/{stripe}/{embryo}_biphasic/logs/visualize.log"
    shell:
        """
        python scripts/03_visualize_results.py \\
            --chain {input.chain} \\
            --degradation-plot {output.degradation} \\
            --halflife-plot {output.halflife} \\
            --overview-plot {output.overview} \\
            --summary {output.summary} \\
            2>&1 | tee {log}
        """


rule validate_mcmc:
    """
    Validate MCMC inference results with posterior predictive checks.
    
    Creates:
    - Posterior predictive check plot comparing observed vs predicted mRNA
    - Residual plot to check for systematic biases
    - Fit statistics (R², RMSE)
    
    This validates that the age-dependent degradation model accurately
    reproduces the observed mRNA data.
    
    Wildcards:
    - {stripe}: Which eve stripe was analyzed
    - {embryo}: Which embryo was analyzed (must not end in _null; use validate_null_mcmc for null model)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        chain="results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/07_validate_MCMC_results.py"
    output:
        validation="results_{max_time}/{stripe}/{embryo}/figures/posterior_predictive_check.pdf"
    params:
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    log:
        "results_{max_time}/{stripe}/{embryo}/logs/validate_mcmc.log"
    shell:
        """
        python scripts/07_validate_MCMC_results.py \
            --chain {input.chain} \
            --transcription {input.transcription} \
            --mrna {input.mrna} \
            --n-ap-bins {params.n_ap_bins} \
            --n-dv-bins {params.n_dv_bins} \
            --output {output.validation} \
            2>&1 | tee {log}
        """


rule counterfactual_spatial:
    """
    Counterfactual analysis for the spatial degradation model.

    Compares three predicted mRNA patterns across all AP × DV bins:
      - Fitted:           each AP bin uses its own inferred D_i
      - Fast-everywhere:  all AP bins use the average of the two edge bin D values
      - Slow-everywhere:  all AP bins use the minimum-D bin (centre) D value

    Produces a two-row figure:
      top row  — four side-by-side AP × DV heatmaps (observed, fitted,
                 fast-everywhere, slow-everywhere) with a shared colour scale
      bottom row — bar chart of the inferred D profile with reference lines
                   for the edge-average and centre rates

    Wildcards:
    - {stripe}: Which eve stripe was analyzed
    - {embryo}: Which embryo was analyzed (spatial model; no model-type suffix)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        chain="results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/10_counterfactual_spatial.py"
    output:
        figure="results_{max_time}/{stripe}/{embryo}/figures/counterfactual_spatial.pdf"
    params:
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    log:
        "results_{max_time}/{stripe}/{embryo}/logs/counterfactual_spatial.log"
    shell:
        """
        python scripts/10_counterfactual_spatial.py \\
            --chain {input.chain} \\
            --transcription {input.transcription} \\
            --mrna {input.mrna} \\
            --n-ap-bins {params.n_ap_bins} \\
            --n-dv-bins {params.n_dv_bins} \\
            --output {output.figure} \\
            2>&1 | tee {log}
        """


rule validate_null_mcmc:
    """
    Validate NULL CONSTANT degradation MCMC results with posterior predictive checks.

    Mirrors validate_mcmc but targets the null model outputs stored in
    {stripe}/{embryo}_null/ directories. Uses mRNA data from the base embryo
    (i.e., strips the _null suffix to find the processed mRNA file).

    Creates:
    - Posterior predictive check plot comparing observed vs predicted mRNA
    - Residual plot to check for systematic biases
    - Fit statistics (R², RMSE)

    Wildcards:
    - {stripe}: Which eve stripe was analyzed
    - {embryo}: Base embryo name WITHOUT the _null suffix (e.g., e8_9um, not e8_9um_null)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        chain="results_{max_time}/{stripe}/{embryo}_null/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/07_validate_MCMC_results.py"
    output:
        validation="results_{max_time}/{stripe}/{embryo}_null/figures/posterior_predictive_check.pdf"
    params:
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    log:
        "results_{max_time}/{stripe}/{embryo}_null/logs/validate_null_mcmc.log"
    shell:
        """
        python scripts/07_validate_MCMC_results.py \
            --chain {input.chain} \
            --transcription {input.transcription} \
            --mrna {input.mrna} \
            --n-ap-bins {params.n_ap_bins} \
            --n-dv-bins {params.n_dv_bins} \
            --output {output.validation} \
            2>&1 | tee {log}
        """


rule validate_age_mcmc:
    """
    Validate delayed-age (GRW) MCMC results with posterior predictive checks.

    Mirrors validate_null_mcmc but targets {stripe}/{embryo}_age/ directories.

    Wildcards:
    - {stripe}: Which eve stripe was analyzed
    - {embryo}: Base embryo name WITHOUT the _age suffix
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        chain="results_{max_time}/{stripe}/{embryo}_age/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/07_validate_MCMC_results.py"
    output:
        validation="results_{max_time}/{stripe}/{embryo}_age/figures/posterior_predictive_check.pdf"
    params:
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    log:
        "results_{max_time}/{stripe}/{embryo}_age/logs/validate_mcmc.log"
    shell:
        """
        python scripts/07_validate_MCMC_results.py \\
            --chain {input.chain} \\
            --transcription {input.transcription} \\
            --mrna {input.mrna} \\
            --n-ap-bins {params.n_ap_bins} \\
            --n-dv-bins {params.n_dv_bins} \\
            --output {output.validation} \\
            2>&1 | tee {log}
        """


rule validate_biphasic_mcmc:
    """
    Validate biphasic MCMC results with posterior predictive checks.

    Mirrors validate_null_mcmc but targets {stripe}/{embryo}_biphasic/ directories.

    Wildcards:
    - {stripe}: Which eve stripe was analyzed
    - {embryo}: Base embryo name WITHOUT the _biphasic suffix
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        chain="results_{max_time}/{stripe}/{embryo}_biphasic/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/07_validate_MCMC_results.py"
    output:
        validation="results_{max_time}/{stripe}/{embryo}_biphasic/figures/posterior_predictive_check.pdf"
    params:
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    log:
        "results_{max_time}/{stripe}/{embryo}_biphasic/logs/validate_mcmc.log"
    shell:
        """
        python scripts/07_validate_MCMC_results.py \\
            --chain {input.chain} \\
            --transcription {input.transcription} \\
            --mrna {input.mrna} \\
            --n-ap-bins {params.n_ap_bins} \\
            --n-dv-bins {params.n_dv_bins} \\
            --output {output.validation} \\
            2>&1 | tee {log}
        """

rule validate_mcmc_combined:
    """
    Combined posterior predictive check for all four degradation models.

    Produces a single scatter plot with all four models (spatial, null, age,
    biphasic) overlaid and colored by model. Each model gets its own line of
    best fit; a red dashed x=y line marks perfect prediction.

    Output lives alongside the model-comparison outputs in comparison/.

    Wildcards:
    - {stripe}: Which eve stripe was analyzed
    - {embryo}: Base embryo name WITHOUT any model suffix
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        spatial_chain="results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv",
        null_chain="results_{max_time}/{stripe}/{embryo}_null/chains/degradation_chain.csv",
        age_chain="results_{max_time}/{stripe}/{embryo}_age/chains/degradation_chain.csv",
        biphasic_chain="results_{max_time}/{stripe}/{embryo}_biphasic/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv",
        script="scripts/07_validate_MCMC_combined.py"
    output:
        plot="results_{max_time}/comparison/{stripe}/{embryo}/posterior_predictive_check_combined.pdf"
    params:
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS
    log:
        "results_{max_time}/comparison/{stripe}/{embryo}/logs/validate_mcmc_combined.log"
    shell:
        """
        python scripts/07_validate_MCMC_combined.py \\
            --spatial-chain  {input.spatial_chain} \\
            --null-chain     {input.null_chain} \\
            --age-chain      {input.age_chain} \\
            --biphasic-chain {input.biphasic_chain} \\
            --transcription  {input.transcription} \\
            --mrna           {input.mrna} \\
            --n-ap-bins      {params.n_ap_bins} \\
            --n-dv-bins      {params.n_dv_bins} \\
            --output         {output.plot} \\
            2>&1 | tee {log}
        """

rule compare_models:
    """
    Compare all four degradation models using LOO-CV and posterior predictive checks.

    Runs scripts/08_compare_models_loo_ppc.py with all four model result directories,
    producing LOO comparison statistics, PPC plots, and a summary report.

    Requires all four model chains to be complete before running.

    Wildcards:
    - {stripe}: Which eve stripe to analyze
    - {embryo}: Which embryo to analyze (base name, without suffix)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        script="scripts/08_compare_models_loo_ppc.py",
        spatial_chain="results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv",
        null_chain="results_{max_time}/{stripe}/{embryo}_null/chains/degradation_chain.csv",
        age_chain="results_{max_time}/{stripe}/{embryo}_age/chains/degradation_chain.csv",
        biphasic_chain="results_{max_time}/{stripe}/{embryo}_biphasic/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv"
    output:
        summary="results_{max_time}/comparison/{stripe}/{embryo}/comparison_summary.txt"
    params:
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS,
        spatial_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}",
        null_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}_null",
        age_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}_age",
        biphasic_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}_biphasic",
        output_dir=lambda wc: f"results_{wc.max_time}/comparison/{wc.stripe}/{wc.embryo}"
    log:
        "results_{max_time}/comparison/{stripe}/{embryo}/logs/compare_models.log"
    shell:
        """
        python scripts/08_compare_models_loo_ppc.py \\
            --delayed-dir {params.age_dir} \\
            --biphasic-dir  {params.biphasic_dir} \\
            --null-dir      {params.null_dir} \\
            --spatial-dir   {params.spatial_dir} \\
            --transcription {input.transcription} \\
            --mrna          {input.mrna} \\
            --output-dir    {params.output_dir} \\
            --n-chains      {params.n_chains} \\
            --n-ap-bins     {params.n_ap_bins} \\
            --n-dv-bins     {params.n_dv_bins} \\
            2>&1 | tee {log}
        """


rule loo_diagnostics:
    """
    Run PSIS-LOO diagnostics for all four degradation models.

    Runs scripts/09_loo_diagnostics.py, producing per-observation Pareto k-hat plots,
    log-likelihood variance plots, posterior D plots, and a diagnostics summary.

    Wildcards:
    - {stripe}: Which eve stripe to analyze
    - {embryo}: Which embryo to analyze (base name, without suffix)
    """
    wildcard_constraints:
        embryo="[^/]*(?<!_null)(?<!_age)(?<!_biphasic)"
    input:
        script="scripts/09_loo_diagnostics.py",
        spatial_chain="results_{max_time}/{stripe}/{embryo}/chains/degradation_chain.csv",
        null_chain="results_{max_time}/{stripe}/{embryo}_null/chains/degradation_chain.csv",
        age_chain="results_{max_time}/{stripe}/{embryo}_age/chains/degradation_chain.csv",
        biphasic_chain="results_{max_time}/{stripe}/{embryo}_biphasic/chains/degradation_chain.csv",
        transcription="data/processed_transcription_data/transcription_traces_no_ids_{stripe}_{max_time}.csv",
        mrna="results_{max_time}/data/processed_mRNA_data_{stripe}/{embryo}_sass_formodel.csv"
    output:
        summary="results_{max_time}/comparison/{stripe}/{embryo}/loo_diagnostics_summary.txt"
    params:
        n_chains=N_MCMC_CHAINS,
        n_ap_bins=N_AP_BINS,
        n_dv_bins=N_DV_BINS,
        spatial_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}",
        null_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}_null",
        age_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}_age",
        biphasic_dir=lambda wc: f"results_{wc.max_time}/{wc.stripe}/{wc.embryo}_biphasic",
        output_dir=lambda wc: f"results_{wc.max_time}/comparison/{wc.stripe}/{wc.embryo}"
    log:
        "results_{max_time}/comparison/{stripe}/{embryo}/logs/loo_diagnostics.log"
    shell:
        """
        python scripts/09_loo_diagnostics.py \\
            --delayed-dir {params.age_dir} \\
            --biphasic-dir  {params.biphasic_dir} \\
            --null-dir      {params.null_dir} \\
            --spatial-dir   {params.spatial_dir} \\
            --transcription {input.transcription} \\
            --mrna          {input.mrna} \\
            --output-dir    {params.output_dir} \\
            --n-chains      {params.n_chains} \\
            --n-ap-bins     {params.n_ap_bins} \\
            --n-dv-bins     {params.n_dv_bins} \\
            2>&1 | tee {log}
        """


rule clean:
    """Remove all generated files."""
    shell:
        """
        rm -rf data/processed_transcription_data/transcription_traces_*.csv
        rm -rf results_*/stripe*/*/chains/*.csv
        rm -rf results_*/stripe*/*/figures/*.pdf
        rm -rf results_*/stripe*/*/summary_statistics*.csv
        rm -rf results_*/stripe*/*/logs/*.log
        rm -rf results_*/stripe*/validation/*.txt
        rm -rf results_*/figures/intermediate/transcription/transcription_heatmap_*.pdf
        rm -rf results_*/figures/intermediate/transcription/stripe_identification.pdf
        rm -rf results_*/figures/intermediate/mrna/*/*.pdf
        """