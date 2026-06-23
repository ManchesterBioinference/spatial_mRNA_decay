#!/bin/bash --login

#SBATCH -p gpuL #multicore #
#SBATCH -G 2
#SBATCH -n 24 #--cpus-per-task=10  # Use this instead of -n for Ray Tune #SBATCH --ntasks=1          # Single task (Ray will handle parallelization)
#SBATCH -t 3:00:00

#####################################################
### run command
#####################################################
source ~/.bashrc
source .venv/bin/activate
#snakemake -c36 --keep-incomplete --keep-going --rerun-incomplete
uv run snakemake --snakefile Snakefile_shift_sensitivity -c24 --rerun-incomplete shiftLeft/results_1200/comparison/stripe2/e7um/comparison_summary.txt shiftLeft/results_1200/comparison/stripe2/e8_9um/comparison_summary.txt shiftLeft/results_1200/comparison/stripe2/e9_10um/comparison_summary.txt shiftRight/results_1200/comparison/stripe2/e7um/comparison_summary.txt shiftRight/results_1200/comparison/stripe2/e8_9um/comparison_summary.txt shiftRight/results_1200/comparison/stripe2/e9_10um/comparison_summary.txt