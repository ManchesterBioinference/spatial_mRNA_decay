#!/bin/bash --login

#SBATCH -p multicore #
#SBATCH -n 36 
#SBATCH -t 2:00:00

#####################################################
### run command
#####################################################
source ~/.bashrc
source .venv/bin/activate
#snakemake --snakefile Snakefile --cores 36 #-f --keep-incomplete --keep-going --rerun-incomplete
uv run snakemake --snakefile Snakefile_shift_sensitivity -c36 --keep-going --rerun-incomplete shiftLeft/results_1200/comparison/stripe2/e7um/comparison_summary.txt shiftLeft/results_1200/comparison/stripe2/e8_9um/comparison_summary.txt shiftLeft/results_1200/comparison/stripe2/e9_10um/comparison_summary.txt shiftRight/results_1200/comparison/stripe2/e7um/comparison_summary.txt shiftRight/results_1200/comparison/stripe2/e8_9um/comparison_summary.txt shiftRight/results_1200/comparison/stripe2/e9_10um/comparison_summary.txt