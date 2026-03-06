#!/bin/bash --login

#SBATCH -p multicore #
#SBATCH -n 36 #--cpus-per-task=10  # Use this instead of -n for Ray Tune #SBATCH --ntasks=1          # Single task (Ray will handle parallelization)
#SBATCH -t 5:00:00

#####################################################
### run command
#####################################################
source ~/.bashrc
source .venv/bin/activate
snakemake -c36 --keep-incomplete --keep-going --rerun-incomplete