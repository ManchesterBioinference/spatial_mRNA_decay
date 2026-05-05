#!/bin/bash --login

#SBATCH -p multicore #
#SBATCH -n 36 
#SBATCH -t 2:00:00

#####################################################
### run command
#####################################################
source ~/.bashrc
source .venv/bin/activate
snakemake --snakefile Snakefile --cores 36 #-f --keep-incomplete --keep-going --rerun-incomplete