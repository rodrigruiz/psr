#!/bin/bash
#SBATCH --partition=work
#SBATCH --ntasks-per-node=4  # run 4 tasks on a single node
#SBATCH --mem=160G
#SBATCH -c 8
#SBATCH -t 23:00:00


WORKFLOW=$1
CONFIG=$2

# Load Nextflow and Singularity environment
#module load nextflow
#module load singularity

# Run Nextflow workflow with Slurm configuration
nextflow run /home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_angular_resolution/km3net_angular_resolution_workflow.nf -c /home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_angular_resolution/profiles.config -profile woody_hannes_local