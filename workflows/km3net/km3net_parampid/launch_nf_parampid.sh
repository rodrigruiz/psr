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
nextflow -C ${CONFIG} run ${WORKFLOW} -profile woody_hannes_orca_data -resume