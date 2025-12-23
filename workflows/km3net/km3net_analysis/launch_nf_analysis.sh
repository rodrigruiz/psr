#!/bin/bash
#SBATCH --partition=work
#SBATCH --ntasks-per-node=4  # run 4 tasks on a single node
#SBATCH --mem=128G
#SBATCH -c 8
#SBATCH -t 12:00:00


WORKFLOW=$1
CONFIG=$2
shift 2
# Load Nextflow and Singularity environment
#module load nextflow
#module load singularity

# Run Nextflow workflow with Slurm configuration
# nextflow plugin install nf-boosts
nextflow -C ${CONFIG} run ${WORKFLOW} -profile woody_hannes_container_analysis "$@" -with-report