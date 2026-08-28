#!/bin/bash
#SBATCH --job-name=batch01
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G           
#SBATCH --time=16:00:00     
#SBATCH --output=%x_%j.log


module load diamond/2.1.0 mafft/7.299 trimal/1.2 iqtree3/3.0.1


export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"

set -euo pipefail

export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"
export TMPDIR="$SCRATCH_JOB/tmp"

cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/batch01
/opt/apps/r/4.4.1-studio/bin/Rscript code/batch01.R

cd /
rm -rf "$SCRATCH_JOB"
