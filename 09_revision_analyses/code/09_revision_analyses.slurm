#!/bin/bash
#SBATCH --job-name=09_revision_analyses
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.log

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -euo pipefail
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/09_revision_analyses
/opt/apps/r/4.4.1-studio/bin/Rscript code/09_revision_analyses.R
cd /; rm -rf "$SCRATCH_JOB"
