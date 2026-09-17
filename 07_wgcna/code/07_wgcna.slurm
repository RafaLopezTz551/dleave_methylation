#!/bin/bash
# Launcher for 07_wgcna/code/07_wgcna.R: sets the R library path and a per-job scratch TMPDIR,
# runs the module script from the module directory, then removes the scratch space.
#SBATCH --job-name=07_wgcna
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.log

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -euo pipefail
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/07_wgcna
/opt/apps/r/4.4.1-studio/bin/Rscript code/07_wgcna.R
cd /; rm -rf "$SCRATCH_JOB"
