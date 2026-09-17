#!/bin/bash
#SBATCH --job-name=05_differential
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.log

# Launcher for 05_differential/code/05_differential.R (DSS differential methylation).
# Submit from main/methylation_pipeline/: sbatch 05_differential/code/05_differential.slurm

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -uo pipefail
SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
if mkdir -p "$SCRATCH_JOB/tmp" 2>/dev/null; then
  export TMPDIR="$SCRATCH_JOB/tmp"
else
  SCRATCH_JOB="/tmp/${USER}_job_${SLURM_JOB_ID}"
  mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
  echo "[05_differential] /scratch unavailable on $(hostname); local TMPDIR=$TMPDIR"
fi
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/05_differential
/opt/apps/r/4.4.1-studio/bin/Rscript code/05_differential.R; rc=$?
cd /; rm -rf "$SCRATCH_JOB"
exit $rc
