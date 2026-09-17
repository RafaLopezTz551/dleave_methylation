#!/bin/bash
# SLURM launcher for 04_TEs: sets R_LIBS and a scratch TMPDIR (local /tmp fallback), then runs
# code/04_TEs.R with Rscript from the module directory. Submit from main/methylation_pipeline/:
#   sbatch 04_TEs/code/04_TEs.slurm
#SBATCH --job-name=04_TEs
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.log

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -uo pipefail   # no -e: a missing /scratch must fall back to /tmp, not kill the job
SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
if mkdir -p "$SCRATCH_JOB/tmp" 2>/dev/null; then
  export TMPDIR="$SCRATCH_JOB/tmp"
else
  SCRATCH_JOB="/tmp/${USER}_job_${SLURM_JOB_ID}"
  mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
  echo "[04_TEs] /scratch unavailable on $(hostname); local TMPDIR=$TMPDIR"
fi
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/04_TEs
/opt/apps/r/4.4.1-studio/bin/Rscript code/04_TEs.R; rc=$?
cd /; rm -rf "$SCRATCH_JOB"
exit $rc
