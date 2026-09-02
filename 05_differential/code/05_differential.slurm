#!/bin/bash
#SBATCH --job-name=05_differential
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06  # local-xfs /scratch decoys, not BeeGFS (node06 found 2026-08-20)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G   # 2026-09-01: the [2b] label-swap DMLtests are uncached on first run (~40G peak each)
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -uo pipefail   # NOT -e: /scratch is node-local (absent on e.g. node09) — a
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
