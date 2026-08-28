#!/bin/bash
#SBATCH --job-name=batch04
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06  # local-xfs /scratch decoys, not BeeGFS (node06 found 2026-08-20)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=12:00:00     # backfillable (added 2026-08-17)
#SBATCH --output=%x_%j.log

# Explicit R library path (audit 2026-08-17): without this, package resolution
# silently depended on the submitting user's ~/.bashrc R_LIBS export.
export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -uo pipefail   # NOT -e: /scratch is node-local (absent/unwritable on e.g. node09)
                   # — a scratch miss must fall back to local /tmp, not kill the job.
# batch04 uses scratch only for TMPDIR (it reads cached objects + the TE table),
# so a local fallback is safe and keeps the job node-agnostic.
SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
if mkdir -p "$SCRATCH_JOB/tmp" 2>/dev/null; then
  export TMPDIR="$SCRATCH_JOB/tmp"
else
  SCRATCH_JOB="/tmp/${USER}_job_${SLURM_JOB_ID}"
  mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
  echo "[batch04] /scratch unavailable on $(hostname); local TMPDIR=$TMPDIR"
fi
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/batch04
/opt/apps/r/4.4.1-studio/bin/Rscript code/batch04.R; rc=$?
cd /; rm -rf "$SCRATCH_JOB"
exit $rc
