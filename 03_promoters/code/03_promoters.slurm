#!/bin/bash
#SBATCH --job-name=03_promoters
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06  # local-xfs /scratch decoys, not BeeGFS (node06 found 2026-08-20)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00     # backfillable (added 2026-08-17)
#SBATCH --output=%x_%j.log

# silently depended on the submitting user's ~/.bashrc R_LIBS export.
export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -euo pipefail
# 03_promoters does no heavy I/O (it reads RDS + classifies), so scratch is optional here.
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
if mkdir -p "$SCRATCH_JOB/tmp" 2>/dev/null; then export TMPDIR="$SCRATCH_JOB/tmp"
else SCRATCH_JOB=""; export TMPDIR="${TMPDIR:-/tmp}"; echo "no scratch on $(hostname); TMPDIR=$TMPDIR"; fi
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/03_promoters
/opt/apps/r/4.4.1-studio/bin/Rscript code/03_promoters.R
cd /; if [ -n "$SCRATCH_JOB" ]; then rm -rf "$SCRATCH_JOB"; fi   # `&& rm` would exit 1 (false SLURM FAILED) when scratch was skipped
