#!/bin/bash
# SLURM launcher for 03_promoters: sets R_LIBS and a scratch TMPDIR, then runs code/03_promoters.R
# with Rscript from the module directory. Submit from main/methylation_pipeline/:
#   sbatch 03_promoters/code/03_promoters.slurm
#SBATCH --job-name=03_promoters
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.log

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -euo pipefail
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
if mkdir -p "$SCRATCH_JOB/tmp" 2>/dev/null; then export TMPDIR="$SCRATCH_JOB/tmp"
else SCRATCH_JOB=""; export TMPDIR="${TMPDIR:-/tmp}"; echo "no scratch on $(hostname); TMPDIR=$TMPDIR"; fi
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/03_promoters
/opt/apps/r/4.4.1-studio/bin/Rscript code/03_promoters.R
cd /; if [ -n "$SCRATCH_JOB" ]; then rm -rf "$SCRATCH_JOB"; fi
