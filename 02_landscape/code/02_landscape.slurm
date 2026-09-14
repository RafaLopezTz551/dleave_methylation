#!/bin/bash
#SBATCH --job-name=02_landscape
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00     # uncached bsseq rebuild + 2x HiFi bed parse; backfillable (2026-08-17)
#SBATCH --exclude=node09,node06  # node09's /scratch is local xfs, not BeeGFS
#SBATCH --output=%x_%j.log

# silently depended on the submitting user's ~/.bashrc R_LIBS export.
export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
# Pure-R module (builds the bsseq object from the 4 CpG reports and reads the HiFi
# pileups). SLURM only sets a scratch TMPDIR and launches the single R script.
set -euo pipefail
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/02_landscape
/opt/apps/r/4.4.1-studio/bin/Rscript code/02_landscape.R
cd /; rm -rf "$SCRATCH_JOB"
