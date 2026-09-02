#!/bin/bash
#SBATCH --job-name=08_motifs
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06  # local-xfs /scratch decoys, not BeeGFS (node06 found 2026-08-20)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log

# silently depended on the submitting user's ~/.bashrc R_LIBS export.
export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -euo pipefail
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/08_motifs
/opt/apps/r/4.4.1-studio/bin/Rscript code/08_motifs.R
cd /; rm -rf "$SCRATCH_JOB"
