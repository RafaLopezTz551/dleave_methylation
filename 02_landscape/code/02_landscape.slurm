#!/bin/bash
#SBATCH --job-name=02_landscape
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --exclude=node09,node06   # node09 /scratch is not BeeGFS
#SBATCH --output=logs/%x_%j.log

# Launcher for 02_landscape: sets a per-job scratch TMPDIR and runs code/02_landscape.R,
# which produces every object, table and figure of the module.

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -euo pipefail
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/02_landscape
/opt/apps/r/4.4.1-studio/bin/Rscript code/02_landscape.R
cd /; rm -rf "$SCRATCH_JOB"
