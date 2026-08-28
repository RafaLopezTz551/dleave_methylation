#!/bin/bash
#SBATCH --job-name=batch03
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00     
#SBATCH --output=%x_%j.log


export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
set -euo pipefail

export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
if mkdir -p "$SCRATCH_JOB/tmp" 2>/dev/null; then export TMPDIR="$SCRATCH_JOB/tmp"
else SCRATCH_JOB=""; export TMPDIR="${TMPDIR:-/tmp}"; echo "no scratch on $(hostname); TMPDIR=$TMPDIR"; fi
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/batch03
/opt/apps/r/4.4.1-studio/bin/Rscript code/batch03.R
cd /; if [ -n "$SCRATCH_JOB" ]; then rm -rf "$SCRATCH_JOB"; fi   
