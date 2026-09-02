#!/bin/bash
#SBATCH --job-name=01_genome_toolkit
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G           # The heavy isoform (DRIMSeq/Salmon) and DEXSeq steps were
#SBATCH --time=16:00:00     # ~1 h original 01_genome_toolkit + ~3 h TF-annotation section (former
#SBATCH --exclude=node09,node06  # node09's /scratch is local xfs, not BeeGFS
#SBATCH --output=%x_%j.log

module load diamond/2.1.0 mafft/7.299 trimal/1.2 iqtree3/3.0.1

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"

set -euo pipefail

export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"
export TMPDIR="$SCRATCH_JOB/tmp"

cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/01_genome_toolkit
/opt/apps/r/4.4.1-studio/bin/Rscript code/01_genome_toolkit.R

cd /
rm -rf "$SCRATCH_JOB"
