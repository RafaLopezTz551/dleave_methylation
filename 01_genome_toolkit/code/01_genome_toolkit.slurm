#!/bin/bash
#SBATCH --job-name=01_genome_toolkit
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=16:00:00
#SBATCH --exclude=node09,node06   # node09 /scratch is not BeeGFS
#SBATCH --output=logs/%x_%j.log

# Launcher for 01_genome_toolkit: loads the alignment/tree tools, sets a per-job scratch TMPDIR
# and runs code/01_genome_toolkit.R, which produces every object, table and figure of the module.

module load diamond/2.1.0 mafft/7.299 trimal/1.2 iqtree3/3.0.1

export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"

set -euo pipefail

export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"
export TMPDIR="$SCRATCH_JOB/tmp"
[ "$(stat -f -c %T /scratch)" = "fhgfs" ] || { echo "FATAL: /scratch is not BeeGFS on $(hostname)"; exit 1; }

cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/01_genome_toolkit
/opt/apps/r/4.4.1-studio/bin/Rscript code/01_genome_toolkit.R

cd /
rm -rf "$SCRATCH_JOB"
