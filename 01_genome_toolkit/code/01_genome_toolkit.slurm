#!/bin/bash
#SBATCH --job-name=01_genome_toolkit
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G           # The heavy isoform (DRIMSeq/Salmon) and DEXSeq steps were
                            #   removed; the heaviest thing left is the genome RDS + gene DE.
                            #   96G is now generous but kept to avoid node13 (~61G, too small
                            #   for a batch that loads the genome RDS).
#SBATCH --time=16:00:00     # ~1 h original 01_genome_toolkit + ~3 h TF-annotation section (former
#SBATCH --exclude=node09,node06  # node09's /scratch is local xfs, not BeeGFS
#SBATCH --output=%x_%j.log

# TF-annotation section tools (former tf_section): DIAMOND, MAFFT, trimAl, IQ-TREE.
# HMMER is project-local (tools/hmmer), driven from inside the R code.
module load diamond/2.1.0 mafft/7.299 trimal/1.2 iqtree3/3.0.1

# silently depended on the submitting user's ~/.bashrc R_LIBS export.
export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
# This SLURM script does NO analysis. It only (1) makes a per-job /scratch
# workspace and (2) launches the single batch code 01_genome_toolkit.R, which produces
# EVERY object/table/figure (incl. the TF bridge + guide).

set -euo pipefail

# per-job scratch workspace (cluster policy: heavy I/O off /mnt/data)
export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"
export TMPDIR="$SCRATCH_JOB/tmp"
[ "$(stat -f -c %T /scratch)" = "fhgfs" ] || { echo "FATAL: /scratch is not BeeGFS on $(hostname)"; exit 1; }

cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/01_genome_toolkit
/opt/apps/r/4.4.1-studio/bin/Rscript code/01_genome_toolkit.R

# remove scratch (BAMs/index were temporary; results are in objects/ data/ figures/)
cd /
rm -rf "$SCRATCH_JOB"
