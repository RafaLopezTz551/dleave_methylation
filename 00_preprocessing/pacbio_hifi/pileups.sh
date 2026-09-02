#!/bin/bash
#SBATCH --job-name=dlaeve_pileup
#SBATCH --partition=defq
#SBATCH --cpus-per-task=14
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=pileup_%A_%a.log
#SBATCH --error=pileup_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
source /home/rlopezt/miniconda3/etc/profile.d/conda.sh
conda activate pbtools

REF=/scratch/groups/alfredvar/rlopezt/33-HiFi/ref.fasta
SAMPLES=(sample1 sample2)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

modkit pileup ${S}.aligned.bam ${S}.cpg.bed \
  --cpg --ref "$REF" \
  --modified-bases 5mC 5hmC \
  --filter-threshold C:0.8 --mod-thresholds m:0.8 --mod-thresholds h:0.8 \
  --threads 14 --log-filepath ${S}.cpg.log

modkit pileup ${S}.aligned.bam ${S}.6mA.bed \
  --motif A 0 --ref "$REF" \
  --modified-bases 6mA \
  --filter-threshold A:0.9 --mod-thresholds a:0.9 \
  --threads 14 --log-filepath ${S}.6mA.log

echo "=== done $S ==="
