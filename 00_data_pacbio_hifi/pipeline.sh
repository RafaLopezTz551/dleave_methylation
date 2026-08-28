#!/bin/bash
#SBATCH --job-name=dlaeve_meth
#SBATCH --partition=defq
#SBATCH --cpus-per-task=14
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=pipeline_%A_%a.log
#SBATCH --error=pipeline_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
source /home/rlopezt/miniconda3/etc/profile.d/conda.sh

REF=/scratch/groups/alfredvar/rlopezt/33-HiFi/ref.fasta
SAMPLES=(sample1 sample2)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# 1. jasmine: call 5mC + 5hmC + 6mA from kinetics 
conda activate jasmine
jasmine ${S}.raw.bam ${S}.allmods.bam \
  -j 14 --log-level INFO

# 2. pbmm2: align reads to the genome (keeps MM/ML tags) 
conda activate pbtools
pbmm2 align "$REF" ${S}.allmods.bam ${S}.aligned.bam \
  --preset HIFI --sort -j 12 -J 2

# 3. modkit CpG methylome: 5mC + 5hmC per CpG 
modkit pileup ${S}.aligned.bam ${S}.cpg.bed \
  --cpg --ref "$REF" \
  --filter-threshold C:0.8 --mod-thresholds m:0.8 --mod-thresholds h:0.8 \
  --threads 14 --log-filepath ${S}.cpg.log

# 4. modkit 6mA: per adenine, genome-wide 
modkit pileup ${S}.aligned.bam ${S}.6mA.bed \
  --motif A 0 --ref "$REF" \
  --filter-threshold A:0.9 --mod-thresholds a:0.9 \
  --threads 14 --log-filepath ${S}.6mA.log

echo "=== done $S ==="
