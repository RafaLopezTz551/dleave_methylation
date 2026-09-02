#!/bin/bash

#SBATCH --job-name=oextract_methylation      # Job name
#SBATCH --output=slurm_logs/extract_met_%j.out  # Standard output log
#SBATCH --error=slurm_logs/extract_met_%j.err   # Standard error log
#SBATCH --nodes=1                   # Run on a single node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=24           # Request 8 CPUs for this task
#SBATCH --mem=200G                  # Request 100GB of memory


# --- Explanation ---
# We request --cpus-per-task=8 because:
# 1. Bismark runs 2 indexing processes in parallel (top/bottom strand).
# 2. You specified --parallel 4, which gives 4 threads *to each* of those processes.
# 3. Total CPUs = 2 processes * 4 threads/process = 8 CPUs.
#
# We request --mem=100G because indexing the human genome (GRCh38) is
# very memory-intensive. 100GB is a safe starting point.
# You may need to adjust this value based on your cluster's limits or
# if the job fails due to an Out-Of-Memory (OOM) error.
# ---------------------

set -e

set -u

echo "Starting Bismark alignment..."
echo "Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"
echo "Allocated Memory: $SLURM_MEM_PER_NODE"

# --- Load Required Modules ---
# Uncomment and update these lines if Bismark or Bowtie2 are
# loaded as modules on your HPC cluster.
#
module load bowtie2/2.5.4 samtools/1.22.1 
export PATH=/mnt/data/alfredvar/jmiranda/80-scripts/81-bin/Bismark-0.25.1:$PATH
export PERL5LIB=$HOME/.perl/lib:$PERL5LIB

deduplicated_bams=$(find 08_deduplication/ -name "*.bam")
GENOME_DIR="00_Genoma"


# --- Run the Command ---
# The --path_to_aligner you provided (/usr/bin/bowtie2/) is used directly.
# If Bowtie2 is loaded as a module, Bismark might find it automatically,
# and you might be able to remove --path_to_aligner.
bismark_methylation_extractor --scaffolds --parallel 8 --ignore_r2 2 --gzip --buffer_size 20G --bedGraph --output_dir 09_methylation_calls --cytosine_report --genome_folder $GENOME_DIR $deduplicated_bams
echo "Bismark alignment finished with exit code $?."
