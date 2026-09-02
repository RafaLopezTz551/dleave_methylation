#!/bin/bash

#SBATCH --job-name=bam2nuc_arr      # Job name
#SBATCH --output=slurm_logs/bam2nuc_%A_%a.out  # Standard output log
#SBATCH --error=slurm_logs/bam2nuc_%A_%a.err   # Standard error log
#SBATCH --array=1-4                   # Run 4 tasks, numbered 1, 2, 3, 4
#SBATCH --ntasks=1                    # Each task is a single task
#SBATCH --cpus-per-task=2           # Request 8 CPUs for this task
#SBATCH --mem=100G                  # Request 100GB of memory

# --- SAFETY FIRST ---
# Stop script on any error
set -e
# Treat unset variables as errors
set -u

echo "Starting Bismark alignment..."
echo "Starting job $SLURM_JOB_ID, array task $SLURM_ARRAY_TASK_ID"
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

SAMPLE_IDS=("A1" "A2" "C1" "C2")

# 2. Get the correct sample ID for this task
#    Task 1 -> Index 0 -> "A1"
#    Task 2 -> Index 1 -> "A2"
#    ...

CURRENT_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
CURRENT_ID=${SAMPLE_IDS[$CURRENT_INDEX]}

echo "Processing Sample ID: $CURRENT_ID"

# 3. Define your main directories
GENOME_DIR="00_Genoma"
INPUT_DIR="06_bismark_alignment"
MAIN_OUTPUT_DIR=" 07_methylation_extraction_and_reports" # Main output folder

# 4. Define the specific files for this task
INPUT_BAM="${INPUT_DIR}/trimmed_${CURRENT_ID}_R1.paired_bismark_bt2_pe.bam"
SAMPLE_OUTPUT_DIR="${MAIN_OUTPUT_DIR}/${CURRENT_ID}_output" # Create a unique output dir for this sample

# --- CREATE DIRECTORIES ---
# Create the main output and slurm log directories
mkdir -p $MAIN_OUTPUT_DIR
mkdir -p $SAMPLE_OUTPUT_DIR
mkdir -p slurm_logs

# --- RUN THE COMMAND ---
echo "Input file: $INPUT_BAM"
echo "Output directory: $SAMPLE_OUTPUT_DIR"

bam2nuc --genome_folder $GENOME_DIR --dir $SAMPLE_OUTPUT_DIR $INPUT_BAM

echo "Finished processing $CURRENT_ID."
