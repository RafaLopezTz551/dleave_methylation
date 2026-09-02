#!/bin/bash

#SBATCH --job-name=ocyt_report      # Job name
#SBATCH --output=slurm_logs/cytReport_%j.out  # Standard output log
#SBATCH --error=slurm_logs/cytReport_%j.err   # Standard error log
#SBATCH --nodes=1                   # Run on a single node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=1           # Request 8 CPUs for this task
#SBATCH --mem=40G                  # Request 100GB of memory

set -euo pipefail

# ---------- Usage ----------
usage() {
  cat <<'USAGE'
Usage:
  sbatch this_script.sh -i <input.cov.gz> -d <out_dir> -o <out_name>

Required:
  -i, --input   Path to input *.bismark.cov.gz file
  -d, --dir     Output directory (will be created if missing)
  -o, --out     Output name (used for -o in coverage2cytosine)

Optional:
  -g, --genome  Genome folder (default: $PWD/00_Genoma)
  -z, --gzip    Keep gzip output (default on)
  -h, --help    Show this help

Examples:
  sbatch this_script.sh \
    -i 09_methylation_calls/C2_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz \
    -d 09_methylation_calls \
    -o C2

  sbatch this_script.sh \
    --input sample.cov.gz --dir 09_methylation_calls --out Sample1 --genome /path/to/genome
USAGE
}

# ---------- Parse CLI args (GNU getopt) ----------
# Shell out early if getopt is missing
command -v getopt >/dev/null 2>&1 || { echo "ERROR: GNU getopt is required."; exit 2; }

GENOME_DIR="$(pwd)/00_Genoma"
GZIP_FLAG="--gzip"

ARGS=$(getopt -o i:d:o:g:zh \
  --long input:,dir:,out:,genome:,gzip,help \
  -n "$(basename "$0")" -- "$@") || { usage; exit 2; }
eval set -- "$ARGS"

INPUT_COV=""
OUT_DIR=""
OUT_NAME=""

while true; do
  case "$1" in
    -i|--input)   INPUT_COV="$2"; shift 2 ;;
    -d|--dir)     OUT_DIR="$2";   shift 2 ;;
    -o|--out)     OUT_NAME="$2";  shift 2 ;;
    -g|--genome)  GENOME_DIR="$2"; shift 2 ;;
    -z|--gzip)    GZIP_FLAG="--gzip"; shift ;;
    -h|--help)    usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Internal parsing error"; exit 2 ;;
  esac
done

# ---------- Validate ----------
[[ -z "${INPUT_COV}" ]] && { echo "ERROR: --input is required."; usage; exit 2; }
[[ -z "${OUT_DIR}"   ]] && { echo "ERROR: --dir is required."; usage; exit 2; }
[[ -z "${OUT_NAME}"  ]] && { echo "ERROR: --out is required."; usage; exit 2; }

if [[ ! -f "${INPUT_COV}" ]]; then
  echo "ERROR: Input file not found: ${INPUT_COV}"
  exit 1
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "ERROR: Genome folder not found: ${GENOME_DIR}"
  exit 1
fi

mkdir -p slurm_logs "${OUT_DIR}"

# ---------- Environment ----------
echo "Starting coverage2cytosine run..."
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-N/A}"
echo "Mem/Node: ${SLURM_MEM_PER_NODE:-N/A}"

# --- Load Required Modules ---
# Uncomment and update these lines if Bismark or Bowtie2 are
# loaded as modules on your HPC cluster.
#
module load bowtie2/2.5.4 samtools/1.22.1 
export PATH=/mnt/data/alfredvar/jmiranda/80-scripts/81-bin/Bismark-0.25.1:$PATH
export PERL5LIB=$HOME/.perl/lib:$PERL5LIB
echo "Genome folder: ${GENOME_DIR}"
echo "Input file   : ${INPUT_COV}"
echo "Output dir   : ${OUT_DIR}"
echo "Output name  : ${OUT_NAME}"


# --- Run the Command ---
# The --path_to_aligner you provided (/usr/bin/bowtie2/) is used directly.
# If Bowtie2 is loaded as a module, Bismark might find it automatically,
# and you might be able to remove --path_to_aligner.
set -x
coverage2cytosine \
  ${GZIP_FLAG} \
  --genome_folder "${GENOME_DIR}" \
  -dir "${OUT_DIR}" \
  -o "${OUT_NAME}" \
  "${INPUT_COV}"
set +

echo "coverage2cytosine finished with exit code $?."
