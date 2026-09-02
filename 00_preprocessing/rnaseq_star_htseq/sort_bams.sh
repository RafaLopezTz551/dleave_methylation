#!/bin/bash
#SBATCH -J sortbam
#SBATCH -o osbam_%A_%a.out
#SBATCH -e osbam_%A_%a.err
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=8          # 8 threads *inside* STAR per sample
#SBATCH --mem=32G                  # tune for your genome & read length
# #SBATCH -p standard

set -euo pipefail
module load samtools/1.22.1
fastq_dir=$(pwd)

#tentacles=($(find $fastq_dir -maxdepth 1 -name "rmoverrep*fwd.fq" | xargs -I{} -n 1 basename -s "_srtd_fwd.fq"  {}))




tentacles=()
while IFS= read -r -d '' f; do
  base=${f##*/}                     # basename
  base=${base%.Aligned.out.bam}         # strip suffix
  tentacles+=("$base")              # e.g. rmoverrep_unclass_C4
done < <(find "$fastq_dir" -maxdepth 1 -type f -name '*.Aligned.out.bam' -print0)

N=${#tentacles[@]}
idx=${SLURM_ARRAY_TASK_ID:-0}

if (( idx >= N )); then
  echo "Array index $idx >= number of samples $N — exiting."
  exit 0
fi


base="${tentacles[$idx]}"

echo "[$(date)] Node: $(hostname)"
echo "[$(date)] Sample index $idx/$((N-1))  base=$base"
echo "[$(date)] Using ${SLURM_CPUS_PER_TASK} threads"

BIND="--cpu-bind=cores"

bamfile="$fastq_dir/${base}.Aligned.out.bam"
# Sanity checks
[[ -s "$bamfile" ]] || { echo "Missing $bamfile"; exit 2; }
srun $BIND samtools sort -o "bams/${base}.sorted.bam" --output-fmt BAM --threads 8 $bamfile


