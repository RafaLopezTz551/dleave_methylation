#!/bin/bash
#SBATCH -J star_align
#SBATCH -o ostar_%A_%a.out
#SBATCH -e ostar_%A_%a.err
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=8          # 8 threads *inside* STAR per sample
#SBATCH --mem=32G                  # tune for your genome & read length
# #SBATCH -p standard

set -euo pipefail
module load star/2.7.11a
fastq_dir=/mnt/data/alfredvar/wgutierrez/30_tentacles/6sortmerna

#tentacles=($(find $fastq_dir -maxdepth 1 -name "rmoverrep*fwd.fq" | xargs -I{} -n 1 basename -s "_srtd_fwd.fq"  {}))




tentacles=()
while IFS= read -r -d '' f; do
  base=${f##*/}                     # basename
  base=${base%_srtd_fwd.fq}         # strip suffix
  tentacles+=("$base")              # e.g. rmoverrep_unclass_C4
done < <(find "$fastq_dir" -maxdepth 1 -type f -name 'rmoverrep*_*fwd.fq' -print0)

N=${#tentacles[@]}
idx=${SLURM_ARRAY_TASK_ID:-0}

if (( idx >= N )); then
  echo "Array index $idx >= number of samples $N — exiting."
  exit 0
fi


base="${tentacles[$idx]}"
sample="${base#rmoverrep_unclass_}"

echo "[$(date)] Node: $(hostname)"
echo "[$(date)] Sample index $idx/$((N-1))  base=$base  sample=$sample"
echo "[$(date)] Using ${SLURM_CPUS_PER_TASK} threads"

BIND="--cpu-bind=cores"

R1="$fastq_dir/${base}_srtd_fwd.fq"
R2="$fastq_dir/${base}_srtd_rev.fq"


# Sanity checks
[[ -s "$R1" ]] || { echo "Missing $R1"; exit 2; }
[[ -s "$R2" ]] || { echo "Missing $R2"; exit 2; }


srun $BIND STAR --genomeDir $genoma/STAR_index --readFilesIn "$R1" "$R2" --outSAMstrandField intronMotif --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --runThreadN "${SLURM_CPUS_PER_TASK}" --outSAMtype BAM Unsorted --outFileNamePrefix "${sample}."

