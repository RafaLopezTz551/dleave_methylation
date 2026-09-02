#!/usr/bin/env bash
#SBATCH --job-name=methbat_s2
#SBATCH --partition=defq
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=14
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --exclude=node09
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -uo pipefail

source /home/rlopezt/miniconda3/etc/profile.d/conda.sh
set +u                      # conda's own scripts reference unset vars
conda activate pbtools
set -u
methbat --version

HERE=/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/00_data_pacbio_hifi
BAM=$HERE/bam/sample2.aligned.bam

export SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
[ "$(stat -f -c %T /scratch)" = "fhgfs" ] || { echo "FATAL: /scratch is not BeeGFS here"; exit 1; }

cd "$SCRATCH_JOB"
echo "[methbat] start $(date)  bam=$BAM"
methbat pileup --threads 14 --input-bam "$BAM" --output-prefix sample2
echo "[methbat] exit=$? $(date)"

ls -la "$SCRATCH_JOB"

for f in sample2.5mC.bed.gz sample2.5mC.bed.gz.tbi \
         sample2.5hmC.bed.gz sample2.5hmC.bed.gz.tbi \
         sample2.6mA.bed.gz sample2.6mA.bed.gz.tbi \
         sample2.summary.json; do
  if [ -s "$SCRATCH_JOB/$f" ]; then
    cp "$SCRATCH_JOB/$f" "$HERE/$f.part" && mv "$HERE/$f.part" "$HERE/$f"
    echo "  kept $f"
  else
    echo "  MISSING $f"
  fi
done

cd /; rm -rf "$SCRATCH_JOB"
echo "[methbat_s2] done"
