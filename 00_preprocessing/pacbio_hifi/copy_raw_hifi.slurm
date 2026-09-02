#!/bin/bash
#SBATCH --job-name=an_hificopy
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --exclude=node09
#SBATCH --output=%x_%j.log

# Import the RAW HiFi movie BAMs (kinetics-bearing: fi/ri/fp/rp tags verified,

set -uo pipefail
SRC=/scratch/groups/alfredvar/jmiranda/33-HiFi
DST=/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/00_data_pacbio_hifi/raw
REF=/mnt/data/alfredvar/30-Genoma/33-HiFi/md5sums.txt
mkdir -p "$DST"
date
for f in m64467e_230308_175634.5mC.hifi_reads.bam \
         m64467e_230308_175634.5mC.hifi_reads.bam.pbi \
         m64467e_230310_031320.5mC.hifi_reads.bam \
         m64467e_230310_031320.5mC.hifi_reads.bam.pbi; do
  if [ -f "$DST/$f" ]; then echo "[copy] $f already present, skipping"; continue; fi
  echo "[copy] $f"
  cp "$SRC/$f" "$DST/$f.part" && mv "$DST/$f.part" "$DST/$f"
done
cp "$REF" "$DST/md5sums.txt"
cd "$DST"
status=0
md5sum -c md5sums.txt || status=$?
if [ "$status" -eq 0 ]; then
  echo "[an_hificopy] VERDICT: NFS copies VERIFIED against md5sums.txt"
else
  echo "[an_hificopy] VERDICT: NFS copy FAILED verification (md5sum exit $status) - do NOT use raw/ until fixed"
fi
date
echo "[an_hificopy] done"
