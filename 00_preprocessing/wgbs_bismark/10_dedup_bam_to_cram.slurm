#!/usr/bin/env bash
#SBATCH --job-name=cram08
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --exclude=node09
#SBATCH --array=0-3
#SBATCH --output=slurm_logs/cram08.%A_%a.out

# 10_dedup_bam_to_cram.slurm — convert the four Bismark deduplicated BAMs in
# 08_deduplication/ to coordinate-sorted CRAM against the canonical genome,
# verify losslessness (quickcheck + identical flagstat, which fully decodes
# every record), then delete the source BAM. 2026-08-20.
#
# NOTE for future re-use: Bismark's methylation extractor needs pair-adjacent
# order, NOT coordinate order. To re-extract from these CRAMs first run:
#   samtools sort -n -T tmp -O bam -o <s>.nsorted.bam --reference $REF <s>.cram

module load samtools/1.22.1
set -euo pipefail

SAMPLES=(A1 A2 C1 C2)
S="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
BASE=/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication
REF=/mnt/data/alfredvar/30-Genoma/Deroceras_laeve_genome_GCA_051403575.fasta
BAM="$BASE/${S}_paired_bismark_bt2_pe.deduplicated.bam"
CRAM="$BASE/${S}_paired_bismark_bt2_pe.deduplicated.cram"

# resume mechanism: task already finished on a previous run
if [[ -s "$CRAM" && -s "$CRAM.crai" && ! -e "$BAM" ]]; then
    echo "$S: CRAM exists and BAM already removed — nothing to do"; exit 0
fi
[[ -s "$BAM" ]] || { echo "$S: missing input $BAM"; exit 1; }

[[ "$(stat -f -c %T /scratch)" == fhgfs ]] || { echo "no BeeGFS on $(hostname)"; exit 1; }
SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_JOB/tmp"
export TMPDIR="$SCRATCH_JOB/tmp"
trap 'cd /; rm -rf "$SCRATCH_JOB"' EXIT

echo "$S: staging to $SCRATCH_JOB"
cp "$REF" "$REF.fai" "$SCRATCH_JOB/"
REF_SC="$SCRATCH_JOB/$(basename "$REF")"
BAM_SC="$SCRATCH_JOB/$(basename "$BAM")"
cp "$BAM" "$BAM_SC"
samtools quickcheck "$BAM_SC"

echo "$S: sorting into CRAM"
CRAM_SC="$SCRATCH_JOB/$(basename "$CRAM")"
samtools sort -@ "$SLURM_CPUS_PER_TASK" -m 2G -T "$TMPDIR/sort_$S" \
    --reference "$REF_SC" -O cram -o "$CRAM_SC" "$BAM_SC"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$CRAM_SC"
samtools quickcheck "$CRAM_SC"

echo "$S: verifying (flagstat fully decodes both files)"
samtools flagstat -@ "$SLURM_CPUS_PER_TASK" "$BAM_SC" > "$SCRATCH_JOB/$S.bam.flagstat.txt"
samtools flagstat -@ "$SLURM_CPUS_PER_TASK" \
    --input-fmt-option reference="$REF_SC" "$CRAM_SC" > "$SCRATCH_JOB/$S.cram.flagstat.txt"
diff "$SCRATCH_JOB/$S.bam.flagstat.txt" "$SCRATCH_JOB/$S.cram.flagstat.txt"
echo "$S: flagstat identical"
cat "$SCRATCH_JOB/$S.cram.flagstat.txt"

echo "$S: copying CRAM back to NFS (atomic)"
cp "$CRAM_SC" "$CRAM.part"
MD5_SC=$(md5sum "$CRAM_SC" | cut -d' ' -f1)
MD5_NFS=$(md5sum "$CRAM.part" | cut -d' ' -f1)
[[ "$MD5_SC" == "$MD5_NFS" ]] || { echo "$S: md5 mismatch after copy-back"; exit 1; }
mv "$CRAM.part" "$CRAM"
cp "$CRAM_SC.crai" "$CRAM.crai.part" && mv "$CRAM.crai.part" "$CRAM.crai"
cp "$SCRATCH_JOB/$S.bam.flagstat.txt" "$SCRATCH_JOB/$S.cram.flagstat.txt" "$BASE/"

echo "$S: all checks passed (md5 $MD5_NFS) — deleting source BAM"
rm -f "$BAM"
ls -lh "$CRAM" "$CRAM.crai"
echo "$S: DONE"
