#!/bin/bash
#SBATCH --job-name=batch05
#SBATCH --partition=defq
#SBATCH --exclude=node09,node06  # local-xfs /scratch decoys, not BeeGFS (node06 found 2026-08-20)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log

# Explicit R library path (audit 2026-08-17): without this, package resolution
# silently depended on the submitting user's ~/.bashrc R_LIBS export.
export R_LIBS="/mnt/data/alfredvar/rlopezt/rlib:/mnt/data/alfredvar/jmiranda/80-scripts/88-R-packages"
# MEMORY (dropped back to 32G on 2026-08-18): the chr+mt DSS recompute completed
#   on 2026-08-17 (job 258923), so objects/dmltest_chrmt*.rds all exist and reruns
#   are cache-load only (measured MaxRSS on the cached path was 13.7G, job 240026
#   — never lower below 32G). RESTORE 128G if the dmltest_chrmt* caches are ever
#   deleted: an uncached whole-genome DMLtest OOMs at 32G.
set -uo pipefail   # NOT -e: /scratch is node-local (absent on e.g. node09) — a
                   # scratch miss must fall back, not kill the job before Rscript.
# batch05 uses scratch only for TMPDIR (its heavy DSS step is cached in objects/),
# so a local fallback is safe and keeps the job node-agnostic.
SCRATCH_JOB="/scratch/groups/alfredvar/$USER/job_$SLURM_JOB_ID"
if mkdir -p "$SCRATCH_JOB/tmp" 2>/dev/null; then
  export TMPDIR="$SCRATCH_JOB/tmp"
else
  SCRATCH_JOB="/tmp/${USER}_job_${SLURM_JOB_ID}"
  mkdir -p "$SCRATCH_JOB/tmp"; export TMPDIR="$SCRATCH_JOB/tmp"
  echo "[batch05] /scratch unavailable on $(hostname); local TMPDIR=$TMPDIR"
fi
cd /mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/batch05
/opt/apps/r/4.4.1-studio/bin/Rscript code/batch05.R; rc=$?
cd /; rm -rf "$SCRATCH_JOB"
exit $rc
