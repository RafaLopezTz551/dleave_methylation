# 00_preprocessing — from raw reads to the inputs the analysis scripts read

Nothing in this folder is run by the analysis scripts. It documents, with the exact job
scripts that were used on the Fénix cluster (UNAM/LAVIS, Slurm), how the raw sequencing
reads became the three inputs that `01_genome_toolkit` … `08_motifs` read by path:

| Input read by the pipeline | Produced by | Scripts |
|---|---|---|
| Bismark CpG reports (`*.CpG_report.txt.gz`, one per WGBS library C1, C2, A1, A2) and their splitting reports | Bismark 0.25.1 / Bowtie 2 2.5.4 | `wgbs_bismark/` |
| HTSeq gene counts (44 libraries of the tissue atlas, EviAnn gene models) | STAR + HTSeq (`--stranded=reverse`) | `rnaseq_star_htseq/` |
| PacBio HiFi per site 5mC pileups (two bodywall animals) | ccs 8.3.0, jasmine, pbmm2 26.2.99, MethBat 1.1.0 | `pacbio_hifi/` |

## WGBS (`wgbs_bismark/`)
Run in this order: `01_prepare_bismarck_genome.sh` (bisulfite genome index),
`02_alignment_bismakr.sh` (paired end Bismark, `--parallel 4`), `03_aln_coverage_report.sh`,
`04_deduplicate_bams.sh` (`deduplicate_bismark`), `05_extract_methylation.sh`
(`bismark_methylation_extractor --ignore_r2 2 --cytosine_report`; the per sample
`*_splitting_report.txt` files give the CHH non conversion floor used in `04_TEs`),
`05b_standalone_coverage2cytosine.sh`, and `10_dedup_bam_to_cram.slurm` (the deduplicated
alignments kept as CRAM against the full assembly; `09_read_patterns` decodes them with
that reference). The scripts are verbatim copies of the lab's run scripts; paths inside
them point at the lab's storage layout and must be edited to yours.

## Bulk RNA-seq (`rnaseq_star_htseq/`)
`align_starment.slurm` / `align_starment_amputated.slurm` (STAR against the EviAnn
genome index), `sort_bams.sh`, then `count_with_HTseq.slurm`
(`htseq-count -f bam -r name --stranded=reverse` on the coding + lncRNA GTF). All
*D. laeve* libraries are stranded (dUTP), hence `--stranded=reverse`.
Note for anyone re-running the alignment: the STAR command carries
`--outSAMstrandField intronMotif`, a flag intended for unstranded data; on stranded
libraries it removes reads spanning non canonical introns. The published counts were
produced with that flag, so they are reproduced exactly by these scripts; a re-alignment
without it would recover those reads and change the counts slightly.

## PacBio HiFi (`pacbio_hifi/`)
`copy_raw_hifi.slurm` stages the movies, `methbat_sample2.slurm` is the MethBat pileup
recipe actually used (per site 5mC/5hmC/6mA, mapq ≥ 1, coverage ≥ 4, 20 bp edge trim);
`pipeline.sh` / `pileups.sh` are the earlier modkit attempt kept as provenance (it was not
used), and `per_6ma_s1.R` is the provenance copy that arrived with the source data.
The honest 6mA / 5hmC treatment (binomial false positive null, fibertools re call) is in
the paper's Methods and in the analysis folder of the working project.
