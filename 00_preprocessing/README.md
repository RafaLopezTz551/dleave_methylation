# 00_preprocessing: From raw reads to the inputs the analysis scripts read

Nothing in this folder is run by the analysis scripts. It documents, with the exact job
scripts that were used on the Fénix cluster (UNAM/LAVIS, Slurm), how the raw sequencing
reads became the three inputs that `01_genome_toolkit` … `08_motifs` read by path:

| Input read by the pipeline | Produced by | Scripts |
|---|---|---|
| Bismark CpG reports (`*.CpG_report.txt.gz`, one per WGBS library C1, C2, A1, A2) and their splitting reports | Bismark 0.25.1 / Bowtie 2 2.5.4 | `wgbs_bismark/` |
| HTSeq gene counts (44 libraries of the tissue atlas, EviAnn gene models) | STAR + HTSeq (`--stranded=reverse`) | `rnaseq_star_htseq/` |
| PacBio HiFi per site 5mC pileups (two bodywall animals) | ccs 8.3.0, jasmine, pbmm2 26.2.99, MethBat 1.1.0 | `pacbio_hifi/` |

## WGBS (`wgbs_bismark/`)
Run in this order: `01_prepare_bismark_genome.sh` (bisulfite genome index),
`02_alignment_bismark.sh` (paired end Bismark, `--parallel 4`), `03_aln_coverage_report.sh`,
`04_deduplicate_bams.sh` (`deduplicate_bismark`), `05_extract_methylation.sh`
(`bismark_methylation_extractor --ignore_r2 2 --cytosine_report`; the per sample
`*_splitting_report.txt` files give the CHH non conversion floor used in `04_TEs`),
`05b_standalone_coverage2cytosine.sh`, and `10_dedup_bam_to_cram.slurm` (the deduplicated
alignments kept as CRAM against the full assembly). The scripts are verbatim copies of the lab's run scripts; paths inside
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
`copy_raw_hifi.slurm` stages the two bodywall movies (BW1, BW2; one of them is the dataset of the genome assembly,
Miranda-Rodríguez et al. 2025) and `methbat_sample2.slurm` is the MethBat pileup recipe used for the per-site 5mC calls
(mapq ≥ 1, coverage ≥ 4, 20 bp edge trim). Circular consensus with kinetics (ccs 8.3.0), jasmine and pbmm2 are run
upstream as described in the manuscript Methods.
