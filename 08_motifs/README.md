# 08_motifs

Transcription factor motif enrichment in the *D. laeve* methylome: MethylSeekR UMR/LMR segmentation, monaLisa binned enrichment along the LMR methylation-change gradient and across Weber promoter classes (slug and human), and HOMER known-motif enrichment of UMR promoters, LMRs and DMRs.

## Inputs

- `01_genome_toolkit/objects/genome_chrmt.rds`, `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/jaspar_ortholog_bridge.tsv`.
- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `03_promoters/dataset/` (GRCh38 reference cache).
- `05_differential/data/dmrs_annotated.tsv`.
- JASPAR2024 CORE SQLite; TE annotation table (`collapsed_te_age_data.tsv`).
- HOMER 5.1 (`findMotifsGenome.pl`) and `rsvg-convert`, called through `system2()`.

Paths are set in the ALL-CAPS constant block at the top of `code/08_motifs.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Setup: seed, minimal packages, paths, cache freshness rule, coverage filter.
2. MethylSeekR UMR/LMR segmentation of the pooled tail methylome (cached).
3. Takai-Jones CGI annotation and the MethylSeekR `calculateFDRs` grid.
4. Attach the motif-enrichment stack and set the parallel backend.
5. JASPAR2024 animal PWMs, sequence-orthology filter, HOMER motif library.
6. Sequence helpers.
7. Per-LMR methylation change (amputated minus control) on the coverage-filtered CpG set.
8. Bin LMRs by methylation change; monaLisa binned motif enrichment.
9. Figure helpers (saver, heatmap wrappers, motif bookkeeping).
10. LMR motif heatmap and LMRs per bin (`fig8_lmr_motif_enrichment`, `figS8_lmr_bin_density`).
11. LMR bin QC: whether the extreme bins are the noisiest LMRs (`figS8_lmr_bin_qc`).
12. Supplementary LMR figures: per-bin composition, full heatmap, segment overview (`figS8_lmr_bindiag_GC`, `figS8_lmr_bindiag_dinuc`, `figS8_lmr_full_enrichment`, `figS8_lmr_overview`).
13. LMR and UMR genomic location, bar and pie (`figS8_lmr_genomic_location`, `figS8_lmr_umr_region_pie`).
14. DMR overlap with LMRs and UMRs against matched random backgrounds.
15. Weber promoter classification, protein-coding genes only.
16. Weber-class promoter motif enrichment with categorical monaLisa bins (`figS8_promoter_weber_motif_enrichment`, `figS8_promoter_weber_gc`).
17. LMR and UMR overlap with Weber promoter classes (`figS8_weber_class_lmr_umr`).
18. HOMER known-motif enrichment for three region sets, UMR promoters, LMRs and DMRs: BED export, `findMotifsGenome.pl`, `knownResults.txt` to TSV, motif selection by effect size and logo figures (`fig8_umr_promoter_homer_motifs`, `figS8_lmr_homer_motifs`, `figS8_dmr_homer_motifs`).
19. Human GRCh38 Weber-class promoter run and cross-species comparison (`figS8_human_promoter_weber_motif_enrichment`, `figS8_hcp_motif_human_vs_dlaeve`).
20. Log notes for analyses handled in other modules; record the package versions (`sessionInfo_08_motifs.txt`).

## Outputs

`data/`

`methylseekr_segments.tsv`, `cgi_takai_jones.bed`, `cgi_takai_jones.tsv`, `cgi_summary.tsv`,
`methylseekr_fdr_table.tsv`, `jaspar_ortholog_bridge.tsv`, `lmr_methylation_change.tsv`,
`lmr_motif_enrichment.tsv`, `lmr_bin_qc.tsv`, `lmr_delta_vs_composition.tsv`,
`lmr_bin_qc_tests.tsv`, `lmr_umr_genomic_location.tsv`, `dmr_lmr_umr_overlap.tsv`,
`promoter_weber_classification_pc.tsv`, `promoter_weber_motif_enrichment.tsv`,
`promoter_weber_motif_summary.tsv`, `weber_class_lmr_umr_overlap.tsv`,
`umr_promoter_homer_known_motifs.tsv`, `lmr_homer_known_motifs.tsv`, `dmr_homer_known_motifs.tsv`,
`human_promoter_weber_motif_enrichment.tsv`, `hcp_motif_enrichment_human_vs_dlaeve.tsv`,
`hcp_motif_human_vs_dlaeve_correlation.tsv`

`objects/`

`qc_alpha_distribution_pmd_check.pdf`, `methylseekr_segments_cov10_chrmt_gr.rds`,
`methylseekr_segmentation_qc.pdf`, `cgi_takai_jones_chrmt.rds`, `qc_methylseekr_fdr_grid.pdf`,
`jaspar_ortholog_homer.motif`, `lmr_deltameth_cov10_chrmt_se.rds`, `promoter_weber_pc_chrmt_se.rds`,
`genome_chrmt.fa`, `human_promoter_weber_hs9606_se.rds`, `homer_<set>_chrmt.bed` (HOMER input
regions for UMR promoters, LMRs and DMRs), `homer_<set>_chrmt/` (HOMER runs with knownResults.txt),
`homer_<set>_chrmt.log` (HOMER logs), `homer_preparsed/`

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

- `fig8_lmr_motif_enrichment`
- `fig8_umr_promoter_homer_motifs`

`figures/supplementary/`

`figS8_lmr_bin_density`, `figS8_lmr_bin_qc`, `figS8_lmr_bindiag_GC`, `figS8_lmr_bindiag_dinuc`,
`figS8_lmr_full_enrichment`, `figS8_lmr_overview`, `figS8_lmr_genomic_location`,
`figS8_lmr_umr_region_pie`, `figS8_promoter_weber_motif_enrichment`, `figS8_promoter_weber_gc`,
`figS8_weber_class_lmr_umr`, `figS8_lmr_homer_motifs`, `figS8_dmr_homer_motifs`,
`figS8_human_promoter_weber_motif_enrichment`, `figS8_hcp_motif_human_vs_dlaeve`

The run also writes `sessionInfo_08_motifs.txt` in the module folder.

## Run

From the pipeline root, after the modules it reads have finished:

```bash
sbatch 08_motifs/code/08_motifs.slurm
```

Resources requested by the launcher: 16 CPUs, 128 GB, 24 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 08_motifs/code/08_motifs.R`. The parallel backend size is read from `SLURM_CPUS_PER_TASK`. HOMER and `rsvg-convert` must be on the paths set at the top of the script.
