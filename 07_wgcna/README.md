# 07_wgcna

Multi-tissue WGCNA coexpression network with a methylation overlay: module-level DMP and DMR enrichment (Fisher test, with a gene-length-stratified sensitivity analysis), GO/KEGG of the enriched modules, tail eigengene scores, and hub-gene analyses.

## Inputs

- HTSeq gene counts of all tissue libraries (the sample map is reconstructed from the file names).
- `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`.
- `02_landscape/data/genebody_methylation_per_gene.tsv`, `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `03_promoters/data/promoter_weber_classification.tsv`.
- `05_differential/data/dmps_annotated.tsv`, `05_differential/data/dmrs_annotated.tsv`.
- STRING v12 enrichment terms; JASPAR2024 SQLite; eggNOG-mapper annotations.

Paths are set in the ALL-CAPS constant block at the top of `code/07_wgcna.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Setup: seed, packages, WGCNA thread count, input paths, output directories, plot theme and figure savers.
2. Sample map reconstructed from the HTSeq count file names.
3. Count matrix, gene universe, low-count filter and VST.
4. Sample clustering and outlier removal (`fig7_s1_sample_clustering`, `fig7_s2_sample_clustering_no_outliers`).
5. Variance pre-filter sweep: the largest cutoff that keeps scale-free topology (`fig7_s2b_variance_cutoff`).
6. Soft-threshold power (`fig7_s3_soft_threshold`).
7. Module detection with `blockwiseModules`, cached as `objects/wgcna.rds` (`fig7_s4_module_dendrogram`).
8. Module DMP-burden enrichment: Fisher test per module (`fig7a_module_dmp_fisher`), plus a length-stratified CMH test kept as a sensitivity table.
9. CMH stratified by gene length and baseline methylation, and by movable CpGs (`figS7_module_dmp_enrichment_forest`).
10. Enrichment split by DMP direction, length-stratified CMH (`figS7_module_dmp_direction`).
11. Expression breadth (tau) across the atlas versus gene-body methylation (`figS7_expression_breadth_vs_genebody_meth`).
12. Raw Fisher versus length-adjusted CMH side by side (`figS7_length_adjustment_effect`).
13. Module DMR-burden enrichment, Fisher (`fig7d_module_dmr_fisher`).
14. HCP (CpG-island promoter) gene enrichment per module and selection of enriched modules (`fig7_s13_module_hcp_enrichment`).
15. GO/KEGG enrichment of the enriched modules (`fig7b_module_go_enrichment`).
16. Module eigengene scores, tail control versus amputated (`fig7c_eigengene_scores_tail`).
17. Module palette, module-trait heatmap over one-hot experiment groups, and module sizes (`fig7_s5_module_trait_heatmap`, `fig7_s6_module_sizes`).
18. Module membership (kME), intramodular connectivity and hub genes (`fig7_s7_hub_counts`, `fig7_s8_module_membership`).
19. Hub genes and DMPs/DMRs: Fisher plus length-tertile CMH (`fig7_s9_hub_dmp_dmr`).
20. Hub DMRs: methylation change versus expression change (`fig7_s10_hub_dmr_meth_vs_expression`).
21. Transcription-factor content of the hubs (`fig7_s12_hub_tf_content`).
22. GO enrichment for all modules, one dot-plot page per module (`fig7_s9_go_all_modules.pdf`).
23. Record the package versions (`sessionInfo_07_wgcna.txt`).

## Outputs

`data/`

`sample_map.tsv`, `variance_filter_sweep.tsv`, `module_assignments.tsv`,
`module_dmp_enrichment.tsv`, `module_dmp_enrichment_length_adjusted.tsv`,
`movable_cpgs_per_gene.tsv`, `module_dmp_enrichment_length_meth_adjusted.tsv`,
`module_dmp_enrichment_by_direction.tsv`, `expression_breadth_vs_genebody_methylation.tsv`,
`expression_breadth_per_gene.tsv`, `expression_breadth_deciles.tsv`, `length_adjustment_effect.tsv`,
`dmp_share_by_length_quintile_network.tsv`, `module_dmr_enrichment.tsv`,
`module_hcp_enrichment.tsv`, `module_go_enrichment.tsv`, `module_eigengene_scores_tail.tsv`,
`module_trait_cor.tsv`, `module_trait_fdr.tsv`, `module_sizes.tsv`, `module_membership_hubs.tsv`,
`hub_genes_dmp_dmr.tsv`, `hub_dmp_dmr_enrichment.tsv`, `hub_dmp_dmr_by_module.tsv`,
`hub_dmr_meth_vs_expression.tsv`, `module_tf_content.tsv`, `hub_transcription_factors.tsv`,
`module_go_all.tsv`

`objects/`

- `wgcna.rds`

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

- `fig7a_module_dmp_fisher`
- `fig7b_module_go_enrichment`
- `fig7c_eigengene_scores_tail`
- `fig7d_module_dmr_fisher`

`figures/supplementary/`

`fig7_s1_sample_clustering`, `fig7_s2_sample_clustering_no_outliers`, `fig7_s2b_variance_cutoff`,
`fig7_s3_soft_threshold`, `fig7_s4_module_dendrogram`, `figS7_module_dmp_enrichment_forest`,
`figS7_module_dmp_direction`, `figS7_expression_breadth_vs_genebody_meth`,
`figS7_length_adjustment_effect`, `fig7_s13_module_hcp_enrichment`, `fig7_s5_module_trait_heatmap`,
`fig7_s6_module_sizes`, `fig7_s7_hub_counts`, `fig7_s8_module_membership`, `fig7_s9_hub_dmp_dmr`,
`fig7_s10_hub_dmr_meth_vs_expression`, `fig7_s12_hub_tf_content`, `fig7_s9_go_all_modules.pdf`

The run also writes `sessionInfo_07_wgcna.txt` in the module folder.

## Run

From the pipeline root, after the modules it reads have finished:

```bash
sbatch 07_wgcna/code/07_wgcna.slurm
```

Resources requested by the launcher: 16 CPUs, 96 GB, 24 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 07_wgcna/code/07_wgcna.R`. The WGCNA thread count is read from `SLURM_CPUS_PER_TASK`.
