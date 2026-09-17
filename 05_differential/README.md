# 05_differential

Differential methylation after tail amputation: DSS DMP and DMR calling, gene assignment, region and GO/KEGG enrichment, a DMR-in-TE permutation test, and control analyses including a label-swap null.

## Inputs

- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`, `01_genome_toolkit/data/argonaute_census.tsv`.
- STRING v12 per-protein GO/KEGG terms.
- TE copy table with Kimura ages.
- OMArk conserved-unknown protein FASTA.

Paths are set in the ALL-CAPS constant block at the top of `code/05_differential.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Seed, packages, paths, palettes, plot theme and figure savers.
2. Load the CpG object and the annotation; derive gene symbols and biotypes.
3. DSS `DMLtest` per chromosome (cached under `objects/`), then `callDML` (DMPs) and `callDMR` (DMRs).
4. Region annotation and multi-gene assignment (2 kb promoter plus gene body).
5. Orphan features overlapping no promoter or gene body.
6. Gene-annotation summary: known versus unknown symbol, by biotype.
7. Feature counts by biotype and on conserved genes without annotation.
8. DMP volcano with a hex-density background of tested CpGs (`fig5a`).
9. DMR volcano, delta beta versus |areaStat| (`fig5b`).
10. Region-distribution pies for DMPs and DMRs (`fig5c`, `fig5d`).
11. Region enrichment against the genomic background, Fisher test (`fig5e`, `fig5e2`).
12. Per-gene DMP and DMR burden tables; top DMP-burden genes (`fig5f`).
13. GO/KEGG over-representation with STRING v12 terms, all and by direction; dot plots (`fig5h_dmp_go_dotplot`; supplementary `fig5h_dmp_hyper_go`, `fig5h_dmp_hypo_go`, `fig5i_dmr_go_dotplot`, `fig5i_dmr_hyper_go`, `fig5i_dmr_hypo_go`).
14. Target-size adjusted GO/KEGG with goseq (Wallenius; bias = CpGs per gene).
15. Three-way overlap of DMP, DMR and DE genes (`fig5j_dmp_dmr_de_venn`).
16. DMP/DMR genes by biotype (`fig5k_gene_coding_pie`).
17. Biotype selectivity: raw versus gene-length-adjusted odds ratio (table only).
18. Known versus unknown genes among DMP/DMR genes (`fig5l_gene_known_unknown_pie`).
19. DMRs in gene-body TEs: permutation enrichment test.
20. Controls for the TE enrichment, direction and function (`fig5m_dmr_te_gene_body`).
21. Label-swap null: the same DSS test on the two mislabelled 2-versus-2 splits (`fig5n_label_swap_null`).
22. Argonaute loci: DMP and DMR counts per locus (table only).
23. Record the package versions (`sessionInfo_05_differential.txt`).

## Outputs

`data/`

`dmps_annotated.tsv`, `dmrs_annotated.tsv`, `dmps_gene_assignments.tsv`,
`dmrs_gene_assignments.tsv`, `dmps_orphan_intergenic.tsv`, `dmrs_orphan_intergenic.tsv`,
`dmp_dmr_gene_annotation_summary.tsv`, `biotype_annotation_feature_counts.tsv`,
`conserved_unknown_dmp_dmr.tsv`, `dmp_region_enrichment.tsv`, `dmr_region_enrichment.tsv`,
`dmp_gene_burden.tsv`, `dmr_gene_burden.tsv`, `dmp_go_kegg_enrichment.tsv`,
`dmp_hyper_go_kegg_enrichment.tsv`, `dmp_hypo_go_kegg_enrichment.tsv`, `dmr_go_kegg_enrichment.tsv`,
`dmr_hyper_go_kegg_enrichment.tsv`, `dmr_hypo_go_kegg_enrichment.tsv`,
`dmp_share_by_length_quintile.tsv`, `goseq_dmp.tsv`, `goseq_dmp_hyper.tsv`, `goseq_dmp_hypo.tsv`,
`goseq_dmr.tsv`, `goseq_summary.tsv`, `dmp_dmr_de_intersection.tsv`,
`biotype_selectivity_length_adjusted.tsv`, `dmr_te_gene_body_enrichment.tsv`,
`dmr_te_gene_body_pairs.tsv`, `dmr_te_direction.tsv`, `dmr_te_go_vs_dmr_universe.tsv`,
`label_swap_counts.tsv`, `argonaute_methylation.tsv`

`objects/`

- `dmltest_chrmt_<chr>.rds`
- `dmltest_chrmt.rds`
- `dmltest_<swap>_chrmt_<chr>.rds`
- `dmltest_<swap>_chrmt.rds`

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

`fig5a_dmp_volcano`, `fig5b_dmr_volcano`, `fig5c_dmp_region_pie`, `fig5d_dmr_region_pie`,
`fig5e2_dmr_region_enrichment`, `fig5e_dmp_region_enrichment`, `fig5f_top_dmp_burden_genes`,
`fig5h_dmp_go_dotplot`

`figures/supplementary/`

`fig5h_dmp_hyper_go`, `fig5h_dmp_hypo_go`, `fig5i_dmr_go_dotplot`, `fig5i_dmr_hyper_go`,
`fig5i_dmr_hypo_go`, `fig5j_dmp_dmr_de_venn`, `fig5k_gene_coding_pie`,
`fig5l_gene_known_unknown_pie`, `fig5m_dmr_te_gene_body`, `fig5n_label_swap_null`

The run also writes `sessionInfo_05_differential.txt` in the module folder.

## Run

From the pipeline root, after the modules it reads have finished:

```bash
sbatch 05_differential/code/05_differential.slurm
```

Resources requested by the launcher: 2 CPUs, 128 GB, 48 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 05_differential/code/05_differential.R`. The DSS tests are cached per chromosome under `objects/`; a rerun after the cache exists is much shorter.
