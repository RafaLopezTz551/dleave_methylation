# 02_landscape

The baseline CpG methylation landscape of the *D. laeve* tail: per sample, per chromosome and per genomic region, its relation to gene expression, and the agreement between PacBio HiFi and WGBS methylation calls.

## Inputs

- Four Bismark CpG reports (C1, C2 control; A1, A2 amputated).
- EviAnn GFF (`01_genome_toolkit/objects/gff_chrmt.rds` when present, otherwise imported from the GFF).
- HTSeq gene counts, tail and bodywall libraries.
- PacBio HiFi per-site 5mC pileups of bodywall (two animals; see `00_preprocessing/pacbio_hifi`).
- `01_genome_toolkit/data/gene_de_tail.tsv`.

Paths are set in the ALL-CAPS constant block at the top of `code/02_landscape.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Seed, packages, paths, chromosome set, palettes, theme and figure saver.
2. Strand-collapsed bsseq object from the four CpG reports, coverage of at least 5 in all samples (cached as `objects/bsseq_cov5_chrmt.rds`).
3. GFF gene models (module 01 object when present, otherwise imported).
4. Tail expression deciles (HTSeq counts to DESeq2 VST) and the per-CpG pooled table.
5. Per-sample methylation in 1 Mb windows along each chromosome (`figS_per_chromosome_methylation`).
6. Global CpG methylation per sample (`fig2a`).
7. Pooled genome-wide methylation in 1 Mb windows (`fig2b`).
8. Region sets: promoter (2 kb upstream of the TSS), exon, intron, genic.
9. Mean methylation per region, regions may overlap (`fig2e`).
10. Share of methylated reads per region as an exclusive partition (`fig2j`).
11. Gene-body methylation per gene by tail expression decile.
12. Methylation versus expression statistics and the top-decile split.
13. Per-gene gene-body methylation by expression decile (`fig2f`).
14. Gene-body metagene by decile (`figS_genebody_metagene_decile`).
15. TSS +/- 5 kb metagene by decile for tail control, tail blastema and bodywall: HiFi bodywall CpG object (cached as `objects/hifi_bodywall_cpg_persample_chrmt.rds`), bodywall expression deciles, the three methylomes on one grid, and the blastema minus control delta (`fig2d`, `figS_tss5kb_delta_metagene`).
16. Mitochondrial and global methylation per sample for the three groups (`figS_mt_methylation_3groups`, `figS_global_methylation_3groups`).
17. HiFi versus WGBS agreement on shared CpGs and on 1 Mb windows.
18. Discrete-region metagene on WGBS tail: longest mRNA per gene, first/internal/last exons and introns, promoter, distal upstream and downstream segments, ten bins per segment (`figS_region_decile_metagene_wgbs_tail`).
19. The same discrete-region metagene on HiFi bodywall methylation (`fig2g`).
20. Sample PCA and correlation on 1 Mb-window beta (`fig2i`, supplementary).
21. Record the package versions (`sessionInfo_02_landscape.txt`).

## Outputs

`data/`

`per_sample_methylation.tsv`, `genomewide_methylation_1mb.tsv`, `region_methylation.tsv`,
`region_methylation_signal.tsv`, `genebody_methylation_per_gene.tsv`,
`genebody_expression_correlation.tsv`, `genebody_decile_summary.tsv`, `decile10_genes.tsv`,
`decile10_unmethylated_vs_methylated.tsv`, `metagene_genebody_decile.tsv`,
`metagene_tss5kb_decile.tsv`, `metagene_tss5kb_delta.tsv`, `mt_methylation_3groups.tsv`,
`global_methylation_3groups.tsv`, `hifi_wgbs_agreement.tsv`, `metagene_region_decile.tsv`,
`metagene_region_decile_bodywall.tsv`, `sample_pca_coords.tsv`

`objects/`

- `bsseq_cov5_chrmt.rds`
- `hifi_bodywall_cpg_persample_chrmt.rds`

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

`fig2a_global_methylation_per_sample`, `fig2b_genomewide_1mb`, `fig2d_tss5kb_metagene_decile`,
`fig2e_region_methylation`, `fig2f_genebody_methylation_decile`,
`fig2g_region_decile_metagene_bodywall`, `fig2j_region_methylation_pie`

`figures/supplementary/`

`figS_per_chromosome_methylation`, `figS_genebody_metagene_decile`, `figS_tss5kb_delta_metagene`,
`figS_mt_methylation_3groups`, `figS_global_methylation_3groups`,
`figS_region_decile_metagene_wgbs_tail`, `fig2i_sample_pca_correlation`

The run also writes `sessionInfo_02_landscape.txt` in the module folder.

## Run

From the pipeline root, after the modules it reads have finished:

```bash
sbatch 02_landscape/code/02_landscape.slurm
```

Resources requested by the launcher: 4 CPUs, 96 GB, 12 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 02_landscape/code/02_landscape.R`.
