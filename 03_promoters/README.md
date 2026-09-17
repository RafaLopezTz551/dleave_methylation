# 03_promoters

Weber promoter classification (HCP, ICP, LCP) of *D. laeve* genes, promoter methylation and expression by class, the same classification applied to human GRCh38 as a positive control, GO/KEGG over-representation by class, and the AP-2 protein domain scan.

## Inputs

- `01_genome_toolkit/objects/genome_chrmt.rds`, `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`.
- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- HTSeq tail gene counts.
- STRING v12 protein enrichment terms for *D. laeve*; cached KEGG hsa pathway tables; `org.Hs.eg.db`.
- GRCh38 RefSeq genome and annotation cache under `dataset/` (a `wget` fallback runs only if the cache is missing).
- EviAnn proteome and HMMER 3.4 with Pfam-A for the AP-2 domain scan.

Paths are set in the ALL-CAPS constant block at the top of `code/03_promoters.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Setup: seed, packages, paths, palette and figure helpers.
2. Load the genome, GFF and bsseq objects from modules 01 and 02.
3. Weber promoter classification from sliding 500 bp windows over TSS -1300/+200; class tables and whole-promoter CpG O/E histograms (`fig3a`, `fig3b`).
4. Pooled promoter methylation per gene over the four WGBS samples; HCP promoters ranked by pooled methylation (`fig3c`).
5. TSS metagene (+/- 2 kb) by promoter class (`fig3d`).
6. Human GRCh38 positive control: the same classification and window on RefSeq protein-coding promoters; two-species histograms with class pies (`fig3e`).
7. Expression versus promoter methylation per Weber class (`fig3f`).
8. GO/KEGG over-representation by Weber class: *D. laeve* with STRING v12 terms, human with `enrichGO` and the cached KEGG tables; dot plots per species (`fig3g_<species>_weber_class_go`, `figS3_<species>_weber_class_go`).
9. Weber-class characterisation: methylation status per class, basal expression per class, and expression of methylated versus unmethylated promoters (`figS3_weber_class_methylation_status`, `figS3_weber_class_basal_expression`, `figS3_weber_class_meth_state_expression`).
10. AP-2 protein domain architecture via `hmmscan` (`figS3_ap2_domain`).
11. Record the package versions (`sessionInfo_03_promoters.txt`).

## Outputs

`data/`

`promoter_weber_classification.tsv`, `weber_class_composition.tsv`,
`promoter_methylation_per_gene.tsv`, `metagene_by_promoter_class.tsv`, `human_promoter_cpg_oe.tsv`,
`promoter_meth_vs_expression.tsv`, `promoter_meth_expr_correlation.tsv`,
`dlaeve_weber_class_go_enrichment.tsv`, `human_weber_class_go_enrichment.tsv`,
`weber_class_methylation_status.tsv`, `weber_class_basal_expression.tsv`,
`weber_class_meth_state_expression.tsv`, `weber_class_meth_state_mannwhitney.tsv`, `ap2_domains.tsv`

`objects/`

- none

`dataset/` (intermediates of the domain scan)

- `ap2_protein.fa`
- `ap2_hmmscan.domtblout`

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

`fig3a_promoter_cpg_oe_histogram`, `fig3b_promoter_oe_distribution_by_class`, `fig3c_hcp_gene_list`,
`fig3d_metagene_by_promoter_class`, `fig3e_promoter_cpg_oe_human_vs_dlaeve`,
`fig3f_promoter_meth_vs_expression`, `fig3g_<species>_weber_class_go`

`figures/supplementary/`

- `figS3_<species>_weber_class_go`
- `figS3_weber_class_methylation_status`
- `figS3_weber_class_basal_expression`
- `figS3_weber_class_meth_state_expression`
- `figS3_ap2_domain`

`<species>` is `dlaeve` or `human`.

The run also writes `sessionInfo_03_promoters.txt` in the module folder.

## Run

From the pipeline root, after the modules it reads have finished:

```bash
sbatch 03_promoters/code/03_promoters.slurm
```

Resources requested by the launcher: 2 CPUs, 64 GB, 12 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 03_promoters/code/03_promoters.R`.
