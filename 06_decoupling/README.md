# 06_decoupling

Methylation change versus expression change: gene-body methylation against expression, differential methylation against differential expression, and DSS track figures for selected DMRs.

## Inputs

- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`.
- HTSeq tail gene counts.
- `05_differential/data/dmrs_annotated.tsv`, `05_differential/data/dmrs_gene_assignments.tsv`, `05_differential/data/dmp_dmr_de_intersection.tsv`.
- `03_promoters/data/promoter_weber_classification.tsv`.

Paths are set in the ALL-CAPS constant block at the top of `code/06_decoupling.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Seed, packages, paths and chromosome filter.
2. Gene-body methylation per gene, pooled and per condition.
3. Tail expression (HTSeq counts to VST) and the upstream DE and DMR tables.
4. Merge and decoupling statistics (Fisher test; length-adjusted Cochran-Mantel-Haenszel).
5. DMR methylation change versus host-gene log2 fold change (Spearman).
6. The decoupling on three axes: direction, occurrence, magnitude (`fig6a_decoupling_three_axes`).
7. Supplementary DSS tracks for the DMP, DMR and DE intersection genes (`figS6_intersection_dmrs.pdf`, `figS6_dmr_<symbol>`).
8. Supplementary: HCP (CpG-rich) promoters carrying a DMR (`figS6_hcp_promoter_dmrs.pdf`, `figS6_hcp_dmr_<symbol>`).
9. Record the package versions (`sessionInfo_06_decoupling.txt`).

## Outputs

`data/`

- `decoupling_per_gene.tsv`
- `decoupling_length_strata.tsv`
- `decoupling_summary.tsv`
- `dmr_de_correlation.tsv`
- `intersection_dmr_tss.tsv`
- `hcp_promoter_dmr_genes.tsv`

`objects/`

- none

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

- `fig6a_decoupling_three_axes`

`figures/supplementary/` (the per-DMR track panels as `.png` and `.svg`, the collections as multi-page `.pdf`)

- `figS6_intersection_dmrs.pdf`
- `figS6_dmr_<symbol>`
- `figS6_hcp_promoter_dmrs.pdf`
- `figS6_hcp_dmr_<symbol>`

The run also writes `sessionInfo_06_decoupling.txt` in the module folder.

## Run

From the pipeline root, after the modules it reads have finished:

```bash
sbatch 06_decoupling/code/06_decoupling.slurm
```

Resources requested by the launcher: 2 CPUs, 48 GB, 24 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 06_decoupling/code/06_decoupling.R`.
