# 06_decoupling

Methylation change versus expression change: gene-body methylation against expression, differential methylation against differential expression, and DSS track figures for selected DMRs.

**Inputs**

- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`.
- HTSeq tail gene counts.
- `05_differential/data/dmrs_annotated.tsv`, `05_differential/data/dmrs_gene_assignments.tsv`, `05_differential/data/dmp_dmr_de_intersection.tsv`.
- `03_promoters/data/promoter_weber_classification.tsv`.

Paths are the ALL-CAPS constants at the top of `code/06_decoupling.R`.

**Outputs**

- `data/` — tables (6 files, named in the script).
- `figures/main/` — `fig6a_decoupling_three_axes`; more in `figures/supplementary/`.
- Paper figures: Fig. S11 (see `docs/FIGURES.md`).

**Run**

`sbatch 06_decoupling/code/06_decoupling.slurm` (2 CPUs, 48 GB, 24 h), after `01_genome_toolkit`, `02_landscape`, `03_promoters`, `05_differential`. Without SLURM: `Rscript 06_decoupling/code/06_decoupling.R`.
