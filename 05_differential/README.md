# 05_differential

Differential methylation after tail amputation: DSS DMP and DMR calling, gene assignment, region and GO/KEGG enrichment, a DMR-in-TE permutation test, and control analyses including a label-swap null.

**Inputs**

- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`, `01_genome_toolkit/data/argonaute_census.tsv`.
- STRING v12 per-protein GO/KEGG terms.
- TE copy table with Kimura ages.
- OMArk conserved-unknown protein FASTA.

Paths are the ALL-CAPS constants at the top of `code/05_differential.R`.

**Outputs**

- `data/` — tables (33 files, named in the script).
- `objects/` — `dmltest_chrmt_<chr>.rds`, `dmltest_chrmt.rds`, `dmltest_<swap>_chrmt_<chr>.rds`, `dmltest_<swap>_chrmt.rds`.
- `figures/main/` — `fig5a_dmp_volcano`, `fig5b_dmr_volcano`, `fig5c_dmp_region_pie`, `fig5d_dmr_region_pie`, `fig5e2_dmr_region_enrichment`, `fig5e_dmp_region_enrichment`, `fig5f_top_dmp_burden_genes`, `fig5h_dmp_go_dotplot`; more in `figures/supplementary/`.
- Paper figures: Fig. 7, Fig. 8, Fig. S10, Table S1 (see `docs/FIGURES.md`).

**Run**

`sbatch 05_differential/code/05_differential.slurm` (2 CPUs, 128 GB, 48 h), after `01_genome_toolkit`, `02_landscape`. Without SLURM: `Rscript 05_differential/code/05_differential.R`.
