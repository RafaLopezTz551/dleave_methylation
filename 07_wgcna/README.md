# 07_wgcna

Multi-tissue WGCNA coexpression network with a methylation overlay: module-level DMP and DMR enrichment (Fisher test, with a gene-length-stratified sensitivity analysis), GO/KEGG of the enriched modules, tail eigengene scores, and hub-gene analyses.

**Inputs**

- HTSeq gene counts of all tissue libraries (the sample map is reconstructed from the file names).
- `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`.
- `02_landscape/data/genebody_methylation_per_gene.tsv`, `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `03_promoters/data/promoter_weber_classification.tsv`.
- `05_differential/data/dmps_annotated.tsv`, `05_differential/data/dmrs_annotated.tsv`.
- STRING v12 enrichment terms; JASPAR2024 SQLite; eggNOG-mapper annotations.

Paths are the ALL-CAPS constants at the top of `code/07_wgcna.R`.

**Outputs**

- `data/` — tables (28 files, named in the script).
- `objects/` — `wgcna.rds`.
- `figures/main/` — `fig7a_module_dmp_fisher`, `fig7b_module_go_enrichment`, `fig7c_eigengene_scores_tail`, `fig7d_module_dmr_fisher`; more in `figures/supplementary/`.
- Paper figures: Fig. 9, Fig. S12 (see `docs/FIGURES.md`).

**Run**

`sbatch 07_wgcna/code/07_wgcna.slurm` (16 CPUs, 96 GB, 24 h), after `01_genome_toolkit`, `02_landscape`, `03_promoters`, `05_differential`. Without SLURM: `Rscript 07_wgcna/code/07_wgcna.R`.
