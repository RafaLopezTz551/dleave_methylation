# 03_promoters

Weber promoter classification (HCP, ICP, LCP) of *D. laeve* genes, promoter methylation and expression by class, the same classification applied to human GRCh38 as a positive control, GO/KEGG over-representation by class, and the AP-2 protein domain scan.

**Inputs**

- `01_genome_toolkit/objects/genome_chrmt.rds`, `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/gene_de_tail.tsv`.
- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- HTSeq tail gene counts.
- STRING v12 protein enrichment terms for *D. laeve*; cached KEGG hsa pathway tables; `org.Hs.eg.db`.
- GRCh38 RefSeq genome and annotation cache under `dataset/` (a `wget` fallback runs only if the cache is missing).
- EviAnn proteome and HMMER 3.4 with Pfam-A for the AP-2 domain scan.

Paths are the ALL-CAPS constants at the top of `code/03_promoters.R`.

**Outputs**

- `data/` — tables (14 files, named in the script).
- `figures/main/` — `fig3a_promoter_cpg_oe_histogram`, `fig3b_promoter_oe_distribution_by_class`, `fig3c_hcp_gene_list`, `fig3d_metagene_by_promoter_class`, `fig3e_promoter_cpg_oe_human_vs_dlaeve`, `fig3f_promoter_meth_vs_expression`, `fig3g_<species>_weber_class_go`; more in `figures/supplementary/`.
- Paper figures: Fig. 4, Fig. 5, Fig. S7 (see `docs/FIGURES.md`).

**Run**

`sbatch 03_promoters/code/03_promoters.slurm` (2 CPUs, 64 GB, 12 h), after `01_genome_toolkit`, `02_landscape`. Without SLURM: `Rscript 03_promoters/code/03_promoters.R`.
