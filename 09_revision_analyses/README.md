# 09_revision_analyses

Analyses added during the internal co-author review of the manuscript (October 2026). The module reads the objects and
tables of the lower-numbered modules and changes none of them.

**What it adds**

- Fig. 1B: genome-wide CpG O/E of five mollusks with the regenerative capacity documented for each species.
- Fig. 1E: the methylation toolkit with absent genes in grey and the eggNOG E-value of each hit beside its bit score.
- Fig. 2B: methylation per feature (promoter, exon, intron, intergenic tile) as violins, with Kruskal-Wallis and
  Benjamini-Hochberg-adjusted pairwise Wilcoxon tests.
- High-methylation calls tested against the bisulfite non-conversion rate (one-sided binomial, BH), pooled and per library.
- HiFi vs WGBS: the upper mode of each platform at shared CpGs, HiFi values by coverage, and the WGBS − HiFi gap per gene
  against tissue-specific expression (tail vs bodywall RNA-seq), to separate a platform from a tissue effect.
- Fig. S3B: genome-wide profile as a rolling mean over five 1 Mb bins.
- TE copies per class and location (Fig. 6 legend).

**Inputs**

- `02_landscape/objects/bsseq_cov5_chrmt.rds`, `02_landscape/objects/hifi_bodywall_cpg_persample_chrmt.rds`,
  `02_landscape/data/genomewide_methylation_1mb{,_bodywall}.tsv`.
- `01_genome_toolkit/data/mollusc_dinucleotide_oe.tsv`, `toolkit_presence.tsv`, `objects/gff_chrmt.rds`; the eggNOG-mapper
  annotation of the proteome.
- `04_TEs/data/non_conversion_rate.tsv`, `te_methylation_per_copy.tsv`, `te_methylation_per_copy_hifi.tsv`.
- `06_decoupling/data/decoupling_per_gene.tsv`; HTSeq counts of the four bodywall libraries.

**Outputs**

- `data/` — `toolkit_presence_with_evalue.tsv`, `region_feature_beta_summary.tsv`, `region_pairwise_wilcoxon.tsv`,
  `high_methylation_false_positive_check.tsv`, `upper_mode_by_platform.tsv`, `hifi_at_wgbs_high_cpgs.tsv`,
  `hifi_at_wgbs_full_cpgs.tsv`, `hifi_at_wgbs_full_by_coverage.tsv`, `gap_by_expression_concordance.tsv`,
  `gene_platform_gap.tsv`, `te_copies_per_class_plotted.tsv`, `revision_summary.tsv`.
- `figures/main/` — `fig1e_mollusc_cpg_oe`, `fig2a_methylation_toolkit_presence`, `fig2e_region_methylation`;
  `figures/supplementary/` — `figS_genomewide_1mb_rolling`.
- Paper figures: Fig. 1BE, Fig. 2B, Fig. S3B (see `docs/FIGURES.md`); these supersede the versions of the same
  stems drawn by `01_genome_toolkit` and `02_landscape`.

**Run**

`sbatch 09_revision_analyses/code/09_revision_analyses.slurm` (2 CPUs, 40 GB, about 20 min), after modules 01–08.
Without SLURM: `Rscript 09_revision_analyses/code/09_revision_analyses.R`.
