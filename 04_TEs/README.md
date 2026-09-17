# 04_TEs

Transposable element methylation by class, age (Kimura divergence) and genomic location, on the PacBio HiFi bodywall methylome (main figures) and on the WGBS tail samples (supplementary), with platform comparisons.

**Inputs**

- TE annotation with Kimura ages (`collapsed_te_age_data.tsv`).
- `01_genome_toolkit/objects/genome_chrmt.rds`, `01_genome_toolkit/objects/gff_chrmt.rds`.
- `02_landscape/objects/bsseq_cov5_chrmt.rds` (WGBS arm) and `02_landscape/objects/hifi_bodywall_cpg_persample_chrmt.rds` (HiFi arm).
- Bismark splitting reports (CHH non-conversion proxy).

Paths are the ALL-CAPS constants at the top of `code/04_TEs.R`.

**Outputs**

- `data/` — tables (21 files, named in the script).
- `figures/main/` — `fig4a_te_methylation_by_class_bodywall`, `fig4b_te_methylation_by_class_location_bodywall`, `fig4c_te_age_by_class_location_bodywall`; more in `figures/supplementary/`.
- Paper figures: Fig. 6, Fig. S9 (see `docs/FIGURES.md`).

**Run**

`sbatch 04_TEs/code/04_TEs.slurm` (2 CPUs, 48 GB, 12 h), after `01_genome_toolkit`, `02_landscape`. Without SLURM: `Rscript 04_TEs/code/04_TEs.R`.
