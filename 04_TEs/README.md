# 04_TEs

Transposable element methylation by class, age (Kimura divergence) and genomic location, on the PacBio HiFi bodywall methylome (main figures) and on the WGBS tail samples (supplementary), with platform comparisons.

## Inputs

- TE annotation with Kimura ages (`collapsed_te_age_data.tsv`).
- `01_genome_toolkit/objects/genome_chrmt.rds`, `01_genome_toolkit/objects/gff_chrmt.rds`.
- `02_landscape/objects/bsseq_cov5_chrmt.rds` (WGBS arm) and `02_landscape/objects/hifi_bodywall_cpg_persample_chrmt.rds` (HiFi arm).
- Bismark splitting reports (CHH non-conversion proxy).

Paths are set in the ALL-CAPS constant block at the top of `code/04_TEs.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Setup: seed, packages, paths, palettes, theme and figure savers.
2. Inputs: TE age table with the chromosome filter, merged repeat fraction, collapse of class/family labels to the five scored classes, WGBS per-condition M and Cov sums, gene set.
3. Per-TE-copy WGBS methylation, pooled per condition.
4. WGBS tail supplementary panels and the non-conversion floor: per-copy methylation ridges by class and condition, by class and location, and by class, location and Kimura quintile (`figS4_te_methylation_by_class_wgbs_tail`, `figS4_te_methylation_by_class_location_wgbs_tail`, `figS4_te_age_by_class_location_wgbs_tail`).
5. PacBio HiFi bodywall arm: per-copy HiFi methylation on the module 02 bodywall CpG set, platform coverage by Kimura bin and per-copy agreement, ridges by class, by class and location, and Kimura quintiles by class and location (`fig4a`, `fig4b`, `fig4c`).
6. Statistics tables for both platforms, copy counts, cross-platform floor concordance, and the WGBS control versus amputated contrast per copy.
7. Record the package versions (`sessionInfo_04_TEs.txt`).

## Outputs

`data/`

`te_genome_fraction.tsv`, `te_unclassified_counts.tsv`, `te_methylation_per_copy.tsv`,
`non_conversion_rate.tsv`, `te_kimura_methylation_correlation.tsv`, `te_age_by_class_location.tsv`,
`te_methylation_per_copy_hifi.tsv`, `te_platform_coverage_by_kimura.tsv`,
`te_platform_coverage_comparison.tsv`, `te_platform_percopy_agreement.tsv`,
`te_kimura_quintile_ranges.tsv`, `te_kimura_methylation_correlation_hifi.tsv`,
`te_age_by_class_location_hifi.tsv`, `te_statistics_bodywall.tsv`, `te_statistics_wgbs_tail.tsv`,
`te_copy_counts.tsv`, `te_platform_floor_concordance.tsv`, `te_condition_contrast_wgbs_tail.tsv`

`objects/`

- none

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

- `fig4a_te_methylation_by_class_bodywall`
- `fig4b_te_methylation_by_class_location_bodywall`
- `fig4c_te_age_by_class_location_bodywall`

`figures/supplementary/`

- `figS4_te_methylation_by_class_wgbs_tail`
- `figS4_te_methylation_by_class_location_wgbs_tail`
- `figS4_te_age_by_class_location_wgbs_tail`

The run also writes `sessionInfo_04_TEs.txt` in the module folder.

## Run

From the pipeline root, after the modules it reads have finished:

```bash
sbatch 04_TEs/code/04_TEs.slurm
```

Resources requested by the launcher: 2 CPUs, 48 GB, 12 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 04_TEs/code/04_TEs.R`.
