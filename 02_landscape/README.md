# 02_landscape

The baseline CpG methylation landscape of the *D. laeve* tail: per sample, per chromosome and per genomic region, its relation to gene expression, and the agreement between PacBio HiFi and WGBS methylation calls.

**Inputs**

- Four Bismark CpG reports (C1, C2 control; A1, A2 amputated).
- EviAnn GFF (`01_genome_toolkit/objects/gff_chrmt.rds` when present, otherwise imported from the GFF).
- HTSeq gene counts, tail and bodywall libraries.
- PacBio HiFi per-site 5mC pileups of bodywall (two animals; see `00_preprocessing/pacbio_hifi`).
- `01_genome_toolkit/data/gene_de_tail.tsv`.

Paths are the ALL-CAPS constants at the top of `code/02_landscape.R`.

**Outputs**

- `data/` — tables (18 files, named in the script).
- `objects/` — `bsseq_cov5_chrmt.rds`, `hifi_bodywall_cpg_persample_chrmt.rds`.
- `figures/main/` — `fig2a_global_methylation_per_sample`, `fig2b_genomewide_1mb`, `fig2d_tss5kb_metagene_decile_blastema`, `fig2e_region_methylation`, `fig2f_genebody_methylation_decile_blastema`, `fig2g_region_decile_metagene_blastema`, `fig2j_region_methylation_pie`; supplementary panels (control and pooled versions, the three-condition TSS metagene, the bodywall region metagene) in `figures/supplementary/`.
- Paper figures: Fig. 2A, Fig. 3ABC, Fig. S3A, Fig. S4, Fig. S5, Fig. S6 (see `docs/FIGURES.md`).

**Run**

`sbatch 02_landscape/code/02_landscape.slurm` (4 CPUs, 96 GB, 12 h), after `00_preprocessing`, `01_genome_toolkit`. Without SLURM: `Rscript 02_landscape/code/02_landscape.R`.
