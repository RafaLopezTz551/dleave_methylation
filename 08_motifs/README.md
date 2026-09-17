# 08_motifs

Transcription factor motif enrichment in the *D. laeve* methylome: MethylSeekR UMR/LMR segmentation, monaLisa binned enrichment along the LMR methylation-change gradient and across Weber promoter classes (slug and human), and HOMER known-motif enrichment of UMR promoters, LMRs and DMRs.

**Inputs**

- `01_genome_toolkit/objects/genome_chrmt.rds`, `01_genome_toolkit/objects/gff_chrmt.rds`, `01_genome_toolkit/data/jaspar_ortholog_bridge.tsv`.
- `02_landscape/objects/bsseq_cov5_chrmt.rds`.
- `03_promoters/dataset/` (GRCh38 reference cache).
- `05_differential/data/dmrs_annotated.tsv`.
- JASPAR2024 CORE SQLite; TE annotation table (`collapsed_te_age_data.tsv`).
- HOMER 5.1 (`findMotifsGenome.pl`) and `rsvg-convert`, called through `system2()`.

Paths are the ALL-CAPS constants at the top of `code/08_motifs.R`.

**Outputs**

- `data/` — tables (22 files, named in the script).
- `objects/` — `methylseekr_segments_cov10_chrmt_gr.rds`, `cgi_takai_jones_chrmt.rds`, `jaspar_ortholog_homer.motif`, `lmr_deltameth_cov10_chrmt_se.rds`, `promoter_weber_pc_chrmt_se.rds`, `genome_chrmt.fa`, `human_promoter_weber_hs9606_se.rds`, `homer_<set>_chrmt/` ….
- `figures/main/` — `fig8_lmr_motif_enrichment`, `fig8_umr_promoter_homer_motifs`; more in `figures/supplementary/`.
- Paper figures: Fig. 5, Fig. 10, Fig. S8, Fig. S13 (see `docs/FIGURES.md`).

**Run**

`sbatch 08_motifs/code/08_motifs.slurm` (16 CPUs, 128 GB, 24 h), after `01_genome_toolkit`, `02_landscape`, `03_promoters`, `05_differential`. Without SLURM: `Rscript 08_motifs/code/08_motifs.R`. Needs HOMER 5.1 and `rsvg-convert` on the paths set in the script.
