# DNA methylation landscape during tail regeneration in the land slug *Deroceras laeve*

Analysis code for the manuscript of the same title. WGBS of control and regenerating tail
(two animals per condition), tail RNA-seq and a PacBio HiFi bodywall methylome are used to ask
how methylation is organised in a genome without DNMT3, what amputation changes, and whether
those changes track transcription.

## Layout

- `01_genome_toolkit/` … `08_motifs/` — one analysis module each: `code/NN_name.R` (the script),
  `code/NN_name.slurm` (the launcher), `README.md`. Outputs go to `data/`, `objects/`,
  `figures/main/`, `figures/supplementary/`.
- `00_preprocessing/` — the job scripts that turned raw reads into the pipeline inputs (documentation only).
- `docs/METHODS.md` — module-by-module methods; `docs/FIGURES.md` — which file is which paper figure.
- `environment/` — R/Bioconductor versions, environment modules, `sessionInfo` of the final run.

## Modules

| Module | What it does | Paper figures |
|---|---|---|
| `01_genome_toolkit` | Genome CpG composition, CpG depletion across molluscs, the methylation toolkit encoded and expressed, tail differential expression, the JASPAR-to-*D. laeve* TF motif library | Fig. 1; S1, S2, S10B, S14, S15 |
| `02_landscape` | The baseline methylome: global level, compartments, gene body methylation vs expression, HiFi vs WGBS agreement | Figs. 2, 3; S3 to S6 |
| `03_promoters` | Weber promoter classes, human control, promoter methylation, expression and GO by class | Fig. 4, Fig. 5B and 5C; S7 |
| `04_TEs` | Transposable element methylation by class, age and location (HiFi bodywall main, WGBS tail supplementary) | Fig. 6; S9 |
| `05_differential` | DSS DMPs and DMRs, their annotation, enrichment and GO, TE permutation test, label-swap null | Figs. 7, 8A; S10A; Table S1 |
| `06_decoupling` | Methylation change vs expression change | S11 |
| `07_wgcna` | Coexpression modules and where methylation change concentrates | Fig. 9; S12 |
| `08_motifs` | UMR/LMR segmentation, motif enrichment by Weber class and along the LMR methylation change, HOMER known motifs | Fig. 5A, Fig. 10; S8, S13 |

Each module reads raw data and modules with a lower number only, and writes only inside its own folder.

## Run

R 4.4.1 / Bioconductor 3.20; package versions in `environment/`. External tools (BLAST+, HMMER,
DIAMOND, MAFFT, trimAl, IQ-TREE, HOMER) are called from the scripts with `system2()`.
Edit the ALL-CAPS path constants at the top of each script and the `cd` line of each launcher, then:

```bash
for m in 01_genome_toolkit 02_landscape 03_promoters 04_TEs 05_differential 06_decoupling 07_wgcna 08_motifs; do
  sbatch $m/code/$m.slurm      # wait for each module before submitting the next
done
```

Without SLURM: `Rscript NN_name/code/NN_name.R` in the same order (threads from `SLURM_CPUS_PER_TASK`).
Slow steps are cached under `objects/`; reruns are short.

## Conventions

- `set.seed(20260426)` is the first statement of every script.
- Every analysis is restricted to `chr1`–`chr31` plus the mitochondrial scaffold (`keep_chr`).
- Every statistical test is marked in the code with `# STAT TEST:`.
- Figures are written as `.pdf`, `.png` and `.svg` with the same stem.

## Data

Raw WGBS, RNA-seq and HiFi reads: NCBI SRA (BioProject in the manuscript). Tables quoted in the
text: `data/` of each module and the Supplementary Data of the paper.
