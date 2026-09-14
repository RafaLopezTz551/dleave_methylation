# DNA methylation landscape during tail regeneration in the land slug *Deroceras laeve*

Analysis code for the manuscript of the same title (López Téllez, Gutiérrez Sarmiento,
Trejo Arellano, Brito Domínguez, Varela Echavarría and Miranda Rodríguez). The terrestrial
slug *D. laeve* regrows its tail after amputation although its genome lacks DNMT3. We
sequenced the methylome of the early tail blastema and of control tail (whole genome
bisulfite sequencing, two animals per condition) together with tail RNA-seq, added a
PacBio HiFi native methylome of bodywall, and asked how methylation is organised, what
amputation changes, and whether those changes track transcription. Methylation is mosaic
and gene body dominated, promoters and transposable elements are sparsely methylated,
amputation remodels the methylome without DNMT3, and the change does not predict
expression gene by gene while concentrating in particular coexpression modules.

## Layout

Eight sequential modules, one R script each. A module reads raw data or the outputs of
lower numbered modules only and writes only into its own folder (`data/` tables,
`objects/` R objects, `figures/main` and `figures/supplementary`).

| Folder | Question | Paper figures |
|---|---|---|
| `00_preprocessing/` | raw reads → Bismark CpG reports, HTSeq counts, HiFi pileups (documentation + the run scripts) | — |
| `01_genome_toolkit/` | genome composition, CpG depletion across molluscs, methylation toolkit, tail differential expression, and the TF motif library by sequence orthology | Fig. 1; S1, S2, S10B, S14, S15 |
| `02_landscape/` | the baseline methylome: global level, compartments, HiFi agreement, gene body methylation vs expression | Figs. 2, 3; S3–S6 |
| `03_promoters/` | Weber promoter classes, human control, promoter methylation and GO by class | Figs. 4, 5B–C; S7 |
| `04_TEs/` | transposable element methylation (HiFi bodywall main, WGBS tail supplementary) | Fig. 6; S9 |
| `05_differential/` | DSS DMPs and DMRs, their annotation and GO, the label swap null | Figs. 7, 8A; S10A; Table S1 |
| `06_decoupling/` | methylation change vs expression change | S11 |
| `07_wgcna/` | coexpression modules and where methylation change concentrates | Fig. 9; S12 |
| `08_motifs/` | UMR/LMR segmentation, motif enrichment by Weber class and at UMRs/LMRs | Figs. 5A, 10; S8, S13 |
| `environment/` | modules, R package versions (sessionInfo per module) | |
| `data/` | small derived tables quoted in the paper; large files are pointed to SRA/GEO | |

Figure numbers refer to the submitted manuscript; each figure file name inside
`NN_name/figures/` starts with its panel stem (for example `fig3b_…`, `figS8_…`).

## How to reproduce

1. Raw reads: NCBI SRA BioProject (accession in the paper). Run `00_preprocessing/` (or
   start from the Bismark CpG reports, the HTSeq counts and the MethBat pileups).
2. Reference data: the *D. laeve* assembly GCA_051403575 restricted to `chr1..chr31` plus
   the mitochondrial scaffold, the EviAnn annotation, the STRING v12 term file for
   *D. laeve*, and the staged external inputs listed in `environment/README.md`.
3. Edit the ALL-CAPS path block at the top of each `NN_name/code/NN_name.R` (and the `cd`
   line of its `.slurm`) to your layout.
4. Run the modules in order, each as one Slurm job:
   `sbatch 01_genome_toolkit/code/01_genome_toolkit.slurm`, then `02_…`, up to `08_…`
   (or `Rscript NN_name/code/NN_name.R` from inside `NN_name/`). Memory and time requests
   are in each launcher; the whole chain takes about ten hours on one node.
5. Every module writes `sessionInfo_<module>.txt`; compare with `environment/`.

Fixed seed `20260426`, chromosome filter, downstream only dependencies and the other
rules that keep the pipeline honest are in `CLAUDE.md`. A module level description of
each analysis, detailed enough to re-implement it, is in `docs/METHODS.md`.

## Data availability

Raw WGBS, RNA-seq and PacBio HiFi reads: NCBI SRA (BioProject accession in the
manuscript). Processed per CpG tables and the tables quoted in the text: `data/` here and
the Supplementary Data of the paper.
