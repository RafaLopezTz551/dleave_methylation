# Changelog

## 1.1.0 — 2026-10-22

Brings the repository in line with the manuscript after the co-author review (manuscript v20.5).

**Code**

- Module scripts `01`–`08` updated to the versions that produced the manuscript figures (figure changes of
  late September and October: per-condition gene-body panels, three-track genome profile, TE histograms,
  module colours in the WGCNA panel, class labels on the Weber heatmaps).
- New module `09_revision_analyses`: Fig. 1B with regenerative capacity, Fig. 1E with E-values, Fig. 2B violins with
  tests, the high-methylation false-positive test, the HiFi vs WGBS platform-vs-tissue test, the rolling-mean genome
  profile and the TE counts per class.
- Comments cleaned of internal review history; the code itself is unchanged apart from three progress messages.

**Repository**

- `supp_data/` added with Data S1–S5, the files the manuscript links to.
- `docs/FIGURES.md` regenerated from the manuscript: every main and supplementary panel mapped to its module and file.
  The LMR ΔMeth-binned motif panels (Fig. 10C and fig. S14 up to manuscript v20.4) are no longer in the paper; `08_motifs`
  still writes them, and the supplementary figures after them are renumbered (S15–S18 became S14–S17).
- `docs/METHODS.md` gains the `09_revision_analyses` section.
- `environment/` holds the sessionInfo of the final run of every module and the versions of the read processing tools
  (pbmm2, HTSeq, eggNOG-mapper, Pfam-A release), as in Table S4.
- `00_preprocessing/`: script names corrected (`01_prepare_bismark_genome.sh`, `02_alignment_bismark.sh`); the unused modkit
  attempt and the 6mA provenance script removed.
- README redesigned; `CITATION.cff`, `CHANGELOG.md` and a fuller `.gitignore` added.

## 1.0.0 — 2026-09-17

First public version.
