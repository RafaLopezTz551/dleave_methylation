# 01_genome_toolkit

Genome CpG composition of *Deroceras laeve*, the DNA methylation toolkit it encodes and expresses, CpG depletion across mollusc assemblies, tail differential expression, and the JASPAR-to-*D. laeve* transcription factor motif library built by sequence orthology.

**Inputs**

- Genome FASTA restricted to chr1-31 plus the mitochondrial scaffold, and the full assembly FASTA (GCA_051403575).
- EviAnn GFF and proteome; eggNOG-mapper annotations of the proteome.
- HTSeq gene counts: tail control and amputated libraries, and the multi-tissue atlas.
- Pfam-A HMM library pressed for HMMER 3.4; staged DNMT and UHRF query proteins; DNMT reference sequences (`dataset/dnmt_refs`).
- Four additional mollusc assemblies (`dataset/mollusc_genomes`).
- JASPAR2024 CORE SQLite; Cis-BP per-family motif-transfer thresholds (`dataset/cisbp`); UniProt bait sequences (staged under `objects/seq/`; fetched with `curl` only when the cache is absent).

Paths are the ALL-CAPS constants at the top of `code/01_genome_toolkit.R`.

**Outputs**

- `data/` — tables (23 files, named in the script).
- `objects/` — `genome_chrmt.rds`, `gff_chrmt.rds`, `seq/`, `diamond/`, `pfam/`, `trees/<family>/`.
- `figures/main/` — `fig1a_dinucleotide_freq`, `fig1b_cpg_density_distribution`, `fig1c_cpg_oe_distribution`, `fig1d_gc_distribution`, `fig1e_mollusc_cpg_oe`, `fig2a_methylation_toolkit_presence`; more in `figures/supplementary/`.
- Paper figures: Fig. 1, Fig. S1, Fig. S2, Fig. S10, Fig. S14, Fig. S15 (see `docs/FIGURES.md`).

**Run**

`sbatch 01_genome_toolkit/code/01_genome_toolkit.slurm` (8 CPUs, 96 GB, 16 h). Without SLURM: `Rscript 01_genome_toolkit/code/01_genome_toolkit.R`. The launcher also loads `diamond/2.1.0 mafft/7.299 trimal/1.2 iqtree3/3.0.1`.
