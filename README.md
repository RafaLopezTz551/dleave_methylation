# DNA methylation landscape during tail regeneration in the land slug *Deroceras laeve*

Analysis code for the manuscript of the same title. The terrestrial slug *D. laeve* regrows
its tail after amputation although its genome lacks DNMT3. We sequenced the methylome of the
early tail blastema and of control tail (whole genome bisulfite sequencing, two animals per
condition) together with tail RNA-seq, added a PacBio HiFi native methylome of bodywall, and
asked how methylation is organised, what amputation changes, and whether those changes track
transcription. Methylation is mosaic and gene body dominated, promoters and transposable
elements are sparsely methylated, amputation remodels the methylome without DNMT3, and the
change does not predict expression gene by gene while concentrating in particular
coexpression modules.

The pipeline is eight sequential R modules, one script each. A module reads raw data or the
outputs of lower numbered modules only and writes only into its own folder. A module level
description of each analysis, detailed enough to re-implement it, is in
[`docs/METHODS.md`](docs/METHODS.md).

## Repository layout

| Folder | Contents |
|---|---|
| `00_preprocessing/` | Run scripts and notes for the steps upstream of the pipeline: Bismark alignment and CpG reports (WGBS), STAR and HTSeq counts (RNA-seq), MethBat per-site pileups (PacBio HiFi). Documentation only; nothing here is called by the analysis scripts. |
| `01_genome_toolkit/` to `08_motifs/` | One analysis module each: `code/NN_name.R` (the script), `code/NN_name.slurm` (the launcher) and `README.md`. Outputs go to `data/` (tables), `objects/` (R objects and caches), `figures/main/` and `figures/supplementary/`. Modules 01 and 03 also carry a `dataset/` folder with staged external inputs. |
| `docs/` | `METHODS.md`, the module by module methods. |
| `environment/` | R and Bioconductor versions, the environment modules loaded by the launchers, the `sessionInfo` files of the final run, and the list of project-local tools. |
| `data/` | Small derived tables quoted in the paper; large files are pointed to the sequence archives. |

## Pipeline overview

| Module | Question it answers | Main inputs | Key outputs | Paper figures |
|---|---|---|---|---|
| `01_genome_toolkit` | Genome CpG composition, CpG depletion across molluscs, the methylation toolkit encoded and expressed, tail differential expression, the JASPAR-to-*D. laeve* TF motif library | Genome FASTA, EviAnn GFF and proteome, eggNOG annotations, HTSeq counts, Pfam-A, mollusc assemblies, JASPAR2024, Cis-BP thresholds | `objects/genome_chrmt.rds`, `objects/gff_chrmt.rds`, `data/gene_de_tail.tsv`, `data/jaspar_ortholog_bridge.tsv` | Fig. 1; S1, S2, S10B, S14, S15 |
| `02_landscape` | The baseline methylome: global level, compartments, gene body methylation vs expression, HiFi vs WGBS agreement | Bismark CpG reports, GFF, HTSeq counts, HiFi pileups, module 01 | `objects/bsseq_cov5_chrmt.rds`, `objects/hifi_bodywall_cpg_persample_chrmt.rds`, `data/genebody_methylation_per_gene.tsv` | Figs. 2, 3; S3 to S6 |
| `03_promoters` | Weber promoter classes, human control, promoter methylation, expression and GO by class | Modules 01 and 02, HTSeq counts, STRING v12 terms, GRCh38 cache, KEGG tables | `data/promoter_weber_classification.tsv`, class methylation and GO tables | Fig. 4, Fig. 5B and 5C; S7 |
| `04_TEs` | Transposable element methylation by class, age and location (HiFi bodywall main, WGBS tail supplementary) | TE table with Kimura ages, modules 01 and 02, Bismark splitting reports | Per-copy methylation, platform comparison and statistics tables | Fig. 6; S9 |
| `05_differential` | DSS DMPs and DMRs, their annotation, enrichment and GO, TE permutation test, label-swap null | Modules 01 and 02, STRING v12 terms, TE table, OMArk protein set | `data/dmps_annotated.tsv`, `data/dmrs_annotated.tsv`, `data/dmr_gene_burden.tsv`, DSS caches | Figs. 7, 8A; S10A; Table S1 |
| `06_decoupling` | Methylation change vs expression change | Modules 01, 02, 03 and 05, HTSeq counts | Decoupling tables, DMR track figures | S11 |
| `07_wgcna` | Coexpression modules and where methylation change concentrates | HTSeq counts (all tissues), modules 01, 02, 03 and 05, STRING v12, JASPAR2024, eggNOG | `objects/wgcna.rds`, `data/module_assignments.tsv`, enrichment and hub tables | Fig. 9; S12 |
| `08_motifs` | UMR/LMR segmentation, motif enrichment by Weber class and along the LMR methylation change, HOMER known motifs | Modules 01, 02, 03 and 05, JASPAR2024, TE table, HOMER | MethylSeekR segments, CGI annotation, motif enrichment tables | Fig. 5A, Fig. 10; S8, S13 |

Figure numbers refer to the submitted manuscript. Each figure file inside `NN_name/figures/`
starts with its stem, and the manuscript includes the `.svg` of each stem. The table below
lists the stems included by the manuscript, by figure.

| Paper figure | Module | File stems (`figures/main/` unless marked supplementary) |
|---|---|---|
| Fig. 1 | `01_genome_toolkit` | `fig1a_dinucleotide_freq`, `fig1e_mollusc_cpg_oe`, `fig1c_cpg_oe_distribution`, `fig1d_gc_distribution`, `fig2a_methylation_toolkit_presence` |
| Fig. 2 | `02_landscape` | `fig2a_global_methylation_per_sample`, `fig2b_genomewide_1mb`, `fig2e_region_methylation` |
| Fig. 3 | `02_landscape` | `fig2f_genebody_methylation_decile`, `fig2g_region_decile_metagene_bodywall`, `fig2d_tss5kb_metagene_decile` |
| Fig. 4 | `03_promoters` | `fig3e_promoter_cpg_oe_human_vs_dlaeve`, `fig3b_promoter_oe_distribution_by_class`, `fig3d_metagene_by_promoter_class` |
| Fig. 5 | `08_motifs`, `03_promoters` | 5A `figS8_promoter_weber_motif_enrichment` (supplementary folder of 08); 5B `fig3g_dlaeve_weber_class_go`; 5C `figS3_human_weber_class_go` (supplementary folder of 03) |
| Fig. 6 | `04_TEs` | `fig4a_te_methylation_by_class_bodywall`, `fig4b_te_methylation_by_class_location_bodywall`, `fig4c_te_age_by_class_location_bodywall` |
| Fig. 7 | `05_differential` | `fig5a_dmp_volcano`, `fig5d_dmr_region_pie`, `fig5e2_dmr_region_enrichment`, `fig5h_dmp_go_dotplot` |
| Fig. 8 | `05_differential` | 8A `fig5j_dmp_dmr_de_venn` (supplementary folder); 8B and 8C are in situ hybridisation images, not produced by the pipeline |
| Fig. 9 | `07_wgcna` | `fig7a_module_dmp_fisher`, `fig7c_eigengene_scores_tail` |
| Fig. 10 | `08_motifs` | `fig8_umr_promoter_homer_motifs`, `fig8_lmr_motif_enrichment` |
| Fig. S1 | `01_genome_toolkit` | `figS_dnmt3_absence` |
| Fig. S2 | `01_genome_toolkit` | `figS_uhrf_domains` |
| Fig. S3 | `02_landscape` | `figS_global_methylation_3groups` |
| Fig. S4 | `02_landscape` | `figS_mt_methylation_3groups` |
| Fig. S5 | `02_landscape` | `figS_region_decile_metagene_wgbs_tail` |
| Fig. S6 | `02_landscape` | `figS_tss5kb_delta_metagene` |
| Fig. S7 | `03_promoters` | `figS3_ap2_domain` |
| Fig. S8 | `08_motifs` | `figS8_promoter_weber_gc` |
| Fig. S9 | `04_TEs` | `figS4_te_methylation_by_class_wgbs_tail`, `figS4_te_methylation_by_class_location_wgbs_tail`, `figS4_te_age_by_class_location_wgbs_tail` |
| Fig. S10 | `05_differential`, `01_genome_toolkit` | S10A `fig5l_gene_known_unknown_pie`; S10B `figS_gene_de_volcano` |
| Fig. S11 | `06_decoupling` | `figS6_dmr_unc_9`, `figS6_dmr_PDLIM3`, `figS6_dmr_Mmp16`, `figS6_dmr_Chrdl2`, `figS6_dmr_LOC_00015983` |
| Fig. S12 | `07_wgcna` | `fig7b_module_go_enrichment` |
| Fig. S13 | `08_motifs` | `figS8_lmr_umr_region_pie`, `figS8_lmr_bin_density` |
| Fig. S14 | `01_genome_toolkit` | `figS_tf_tf_annotation_workflow` |
| Fig. S15 | `01_genome_toolkit` | `figS_tf_tf_family_complement` |
| Table S1 | `05_differential` | `data/dmr_gene_burden.tsv` |

Supplementary stems live in `figures/supplementary/` of their module. Each module writes
more figures than the manuscript includes; the per-module `README.md` lists all of them.

## How to reproduce

### Requirements

- R 4.4.1 with Bioconductor 3.20. Package versions are recorded in
  `environment/sessionInfo_NN_name.txt`, one file per module, written at the end of each run.
- Bioconductor packages used across modules: `GenomicRanges`, `IRanges`, `Biostrings`,
  `rtracklayer`, `bsseq`, `DESeq2`, `SummarizedExperiment`, `BiocParallel`. Module specific:
  `DSS` (05, 06), `goseq` (05), `clusterProfiler` and `org.Hs.eg.db` (03, 05), `WGCNA`
  (07), `MethylSeekR`, `monaLisa`, `TFBSTools`, `universalmotif`, `BSgenome` and
  `ComplexHeatmap` (08).
- CRAN packages: `data.table`, `ggplot2` (4.0.2), `patchwork`, `scales`, `ggrepel`,
  `svglite`, `ape`, `phangorn`, `RSQLite`, `EnhancedVolcano` (01), `ggridges` (04),
  `hexbin`, `VennDiagram` (05), `fastcluster`, `dynamicTreeCut` (07).
- Reference data: the *D. laeve* assembly GCA_051403575 restricted to `chr1` to `chr31` plus
  the mitochondrial scaffold, the EviAnn annotation and proteome, the STRING v12 term file
  for *D. laeve*, and the staged external inputs listed in `environment/README.md`
  (JASPAR2024 SQLite, Pfam-A, cached KEGG tables, GRCh38 RefSeq cache, mollusc assemblies,
  Cis-BP thresholds). The scripts assume compute nodes without network access: every external
  resource is staged once and read by path (module 01 fetches UniProt bait sequences and
  module 03 the GRCh38 cache only when the staged copy is missing).

External tools are called from inside the R scripts with `system2()`, so the exact command
lines are in the scripts:

| Tool | Used by | Purpose |
|---|---|---|
| BLAST+ 2.13.0 (`tblastn`) | 01 | DNMT queries against the genome |
| HMMER 3.4 (`hmmscan`, `hmmsearch`) with Pfam-A | 01, 03 | Protein domain scans |
| DIAMOND 2.1.0 | 01 | Reciprocal best hits for the TF orthology |
| MAFFT 7.299, trimAl 1.2, IQ-TREE 3.0.1 | 01 | Domain alignments and gene trees |
| HOMER 5.1 (`findMotifsGenome.pl`) | 08 | Known-motif enrichment |
| `rsvg-convert` | 08 | Motif logo composition |
| `curl`, `wget` | 01, 03 | Only as a fallback when a staged external file is missing |

### Inputs from `00_preprocessing`

The analysis starts from the Bismark CpG reports and splitting reports (one per WGBS library),
the HTSeq gene counts (tail and tissue atlas) and the MethBat per-site 5mC pileups of the two
bodywall HiFi samples. `00_preprocessing/README.md` documents how these were produced from the
raw reads.

### Path constants

Each script has one block of ALL-CAPS path constants at its top (genome, GFF, count folders,
methylation calls, the pipeline root `PIPE` and the module folder `BATCH`). Edit those lines
to the location of your data and nothing else. Each `.slurm` launcher also has a `cd` line
with the absolute module path and an `R_LIBS` export; edit both to your layout.

### Running with SLURM

From the pipeline root, submit the modules in order; each one reads only modules with a
lower number, so module N must finish before module N+1 starts.

```bash
sbatch 01_genome_toolkit/code/01_genome_toolkit.slurm
sbatch 02_landscape/code/02_landscape.slurm
sbatch 03_promoters/code/03_promoters.slurm
sbatch 04_TEs/code/04_TEs.slurm
sbatch 05_differential/code/05_differential.slurm
sbatch 06_decoupling/code/06_decoupling.slurm
sbatch 07_wgcna/code/07_wgcna.slurm
sbatch 08_motifs/code/08_motifs.slurm
```

Each launcher sets `R_LIBS`, creates a per-job scratch `TMPDIR`, runs the single `Rscript`
call and removes the scratch space. Module 01 also loads the environment modules
`diamond/2.1.0 mafft/7.299 trimal/1.2 iqtree3/3.0.1`. Resources requested by the launchers:

| Module | CPUs | Memory | Time |
|---|---|---|---|
| `01_genome_toolkit` | 8 | 96 GB | 16 h |
| `02_landscape` | 4 | 96 GB | 12 h |
| `03_promoters` | 2 | 64 GB | 12 h |
| `04_TEs` | 2 | 48 GB | 12 h |
| `05_differential` | 2 | 128 GB | 48 h |
| `06_decoupling` | 2 | 48 GB | 24 h |
| `07_wgcna` | 16 | 96 GB | 24 h |
| `08_motifs` | 16 | 128 GB | 24 h |

Several modules cache their slowest step under `objects/` (the bsseq object in 02, the DSS
tests in 05, the network in 07, the MethylSeekR segmentation in 08); a rerun after the cache
exists is much shorter than the requested time.

### Running without SLURM

After editing the path constants, run each script with `Rscript` from the pipeline root, in
the same order:

```bash
Rscript 01_genome_toolkit/code/01_genome_toolkit.R
Rscript 02_landscape/code/02_landscape.R
# ... up to 08_motifs
```

The scripts resolve every path from their constants, so the working directory does not
matter. Thread counts are read from the environment variable `SLURM_CPUS_PER_TASK` and fall
back to a small default when it is unset; export it to use more cores. The external tools
listed above must be reachable at the paths set in the scripts.

## Conventions

- One script per module. `NN_name/code/NN_name.R` produces every table, object and figure of
  its module; the launcher only sets the environment and calls that script.
- Downstream only. Module N reads modules 1 to N-1 and raw data, never a higher module, and
  writes only inside its own folder.
- Fixed seed. `set.seed(20260426)` is the first statement of every script.
- Chromosome set. Every data type is restricted to `chr1` to `chr31` plus the mitochondrial
  scaffold (`keep_chr` at the top of each script). Module 07 works on gene tables that
  already carry this filter.
- Figures are saved as `.pdf`, `.png` and `.svg` with the same stem.
- Every statistical test is marked in the code with a `# STAT TEST:` comment.
- Caches. The slowest step of a module is saved under `objects/` and reused on a rerun; module
  05 stops with a message if its DSS cache is older than the bsseq object it was computed from.
- Atomic writes. In the TF section of module 01, tables and sequence files are written to a
  `.part` file and renamed on completion.
- Every module ends by writing `sessionInfo_NN_name.txt` (collected in `environment/`).

## Data availability

Raw WGBS, RNA-seq and PacBio HiFi reads: NCBI SRA (BioProject accession in the manuscript).
Processed per CpG tables and the tables quoted in the text: `data/` here and the
Supplementary Data of the paper.
