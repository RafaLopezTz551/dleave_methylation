# 01_genome_toolkit

Genome CpG composition of *Deroceras laeve*, the DNA methylation toolkit it encodes and expresses, CpG depletion across mollusc assemblies, tail differential expression, and the JASPAR-to-*D. laeve* transcription factor motif library built by sequence orthology.

## Inputs

- Genome FASTA restricted to chr1-31 plus the mitochondrial scaffold, and the full assembly FASTA (GCA_051403575).
- EviAnn GFF and proteome; eggNOG-mapper annotations of the proteome.
- HTSeq gene counts: tail control and amputated libraries, and the multi-tissue atlas.
- Pfam-A HMM library pressed for HMMER 3.4; staged DNMT and UHRF query proteins; DNMT reference sequences (`dataset/dnmt_refs`).
- Four additional mollusc assemblies (`dataset/mollusc_genomes`).
- JASPAR2024 CORE SQLite; Cis-BP per-family motif-transfer thresholds (`dataset/cisbp`); UniProt bait sequences (staged under `objects/seq/`; fetched with `curl` only when the cache is absent).
- No earlier pipeline module is read.

Paths are set in the ALL-CAPS constant block at the top of `code/01_genome_toolkit.R`.

## Steps

The script is divided by `# Step N - ...` banners in this order.

1. Seed, packages, paths, chromosome set, palette, theme and figure savers.
2. Load the genome and the GFF (chr1-31 plus mitochondrial scaffold); cache them as `objects/genome_chrmt.rds` and `objects/gff_chrmt.rds`.
3. Genome dinucleotide frequencies, observed versus expected (`fig1a`).
4. Per-region, per-chromosome CpG density, CpG O/E and GC content.
5. Density plots of the per-chromosome region statistics (`fig1b`, `fig1c`, `fig1d`).
6. Methylation toolkit presence from eggNOG-mapper orthology (`fig2a`).
7. DNMT3 absence: tBLASTn of DNMT query proteins against the genome (DNMT1, DNMT2 and TET3 as controls) and HMMER/Pfam scans for catalytic C5-MTase and DNMT3-specific ADD domains; best-hit and domain-architecture panels (`figS_dnmt3_absence`).
8. Toolkit mRNA in tail, control versus amputated (`figS_toolkit_mrna_tail`).
9. Tail gene-level differential expression (DESeq2 Wald test, apeglm shrinkage), Argonaute superfamily census from Pfam Piwi domains, volcano and top-gene figures (`figS_gene_de_volcano`, `figS_top_de_genes`).
10. UHRF and DNMT1 domain architecture against human UHRF1/UHRF2 (`figS_uhrf_domains`).
11. Genome-wide CpG O/E across five mollusc assemblies (`fig1e`, `figS_dinuc_<species>`).
12. Toolkit expression across the tissue atlas (`figS_toolkit_mrna_atlas`).
13. Maximum-likelihood tree of C5-methyltransferase (PF00145) domains with MAFFT, trimAl and IQ-TREE (`figS_dnmt_c5_tree`).
14. TF annotation, JASPAR motifs to *D. laeve* orthologs: JASPAR2024 CORE matrices to UniProt accessions; bait sequences; DIAMOND blastp in both directions with isoforms collapsed to locus and reciprocal best hits; Pfam DNA-binding-domain licensing of bait and locus; pairwise DBD identity (MAFFT) against the Cis-BP per-family thresholds; per-family DBD gene trees (MAFFT, trimAl, IQ-TREE) and co-ortholog clades; the motif-to-gene bridge, the TF guide, funnel, conservation and decision-grid tables (`figS_tf_tf_family_complement`, `figS_tf_tf_annotation_workflow`).
15. Record the package versions (`sessionInfo_01_genome_toolkit.txt`).

## Outputs

`data/`

`dinucleotide_frequencies.tsv`, `region_chr_cpg_stats.tsv`, `toolkit_presence.tsv`,
`dnmt3_tblastn_besthit.tsv`, `dnmt3_domain_architecture.tsv`, `toolkit_mrna_tail.tsv`,
`gene_de_tail.tsv`, `argonaute_census.tsv`, `uhrf_dnmt1_domain_architecture.tsv`,
`mollusc_dinucleotide_oe.tsv`, `toolkit_mrna_atlas_per_library.tsv`, `toolkit_mrna_atlas.tsv`,
`dnmt_c5_tree_members.tsv`, `dnmt_c5_tree_unrooted.nwk`, `dnmt_c5_tree_midpoint.nwk`,
`dnmt_c5_tree_iqtree_report.txt`, `matrix_to_acc.tsv`, `dbd_identity.tsv`,
`tree_coortholog_clades.tsv`, `motif_to_dlaeve.tsv`, `jaspar_ortholog_bridge.tsv`,
`tf_dlaeve_guide.tsv`, `tf_dlaeve_orthologs.faa`, `tf_orthology_funnel.tsv`,
`dbd_conservation_by_family.tsv`, `pwm_library_decision_grid.tsv`, `tf_family_complement.tsv`

`objects/`

- `genome_chrmt.rds`
- `gff_chrmt.rds`
- `seq/` (bait and locus protein sequences)
- `diamond/` (forward and reverse hits)
- `pfam/` (domain tables and DBD profiles)
- `trees/<family>/` (per-family gene trees)

`figures/main/` (each stem as `.pdf`, `.png` and `.svg`)

- `fig1a_dinucleotide_freq`
- `fig1b_cpg_density_distribution`
- `fig1c_cpg_oe_distribution`
- `fig1d_gc_distribution`
- `fig1e_mollusc_cpg_oe`
- `fig2a_methylation_toolkit_presence`

`figures/supplementary/`

`figS_dnmt3_absence`, `figS_toolkit_mrna_tail`, `figS_gene_de_volcano`, `figS_top_de_genes`,
`figS_uhrf_domains`, `figS_dinuc_<species>`, `figS_toolkit_mrna_atlas`, `figS_dnmt_c5_tree`,
`figS_tf_tf_family_complement`, `figS_tf_tf_annotation_workflow`

The run also writes `sessionInfo_01_genome_toolkit.txt` in the module folder.

## Run

From the pipeline root:

```bash
sbatch 01_genome_toolkit/code/01_genome_toolkit.slurm
```

Resources requested by the launcher: 8 CPUs, 96 GB, 16 h. Edit the `cd` line of the launcher and the path constants at the top of the script before running. Without SLURM: `Rscript 01_genome_toolkit/code/01_genome_toolkit.R`. The launcher also loads the environment modules `diamond/2.1.0 mafft/7.299 trimal/1.2 iqtree3/3.0.1`.
