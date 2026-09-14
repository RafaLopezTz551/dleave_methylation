# Methods, module by module

What each module computes, in the order it runs, with the inputs it reads and the outputs
the paper quotes. Statistical tests are named where they are made; every test in the code
carries a `# STAT TEST:` comment in the working copy of the scripts. Numbers are not
repeated here: they are in the module `data/` tables and in the manuscript.

Conventions shared by all modules: R 4.4.1 / Bioconductor 3.20; `set.seed(20260426)`;
every data type restricted to `chr1..chr31` plus the mitochondrial scaffold
`HiC_scaffold_1563`; gene sets from the genome annotation (EviAnn; protein coding and
lncRNA); a module reads only lower numbered modules; every figure written as PDF, PNG and
SVG; every p value family Benjamini–Hochberg adjusted and reported as FDR next to an effect
size.

## 01_genome_toolkit — genome composition, toolkit, tail differential expression, TF library
* **Inputs**: the assembly (32 sequence subset), the EviAnn GFF, the TE table, HTSeq tail
  counts (4 control, 3 amputated libraries), eggNOG-mapper annotations of the proteome,
  Pfam-A, four staged mollusc assemblies (NCBI), the JASPAR2024 CORE SQLite, UniProt
  sequences of each JASPAR TF (staged), Cis-BP 3.10 per family thresholds (staged, see
  `dataset/cisbp/`).
* **Dinucleotide composition and CpG O/E** by region (promoter, exon, intron, TE,
  intergenic; one value per chromosome) and genome wide for five molluscs on their whole
  assemblies; O/E = (CpG count × length) / (C count × G count).
* **Methylation toolkit** presence: eggNOG-mapper orthology at bit score ≥ 60, scored at the
  family level for genes that are paralog sets only in vertebrates (DNMT3, TET, MBD2/3,
  UHRF, EHMT, EZH, SUV39H, GADD45). DNMT3 absence re-checked by tBLASTn of human DNMT3A/B
  against the assembly plus Pfam domain architecture (hmmscan) of every hit; UHRF/DNMT1
  domain architecture likewise.
* **Tail differential expression**: HTSeq counts → DESeq2 (Wald), apeglm shrunken log2 fold
  changes; the paper's DE set is `padj < 0.05 & |log2FC| ≥ 1`.
* **TF motif library by sequence orthology** (the `[TF]` section): JASPAR matrix → UniProt
  protein → DIAMOND blastp in both directions against the *D. laeve* proteome (isoforms
  collapsed to locus) → reciprocal best hits → both proteins must carry the Pfam DNA binding
  domain the matrix's JASPAR class expects (hmmscan, gathering thresholds) → DNA binding
  domain identity (identities over the longer domain) against the Cis-BP per family threshold
  → per family gene trees on the domain envelope (MAFFT, trimAl, IQ-TREE 3 with ModelFinder
  and UFBoot2, midpoint rooted) as a second line of evidence. The deliverable is
  `data/jaspar_ortholog_bridge.tsv` (one row per JASPAR matrix, `has_ortholog` and the locus)
  and `data/tf_dlaeve_guide.tsv`; `08_motifs` reads the bridge by path.
* **Argonaute census**: every chromosome filtered locus whose eggNOG annotation carries the
  Pfam Piwi domain, best scoring isoform per locus, kept at the toolkit rule (a Preferred
  name and bit score ≥ 60), assigned to the PIWI or AGO clade by name, with tail baseMean,
  log2 fold change and adjusted P from the DE table (`argonaute_census.tsv`).

## 02_landscape — the baseline methylome
* **Inputs**: four Bismark CpG reports; the GFF object from module 01; HTSeq counts (tail and
  bodywall); the two MethBat HiFi 5mC pileups.
* **CpG object**: strand collapsed CpG dyads (minus strand counts added to the plus strand),
  kept at coverage ≥ 5 in all four libraries (`objects/bsseq_cov5_chrmt.rds`, reused by
  every later module).
* Global methylation per sample (mean β over retained CpGs), 1 Mb profiles, compartment
  means (methylated reads / total reads over the compartment's CpGs; promoter = 2 kb
  upstream; compartments overlap by construction), gene body β per gene against expression
  deciles (Spearman), TSS ± 5 kb and discrete region metagenes by decile (tail control, tail
  blastema, HiFi bodywall with bodywall deciles), the blastema − control TSS difference
  profile, mitochondrial methylation, and HiFi versus WGBS agreement (Pearson per CpG and per
  1 Mb window on the shared CpGs).

## 03_promoters — Weber promoter classes
* **Classifier**: Weber et al. 2007 on the −1300/+200 region of every protein coding TSS,
  sliding 500 bp windows at 5 bp offset; HCP = any window with CpG O/E > 0.75 and G+C > 55 %,
  LCP = no window with O/E > 0.48, else ICP. Classification uses the maximum window O/E;
  figures plot the whole promoter O/E.
* **Human control**: identical classifier on GRCh38 RefSeq protein coding promoters (primary
  chromosomes); the Weber 2007 split is reproduced.
* Promoter methylation by class (pooled β over the -1,300/+200 classifier window), TSS metagenes by
  class, expression versus promoter methylation by class (Spearman), Gene Ontology and KEGG
  overrepresentation per class (clusterProfiler `enricher` with STRING v12 terms for
  *D. laeve*, `enrichGO` + cached KEGG for human; universe = classified genes carrying an
  annotation in the category tested), and the AP-2 domain architecture.

## 04_TEs — transposable element methylation
* Per copy methylation of the five classified repeat classes (LINE, LTR, DNA, RC, SINE) on
  the HiFi bodywall methylome (main figures) and on WGBS tail (supplementary), by class,
  location (overlapping a gene body versus intergenic) and Kimura divergence (CpG adjusted
  column) quintile. Tests: Mann–Whitney with rank biserial effect size, Fisher's exact test
  on the β > 0.5 fraction, Spearman of β against divergence, oldest versus youngest quintile
  per class and location. Bisulfite non conversion floor from the CHH context of the Bismark
  splitting reports.

## 05_differential — DMPs and DMRs
* **DSS**: `DMLtest` (smoothed) per chromosome on the coverage ≥ 5 CpG object, control
  (C1, C2) versus amputated (A1, A2); `callDML` with p < 0.05 and |Δβ| ≥ 0.10 (DMPs);
  `callDMR` with p < 0.05, Δβ ≥ 0.10, minimum 50 bp and 3 CpGs (DMRs).
* Annotation to promoter (2 kb upstream) / exon / intron / intergenic with one primary gene
  per feature; region enrichment against a geometry matched background (single CpGs for
  DMPs, width matched random intervals for DMRs; Fisher); gene biotype enrichment with a
  length stratified Cochran–Mantel–Haenszel version; GO/KEGG overrepresentation of DMP and
  DMR genes (all, hyper, hypo) with the universe restricted to annotated genes of the 32
  sequence universe; DMP/DMR/DE overlap; DMRs in gene body TEs against 1,000 random
  placements.
* **Label swap null (reviewer control)**: the identical DSS procedure on the two mislabelled
  two versus two splits (C1+A1 vs C2+A2; C1+A2 vs C2+A1); counts and the fraction of true
  label DMPs recovered are in `data/label_swap_counts.tsv`.
* **Target size adjusted re-test (GOmeth)**: every DMP/DMR gene set is re-tested with the
  GOmeth procedure (missMethyl; the methylation form of the GOseq selection bias correction),
  implemented with goseq 1.58.0: a Wallenius non central hypergeometric test whose per gene
  bias is the number of analysed CpGs over the same target the features were assigned to
  (gene body plus 2 kb upstream), universe equal to the enricher universe, BH over the terms
  with at least one hit gene; gene length is run as a sensitivity bias. `data/goseq_*.tsv` and
  `goseq_summary.tsv` report, per set, the enricher significant terms and how many survive.
  `dmp_share_by_length_quintile.tsv` gives the share of genes carrying a DMP per length quintile.
* **Argonaute loci**: the DMP and DMR counts of the PIWI/AGO loci listed by module 01
  (`argonaute_methylation.tsv`).

## 06_decoupling — methylation change versus expression change
* Per gene: pooled and per condition gene body β (≥ 5 CpGs), VST expression, DESeq2 log2 fold
  change, DMR status, gene length. Baseline R² (β vs expression) and differential R²
  (Δβ vs log2FC); Fisher's exact test of DMR bearing versus differentially expressed and its
  gene length tertile stratified CMH; Spearman of per DMR Δβ against the gene's log2FC
  overall and per region. The three axes are drawn together in the main figure; DSS tracks of
  the DMP∩DMR∩DE genes and of HCP promoter DMRs are supplementary.

## 07_wgcna — coexpression modules
* Signed WGCNA on the 40 library tissue atlas (44 minus four outliers; intestine excluded),
  VST expression, variance quantile chosen by a sweep against the scale free criterion,
  soft power 20, eigengenes; experiment groups kept separate.
* DMP enrichment per module: Fisher, then Cochran–Mantel–Haenszel stratified by gene length
  quintile (the interpreted test), and (reviewer control) by length × pooled gene body β and
  by the number of movable CpGs (pooled β within 0.10–0.90 over gene body + 2 kb upstream);
  GO/KEGG of the enriched modules; tail eigengene shifts; hubs (top decile intramodular
  connectivity and |kME| ≥ 0.80) against DMP/DMR status.
* **Direction split**: the same Fisher/CMH per module for hypermethylated and hypomethylated
  DMP genes separately (`module_dmp_enrichment_by_direction.tsv`).
* **Expression breadth**: tau (Yanai) over the eleven experiment groups of the atlas, computed
  on every low count passing gene BEFORE the variance pre-filter (the filter removes the evenly
  expressed genes the test is about); Spearman of tau and of the coefficient of variation
  against pooled gene body β, Mann–Whitney of tau in DMP bearing versus other genes
  (`expression_breadth_*.tsv`).
* **With and without the length adjustment**: the unadjusted Fisher and the length
  stratified CMH odds ratio per module side by side, with each module's gene length
  distribution and the DMP share per length quintile (`length_adjustment_effect.tsv`,
  `dmp_share_by_length_quintile_network.tsv`).
* The network object is cached with a content fingerprint (gene identity, matrix dimensions,
  matrix sum, soft power); any mismatch stops the run so a changed input is never scored on a
  stale network.

## 08_motifs — UMRs, LMRs and TF motifs
* **Segmentation**: MethylSeekR on the pooled tail methylome restricted to CpGs at ≥ 10× in
  every library (methylation cutoff 0.5, ≥ 3 CpGs), with the Takai–Jones CGI track and the
  `calculateFDRs` grid as the cutoff justification.
* **Motif library**: JASPAR2024 CORE (vertebrates, insects, nematodes, urochordates) filtered
  by the module 01 bridge; HOMER thresholds converted to natural log odds.
* **LMR gradient**: per LMR Δβ (amputated − control, ≥ 3 CpGs), seven equal size bins,
  monaLisa binned enrichment against the other bins; **bin QC (reviewer control)**: width,
  CpG count, depth, control β and GC per bin with Kruskal–Wallis and Mann–Whitney tests
  (rank biserial r, positive when the extreme bin is larger), and Spearman of |Δβ| against each
  metric over all binned LMRs (`lmr_bin_qc*.tsv`, `lmr_delta_vs_composition.tsv`).
* **DMRs in LMRs and UMRs**: each DMR scored for overlap with any LMR and any UMR against
  1,000 random interval sets anchored on random analysed CpGs, once width matched and once
  matched on the number of analysed CpGs (both region classes are defined by CpG content);
  fold, two-sided empirical P and a Fisher test against the pooled random intervals
  (`dmr_lmr_umr_overlap.tsv`).
* **Promoter classes**: monaLisa enrichment across the Weber classes of *D. laeve* and of
  human (GC panel shown beside every class heatmap); the per class log2 enrichment and
  BH adjusted P of every motif are written (`promoter_weber_motif_enrichment.tsv`, human
  likewise) with the count of motifs significant in at least one class
  (`promoter_weber_motif_summary.tsv`); the human versus *D. laeve* HCP enrichment
  correlation is a two-sided Pearson test with its 95% CI, plus Spearman
  (`hcp_motif_human_vs_dlaeve_correlation.tsv`).
* **Region sets**: HOMER `findMotifsGenome.pl -size given -nomotif -mknown` on UMR promoters
  (CpG matched background), LMRs and DMRs (GC matched background); results across region
  sets are not compared because their null models differ.
