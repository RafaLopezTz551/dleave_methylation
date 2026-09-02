#!/usr/bin/env Rscript
set.seed(20260426)  # reproducibility

suppressPackageStartupMessages({
  library(Biostrings)      # genome sequence, di/mononucleotide counting
  library(rtracklayer)     # import the GFF
  library(GenomicRanges)   # region sets and overlaps
  library(data.table)      # fast table read/write
  library(ggplot2)         # all figures
  library(DESeq2)          # normalized counts + gene DE (toolkit heatmap, volcano)
})

# Pinned package versions (provenance; validated under R 4.4.1 / Bioc 3.20, Fenix):
#   Biostrings 2.74.1     rtracklayer 1.66.0    GenomicRanges 1.58.0   GenomeInfoDb 1.42.3
#   data.table 1.18.2.1   ggplot2 4.0.2         DESeq2 1.46.0          apeglm 1.28.0
#   EnhancedVolcano 1.24.0  patchwork 1.3.2
# Machine-checkable record: sessionInfo_01_genome_toolkit.txt, written at the end of the run.

# ---- [0] Paths, palette, theme, figure savers --------------------------------
# ID-list subset of the canonical GCA_051403575 assembly, coordinates unchanged.
GENOME <- "/mnt/data/alfredvar/rlopezt/genoma/dlaeve/genoma_chr_mt/Dlaeve_chr_mt.fasta"
# Whole canonical assembly, used ONLY in section 8 (cross-species CpG O/E), where
# every species must be measured on the same kind of sequence (whole assembly).
GENOME_FULL <- "/mnt/data/alfredvar/30-Genoma/Deroceras_laeve_genome_GCA_051403575.fasta"
GFF    <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff"
TE     <- "/mnt/data/alfredvar/30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv"
HTSEQ  <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"

BATCH  <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline/01_genome_toolkit"
OBJ    <- file.path(BATCH, "objects")
DAT    <- file.path(BATCH, "data")
FIG    <- file.path(BATCH, "figures/main")
SUPP   <- file.path(BATCH, "figures/supplementary")
for (d in c(OBJ, DAT, FIG, SUPP)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")

# Project palette (colour-blind safe, consistent across the paper)
COL_REGION <- c(Promoter = "#56B4E9", `Gene body` = "#117733",
                Exon = "#009E73", Intron = "#E6AB02", TE = "#B15928",
                Intergenic = "#CC79A7")
COL_COND   <- c(Control = "#0072B2", Amputated = "#D55E00")

# Minimal Nature/Science-style theme
theme_pub <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          legend.title = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"))
}

# Save a figure as both PDF (manuscript) and PNG (quick view)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIG, paste0(name, ".svg")), p, width = w, height = h)   # vector (svg)
  cat(sprintf("  saved %s (%g x %g in)\n", name, w, h))
}
# Same, but into figures/supplementary/
save_supp <- function(p, name, w, h) {
  ggsave(file.path(SUPP, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(SUPP, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(SUPP, paste0(name, ".svg")), p, width = w, height = h)   # vector (svg)
  cat(sprintf("  saved (supp) %s (%g x %g in)\n", name, w, h))
}

# ---- [1] Load genome + GFF (built once into objects/, reused on rerun) -------
# Writes genome_chrmt.rds + gff_chrmt.rds, the foundational objects later batches read.
cat("[1] loading genome + GFF (chr1..chr31 + mito scaffold)\n")
genome_rds <- file.path(OBJ, "genome_chrmt.rds")
gff_rds    <- file.path(OBJ, "gff_chrmt.rds")

stale_cache <- function(cache, ...) !file.exists(cache) || any(file.mtime(c(...)) > file.mtime(cache))
if (!stale_cache(genome_rds, GENOME)) {
  genome <- readRDS(genome_rds)
} else {
  genome <- readDNAStringSet(GENOME)
  names(genome) <- sub("\\s.*", "", names(genome))   # ">chr1 ..." -> "chr1"
  genome <- genome[names(genome) %in% keep_chr]
  genome <- genome[keep_chr]
  saveRDS(genome, genome_rds)
}

if (!stale_cache(gff_rds, GFF)) {
  gff <- readRDS(gff_rds)
} else {
  gff <- import(GFF)
  gff <- gff[as.character(seqnames(gff)) %in% keep_chr]
  # 32-level universe (a pre-normalisation cache once broke 03_promoters's keepSeqlevels).
  gff <- GenomeInfoDb::keepSeqlevels(gff, intersect(keep_chr, GenomeInfoDb::seqlevels(gff)),
                                     pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(gff) <- keep_chr
  saveRDS(gff, gff_rds)
}
# Re-normalise on BOTH paths (fresh or cached): EviAnn has ZERO features on the mito
# scaffold, so that seqlevel is absent from the imported GFF — keep the intersection,
# then set the full 32-level keep_chr so the GFF seqlevels match the genome exactly.
gff <- GenomeInfoDb::keepSeqlevels(gff, intersect(keep_chr, GenomeInfoDb::seqlevels(gff)),
                                   pruning.mode = "coarse")
GenomeInfoDb::seqlevels(gff) <- keep_chr
chr_len <- setNames(width(genome), names(genome))
cat(sprintf("  genome: %d chromosomes, %.1f Mb\n",
            length(genome), sum(as.numeric(chr_len)) / 1e6))

# ---- [2] fig1a: dinucleotide composition (observed vs expected) --------------
# Computes all 16 dinucleotide frequencies over the chr+mt genome; writes
# dinucleotide_frequencies.tsv and fig1a (CG bar highlighted, diamonds = expected).
cat("[2] fig1a dinucleotide frequencies\n")
di_counts <- colSums(dinucleotideFrequency(genome))     # all 16 dinucleotides
observed  <- di_counts / sum(di_counts)
mono      <- colSums(alphabetFrequency(genome, baseOnly = TRUE))[c("A","C","G","T")]
mono_freq <- mono / sum(mono)
# expected dinuc freq = product of the two single-base frequencies
expected  <- sapply(names(observed), function(d) {
  b <- strsplit(d, "")[[1]]; mono_freq[b[1]] * mono_freq[b[2]]
})

dn <- data.table(dinucleotide = names(observed),
                 observed = as.numeric(observed),
                 expected = as.numeric(expected))
dn[, ratio_obs_exp := observed / expected]
fwrite(dn[order(-observed)], file.path(DAT, "dinucleotide_frequencies.tsv"), sep = "\t")

dn[, dinucleotide := factor(dinucleotide, levels = dinucleotide[order(-observed)])]
dn[, is_cg := ifelse(dinucleotide == "CG", "CG", "other")]

pa <- ggplot(dn, aes(dinucleotide, 100 * observed, fill = is_cg)) +
  geom_col(width = 0.75) +
  geom_point(aes(y = 100 * expected), shape = 18, size = 2.2, colour = "black") +
  scale_fill_manual(values = c(CG = "#D55E00", other = "#9FB1BC"), guide = "none") +
  labs(x = "Dinucleotide", y = "Frequency (%)",
       title = "Genomic dinucleotide composition") +
  theme_pub()
save_fig(pa, "fig1a_dinucleotide_freq", 3.4, 2.3)
cat(sprintf("  CpG O/E (genome) = %.3f\n", dn[dinucleotide == "CG", ratio_obs_exp]))

# ---- [3] Region x chromosome CpG density / O/E / GC (feeds fig1b/c/d) --------
# Region sets from the GFF (promoter = 2 kb upstream ending AT the TSS; intergenic =
# genome minus gene body + promoter) plus the approved TE table; one row per
# region x chromosome; writes region_chr_cpg_stats.tsv.
cat("[3] per-region CpG density / O/E / GC\n")
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]
seqlengths(genes) <- chr_len[seqlevels(genes)]

prom   <- trim(promoters(genes, upstream = 2000, downstream = 0))  # 2 kb upstream
# strand-neutralise before reduce() so two genes overlapping on OPPOSITE strands merge
# (a stranded reduce keeps both, double-counting shared bases in the length/CpG stats).
gu <- genes; strand(gu) <- "*"; eu <- exons; strand(eu) <- "*"
body   <- reduce(gu)                                              # gene bodies
exon_r <- reduce(eu)
intron <- GenomicRanges::setdiff(body, exon_r)                     # body minus exons (both unstranded)
# Intergenic = genome minus (gene body + promoter)
genic  <- reduce(c(granges(body), granges(prom)))
strand(genic) <- "*"; genic <- reduce(genic)
seqlengths(genic) <- chr_len[seqlevels(genic)]
inter  <- gaps(genic); inter <- inter[strand(inter) == "*"]

# Transposable elements (approved TE table; chr1..31 only, like everything else)
te <- fread(TE)
te <- te[chrom %in% keep_chr]
te_gr <- GRanges(te$chrom, IRanges(te$start, te$end))

regions <- list(Promoter = prom, `Gene body` = body, Exon = exon_r,
                Intron = intron, TE = te_gr, Intergenic = inter)

# For a set of ranges on one chromosome, sum CpG, C, G over the actual sequence.
# reduce() merges overlapping ranges first so shared bases aren't counted twice.
region_chr_stat <- function(reg, chr) {
  r <- reg[as.character(seqnames(reg)) == chr]
  if (length(r) == 0) return(NULL)
  strand(r) <- "*"; r <- reduce(trim(r)); r <- r[width(r) > 0]   # neutralise strand before reduce (no double-count)
  len <- sum(as.numeric(width(r)))
  if (len < 10000) return(NULL)                # skip tiny region/chr combos
  v   <- Views(genome[[chr]], start = start(r), end = end(r))
  ncpg <- sum(vcountPattern("CG", v))
  cc   <- sum(letterFrequency(v, "C")); gg <- sum(letterFrequency(v, "G"))
  data.table(chr = chr,
             cpg_per_kb = 1000 * ncpg / len,
             cpg_oe     = if (cc > 0 && gg > 0) ncpg / (as.numeric(cc) * gg / len) else NA_real_,
             gc_pct     = 100 * (cc + gg) / len)
}

# mito scaffold would enter ONLY the Intergenic series, where its 28.7% GC (chromosomal
# median 44.2%) is a 15-point outlier inflating that series' SD 0.60 -> 2.82. Organelle
# DNA is not "intergenic" in the nuclear sense (same reason the LMR segmentation drops it).
chr_nuc <- paste0("chr", 1:31)
stats <- rbindlist(lapply(names(regions), function(nm) {
  rr <- rbindlist(lapply(chr_nuc, function(chr) region_chr_stat(regions[[nm]], chr)))
  if (nrow(rr)) rr[, region := nm]
  rr
}))
stats <- stats[is.finite(cpg_oe)]   # drop rare non-finite O/E rows (a region/chr with no C or G)
stats[, region := factor(region, levels = names(regions))]
fwrite(stats, file.path(DAT, "region_chr_cpg_stats.tsv"), sep = "\t")

# Genome-wide reference values (dotted lines)
g_cpg_per_kb <- 1000 * sum(as.numeric(di_counts["CG"])) / sum(as.numeric(chr_len))
g_gc_pct     <- 100 * (mono["C"] + mono["G"]) / sum(mono)
g_cpg_oe     <- dn[dinucleotide == "CG", ratio_obs_exp]

# ---- [4] fig1b/c/d: per-region distributions ---------------------------------
# Density curves of §3's per-chromosome stats; dotted vline = genome-wide value.
cat("[4] fig1b/c/d region distributions\n")
pb <- ggplot(stats, aes(cpg_per_kb, fill = region, colour = region)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = g_cpg_per_kb, linetype = "dotted") +
  scale_fill_manual(values = COL_REGION) + scale_colour_manual(values = COL_REGION) +
  labs(x = "CpG per kb", y = "Density", title = "CpG density by region", subtitle = "one point per chromosome (n = 31)") + theme_pub()
save_fig(pb, "fig1b_cpg_density_distribution", 3.4, 2.2)

# fig1c axis decision: all regions sit at CpG O/E ~0.45-0.70, so the x=1 "parity"
# line is deliberately NOT drawn (it left ~40% of the panel empty); the dotted
# genome-wide O/E line (g_cpg_oe) is the anchor instead.
pc <- ggplot(stats, aes(cpg_oe, fill = region, colour = region)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = g_cpg_oe, linetype = "dotted") +
  scale_fill_manual(values = COL_REGION) + scale_colour_manual(values = COL_REGION) +
  labs(x = "CpG observed / expected", y = "Density", title = "CpG O/E by region", subtitle = "one point per chromosome (n = 31)") + theme_pub()
save_fig(pc, "fig1c_cpg_oe_distribution", 3.4, 2.2)

pd <- ggplot(stats, aes(gc_pct, fill = region, colour = region)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = g_gc_pct, linetype = "dotted") +
  scale_fill_manual(values = COL_REGION) + scale_colour_manual(values = COL_REGION) +
  labs(x = "GC content (%)", y = "Density", title = "GC content by region", subtitle = "one point per chromosome (n = 31)") + theme_pub()
save_fig(pd, "fig1d_gc_distribution", 3.4, 2.2)

# ---- [5] fig2a: methylation toolkit presence (eggNOG orthology) --------------
# 20 canonical animal toolkit genes/families; writes toolkit_presence.tsv + fig2a.
# MBD2/3, UHRF1/2, EHMT1/2, EZH1/2, SUV39H1/2, GADD45A/B/G are vertebrate/mammalian
# duplications, so an invertebrate carries ONE ancestral member named by eggNOG after
# whichever paralog it best matches. Per-paralog scoring called G9A/EZH2/SUV39H1
# "absent" while the same loci were present as EHMT1/EZH1/SUV39H2 (a naming artefact,
# not a finding). Each family = ONE row matched against every paralog name (tk_syn);
# the bar carries the name eggNOG assigned. MBD1 and MeCP2 stay as rows because they
# are vertebrate-specific and their expected absence should be visible.
cat("[5] fig2a toolkit presence\n")
toolkit <- data.table(
  gene_symbol = c("DNMT1","DNMT3 (A/B/L)","DNMT2",
                  "TET (1/2/3)","MBD1","MBD2/3","MBD4","MECP2",
                  "UHRF (1/2)","HELLS","DMAP1","PCNA",
                  "EHMT (G9A/GLP)","EZH (1/2)","SETDB1","SUV39H (1/2)",
                  "GADD45 (A/B/G)","AICDA","SMUG1","TDG"),
  category = c("Writer (maintenance)","Writer (de novo)","Writer (tRNA)",
               "Eraser (TET)","Reader (MBD)","Reader (MBD)","Reader (MBD)","Reader (MBD)",
               "Support","Support","Cofactor","Cofactor",
               "Cofactor","Cofactor","Cofactor","Cofactor",
               "Demethylation pathway","Demethylation pathway",
               "Demethylation pathway","Demethylation pathway"))

# assigned that name by eggNOG-mapper (the ortholog database), and the CONFIDENCE
# of the ortholog call is the seed-ortholog bit score. This replaces the older
# GFF-"Similar to" text match with an actual orthology assignment.
EMAPPER <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/eggnog_mapper/dlasi_proteome.emapper.annotations"
gene_id <- sub(";.*", "", as.character(mcols(genes)$ID))          # chr1-31 gene loci
# GFF "Similar to SYM:" note per gene — still used downstream (§7 DE-volcano labels)
gene_note <- sapply(mcols(genes)$Note, function(x)
  if (length(x) == 0) NA_character_ else as.character(x)[1])
egg <- fread(EMAPPER, sep = "\t", quote = "", header = TRUE, skip = "#query",
             na.strings = c("-", "", "NA"), fill = TRUE)
setnames(egg, 1, "query"); egg <- egg[!startsWith(query, "##")]
egg[, locus := sub("-mRNA-.*$", "", query)]
egg <- egg[locus %in% gene_id & !is.na(Preferred_name)]          # chr1-31 only
egg[, PN := toupper(Preferred_name)]
# Every name a family row accepts (eggNOG Preferred_name, case-insensitive). DNMT2 is
# annotated under its HGNC synonym TRDMT1; G9A is HGNC EHMT2 and GLP is EHMT1.
tk_syn <- list(DNMT2 = "TRDMT1",
               `DNMT3 (A/B/L)`  = c("DNMT3A", "DNMT3B", "DNMT3L"),
               `TET (1/2/3)`    = c("TET1", "TET2", "TET3"),
               `MBD2/3`         = c("MBD2", "MBD3"),
               `UHRF (1/2)`     = c("UHRF1", "UHRF2"),
               `EHMT (G9A/GLP)` = c("EHMT1", "EHMT2", "G9A", "GLP"),
               `EZH (1/2)`      = c("EZH1", "EZH2"),
               `SUV39H (1/2)`   = c("SUV39H1", "SUV39H2"),
               `GADD45 (A/B/G)` = c("GADD45A", "GADD45B", "GADD45G"))
lookup <- function(sym) {                                        # best eggNOG hit by bit score
  names_try <- toupper(c(sym, tk_syn[[sym]]))
  h <- egg[PN %in% names_try][order(-score)]
  if (nrow(h)) list(gene_id = h$locus[1], bit = as.numeric(h$score[1]),
                    ortholog = h$Preferred_name[1], seed = h$seed_ortholog[1], n = uniqueN(h$locus))
  else list(gene_id = NA_character_, bit = 0, ortholog = NA_character_, seed = NA_character_, n = 0L)
}
mm <- lapply(toolkit$gene_symbol, lookup)
toolkit[, gene_id    := sapply(mm, `[[`, "gene_id")]
toolkit[, bit        := as.numeric(sapply(mm, `[[`, "bit"))]
toolkit[, ortholog   := sapply(mm, `[[`, "ortholog")]
toolkit[, seed       := sapply(mm, `[[`, "seed")]
toolkit[, n_paralogs := as.integer(sapply(mm, `[[`, "n"))]
# PRESENT = best eggNOG hit clears the mapper's default seed-ortholog acceptance
# threshold (bit >= 60); "absent" = NOT DETECTED above threshold, weaker than "lost"
# DNMT3 absence is additionally confirmed by tBLASTn + HMMER in §5b.
BIT_MIN <- 60
toolkit[, present    := !is.na(gene_id) & bit >= BIT_MIN]
fwrite(toolkit, file.path(DAT, "toolkit_presence.tsv"), sep = "\t")
cat(sprintf("  present %d / %d by eggNOG ortholog (DNMT1=%s, DNMT3A=%s)\n",
            sum(toolkit$present), nrow(toolkit),
            toolkit[gene_symbol == "DNMT1", present],
            toolkit[gene_symbol == "DNMT3 (A/B/L)", present]))

# Category levels: detailed order for the §6 heatmap, broad name for colour.
cat_levels <- c("Writer (maintenance)","Writer (de novo)","Writer (tRNA)",
                "Eraser (TET)","Reader (MBD)","Support","Cofactor",
                "Demethylation pathway")
broad_levels <- c("Writer","Eraser","Reader","Support","Cofactor",
                  "Demethylation pathway")
broad_cat <- function(x) sub(" *\\(.*", "", x)   # "Writer (de novo)" -> "Writer"

# fig2a: ORTHOLOG CONFIDENCE bar chart. Bar length = eggNOG bit score, labelled
# with the ortholog name eggNOG assigned (e.g. the EHMT family row reads "EHMT1");
# families with no ortholog above the cutoff (e.g. DNMT3 (A/B/L)) read "absent".
tk <- copy(toolkit)
tk[, cat_broad := factor(broad_cat(as.character(category)), levels = broad_levels)]
tk[, label := ifelse(present, ortholog, gene_symbol)]            # ortholog name when present
setorder(tk, -present, -bit, gene_symbol)                        # strongest first, absent last
tk[, label := factor(make.unique(label), levels = rev(make.unique(label)))]
COL_TKCAT <- c(Writer = "#1B9E77", Eraser = "#D95F02", Reader = "#7570B3",   # ColorBrewer Dark2
               Support = "#E7298A", Cofactor = "#66A61E", `Demethylation pathway` = "#E6AB02")
pe <- ggplot(tk, aes(bit, label, fill = cat_broad)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = ifelse(present, format(round(bit), big.mark = ","), "absent")),
            hjust = -0.1, size = 2.5, colour = "grey20") +
  facet_grid(cat_broad ~ ., scales = "free_y", space = "free_y", switch = "y") +  # divide by class
  scale_fill_manual(values = COL_TKCAT, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = "eggNOG ortholog bit score", y = NULL,
       title = "DNA methylation toolkit orthologs") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7), panel.grid.major.y = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, colour = "grey25"),
        plot.title = element_text(size = rel(1)), plot.title.position = "plot")
        # "panel" anchoring pushes the title right of the long y-label block
save_fig(pe, "fig2a_methylation_toolkit_presence", 4.4, 5.1)

# ---- [5b] figS_dnmt3_absence: tBLASTn + Pfam domain architecture -------------
# Hardens §5's "DNMT3 absent" (annotation-derived, not proven loss) two ways:
#   (A) tBLASTn of DNMT3 proteins vs the 6-frame genome (annotation-INDEPENDENT),
#       with DNMT1/DNMT2/TET3 as positive controls proving the search sensitivity.
#   (B) HMMER/Pfam: enumerate every catalytic C5-MTase (PF00145) in the proteome and
#       search the DNMT3-specific ADD domain (PF17980) genome-wide (expect 0).
# (shared with the unrelated PWWP2A gene); a naive best-hit reads "DNMT3 present".
# only small result tables come back. Writes dnmt3_*.tsv + figS_dnmt3_absence.
cat("[5b] figS_dnmt3_absence: tBLASTn + Pfam domain architecture\n")
suppressPackageStartupMessages(library(patchwork))
PROJ      <- "/mnt/data/alfredvar/rlopezt/meth_paper"
BLAST_BIN <- "/opt/apps/blast+/2.13.0/bin"
HMMER_BIN <- file.path(PROJ, "tools/hmmer/bin")
PFAM      <- file.path(PROJ, "tools/pfam/Pfam-A.hmm")           # gunzipped + hmmpress'd
DNMT_Q    <- file.path(PROJ, "tools/dnmt/dnmt_queries.fasta")   # staged query proteins
PROTE     <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.proteins.fasta"

job <- Sys.getenv("SLURM_JOB_ID", "")     # heavy I/O off /mnt/data, into /scratch
scr <- if (nzchar(job)) file.path("/scratch/groups/alfredvar", Sys.getenv("USER"),
                                  paste0("job_", job), "dnmt3") else file.path(tempdir(), "dnmt3")
dir.create(scr, recursive = TRUE, showWarnings = FALSE)

# (A) tBLASTn: write the chr1..31+mito genome, build the nucleotide DB, search all queries.
gfa <- file.path(scr, "genome_chrmt.fa"); writeXStringSet(genome, gfa)
system2(file.path(BLAST_BIN, "makeblastdb"),
        c("-in", gfa, "-dbtype", "nucl", "-out", file.path(scr, "db")), stdout = FALSE)
bt <- file.path(scr, "tblastn.tsv")
system2(file.path(BLAST_BIN, "tblastn"),
        c("-query", DNMT_Q, "-db", file.path(scr, "db"), "-evalue", "10",
          "-num_threads", "8", "-max_target_seqs", "20", "-seg", "yes",
          # must be one shell token or tblastn errors out and writes nothing.
          "-outfmt", shQuote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"),
          "-out", bt))
bl <- fread(bt, header = FALSE,
            col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen",
                          "qstart","qend","sstart","send","evalue","bitscore","qcovs"))
best <- bl[order(evalue)][, .SD[1], by = qseqid]      # best hit per query

# Which gene does the best DNMT3 hit fall in? (expect PWWP2A, NOT a methyltransferase)
d3 <- bl[grepl("DNMT3", qseqid) & evalue < 1e-3][order(evalue)][1]
genes <- gff[gff$type == "gene"]
hit_gr   <- GRanges(d3$sseqid, IRanges(min(d3$sstart, d3$send), max(d3$sstart, d3$send)))
hit_gene <- genes[subjectHits(findOverlaps(hit_gr, genes))[1]]
sym_of_gene <- function(g) { s <- sub(".*Similar to ([^:]+):.*", "\\1", as.character(g$Note)[1])
                        if (is.na(s) || !nzchar(s)) as.character(g$ID)[1] else s }
hit_loc <- as.character(hit_gene$ID)[1]; hit_sym <- sym_of_gene(hit_gene)

# (B) HMMER/Pfam. Fixed leading columns of the tables are parsed in pure R (the
# trailing description column contains spaces). hmmfetch needs an SSI index; tools/pfam
# is read-only in practice, so it is built only if BOTH possible index names are absent.
if (!file.exists(paste0(PFAM, ".ssi")) && !file.exists(paste0(PFAM, ".h3m.ssi")))
  system2(file.path(HMMER_BIN, "hmmfetch"), c("--index", PFAM), stdout = FALSE)
qcat <- file.path(scr, "q_cat.hmm"); qadd <- file.path(scr, "q_add.hmm")
system2(file.path(HMMER_BIN, "hmmfetch"), c(PFAM, "DNA_methylase"), stdout = qcat)
system2(file.path(HMMER_BIN, "hmmfetch"), c(PFAM, "ADD_DNMT3"),     stdout = qadd)
cat_tbl <- file.path(scr, "cat.tbl"); add_tbl <- file.path(scr, "add.tbl")
system2(file.path(HMMER_BIN, "hmmsearch"), c("--cut_ga", "--tblout", cat_tbl, qcat, PROTE), stdout = FALSE)
system2(file.path(HMMER_BIN, "hmmsearch"), c("--cut_ga", "--tblout", add_tbl, qadd, PROTE), stdout = FALSE)
tbl_targets <- function(f) { l <- grep("^#", readLines(f), invert = TRUE, value = TRUE)
  if (!length(l)) character(0) else sub("\\s.*", "", trimws(l)) }
cat_ids <- unique(sub("-mRNA.*", "", tbl_targets(cat_tbl)))   # genes with a catalytic domain
n_add   <- length(tbl_targets(add_tbl))                       # DNMT3 ADD hits (expect 0)

# Domain cartoon set: the catalytic MTases + the PWWP2A hit + human DNMT3A reference.
prot <- readAAStringSet(PROTE); names(prot) <- sub("\\s.*", "", names(prot))
pick <- unique(c(cat_ids, hit_loc))
sel  <- prot[names(prot) %in% paste0(pick, "-mRNA-1")]   # assumes isoform -mRNA-1 exists; a gene lacking it drops from the cartoon silently
ref  <- readAAStringSet(DNMT_Q); ref <- ref[grepl("DNMT3A", names(ref))]; names(ref) <- "Human_DNMT3A"
sfa  <- file.path(scr, "scan.faa"); writeXStringSet(c(sel, ref), sfa)
dtbl <- file.path(scr, "dom.domtbl")
system2(file.path(HMMER_BIN, "hmmscan"), c("--domtblout", dtbl, "--cut_ga", PFAM, sfa), stdout = FALSE)
dl <- grep("^#", readLines(dtbl), invert = TRUE, value = TRUE)
pr <- strsplit(trimws(dl), "\\s+")
dom <- data.table(protein = vapply(pr, `[`, "", 4), domain = vapply(pr, `[`, "", 1),
                  start = as.integer(vapply(pr, `[`, "", 20)), end = as.integer(vapply(pr, `[`, "", 21)),
                  plen  = as.integer(vapply(pr, `[`, "", 6)))
dom <- dom[grepl("-mRNA-1$|Human_DNMT3A", protein)]

# friendly labels: symbol from the GFF Note, PWWP2A hit flagged, human ref
lab_of <- function(p) {
  if (p == "Human_DNMT3A") return("Human DNMT3A (reference)")
  loc <- sub("-mRNA.*", "", p); s <- sym_of_gene(genes[genes$ID == loc])
  if (loc == hit_loc) sprintf("D. laeve %s\n(the BLAST hit)", s) else sprintf("D. laeve %s", s)
}
dom[, plab := vapply(protein, lab_of, "")]
# order the cartoon rows: human DNMT3A reference on top, slug proteins below
dom[, plab := factor(plab, levels = rev(c("Human DNMT3A (reference)",
                     sort(setdiff(unique(plab), "Human DNMT3A (reference)")))))]

fwrite(best, file.path(DAT, "dnmt3_tblastn_besthit.tsv"), sep = "\t")
fwrite(dom,  file.path(DAT, "dnmt3_domain_architecture.tsv"), sep = "\t")
cat(sprintf("  catalytic MTases: %s | ADD_DNMT3 hits: %d | best DNMT3 hit -> %s (%s)\n",
            paste(cat_ids, collapse = ", "), n_add, hit_loc, hit_sym))

# ---- panel A: tBLASTn best hit per query vs the genome ----
COL_PA <- c("control (known present)" = "#0072B2", "DNMT3 query" = "#D55E00")
a <- best[, .(query = qseqid, evalue, bitscore)]
a[, grp := ifelse(grepl("DNMT3", query), "DNMT3 query", "control (known present)")]
a[, neglogE := -log10(pmax(evalue, 1e-200))]
# annotate each row with the gene its best hit lands in
best_gr <- GRanges(best$sseqid, IRanges(pmin(best$sstart, best$send), pmax(best$sstart, best$send)))
hg <- findOverlaps(best_gr, genes)
a[, target := "no significant hit"]
a[queryHits(hg), target := sprintf("%s: %s", best$sseqid[queryHits(hg)],
                                   vapply(subjectHits(hg), function(i) sym_of_gene(genes[i]), ""))]
a[evalue > 1e-3, target := "no significant hit"]
lbl <- c(DNMT1_human="DNMT1 (human)", DNMT2_human="DNMT2 (human)", TET3_human="TET3 (human)",
         DNMT3A_human="DNMT3A (human)", DNMT3B_human="DNMT3B (human)",
         DNMT3_insect="DNMT3 (insect)", DNMT3_mollusc="DNMT3 (mollusc)")
a[, lab := ifelse(query %in% names(lbl), lbl[query], query)]
a <- a[order(grp, neglogE)]; a[, lab := factor(lab, levels = lab)]
pA <- ggplot(a, aes(neglogE, lab, colour = grp)) +
  geom_segment(aes(x = 0, xend = neglogE, yend = lab), linewidth = 0.5) +
  geom_point(size = 2.4) +
  geom_vline(xintercept = -log10(0.001), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_text(aes(label = target), hjust = -0.1, size = 2.4, colour = "grey20") +
  scale_colour_manual(values = COL_PA, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = expression(-log[10]~"(best tBLASTn E-value vs genome)"), y = NULL,
       title = "A  DNMT3 queries hit only PWWP2A, not a methyltransferase") +
  theme_pub() + theme(legend.position = "top", legend.justification = "left")

# ---- panel B: domain architecture (DNMT3 = PWWP + ADD + catalytic) ----
bb <- unique(dom[, .(plab, plen)])
dom[, dclass := fifelse(domain == "DNA_methylase", "catalytic (C5-MTase)",
              fifelse(domain %in% c("ADD_DNMT3","DNMT3_ADD_GATA1-like"), "ADD (DNMT3-specific)",
              fifelse(domain == "PWWP", "PWWP", "other domain")))]
dom[, dclass := factor(dclass, levels = c("catalytic (C5-MTase)","ADD (DNMT3-specific)","PWWP","other domain"))]
COL_PB <- c("catalytic (C5-MTase)" = "#D55E00", "ADD (DNMT3-specific)" = "#009E73",
            "PWWP" = "#56B4E9", "other domain" = "#BBBBBB")
pB <- ggplot() +
  geom_segment(data = bb, aes(x = 1, xend = plen, y = plab, yend = plab), colour = "grey65", linewidth = 0.6) +
  geom_rect(data = dom, aes(xmin = start, xmax = end, ymin = as.numeric(plab) - 0.3,
            ymax = as.numeric(plab) + 0.3, fill = dclass), colour = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = COL_PB, name = NULL) +
  # main.tex caption, not burned into the plot.
  labs(x = "amino-acid position", y = NULL,
       title = "B  No D. laeve protein has the DNMT3 architecture") +
  theme_pub() + theme(legend.position = "top", legend.justification = "left")

save_supp(pA / pB + patchwork::plot_layout(heights = c(1, 0.95)), "figS_dnmt3_absence", 8.5, 7.2)
unlink(scr, recursive = TRUE)     # drop the heavy BLAST DB from scratch
rm(bl, best, prot, sel, ref, dom, a); gc(verbose = FALSE)

# ---- [6] figS_toolkit_mrna_tail: toolkit mRNA, tail control vs amputated -----
# fig2a is the main toolkit panel.) Also builds the dds reused by §7 for gene DE.
# Sample map rebuilt from filenames (NO external metadata):
# tail control = C1S1..C4S4 ; tail amputated = T2S6,T3S7,T4S8 (T1S5 excluded).
cat("[6] fig2c toolkit mRNA, tail control vs amputated\n")
samples <- data.table(
  sample    = c("C1S1","C2S2","C3S3","C4S4","T2S6","T3S7","T4S8"),
  condition = c(rep("Control", 4), rep("Amputated", 3)))
samples[, file := file.path(HTSEQ, paste0(sample, "_htseq_gene_counts.txt"))]
stopifnot(all(file.exists(samples$file)))

# Read each HTSeq file (gene_id, count); drop the trailing "__" summary rows.
counts_list <- lapply(samples$file, function(f) {
  x <- fread(f, header = FALSE, col.names = c("gene_id", "count"))
  x[!startsWith(gene_id, "__")]
})
genes_in_counts <- counts_list[[1]]$gene_id
count_mat <- sapply(counts_list, function(x) x$count)
rownames(count_mat) <- genes_in_counts
colnames(count_mat) <- samples$sample

# Keep only genes annotated on the keep_chr universe (gene_id from §5) — the
# transcriptomic side of the chromosome filter.
count_mat <- count_mat[rownames(count_mat) %in% gene_id, , drop = FALSE]
count_mat <- count_mat[rowSums(count_mat >= 5) >= 2, , drop = FALSE]  # expressed
cat(sprintf("  %d genes x %d tail libraries after filtering\n",
            nrow(count_mat), ncol(count_mat)))

coldata <- data.frame(condition = factor(samples$condition,
                                          levels = c("Control","Amputated")),
                      row.names = samples$sample)
dds <- DESeqDataSetFromMatrix(count_mat, coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
norm_mat <- counts(dds, normalized = TRUE)   # DESeq2 library-size normalized counts (not VST)

tk_present <- toolkit[present == TRUE & gene_id %in% rownames(norm_mat)]   # present AND expressed: a toolkit gene failing the expression filter is not drawn
mean_ctrl <- rowMeans(norm_mat[tk_present$gene_id, samples[condition == "Control", sample], drop = FALSE])
mean_amp  <- rowMeans(norm_mat[tk_present$gene_id, samples[condition == "Amputated", sample], drop = FALSE])

mrna <- data.table(gene_symbol = tk_present$ortholog,           # the name eggNOG assigned (family rows carry e.g. EHMT1)
                   category    = tk_present$category,
                   Control = mean_ctrl, Amputated = mean_amp)
fwrite(mrna, file.path(DAT, "toolkit_mrna_tail.tsv"), sep = "\t")

mrna[, category  := factor(category, levels = cat_levels)]
mrna[, cat_broad := factor(broad_cat(as.character(category)), levels = broad_levels)]
setorder(mrna, category, gene_symbol)
mrna[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
long <- data.table::melt(mrna, id.vars = c("gene_symbol","cat_broad"),   # namespaced (masked-generic rule)
             measure.vars = c("Control","Amputated"),
             variable.name = "condition", value.name = "expr")

pf <- ggplot(long, aes(condition, gene_symbol, fill = expr)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Norm.\ncounts") +
  facet_grid(cat_broad ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "Methylation toolkit - Tail blastema mRNA levels") +
  theme_pub() +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, colour = "grey25"),
        panel.grid = element_blank())
save_supp(pf, "figS_toolkit_mrna_tail", 5.0, 4.5)   # supplementary (fig2a is the main toolkit panel)

# ---- [7] Supplementary gene DE (tail): DESeq2/apeglm volcano + top-20 --------
# HTSeq gene counts -> DESeq2/apeglm -> gene_de_tail.tsv (the project's canonical DE
# table), figS_gene_de_volcano, figS_top_de_genes.
cat("[7] supplementary: gene DE\n")
suppressPackageStartupMessages({ library(EnhancedVolcano) })

# Gene-symbol map from the GFF Note (for readable volcano / table labels)
has_sym <- grepl("^Similar to [^:]+:", gene_note)
sym_of  <- ifelse(has_sym, sub("^Similar to ([^:]+):.*$", "\\1", gene_note), NA_character_)
names(sym_of) <- gene_id
lab_for <- function(ids) ifelse(is.na(sym_of[ids]), ids, sym_of[ids])   # tables: symbol or LOC id
# "Non annotated gene" placeholder is retired.
disp_name <- function(ids) unname(ifelse(is.na(sym_of[ids]), ids, sym_of[ids]))

# ---- [7a] Gene-level DE (tail HTSeq, 4 control vs 3 amputated) ---------------
# Reuses the dds built in §6 (design ~ condition).
dds <- DESeq(dds)
# STAT TEST: DESeq2 differential expression = Wald test on the per-gene negative-binomial
# GLM (p-values/padj from that test, BH-adjusted); apeglm shrinkage applied to the LFC (Zhu et al.).
res_g <- lfcShrink(dds, coef = "condition_Amputated_vs_Control", type = "apeglm")
deg <- as.data.table(as.data.frame(res_g), keep.rownames = "gene_id")
deg[, symbol := lab_for(gene_id)]
setorder(deg, padj, na.last = TRUE)
fwrite(deg, file.path(DAT, "gene_de_tail.tsv"), sep = "\t")
nsig_g <- deg[!is.na(padj) & padj < 0.05, .N]
cat(sprintf("  gene DE: %d genes FDR<0.05\n", nsig_g))

# DE volcano — EnhancedVolcano. Labels are FORCED to the top 5 up + top 5 down
# NAMED genes by FDR (selectLab), so the down side is never silently dropped.
# Cool "winter" palette (blue-green).
de_df  <- as.data.frame(res_g); de_df <- de_df[!is.na(de_df$padj), ]
de_dt  <- as.data.table(de_df, keep.rownames = "gene_id")[, symbol := disp_name(gene_id)]
# gets drawn on all of them, turning 10 intended labels into 15. Repeated symbols are
# disambiguated with their LOC id. (Same duplicate-symbol trap as 03_promoters fig3c.)
de_dt[, lab_uniq := fifelse(duplicated(symbol) | duplicated(symbol, fromLast = TRUE),
                            paste0(symbol, " (", gene_id, ")"), symbol)]
de_lab <- de_dt$lab_uniq                                 # row-aligned with de_df
n_up   <- de_dt[padj < 0.05 & log2FoldChange >=  1, .N]
n_down <- de_dt[padj < 0.05 & log2FoldChange <= -1, .N]
sel_de <- c(de_dt[padj < 0.05 & log2FoldChange >=  1][order(padj)][seq_len(min(5, .N)), lab_uniq],
            de_dt[padj < 0.05 & log2FoldChange <= -1][order(padj)][seq_len(min(5, .N)), lab_uniq])
stopifnot(!anyDuplicated(sel_de))
pvg <- EnhancedVolcano(de_df, lab = de_lab, x = "log2FoldChange", y = "padj",
  # the y column is ADJUSTED P, so the default "-Log10 P" axis would misstate it
  xlab = bquote(Log[2]~"fold change (regenerated vs control)"),
  ylab = bquote(-Log[10]~"adjusted"~italic(P)),
  selectLab = sel_de, pCutoff = 0.05, FCcutoff = 1,
  col = c("grey80", "#66C2A5", "#3288BD", "#5E4FA2"), colAlpha = 0.65,
  pointSize = 1.4, labSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.3,
  maxoverlapsConnectors = Inf, boxedLabels = FALSE,
  legendLabels = c("NS", "log2FC", "FDR", "FDR & log2FC"), legendPosition = "top",
  title = "Tail regeneration differential expression", caption = NULL,
  subtitle = sprintf("FDR < 0.05 & |log2FC| >= 1: %d up, %d down", n_up, n_down))
save_supp(pvg, "figS_gene_de_volcano", 7.0, 6.5)

top_g <- head(deg[!is.na(padj)][order(padj)], 20)       # 20 most significant genes
setorder(top_g, log2FoldChange)                          # order bars high -> low (up top, down bottom)
top_g[, dir := ifelse(log2FoldChange > 0, "Up", "Down")]
top_g[, dlab := disp_name(gene_id)]                      # symbol, or LOC id when unnamed
top_g[, dlab := factor(make.unique(dlab), levels = make.unique(dlab))]
ptg <- ggplot(top_g, aes(log2FoldChange, dlab, fill = dir)) +
  geom_col() +
  scale_fill_manual(values = c(Up = unname(COL_COND["Amputated"]),
                               Down = unname(COL_COND["Control"])), name = NULL) +
  labs(x = "log2 fold-change (Amp / Ctrl)", y = NULL,
       title = "Top 20 DE genes tail blastema") + theme_pub() +
  theme(axis.text.y = element_text(size = 7))
save_supp(ptg, "figS_top_de_genes", 5, 5)

# NOTE: an isoform-level transcript analysis that used to live here (formerly §7c)
# objects/salmon/ index+quants it produced are now unused.

# ---- [5c] figS_uhrf_domains: can DNMT1 maintain methylation with no UHRF1? ---
# (Placed after §7 by file history; reuses §5b's PROJ/PFAM/HMMER_BIN paths.)
# eggNOG calls UHRF1 absent but UHRF2 PRESENT (LOC_00011698). The UHRF1/2 pair is a
# VERTEBRATE duplication, so "UHRF1 lost" is a paralog-level claim the data cannot
# support: D. laeve has ONE UHRF-family gene, and function follows domains, not names.
# Asks whether that UHRF keeps SRA (PF02182) + RING (PF13445/PF13639/PF00097) and
# whether DNMT1 keeps RFTS (PF12047), against human UHRF1/UHRF2 references.
cat("[5c] figS_uhrf_domains: UHRF/DNMT1 domain architecture\n")
UHRF_REF <- file.path(PROJ, "tools/dnmt/uhrf_queries.fasta")
PROTEOME <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.proteins.fasta"
if (file.exists(UHRF_REF) && file.exists(PROTEOME) && dir.exists(HMMER_BIN)) {
  prot <- readAAStringSet(PROTEOME)
  loc_of <- sub("-mRNA.*", "", sub(" .*", "", names(prot)))
  pick_longest <- function(loc, lab) {            # one representative protein per gene
    i <- which(loc_of == loc); if (!length(i)) return(NULL)
    s <- prot[i[which.max(width(prot[i]))]]; names(s) <- lab; s
  }
  uhrf_locus <- toolkit[gene_symbol == "UHRF (1/2)", gene_id]; dnmt1_locus <- toolkit[gene_symbol == "DNMT1", gene_id]
  stopifnot(length(uhrf_locus) == 1L, nzchar(uhrf_locus), length(dnmt1_locus) == 1L, nzchar(dnmt1_locus))   # from §5's eggNOG assignment, never hardcoded
  qry <- c(pick_longest(uhrf_locus, "UHRF2_D.laeve"),
           pick_longest(dnmt1_locus, "DNMT1_D.laeve"),
           readAAStringSet(UHRF_REF))
  scr2 <- file.path(tempdir(), "uhrf"); dir.create(scr2, showWarnings = FALSE)
  ufa <- file.path(scr2, "uhrf.faa"); writeXStringSet(qry, ufa)
  udt <- file.path(scr2, "uhrf.domtbl")
  system2(file.path(HMMER_BIN, "hmmscan"), c("--domtblout", udt, "--cut_ga", PFAM, ufa), stdout = FALSE)
  ul <- grep("^#", readLines(udt), value = TRUE, invert = TRUE)
  up <- strsplit(trimws(ul), "\\s+")
  udom <- data.table(protein = vapply(up, `[`, "", 4), domain = vapply(up, `[`, "", 1),
                     acc = sub("\\..*", "", vapply(up, `[`, "", 2)),
                     start = as.integer(vapply(up, `[`, "", 18)),
                     end   = as.integer(vapply(up, `[`, "", 19)))
  plen <- data.table(protein = names(qry), len = width(qry))
  # keep the domains that carry the maintenance-methylation logic
  KEEP <- c(PF00240 = "UBL", PF12148 = "TTD", PF00628 = "PHD", PF02182 = "SRA",
            PF13445 = "RING", PF13639 = "RING", PF00097 = "RING", PF12047 = "RFTS",
            PF02008 = "CXXC", PF01426 = "BAH", PF00145 = "MTase")
  udom <- udom[acc %in% names(KEEP)][, dom := KEEP[acc]]
  udom <- unique(udom[, .(protein, dom, start, end)])
  fwrite(udom, file.path(DAT, "uhrf_dnmt1_domain_architecture.tsv"), sep = "\t")
  ord <- c("UHRF1_human", "UHRF2_human", "UHRF2_D.laeve", "DNMT1_D.laeve")
  plen[, protein := factor(protein, levels = rev(ord))]; udom[, protein := factor(protein, levels = rev(ord))]
  COL_DOM <- c(UBL = "#999999", TTD = "#56B4E9", PHD = "#009E73", SRA = "#D55E00",
               RING = "#CC79A7", RFTS = "#0072B2", CXXC = "#E69F00", BAH = "#7FB3D5", MTase = "#117733")
  pu <- ggplot() +
    geom_segment(data = plen, aes(x = 1, xend = len, y = protein, yend = protein),
                 colour = "grey75", linewidth = 0.7) +
    geom_rect(data = udom, aes(xmin = start, xmax = end,
                               ymin = as.numeric(protein) - 0.28, ymax = as.numeric(protein) + 0.28,
                               fill = dom), colour = "grey25", linewidth = 0.2) +
    geom_text(data = udom, aes(x = (start + end)/2, y = as.numeric(protein) + 0.42, label = dom),
              size = 2.2, colour = "grey20") +
    scale_fill_manual(values = COL_DOM, name = NULL) +
    scale_y_discrete(limits = rev(ord)) +
    # RING/H3-ubiquitination step this scan does not recover — it overstated the data.
    labs(x = "Amino-acid position", y = NULL,
         title = "Maintenance-methylation domain architecture") +
    theme_pub() + theme(legend.position = "bottom")
  save_supp(pu, "figS_uhrf_domains", 7.5, 3.6)
  for (p in ord) cat(sprintf("  %-16s domains: %s\n", p,
    paste(sort(unique(udom[protein == p]$dom)), collapse = ", ")))
} else cat("  SKIPPED - UHRF refs / proteome / HMMER not available\n")

# ---- [8] fig1e: cross-species mollusc CpG O/E (log tag "[N]") ----------------
# Germline deamination (CpG -> TpG/CpA) depletes CpG in methylated genomes; computes
# genome-wide dinucleotide O/E for D. laeve + 4 molluscs (genomes pre-staged in
# dataset/mollusc_genomes/, NCBI accessions below). Writes mollusc_dinucleotide_oe.tsv,
# fig1e and one figS_dinuc_<species> per species.
# not somatic methylation, which was measured only in D. laeve.
cat("[N] cross-species mollusc CpG O/E (deamination signature)\n")
MOLLDIR <- file.path(BATCH, "dataset", "mollusc_genomes")
molluscs <- data.table(
  species = c("Deroceras laeve", "Pomacea canaliculata", "Elysia atroviridis",
              "Aplysia californica", "Octopus bimaculoides"),
  group   = c("Gastropoda", "Gastropoda", "Gastropoda", "Gastropoda", "Cephalopoda"),
  # are whole NCBI assemblies (unplaced, repeat-rich scaffolds included) and repeat
  # content shifts CpG O/E, so a chromosome-only subset would not be like with like.
  fa = c(GENOME_FULL,                                                                  # ours (whole assembly)
         file.path(MOLLDIR, "GCF_003073045.1_ASM307304v1_genomic.fna.gz"),            # Pomacea
         file.path(MOLLDIR, "GCA_059052615.1_ASM5905261v1_genomic.fna.gz"),           # Elysia atroviridis
         file.path(MOLLDIR, "GCF_000002075.1_AplCal3.0_genomic.fna.gz"),              # Aplysia
         file.path(MOLLDIR, "GCF_001194135.2_ASM119413v2_genomic.fna.gz")))           # Octopus bimaculoides
COL_GROUP <- c(Gastropoda = "#0072B2", Cephalopoda = "#D55E00")
# genome-wide dinucleotide O/E for one FASTA (one genome loaded at a time, then freed)
dinuc_one <- function(fa) {
  g   <- readDNAStringSet(fa)
  di  <- colSums(dinucleotideFrequency(g)); obs <- di / sum(di)
  mono <- colSums(alphabetFrequency(g, baseOnly = TRUE))[c("A","C","G","T")]; mf <- mono / sum(mono)
  exp <- vapply(names(obs), function(d){ b <- strsplit(d, "")[[1]]; unname(mf[b[1]] * mf[b[2]]) }, numeric(1))
  rm(g); invisible(gc())
  data.table(dinucleotide = names(obs), observed = as.numeric(obs),
             expected = as.numeric(exp),
             gc = unname(mf["C"] + mf["G"]))[, ratio_obs_exp := observed / expected][]
}
oe_all <- rbindlist(lapply(seq_len(nrow(molluscs)), function(i) {
  if (!file.exists(molluscs$fa[i])) { cat(sprintf("  MISSING genome, skipped: %s\n", molluscs$fa[i])); return(NULL) }
  dt <- dinuc_one(molluscs$fa[i]); dt[, `:=`(species = molluscs$species[i], group = molluscs$group[i])]
  # per-species dinucleotide-frequency SUPPLEMENTARY figure (fig1a style: bars = observed, diamonds = expected)
  d <- copy(dt); d[, dinucleotide := factor(dinucleotide, levels = dinucleotide[order(-observed)])]
  d[, is_cg := ifelse(dinucleotide == "CG", "CG", "other")]
  ps <- ggplot(d, aes(dinucleotide, 100 * observed, fill = is_cg)) +
    geom_col(width = 0.75) + geom_point(aes(y = 100 * expected), shape = 18, size = 2.2, colour = "black") +
    scale_fill_manual(values = c(CG = "#D55E00", other = "#9FB1BC"), guide = "none") +
    labs(x = "Dinucleotide", y = "Frequency (%)",
         title = bquote(italic(.(molluscs$species[i]))~"dinucleotide composition")) + theme_pub()
  save_supp(ps, sprintf("figS_dinuc_%s", gsub(" ", "_", tolower(molluscs$species[i]))), 4.5, 3.0)
  cat(sprintf("  %-22s CpG O/E = %.3f\n", molluscs$species[i], dt[dinucleotide == "CG", ratio_obs_exp]))
  dt
}))
fwrite(oe_all, file.path(DAT, "mollusc_dinucleotide_oe.tsv"), sep = "\t")
# COMBINED figure: genome-wide CpG O/E across all five molluscs. Every species below O/E = 1
# is the shared deamination signature; colour separates gastropods from the cephalopod.
cg <- oe_all[dinucleotide == "CG"][order(ratio_obs_exp)]
# Abbreviate the genus on the axis ("Deroceras laeve" -> "D. laeve"). At the canvas width
# this panel is actually rendered at, the full binomials run off the left edge.
cg[, label := sub("^([A-Z])[a-z]+ ", "\\1. ", species)]
cg[, label := factor(label, levels = label)]
cg[, species := factor(species, levels = species)]
pe <- ggplot(cg, aes(label, ratio_obs_exp, fill = group)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", ratio_obs_exp)), vjust = -0.4, size = 3) +
  # GC under each bar: CpG O/E covaries with G+C content as a property of the measure itself
  # (Duret & Galtier 2000), so the reader must be able to see the GC of every bar being compared.
  geom_text(aes(y = 0.04, label = sprintf("%.0f%%", 100 * gc)), size = 2.0, colour = "white") +
  scale_fill_manual(values = COL_GROUP, name = NULL) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(x = NULL, y = "Genome-wide CpG O/E",
       title = "CpG depletion across molluscs") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "italic"))
save_fig(pe, "fig1e_mollusc_cpg_oe", 3.4, 2.1)

# ---- [8b] figS_toolkit_mrna_atlas: toolkit expression across the tissue atlas ---
# The first question after "no DNMT3, yet gains" is whether the writer is induced. [6]
# answers it for the tail contrast; this panel shows every present toolkit gene across
# the 40-library atlas (44 minus the four outliers of 07_wgcna: T1S5, dcrep4, R6, irrep7)
# as the mean DESeq2 normalised count per experiment group. Groups are never pooled
# across experiments (label rule): each column is one group and each experiment keeps
# its own control. Decode = the lab's sample sheet (META, read-only).
cat("[8b] toolkit mRNA across the 40-library atlas\n")
META <- "/mnt/data/alfredvar/rlopezt/WGCNA/metadata.tsv"
meta_atlas <- fread(META)[!orig_sample %in% c("T1S5", "dcrep4", "R6", "irrep7")]
ht_files <- list.files(HTSEQ, pattern = "_htseq_gene_counts\\.txt$")
ht_base  <- sub("\\.Aligned\\.out\\.bam_htseq_gene_counts\\.txt$|_htseq_gene_counts\\.txt$", "", ht_files)
meta_atlas[, path := file.path(HTSEQ, ht_files[match(orig_sample, ht_base)])]
if (anyNA(meta_atlas$path)) {
  cat(sprintf("  WARNING: %d atlas libraries have no HTSeq file and are dropped: %s\n",
              sum(is.na(meta_atlas$path)), paste(meta_atlas[is.na(path), orig_sample], collapse = ", ")))
  meta_atlas <- meta_atlas[!is.na(path)]
}
cl_atlas <- lapply(meta_atlas$path, function(f)
  fread(f, header = FALSE, col.names = c("gene_id", "count"))[!startsWith(gene_id, "__")])
stopifnot(all(vapply(cl_atlas, function(x) identical(x$gene_id, cl_atlas[[1]]$gene_id), logical(1))))
cm_atlas <- sapply(cl_atlas, function(x) x$count)
rownames(cm_atlas) <- cl_atlas[[1]]$gene_id; colnames(cm_atlas) <- meta_atlas$sample
cm_atlas <- cm_atlas[rownames(cm_atlas) %in% gene_id, , drop = FALSE]       # chromosome universe
cm_atlas <- cm_atlas[rowSums(cm_atlas >= 5) >= 2, , drop = FALSE]           # expressed
dds_atlas <- estimateSizeFactors(DESeqDataSetFromMatrix(
  cm_atlas, data.frame(group = meta_atlas$group, row.names = meta_atlas$sample), design = ~ 1))
norm_atlas <- counts(dds_atlas, normalized = TRUE)
cat(sprintf("  %d genes x %d atlas libraries, %d experiment groups\n",
            nrow(norm_atlas), ncol(norm_atlas), uniqueN(meta_atlas$group)))
tk_atlas <- toolkit[present == TRUE & gene_id %in% rownames(norm_atlas)]
grp_levels <- c("TailControl", "TailAmputated", "EyeControl", "EyeAmputated",
                "IrradiatedControl", "IrradiatedTreated", "FungicideControl", "FungicideTreated",
                "Head", "Juvenile", "Ovotestis")
atlas <- rbindlist(lapply(seq_len(nrow(tk_atlas)), function(i)
  data.table(gene_symbol = tk_atlas$ortholog[i], category = tk_atlas$category[i],
             gene_id = tk_atlas$gene_id[i], group = meta_atlas$group, sample = meta_atlas$sample,
             norm_count = as.numeric(norm_atlas[tk_atlas$gene_id[i], ]))))
fwrite(atlas, file.path(DAT, "toolkit_mrna_atlas_per_library.tsv"), sep = "\t")
atlas_mean <- atlas[, .(mean_norm = mean(norm_count), sd_norm = sd(norm_count), n = .N),
                    by = .(gene_symbol, category, gene_id, group)]
fwrite(atlas_mean, file.path(DAT, "toolkit_mrna_atlas.tsv"), sep = "\t")
atlas_mean[, `:=`(group = factor(group, levels = intersect(grp_levels, unique(group))),
                  cat_broad = factor(broad_cat(category), levels = broad_levels))]
setorder(atlas_mean, cat_broad, gene_symbol)
atlas_mean[, gene_symbol := factor(gene_symbol, levels = rev(unique(gene_symbol)))]
pa8 <- ggplot(atlas_mean, aes(group, gene_symbol, fill = log10(mean_norm + 1))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "log10 mean\nnorm. counts") +
  facet_grid(cat_broad ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "Methylation toolkit across the tissue atlas (40 libraries)") +
  theme_pub() +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, colour = "grey25"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7), panel.grid.major.y = element_blank())
save_supp(pa8, "figS_toolkit_mrna_atlas", 6.4, 5.4)

# ---- [8c] figS_dnmt_c5_tree: C5-cytosine methyltransferase gene tree ---------------
# DNMT3 absence should be phylogenetic, not a BLAST miss. The DNA_methylase (Pfam
# PF00145) domains of the D. laeve DNMT1 and DNMT2/TRDMT1 loci are placed among the
# reviewed vertebrate DNMT1, DNMT3A, DNMT3B, DNMT3L and TRDMT1 proteins and every UniProt
# "DNA (cytosine-5)-methyltransferase" of Mollusca (dataset/dnmt_refs, staged once from
# node12 by fetch_dnmt_refs.sh), with the bacterial M.HhaI as outgroup. Domains are cut to
# the Pfam envelope (hmmsearch --cut_ga), aligned (MAFFT --auto), trimmed (trimAl
# -automated1) and a maximum likelihood tree is built (IQ-TREE 3, ModelFinder, 1,000
# clade would contradict the toolkit call; none is expected.
cat("[8c] C5-methyltransferase gene tree\n")
suppressPackageStartupMessages({ library(ape); library(phangorn) })
DNMT_REF <- file.path(BATCH, "dataset", "dnmt_refs")
ref_files <- file.path(DNMT_REF, c("vertebrate_reviewed.fasta", "mollusca_uniprot.fasta", "outgroup_MHhaI.fasta"))
stopifnot(all(file.exists(ref_files)), all(file.size(ref_files) > 100))
refs <- do.call(c, lapply(ref_files, readAAStringSet))
names(refs) <- sub(" .*", "", names(refs))                        # sp|P26358|DNMT1_HUMAN
ref_meta <- fread(file.path(DNMT_REF, "mollusca_uniprot.tsv"))
setnames(ref_meta, c("acc", "entry", "gene_names", "organism", "length", "reviewed", "protein_names"))
prot_all <- readAAStringSet(PROTE); names(prot_all) <- sub(" .*", "", names(prot_all))
dl_loci <- toolkit[gene_symbol %in% c("DNMT1", "DNMT2") & present == TRUE, gene_id]
dl_seq <- do.call(c, lapply(dl_loci, function(L) {                 # longest isoform per locus
  s <- prot_all[startsWith(names(prot_all), paste0(L, "-"))]; s[which.max(width(s))] }))
names(dl_seq) <- paste0("DLAEVE|", dl_loci, "|", toolkit[match(dl_loci, gene_id), gene_symbol])
all_aa <- c(dl_seq, refs)
all_aa <- all_aa[!duplicated(names(all_aa))]
# simple tip ids: MAFFT/IQ-TREE keep names, but "|" and spaces in Newick are fragile
tipmap <- data.table(tip = sprintf("t%03d", seq_along(all_aa)), name = names(all_aa))
names(all_aa) <- tipmap$tip
acc_of <- function(nm) ifelse(startsWith(nm, "DLAEVE"), sub("^DLAEVE\\|([^|]+)\\|.*$", "\\1", nm),
                              sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", nm))
tipmap[, acc := acc_of(name)]
tipmap[, entry := ifelse(startsWith(name, "DLAEVE"), sub(".*\\|", "", name), sub("^(sp|tr)\\|[^|]+\\|", "", name))]
tipmap <- merge(tipmap, ref_meta[, .(acc, gene_names, organism)], by = "acc", all.x = TRUE)
tipmap[, species := ifelse(startsWith(name, "DLAEVE"), "Deroceras laeve",
                    ifelse(!is.na(organism), sub("^(\\S+ \\S+).*$", "\\1", organism),
                    ifelse(grepl("_HUMAN$", entry), "Homo sapiens", ifelse(grepl("_MOUSE$", entry), "Mus musculus",
                    ifelse(grepl("_HAEHA$|MTH1", entry), "Haemophilus haemolyticus (M.HhaI)", "unknown")))))]
gene_lab <- toupper(ifelse(startsWith(tipmap$name, "DLAEVE"), tipmap$entry,
                    ifelse(!is.na(tipmap$gene_names), sub(" .*", "", tipmap$gene_names), sub("_.*$", "", tipmap$entry))))
tipmap[, gene := gene_lab]
tipmap[, class := fifelse(startsWith(name, "DLAEVE"), "D. laeve",
                  fifelse(grepl("^DNMT3|^DNMT3A|^DNMT3B|^DNMT3L", gene), "DNMT3 family",
                  fifelse(grepl("^DNMT1", gene), "DNMT1",
                  fifelse(grepl("^TRDMT1|^DNMT2", gene), "TRDMT1 / DNMT2",
                  fifelse(species == "Haemophilus haemolyticus (M.HhaI)", "Outgroup (M.HhaI)", "Mollusc, unnamed")))))]
tipmap[, label := sprintf("%s  %s  (%s)", species, gene, acc)]
scr8 <- file.path(if (nzchar(Sys.getenv("SLURM_JOB_ID"))) file.path("/scratch/groups/alfredvar", Sys.getenv("USER"),
                                                                     paste0("job_", Sys.getenv("SLURM_JOB_ID"))) else tempdir(), "dnmt_tree")
dir.create(scr8, recursive = TRUE, showWarnings = FALSE)
fa8 <- file.path(scr8, "all.faa"); writeXStringSet(all_aa, fa8)
hmm8 <- file.path(scr8, "DNA_methylase.hmm")
system2(file.path(HMMER_BIN, "hmmfetch"), c(PFAM, "DNA_methylase"), stdout = hmm8)
dom8 <- file.path(scr8, "dna_methylase.domtblout")
system2(file.path(HMMER_BIN, "hmmsearch"), c("--cut_ga", "--domtblout", dom8, "-o", "/dev/null", hmm8, fa8))
dt8 <- fread(cmd = sprintf("awk '!/^#/{print $1\"\\t\"$14\"\\t\"$20\"\\t\"$21}' %s", shQuote(dom8)),
             header = FALSE, col.names = c("tip", "score", "env_from", "env_to"))
dt8 <- dt8[order(tip, -score)][, .SD[1], by = tip][env_to - env_from + 1L >= 150L]   # best, near-complete domain per protein
tipmap[, in_tree := tip %in% dt8$tip]
fwrite(tipmap[, .(tip, acc, species, gene, class, in_tree, name)], file.path(DAT, "dnmt_c5_tree_members.tsv"), sep = "\t")
cat(sprintf("  %d of %d proteins carry a DNA_methylase domain of >= 150 aa; D. laeve tips in the tree: %s\n",
            nrow(dt8), length(all_aa), paste(tipmap[in_tree == TRUE & class == "D. laeve", gene], collapse = ", ")))
if (nrow(dt8) >= 4) {
  dom_aa <- AAStringSet(vapply(seq_len(nrow(dt8)), function(i)
    as.character(subseq(all_aa[[dt8$tip[i]]], dt8$env_from[i], dt8$env_to[i])), character(1)))
  names(dom_aa) <- dt8$tip
  writeXStringSet(dom_aa, file.path(scr8, "dom.faa"))
  NP8 <- Sys.getenv("SLURM_CPUS_PER_TASK", "4")
  st <- system2("mafft", c("--auto", "--quiet", "--anysymbol", "--thread", NP8, shQuote(file.path(scr8, "dom.faa"))),
                stdout = file.path(scr8, "aln.faa")); stopifnot(st == 0)
  st <- system2("trimal", c("-in", shQuote(file.path(scr8, "aln.faa")), "-out", shQuote(file.path(scr8, "aln.trim.faa")),
                            "-automated1")); stopifnot(st == 0)
  st <- system2("iqtree3", c("-s", shQuote(file.path(scr8, "aln.trim.faa")), "-m", "MFP", "-bb", "1000", "-nt", NP8,
                             "-seed", "20260426", "-pre", shQuote(file.path(scr8, "dnmt")), "-quiet", "-redo")); stopifnot(st == 0)
  tr8 <- phangorn::midpoint(read.tree(file.path(scr8, "dnmt.treefile")), node.labels = "support")
  file.copy(file.path(scr8, "dnmt.treefile"), file.path(DAT, "dnmt_c5_tree_unrooted.nwk"), overwrite = TRUE)
  write.tree(tr8, file.path(DAT, "dnmt_c5_tree_midpoint.nwk"))
  file.copy(file.path(scr8, "dnmt.iqtree"), file.path(DAT, "dnmt_c5_tree_iqtree_report.txt"), overwrite = TRUE)
  tips8 <- tipmap[match(tr8$tip.label, tip)]
  col8 <- c(`D. laeve` = "#D55E00", DNMT1 = "#0072B2", `DNMT3 family` = "#009E73",
            `TRDMT1 / DNMT2` = "#CC79A7", `Mollusc, unnamed` = "grey45", `Outgroup (M.HhaI)` = "black")
  draw_tree <- function() {
    op <- par(mar = c(3, 0.5, 2, 0.5)); on.exit(par(op))
    plot(tr8, show.tip.label = TRUE, tip.color = col8[tips8$class], cex = 0.42, label.offset = 0.01,
         font = ifelse(tips8$class == "D. laeve", 2, 1), no.margin = FALSE,
         main = "C5-cytosine methyltransferase domains (PF00145), midpoint rooted", cex.main = 0.8)
    tiplabels(pch = 16, col = col8[tips8$class], cex = 0.5, adj = 0.5)
    sup8 <- suppressWarnings(as.numeric(tr8$node.label)); ok8 <- !is.na(sup8) & sup8 >= 70
    nodelabels(text = ifelse(ok8, sup8, ""), frame = "none", cex = 0.35, adj = c(1.2, -0.3), col = "grey30")
    add.scale.bar(cex = 0.5); legend("bottomleft", legend = names(col8), col = col8, pch = 16, bty = "n", cex = 0.55)
  }
  tr8$tip.label <- tips8$label
  for (ext in c("pdf", "png", "svg")) {
    f8 <- file.path(SUPP, paste0("figS_dnmt_c5_tree.", ext))
    if (ext == "pdf") grDevices::cairo_pdf(f8, width = 7.5, height = 0.14 * length(tr8$tip.label) + 1.5)
    else if (ext == "png") png(f8, width = 7.5, height = 0.14 * length(tr8$tip.label) + 1.5, units = "in", res = 150)
    else svglite::svglite(f8, width = 7.5, height = 0.14 * length(tr8$tip.label) + 1.5)
    draw_tree(); dev.off()
  }
  cat("  saved (supp) figS_dnmt_c5_tree\n")
} else cat("  fewer than 4 domain sequences: tree skipped\n")

# ---- [TF] TRANSCRIPTION-FACTOR ANNOTATION ----
# codes". The complete sequence-orthology pipeline (JASPAR -> UniProt -> DIAMOND RBH
# -> Pfam DBD licensing -> Cis-BP thresholds -> gene trees -> the bridge + the
# tf_dlaeve_guide) now runs here; 08_motifs (motif enrichment) reads
# 01_genome_toolkit/data/jaspar_ortholog_bridge.tsv. Wrapped in local() so its helpers and
# constants (run/need/fw, theme_pub, BATCH/OBJ/DAT/FIGS, palettes) cannot collide
# with 01_genome_toolkit's; all outputs land under 01_genome_toolkit/{data,objects,figures}. Staged
# network inputs (bait FASTA, Cis-BP thresholds) live in 01_genome_toolkit/objects now; the
# BATCH02B_FETCH_ONLY=1 node12 staging mode still works and exits the whole script.
local({
set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(Biostrings); library(RSQLite)
  library(ggplot2); library(patchwork); library(svglite)
  library(ape); library(phangorn); library(parallel); library(GenomicRanges)
})

# ---- [0] Configuration -------------------------------------------------------
PIPE   <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
ROOT   <- "/mnt/data/alfredvar/rlopezt/meth_paper"
BATCH  <- file.path(PIPE, "01_genome_toolkit")
OBJ    <- file.path(BATCH, "objects"); DAT <- file.path(BATCH, "data")
FIGS   <- file.path(BATCH, "figures/supplementary")
SEQD   <- file.path(OBJ, "seq"); DIAD <- file.path(OBJ, "diamond")
PFAMD  <- file.path(OBJ, "pfam"); TREED <- file.path(OBJ, "trees")
for (d in c(OBJ, DAT, FIGS, SEQD, DIAD, PFAMD, TREED))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

JASPAR_DB <- file.path(ROOT, "tools/jaspar/JASPAR2024.sqlite")
PFAM_HMM  <- file.path(ROOT, "tools/pfam/Pfam-A.hmm")     # already hmmpress'ed, never redo
HMMER_BIN <- file.path(ROOT, "tools/hmmer/bin")           # 3.4: the build that pressed Pfam-A
# Canonical read-only EviAnn proteome in the shared lab tree (the old implementation
PROTEOME  <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.proteins.fasta"
TAXGROUPS <- c("vertebrates", "insects", "nematodes", "urochordates")  # mollusks are not a JASPAR group
NPROC     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

# motif library for the WHOLE paper, so they are named here, never buried in a
# filter; the pwm_library_decision_grid table records what every combination gives.
#   ORTH_HOMOLOGY  "rbh" = reciprocal best hit only (strict) | "rbh_or_tree" = RBH or resolved tree placement
ORTH_HOMOLOGY <- "rbh"
ORTH_IDENTITY <- "cisbp"
ORTH_DIMER    <- "any"     # heterodimer matrices (A::B): "any" = kept when at least one subunit passes the gate,
                           # "both" = every subunit must pass. Recorded per matrix as n_subunits_kept / n_subunits.
TREE_MAX_CLADE   <- 10L   # a placement into a clade larger than this resolved nothing
TREE_MIN_SUPPORT <- 70    # ultrafast bootstrap at the clade-defining node

stopifnot(file.exists(JASPAR_DB), file.exists(PFAM_HMM), file.exists(PROTEOME),
          dir.exists(HMMER_BIN))

# ---- [0b] Helpers ------------------------------------------------------------
# two once wrote diagnostics into dbd.hmm -> hmmpress -> garbage domain table.
run <- function(cmd, args, out = NULL, quiet = TRUE) {
  t0 <- Sys.time()
  err <- tempfile()
  st <- system2(cmd, args, stdout = if (is.null(out)) "" else out, stderr = err)
  cat(sprintf("    %s -> exit %d (%.1f min)\n", basename(cmd), st,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (st != 0) {
    msg <- paste(utils::head(readLines(err, warn = FALSE), 15), collapse = "\n")
    unlink(err); stop(sprintf("%s failed (exit %d):\n%s", basename(cmd), st, msg), call. = FALSE)
  }
  unlink(err); invisible(st)
}
# Trust outputs, not exit codes: jobs here have exited 0 with truncated results.
need <- function(path, what, min_rows = 1L, comment = "#") {
  if (!file.exists(path) || file.size(path) == 0)
    stop(sprintf("%s produced no output at %s", what, path), call. = FALSE)
  # 374 real rows that passed an existence-only check; hence the row-count floor.
  n <- length(grep(sprintf("^%s", comment), readLines(path, warn = FALSE), invert = TRUE, value = TRUE))
  if (n < min_rows)
    stop(sprintf("%s produced only %d data rows at %s (expected >= %d)", what, n, path, min_rows), call. = FALSE)
  path
}
fw <- function(x, path) {                       # atomic: never leave a truncated file
  p <- paste0(path, ".part"); fwrite(x, p, sep = "\t"); stopifnot(file.rename(p, path)); path
}
read_domtbl <- function(path, fam_col, q_col) {
  fread(cmd = sprintf("awk '!/^#/{split($%d,a,\"|\"); print $%d\"\\t\"(a[2]!=\"\"?a[2]:$%d)\"\\t\"$14\"\\t\"$20\"\\t\"$21}' %s",
                      q_col, fam_col, q_col, shQuote(path)),
        header = FALSE, col.names = c("pfam", "query", "score", "env_from", "env_to"))
}
COL_A <- "#009E73"; COL_B <- "#E69F00"; COL_MARK <- "#333333"   # Okabe-Ito; #0072B2/#D55E00 reserved
theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          legend.key.size = unit(0.8, "lines"), legend.background = element_blank(),
          legend.key = element_blank(), plot.title = element_text(face = "plain", size = base_size + 1))
}
save_supp <- function(p, stem, width, height) {
  ggsave(file.path(FIGS, paste0(stem, ".pdf")), p, width = width, height = height, device = grDevices::cairo_pdf)
  ggsave(file.path(FIGS, paste0(stem, ".png")), p, width = width, height = height, dpi = 400)
  ggsave(file.path(FIGS, paste0(stem, ".svg")), p, width = width, height = height, device = svglite::svglite)
  cat(sprintf("  saved %s (%.1f x %.1f in)\n", stem, width, height))
}

# ---- [1] JASPAR matrix -> UniProt accession ----------------------------------
# Latest VERSION per BASE_ID in CORE, four animal tax groups, joined to MATRIX_PROTEIN.
# Long format: a heterodimer (ELK1::HOXA1) expands to one row per partner, never
# dropped. Writes matrix_to_acc.tsv.
cat("[TF 1] JASPAR2024 CORE -> UniProt accessions\n")
con <- dbConnect(SQLite(), JASPAR_DB)
m2a <- as.data.table(dbGetQuery(con, sprintf("
  WITH latest AS (SELECT BASE_ID, MAX(VERSION) AS V FROM MATRIX WHERE COLLECTION='CORE' GROUP BY BASE_ID),
       m AS (SELECT M.ID, M.BASE_ID || '.' || M.VERSION AS matrix_id, M.NAME AS tf_name
             FROM MATRIX M JOIN latest L ON M.BASE_ID=L.BASE_ID AND M.VERSION=L.V WHERE M.COLLECTION='CORE'),
       mt AS (SELECT m.ID, m.matrix_id, m.tf_name, A.VAL AS tax_group
              FROM m JOIN MATRIX_ANNOTATION A ON A.ID=m.ID AND A.TAG='tax_group'
              WHERE A.VAL IN (%s))
  SELECT mt.matrix_id, mt.tf_name, mt.tax_group, COALESCE(P.ACC,'') AS acc
  FROM mt LEFT JOIN MATRIX_PROTEIN P ON P.ID=mt.ID ORDER BY mt.matrix_id;",
  paste(sprintf("'%s'", TAXGROUPS), collapse = ","))))
mclass <- as.data.table(dbGetQuery(con, "
  WITH latest AS (SELECT BASE_ID, MAX(VERSION) AS V FROM MATRIX WHERE COLLECTION='CORE' GROUP BY BASE_ID)
  SELECT M.BASE_ID || '.' || M.VERSION AS matrix_id, A.VAL AS class
  FROM MATRIX M JOIN latest L ON M.BASE_ID=L.BASE_ID AND M.VERSION=L.V
  JOIN MATRIX_ANNOTATION A ON A.ID=M.ID AND A.TAG='class' WHERE M.COLLECTION='CORE';"))
dbDisconnect(con)
N_MATRIX <- uniqueN(m2a$matrix_id)
accs <- sort(unique(m2a[nzchar(acc), acc]))
cat(sprintf("  %d matrices, %d with an accession, %d unique accessions\n",
            N_MATRIX, uniqueN(m2a[nzchar(acc), matrix_id]), length(accs)))
fw(m2a, file.path(DAT, "matrix_to_acc.tsv"))

# ---- [2] Bait sequences from UniProt (needs network: node12, not defq) -------
# /uniprotkb/accessions, comma-separated, batches of 150; canonical sequence only.
# returns an EMPTY body with HTTP 200.
BAITS <- file.path(SEQD, "bait_tfs.faa")
if (!file.exists(BAITS)) {
  cat("[TF 2] fetching bait sequences from UniProt (needs network: node12, not defq)\n")
  chunks <- split(accs, ceiling(seq_along(accs) / 150))
  # so only the final batch survives (30 of 1,230 sequences, every curl exit 0).
  got <- character(0)
  for (i in seq_along(chunks)) {
    out <- system2("curl", c("-sS", "--retry", "3", "--retry-delay", "2", "--max-time", "180",
                             "-G", "https://rest.uniprot.org/uniprotkb/accessions",
                             "--data-urlencode", paste0("accessions=", paste(chunks[[i]], collapse = ",")),
                             "--data", "format=fasta"), stdout = TRUE, stderr = TRUE)
    n <- sum(startsWith(out, ">"))
    if (n == 0) stop(sprintf("UniProt returned 0 sequences for batch %d. On a defq compute node? Run once on node12.", i), call. = FALSE)
    got <- c(got, out)
    cat(sprintf("  batch %d/%d: %d sequences (cumulative %d)\n", i, length(chunks), n, sum(startsWith(got, ">"))))
  }
  tmp <- paste0(BAITS, ".part"); writeLines(got, tmp); file.rename(tmp, BAITS)
} else cat("[TF 2] bait sequences already staged\n")
bait_aa <- readAAStringSet(need(BAITS, "UniProt fetch"))
names(bait_aa) <- vapply(strsplit(names(bait_aa), "\\|"), function(x) if (length(x) > 1) x[2] else x[1], character(1))
bait_aa <- bait_aa[!duplicated(names(bait_aa))]
missing_acc <- setdiff(accs, names(bait_aa))
cat(sprintf("  %d sequences; %d accessions missing\n", length(bait_aa), length(missing_acc)))
# whole library; >10% missing baits means restage on node12, not proceed
stopifnot(length(missing_acc) < 0.10 * length(accs))

# Fetch-only mode (BATCH02B_FETCH_ONLY=1, once on node12): stage [2] and verify the
# [7] thresholds file, then stop — defq compute nodes have no network (see header).
if (Sys.getenv("BATCH02B_FETCH_ONLY") == "1") {
  stopifnot(file.exists(THRESH_CHECK <- file.path(OBJ, "cisbp_family_thresholds.tsv")))
  cat("[TF fetch-only] bait sequences fetched; the Cis-BP thresholds file (external, dataset/cisbp) is present; stopping before the compute\n")
  quit(save = "no", status = 0)
}

# ---- [3] DIAMOND both directions, isoforms collapsed to locus, RBH -----------
cat("[TF 3] DIAMOND blastp both directions vs the EviAnn proteome\n")
DB   <- file.path(DIAD, "proteome")
FWD  <- file.path(DIAD, "fwd.tsv"); REV <- file.path(DIAD, "rev.tsv")
OFMT <- c("6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore")
run("diamond", c("makedb", "--in", shQuote(PROTEOME), "--db", shQuote(DB), "--quiet"))
writeXStringSet(bait_aa, file.path(SEQD, "baits.clean.faa"))
run("diamond", c("blastp", "--ultra-sensitive", "--quiet", "-p", NPROC, "-k", "25",
                 "-q", shQuote(file.path(SEQD, "baits.clean.faa")), "-d", shQuote(paste0(DB, ".dmnd")),
                 "--outfmt", OFMT, "-o", shQuote(FWD)))
fwd <- fread(need(FWD, "DIAMOND forward"), header = FALSE,
             col.names = c("acc", "protein", "pident", "len", "evalue", "bits"))
fwd[, locus := sub("-mRNA-[0-9]+$", "", protein)]
# collapse isoforms: one representative per locus, the highest-scoring isoform
best_iso <- fwd[order(-bits), .SD[1], by = .(acc, locus)]
# forward best locus per bait
fwd_best <- best_iso[order(-bits), .SD[1], by = acc][, .(acc, locus, fwd_bits = bits)]
# reverse: each candidate locus representative back against the baits
rep_iso <- unique(best_iso[order(-bits), .SD[1], by = locus][, .(locus, protein)])
prot_all <- readAAStringSet(PROTEOME)
names(prot_all) <- sub("\\s.*$", "", names(prot_all))
writeXStringSet(prot_all[rep_iso$protein], file.path(SEQD, "locus_reps.faa"))
run("diamond", c("makedb", "--in", shQuote(file.path(SEQD, "baits.clean.faa")),
                 "--db", shQuote(file.path(DIAD, "baits")), "--quiet"))
run("diamond", c("blastp", "--ultra-sensitive", "--quiet", "-p", NPROC, "-k", "25",
                 "-q", shQuote(file.path(SEQD, "locus_reps.faa")), "-d", shQuote(file.path(DIAD, "baits.dmnd")),
                 "--outfmt", OFMT, "-o", shQuote(REV)))
rev_ <- fread(need(REV, "DIAMOND reverse"), header = FALSE,
              col.names = c("protein", "acc", "pident", "len", "evalue", "bits"))
rev_ <- merge(rev_, rep_iso, by = "protein")
rev_best <- rev_[order(-bits), .SD[1], by = locus][, .(locus, acc_back = acc)]
cand <- merge(best_iso[, .(acc, locus, protein, fwd_bits = bits)], rev_best, by = "locus", all.x = TRUE)
cand <- merge(cand, fwd_best[, .(acc, best_locus = locus)], by = "acc", all.x = TRUE)
# secondary hit whose locus points back and roughly doubles the RBH count.
cand[, is_rbh := fifelse(!is.na(acc_back) & acc_back == acc &
                         !is.na(best_locus) & locus == best_locus, "yes", "no")]
cand <- cand[acc %in% fwd_best$acc | is_rbh == "yes"]
cat(sprintf("  %d candidate pairs, %d reciprocal best hits\n", nrow(cand), sum(cand$is_rbh == "yes")))

# ---- [4] Pfam domains --------------------------------------------------------
cat("[TF 4] Pfam domains (hmmscan on baits; hmmsearch of DBD profiles over the proteome)\n")
BAIT_DOM <- file.path(PFAMD, "baits.domtbl")
run(file.path(HMMER_BIN, "hmmscan"),
    c("--cut_ga", "--cpu", NPROC, "--domtblout", shQuote(BAIT_DOM), "-o", "/dev/null",
      shQuote(PFAM_HMM), shQuote(file.path(SEQD, "baits.clean.faa"))))
bd <- read_domtbl(need(BAIT_DOM, "hmmscan on baits", min_rows = 500L), 1, 4)

# The JASPAR structural class -> Pfam DBD map, reused verbatim from
# analysis/tf_dbd/01_dbd_verification.R so the two agree by construction.
CLASS2PFAM <- list(
  "Homeo domain factors" = c("Homeodomain","Homeobox","Homeobox_KN","Pou","CUT","HPD","PBC","SIX1_SD"),
  "C2H2 zinc finger factors" = c("zf-C2H2","zf-C2H2_2","zf-C2H2_3","zf-C2H2_4","zf-C2H2_5","zf-C2H2_6",
    "zf-C2H2_7","zf-C2H2_8","zf-C2H2_9","zf-C2H2_10","zf-C2H2_11","zf-C2H2_jaz","zf-H2C2","zf-H2C2_2",
    "zf-H2C2_5","zf-met","zf-met2","zf-BED"),
  "Basic helix-loop-helix factors (bHLH)" = "HLH",
  "Basic leucine zipper factors (bZIP)" = c("bZIP_1","bZIP_2","bZIP_Maf"),
  "Nuclear receptors with C4 zinc fingers" = c("zf-C4"),
  "Tryptophan cluster factors" = c("Ets","IRF","Myb_DNA-binding","Myb_DNA-bind_6"),
  "Fork head/winged helix factors" = c("Forkhead","E2F_TDP","WHD_E2F_TDP","RFX_DNA_binding"),
  "High-mobility group (HMG) domain factors" = c("HMG_box","HMG_box_2"),
  "Rel homology region (RHR) factors" = c("RHD_DNA_bind","COE1_DBD","BTD","LAG1-DNAbind"),
  "SMAD/NF-1 DNA-binding domain factors" = c("MH1","CTF_NFI"),
  "Paired box factors" = "PAX", "T-Box factors" = "T-box",
  "Other C4 zinc finger-type factors" = c("GATA","zf-C4"),
  "Basic helix-span-helix factors (bHSH)" = "TF_AP-2", "TEA domain factors" = "TEA",
  "STAT domain factors" = "STAT_bind", "MADS box factors" = "SRF-TF",
  "Heteromeric CCAAT-binding factors" = c("CBFB_NFYA","CBFD_NFYB_HMF"),
  "Runt domain factors" = "Runt", "Grainyhead domain factors" = "CP2",
  "p53 domain factors" = "P53", "TATA-binding proteins" = "TBP", "SAND domain factors" = "SAND",
  "Heat shock factors" = "HSF_DNA-bind", "GCM domain factors" = "GCM",
  "DM-type intertwined zinc finger factors" = "DM", "CRC domain" = "TCR",
  "C2CH THAP-type zinc finger factors" = "THAP", "ARID" = "ARID")
# Two matrices carry no class row in JASPAR2024 but have an unambiguous DBD.
NOCLASS2PFAM <- list("MA0506.3" = "Nrf1_DNA-bind", "MA1618.2" = "HLH")
DBD_NAMES <- sort(unique(unlist(CLASS2PFAM, use.names = FALSE)))

# Genome-wide DBD presence: hmmsearch the DBD profiles only (far faster than hmmscan
# over all 85,163 proteins). NOT eggNOG: its vocabulary says Homeobox, so a literal
# lookup of Homeodomain returns zero loci.
DBD_HMM <- file.path(PFAMD, "dbd.hmm")
# missing name AFTER writing earlier profiles, leaving a truncated, plausible dbd.hmm.
# This release lacks E2F_TDP and Homeobox (redundant aliases of WHD_E2F_TDP and
# Homeodomain), so dropping them costs nothing; the drop is REPORTED, never silent.
pfam_on_disk <- unique(sub("^NAME\\s+", "", grep("^NAME", readLines(PFAM_HMM, warn = FALSE), value = TRUE)))
absent <- setdiff(DBD_NAMES, pfam_on_disk)
if (length(absent))
  cat(sprintf("  %d mapped DBD name(s) absent from Pfam-A, dropped: %s\n",
              length(absent), paste(absent, collapse = ", ")))
DBD_NAMES <- intersect(DBD_NAMES, pfam_on_disk)
stopifnot(length(DBD_NAMES) >= 50)
writeLines(DBD_NAMES, file.path(PFAMD, "dbd_names.txt"))
run(file.path(HMMER_BIN, "hmmfetch"), c("-f", shQuote(PFAM_HMM), shQuote(file.path(PFAMD, "dbd_names.txt"))),
    out = DBD_HMM)
need(DBD_HMM, "hmmfetch", min_rows = 100L, comment = "$^")
stopifnot(sum(startsWith(readLines(DBD_HMM, warn = FALSE), "NAME")) >= length(DBD_NAMES) * 0.8)
run(file.path(HMMER_BIN, "hmmpress"), c("-f", shQuote(DBD_HMM)))
PROT_DOM <- file.path(PFAMD, "proteome_dbd.domtbl")
run(file.path(HMMER_BIN, "hmmsearch"),
    c("--cut_ga", "--cpu", NPROC, "--domtblout", shQuote(PROT_DOM), "-o", "/dev/null",
      shQuote(DBD_HMM), shQuote(PROTEOME)))
pd <- read_domtbl(need(PROT_DOM, "hmmsearch on proteome", min_rows = 1000L), 4, 1)   # hmmsearch swaps the columns
pd[, locus := sub("-mRNA-[0-9]+$", "", query)]

# ---- [5] DBD licensing -------------------------------------------------------
cat("[TF 5] DBD licensing: both proteins must carry the class-expected domain\n")
acc2class <- merge(m2a[nzchar(acc), .(matrix_id, acc)], mclass, by = "matrix_id", all.x = TRUE, allow.cartesian = TRUE)
expected <- acc2class[, .(pfam = unique(unlist(c(CLASS2PFAM[class],
                                                 NOCLASS2PFAM[intersect(matrix_id, names(NOCLASS2PFAM))])))),
                      by = acc]
expected <- expected[!is.na(pfam)]
bait_fam <- unique(bd[, .(acc = query, pfam, score, env_from, env_to)])
locus_fam <- unique(pd[, .(protein = query, pfam, score, env_from, env_to)])
lic <- merge(cand[, .(acc, locus, protein, fwd_bits, is_rbh)],
             merge(bait_fam, expected, by = c("acc", "pfam")), by = "acc", allow.cartesian = TRUE)
lic <- merge(lic, locus_fam, by = c("protein", "pfam"), suffixes = c("_b", "_l"))
# one row per (acc, locus): the best-scoring shared expected family
lic <- lic[order(-score_b)][, .SD[1], by = .(acc, locus)]
cat(sprintf("  %d licensed pairs over %d loci\n", nrow(lic), uniqueN(lic$locus)))

# ---- [6] DBD identity, on both denominators ----------------------------------
# PctID_O = identities / both-residue columns (ours); PctID_L = identities / length
# on PctID_L, and PctID_O runs larger whenever the alignment gaps: both are written,
# THE GATE USES PctID_L.
cat(sprintf("[6] pairwise DBD alignment for %d licensed pairs (MAFFT, %d cores)\n", nrow(lic), NPROC))
prot_rep <- prot_all[unique(lic$protein)]
# shared storage once failed 9,175 of 9,175 alignments on the file race.
pair_id <- function(i) {
  b <- subseq(bait_aa[[lic$acc[i]]], lic$env_from_b[i], lic$env_to_b[i])
  l <- subseq(prot_rep[[lic$protein[i]]], lic$env_from_l[i], lic$env_to_l[i])
  out <- suppressWarnings(system2("mafft", c("--auto", "--quiet", "--anysymbol", "-"),
                                  stdout = TRUE,
                                  input = c(">a", as.character(b), ">b", as.character(l))))
  if (length(out) < 3) return(c(NA_real_, NA_real_, NA_real_))
  s <- strsplit(paste(out, collapse = "\n"), ">")[[1]]; s <- s[nzchar(s)]
  aln <- vapply(s, function(x) gsub("[^A-Za-z-]", "", sub("^[^\n]*\n", "", x)), character(1))
  if (length(aln) != 2) return(c(NA_real_, NA_real_, NA_real_))
  a1 <- strsplit(aln[1], "")[[1]]; a2 <- strsplit(aln[2], "")[[1]]
  k <- a1 != "-" & a2 != "-"
  if (!any(k)) return(c(NA_real_, NA_real_, NA_real_))
  ident <- sum(a1[k] == a2[k])
  c(100 * ident / sum(k),                                             # PctID_O
    100 * ident / max(nchar(b), nchar(l)),                            # PctID_L
    sum(k))
}
# unguarded do.call(rbind, ...) then silently stops being a matrix, and a swallow-
# everything run once died at [8] with an empty table, far from the real cause.
probe <- tryCatch(pair_id(1L), error = function(e) stop(sprintf(
  "the first DBD alignment failed, so the rest will too: %s", conditionMessage(e)), call. = FALSE))
if (is.na(probe[2])) stop("the first DBD alignment returned no result; is mafft on PATH?", call. = FALSE)
pid_safe <- function(i) tryCatch(pair_id(i), error = function(e) c(NA_real_, NA_real_, NA_real_))
pid <- do.call(rbind, mclapply(seq_len(nrow(lic)), pid_safe, mc.cores = NPROC))
stopifnot(is.matrix(pid), nrow(pid) == nrow(lic))
ok <- sum(!is.na(pid[, 2]))
cat(sprintf("  %d of %d pairs aligned (%.1f%%)\n", ok, nrow(pid), 100 * ok / nrow(pid)))
# A few unalignable pairs are normal; a mass failure is a broken run, not a result.
if (ok < 0.5 * nrow(pid))
  stop(sprintf("only %d of %d DBD alignments succeeded: refusing to continue on a broken step",
               ok, nrow(pid)), call. = FALSE)
lic[, `:=`(pid_dbd = pid[, 1], pid_dbd_L = pid[, 2], dbd_aln_cols = pid[, 3])]
lic <- lic[!is.na(pid_dbd_L)]
fw(lic[, .(acc, locus, protein, pfam_dbd = pfam, pid_dbd = round(pid_dbd, 2),
           pid_dbd_L = round(pid_dbd_L, 2), dbd_aln_cols, is_rbh, fwd_bits)],
   file.path(DAT, "dbd_identity.tsv"))
cat(sprintf("  %d pairs with an identity; median PctID_L %.1f\n", nrow(lic), median(lic$pid_dbd_L)))

# ---- [7] Cis-BP / Weirauch motif-transfer thresholds (staged file) -----------
# Per-DBD-family identity above which a measured motif transfers. Weirauch et al.
# 2014 Cell (PMID 25215497) set the per-family thresholds; Lambert et al. 2019 Nat
# Genet (PMID 31133749) showed motif divergence is pervasive, worst in C2H2 zinc
# fingers, which is why those families end up near-untransferable.
THRESH <- file.path(OBJ, "cisbp_family_thresholds.tsv")
if (!file.exists(THRESH))
  stop(sprintf(paste("Cis-BP thresholds not staged at %s.\n",
       "This file is an EXTERNAL staged input: dataset/cisbp/00_fetch.sh downloads Cis-BP 3.10 and",
       "dataset/cisbp/01_build_thresholds.py writes cisbp_family_thresholds.tsv (see dataset/cisbp/PROVENANCE.md);",
       "copy it into objects/ before running. Compute nodes have no network."), THRESH), call. = FALSE)
thr <- fread(THRESH)
thr[, `:=`(t_cis = suppressWarnings(as.numeric(threshold_pct)),
           t_w14 = suppressWarnings(as.numeric(weirauch2014_threshold_pct)))]
lic <- merge(lic, thr[, .(pfam_dbd, t_cis, t_w14)], by.x = "pfam", by.y = "pfam_dbd", all.x = TRUE)
# gate then EXCLUDES those families; that exclusion was silent — now it is counted
# (and IRF's literal 100 is a known Cis-BP failure code, applied as never-pass)
cat(sprintf("  Cis-BP cutoffs: %d licensed pairs across %d families carry a non-numeric threshold (NA -> excluded by the cisbp gate)\n",
            lic[is.na(t_cis), .N], lic[is.na(t_cis), uniqueN(pfam)]))

# ---- [8] Per-family gene trees -----------------------------------------------
# -automated1 strips a full-length alignment to almost nothing (HLH: 5 columns over
# 194 sequences), giving an uninformative tree that returns the whole family.
cat("[TF 8] per-family DBD gene trees\n")
fam_tab <- lic[, .(n_baits = uniqueN(acc), n_loci = uniqueN(locus)), by = pfam][n_baits >= 2 & n_loci >= 1]
# rejects bootstrap-less placements anyway, and their loci still reach the bridge
# through the RBH arm, which needs no tree.
fam_small <- fam_tab[n_baits + n_loci < 4L]
if (nrow(fam_small))
  cat(sprintf("  %d famil%s too small to bootstrap (< 4 sequences), no tree: %s\n",
              nrow(fam_small), if (nrow(fam_small) == 1L) "y is" else "ies are",
              paste(sort(fam_small$pfam), collapse = ", ")))
fam_tab <- fam_tab[n_baits + n_loci >= 4L]
tree_rows <- rbindlist(lapply(fam_tab$pfam, function(f) {
  sub <- lic[pfam == f]
  fa <- c(setNames(vapply(unique(sub$acc), function(a) {
            r <- sub[acc == a][1]; as.character(subseq(bait_aa[[a]], r$env_from_b, r$env_to_b)) },
            character(1)), paste0("BAIT_", unique(sub$acc))),
          setNames(vapply(unique(sub$locus), function(L) {
            r <- sub[locus == L][1]; as.character(subseq(prot_rep[[r$protein]], r$env_from_l, r$env_to_l)) },
            character(1)), paste0("DLAE_", unique(sub$locus))))
  d <- file.path(TREED, f); dir.create(d, showWarnings = FALSE)
  writeXStringSet(AAStringSet(fa), file.path(d, "fam.faa"))
  run("mafft", c("--auto", "--quiet", "--anysymbol", "--thread", NPROC, shQuote(file.path(d, "fam.faa"))),
      out = file.path(d, "aln.faa"))
  run("trimal", c("-in", shQuote(file.path(d, "aln.faa")), "-out", shQuote(file.path(d, "aln.trim.faa")), "-automated1"))
  ncol_trim <- tryCatch(width(readAAStringSet(file.path(d, "aln.trim.faa")))[1], error = function(e) 0L)
  run("iqtree3", c("-s", shQuote(file.path(d, "aln.trim.faa")), "-m", "MFP", "-bb", "1000",
                   "-nt", NPROC, "-seed", "20260426", "-pre", shQuote(file.path(d, f)), "-quiet", "-redo"))
  tf <- file.path(d, paste0(f, ".treefile")); if (!file.exists(tf)) return(NULL)
  # ancestor is only defined on a ROOTED tree. Midpoint root first.
  tr <- tryCatch(phangorn::midpoint(read.tree(tf), node.labels = "support"), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  nt <- length(tr$tip.label)
  sup <- suppressWarnings(as.numeric(tr$node.label))   # absent label -> NA, never 1
  parent_of <- function(n) { e <- tr$edge[tr$edge[, 2] == n, 1]; if (length(e)) e else NA_integer_ }
  rbindlist(lapply(which(startsWith(tr$tip.label, "DLAE_")), function(ti) {
    node <- ti; clade <- character(0); nd <- NA_integer_
    repeat {
      node <- parent_of(node); if (is.na(node)) break
      tips <- extract.clade(tr, node)$tip.label
      b <- sub("^BAIT_", "", tips[startsWith(tips, "BAIT_")])
      if (length(b)) { clade <- b; nd <- node; break }
    }
    if (!length(clade)) return(NULL)
    data.table(pfam = f, trim_cols = ncol_trim, locus = sub("^DLAE_", "", tr$tip.label[ti]),
               n_baits_in_clade = length(clade), support = sup[nd - nt], baits = paste(clade, collapse = ";"))
  }))
}), fill = TRUE)
fw(tree_rows, file.path(DAT, "tree_coortholog_clades.tsv"))
tree_ok <- tree_rows[n_baits_in_clade <= TREE_MAX_CLADE & !is.na(support) & support >= TREE_MIN_SUPPORT]
tree_pairs <- unique(tree_ok[, .(acc = unlist(strsplit(baits, ";"))), by = locus][, .(acc, locus)])
cat(sprintf("  %d families, %d placements, %d passing the clade and support guard\n",
            uniqueN(tree_rows$pfam), nrow(tree_rows), nrow(tree_pairs)))

# ---- [9] The gate -> the bridge ----------------------------------------------
cat(sprintf("[9] applying the gate: homology=%s identity=%s\n", ORTH_HOMOLOGY, ORTH_IDENTITY))
lic[, tree_ok := paste(acc, locus) %in% paste(tree_pairs$acc, tree_pairs$locus)]
gate <- function(d, hom, idn) {
  if (idn != "none") { tc <- if (idn == "cisbp") "t_cis" else "t_w14"
                       d <- d[!is.na(get(tc)) & pid_dbd_L >= get(tc)] }
  switch(hom, any = d, rbh_or_tree = d[is_rbh == "yes" | tree_ok], rbh = d[is_rbh == "yes"])
}
kept <- gate(lic, ORTH_HOMOLOGY, ORTH_IDENTITY)
kept[, evidence := fifelse(is_rbh == "yes", "rbh", "tree")]
m_keep <- merge(m2a[nzchar(acc)], kept[, .(acc, locus, evidence, pid_dbd_L, pfam_dbd = pfam)],
                by = "acc", allow.cartesian = TRUE)
fw(m_keep[order(matrix_id, -pid_dbd_L)], file.path(DAT, "motif_to_dlaeve.tsv"))
loci_by_m <- m_keep[, .(dlaeve_locus = paste(sort(unique(locus)), collapse = ";"),
                        n_subunits_kept = uniqueN(acc)), by = matrix_id]
n_sub <- m2a[nzchar(acc), .(n_subunits = uniqueN(acc)), by = matrix_id]
bridge <- merge(unique(m2a[, .(motif_id = matrix_id, tf_name)]), loci_by_m,
                by.x = "motif_id", by.y = "matrix_id", all.x = TRUE)
bridge <- merge(bridge, n_sub, by.x = "motif_id", by.y = "matrix_id", all.x = TRUE)
bridge[is.na(dlaeve_locus), dlaeve_locus := ""]
bridge[is.na(n_subunits_kept), n_subunits_kept := 0L]
# ORTH_DIMER: "any" keeps a multi-subunit matrix when >= 1 subunit passed; "both" needs all
if (ORTH_DIMER == "both") bridge[n_subunits_kept < n_subunits, dlaeve_locus := ""]
bridge[, has_ortholog := fifelse(nzchar(dlaeve_locus), "TRUE", "FALSE")]
n_dimer <- bridge[n_subunits > 1 & has_ortholog == "TRUE"]
cat(sprintf("  heterodimer rule ORTH_DIMER = %s: %d kept multi-subunit matrices, %d of them with only one subunit passing\n",
            ORTH_DIMER, nrow(n_dimer), sum(n_dimer$n_subunits_kept < n_dimer$n_subunits)))
fw(bridge[order(motif_id), .(motif_id, tf_name, has_ortholog, dlaeve_locus, n_subunits, n_subunits_kept)],
   file.path(DAT, "jaspar_ortholog_bridge.tsv"))
N_KEPT <- bridge[has_ortholog == "TRUE", .N]
cat(sprintf("  BRIDGE: %d of %d matrices keep a D. laeve ortholog\n", N_KEPT, N_MATRIX))
stopifnot(N_KEPT >= 20, nrow(bridge) == N_MATRIX)

# ---- [9b] The D. laeve TF guide ----------------------------------------------
# The bridge says WHICH motifs survive; the guide records WHY: one row per kept
# motif -> locus (coordinates, symbol, DBD sequence, every filter's passing value).
# chromosome filter is INHERITED, never re-implemented. The proteome covers the
# whole assembly: a kept locus off chr1-31+mt stays, on_assembly=FALSE, never
# silently dropped.
cat("[TF 9b] the D. laeve TF guide\n")
gff <- readRDS(file.path(PIPE, "01_genome_toolkit/objects/gff_chrmt.rds"))
gene_gr <- gff[gff$type == "gene"]
note <- vapply(as.list(mcols(gene_gr)$Note), function(x)
  if (length(x)) as.character(x)[1] else NA_character_, character(1))
gtab <- data.table(locus  = sub(";.*", "", as.character(mcols(gene_gr)$ID)),
                   symbol = fifelse(grepl("^Similar to [^:]+:", note),
                                    sub("^Similar to ([^:]+):.*$", "\\1", note), NA_character_),
                   chrom  = as.character(seqnames(gene_gr)),
                   gene_start = start(gene_gr), gene_end = end(gene_gr),
                   gene_strand = as.character(strand(gene_gr)))
kd <- kept[, .(acc, locus, dlaeve_protein = protein, pfam_dbd = pfam,
               is_rbh, tree_ok, evidence, diamond_bits = fwd_bits,
               pid_dbd_overlap = round(pid_dbd, 2), pid_dbd_L = round(pid_dbd_L, 2),
               dbd_aln_cols, cisbp_threshold = t_cis, weirauch_threshold = t_w14,
               dbd_from = env_from_l, dbd_to = env_to_l)]
kd[, passes_cisbp := pid_dbd_L >= cisbp_threshold]
kd[, dbd_seq_dlaeve := vapply(seq_len(.N), function(i)
     as.character(subseq(prot_all[[dlaeve_protein[i]]], dbd_from[i], dbd_to[i])), character(1))]
jclass <- mclass[, .(jaspar_class = paste(sort(unique(class)), collapse = ";")), by = matrix_id]
guide <- merge(m2a[nzchar(acc), .(motif_id = matrix_id, tf_name, tax_group, acc)],
               kd, by = "acc", allow.cartesian = TRUE)
guide <- merge(guide, jclass, by.x = "motif_id", by.y = "matrix_id", all.x = TRUE)
guide <- merge(guide, gtab, by = "locus", all.x = TRUE)
guide[, on_assembly := !is.na(chrom)]
setcolorder(guide, c("motif_id", "tf_name", "jaspar_class", "tax_group", "acc",
                     "locus", "symbol", "chrom", "gene_start", "gene_end", "gene_strand",
                     "on_assembly", "dlaeve_protein", "evidence", "is_rbh", "tree_ok",
                     "diamond_bits", "pfam_dbd", "dbd_from", "dbd_to", "dbd_aln_cols",
                     "pid_dbd_overlap", "pid_dbd_L", "cisbp_threshold",
                     "weirauch_threshold", "passes_cisbp", "dbd_seq_dlaeve"))
setorder(guide, tf_name, motif_id, locus)
fw(guide, file.path(DAT, "tf_dlaeve_guide.tsv"))
# companion FASTA: the full protein sequence of every kept D. laeve TF
hdr <- guide[, .(lab = sprintf("locus=%s symbol=%s jaspar=%s", locus[1],
                 fifelse(is.na(symbol[1]), locus[1], symbol[1]),
                 paste(sort(unique(tf_name)), collapse = ";"))), by = dlaeve_protein]
fa_out <- prot_all[hdr$dlaeve_protein]
names(fa_out) <- paste(hdr$dlaeve_protein, hdr$lab)
fa_tmp <- file.path(DAT, "tf_dlaeve_orthologs.faa.part")
writeXStringSet(fa_out, fa_tmp)
stopifnot(file.rename(fa_tmp, file.path(DAT, "tf_dlaeve_orthologs.faa")))
cat(sprintf("  GUIDE: %d motif-locus rows, %d loci, %d motifs; %d locus/loci off chr1-31+mt\n",
            nrow(guide), uniqueN(guide$locus), uniqueN(guide$motif_id),
            uniqueN(guide[on_assembly == FALSE, locus])))

# ---- [10] Figures ------------------------------------------------------------
# funnel, per-family conservation and decision-grid FIGURES are retired, but their
# TABLES are still written (counts feed the diagram; the grid records the decisions).
cat("[TF 10] figures\n")
unlink(list.files(FIGS, pattern = "^(figS2b_|figS_tf_)", full.names = TRUE))   # TF stems ONLY: 01_genome_toolkit's other supplementary figures must survive
mot_with <- function(a) uniqueN(m2a[acc %in% a, matrix_id])
fun <- data.table(
  step = 1:8,
  gate = c("All animal JASPAR CORE motifs", "JASPAR names the protein (UniProt)",
           "Protein sequence retrieved", "Protein carries a Pfam domain",
           "A similar D. laeve protein exists", "Both carry the SAME expected DBD",
           sprintf("DBD similar enough (%s)", if (ORTH_IDENTITY == "cisbp") "Cis-BP 3.10" else "Weirauch 2014"),
           sprintf("...and %s", if (ORTH_HOMOLOGY == "rbh") "reciprocal best hit" else "RBH or gene tree")),
  n = c(N_MATRIX, uniqueN(m2a[nzchar(acc), matrix_id]), mot_with(names(bait_aa)),
        mot_with(unique(bd$query)), mot_with(unique(cand$acc)), mot_with(unique(lic$acc)),
        mot_with(unique(gate(lic, "any", ORTH_IDENTITY)$acc)), N_KEPT))
fun[, gate := factor(gate, levels = rev(gate))]
fw(fun, file.path(DAT, "tf_orthology_funnel.tsv"))

famf <- lic[, .(n_pairs = .N, med = median(pid_dbd_L), lo = quantile(pid_dbd_L, .25),
                hi = quantile(pid_dbd_L, .75)), by = .(pfam_dbd = pfam)]
famf <- merge(famf, thr[, .(pfam_dbd, t_cis, t_w14)], by = "pfam_dbd", all.x = TRUE)
famf <- merge(famf, pd[, .(n_dlaeve_loci = uniqueN(locus)), by = .(pfam_dbd = pfam)], by = "pfam_dbd", all.x = TRUE)
famf <- merge(famf, bd[, .(n_jaspar = uniqueN(query)), by = .(pfam_dbd = pfam)], by = "pfam_dbd", all.x = TRUE)
setorder(famf, med); famf[, pfam_dbd := factor(pfam_dbd, levels = pfam_dbd)]
fw(famf, file.path(DAT, "dbd_conservation_by_family.tsv"))

HOM <- c(any = "Any similar protein\n(domain licensed)",
         rbh_or_tree = sprintf("Reciprocal best hit OR gene tree\n(clade <= %d, support >= %d)", TREE_MAX_CLADE, TREE_MIN_SUPPORT),
         rbh = "Reciprocal best hit only")
RUL <- c(none = "No identity gate", w2014 = "Weirauch 2014", cisbp = "Cis-BP 3.10")
gr <- CJ(hom = names(HOM), rule = names(RUL), sorted = FALSE)
gr[, n := mapply(function(h, r) mot_with(unique(gate(lic, h, r)$acc)), hom, rule)]
gr[, `:=`(hom_l = factor(HOM[hom], levels = rev(HOM)), rule_l = factor(RUL[rule], levels = RUL),
          defensible = hom != "any" & rule != "none",
          chosen = hom == ORTH_HOMOLOGY & rule == ORTH_IDENTITY)]
fw(gr[, .(homology = hom, identity_rule = rule, motifs = n, defensible, chosen)],
   file.path(DAT, "pwm_library_decision_grid.tsv"))

# ---- [10a] Figure: the TF DBD family complement ------------------------------
# accessory bait domains were never scanned and stay in tf_family_complement.tsv
# (the retired first version drew them as if they were absent TF families). Count
# LOCI, not proteins: a protein count inflates every family by isoform multiplicity.
comp <- pd[, .(n_loci = uniqueN(locus)), by = .(pfam_dbd = pfam)]
comp <- merge(comp, bd[, .(n_jaspar = uniqueN(query)), by = .(pfam_dbd = pfam)], by = "pfam_dbd", all = TRUE)
comp[is.na(n_loci), n_loci := 0L][is.na(n_jaspar), n_jaspar := 0L]
comp[, class_expected_dbd := pfam_dbd %in% DBD_NAMES]
fw(comp[order(-n_loci)], file.path(DAT, "tf_family_complement.tsv"))
dbd_comp <- data.table(pfam_dbd = DBD_NAMES)
dbd_comp <- merge(dbd_comp, comp[, .(pfam_dbd, n_loci, n_jaspar)], by = "pfam_dbd", all.x = TRUE)
dbd_comp[is.na(n_loci), n_loci := 0L][is.na(n_jaspar), n_jaspar := 0L]
dbd_comp[, status := fifelse(n_loci == 0L, "Not detected in the proteome",
                     fifelse(n_jaspar > 0L, "Present, JASPAR motif available",
                                            "Present, no JASPAR motif"))]
setorder(dbd_comp, n_loci)
dbd_comp[, pfam_dbd := factor(pfam_dbd, levels = pfam_dbd)]
n_present <- dbd_comp[n_loci > 0, .N]
# A log axis has no zero, so absent families sit at a marked position left of 1,
# drawn as open circles and labelled 0.
dbd_comp[, x_plot := fifelse(n_loci == 0L, 0.72, as.numeric(n_loci))]
figD <- ggplot(dbd_comp, aes(y = pfam_dbd)) +
  geom_segment(data = dbd_comp[n_loci > 0], aes(x = 1, xend = n_loci, yend = pfam_dbd),
               linewidth = 0.35, colour = "grey85") +
  geom_point(data = dbd_comp[n_loci > 0], aes(x = n_loci, colour = status), size = 1.8) +
  geom_point(data = dbd_comp[n_loci == 0L], aes(x = x_plot, colour = status),
             size = 1.8, shape = 1, stroke = 0.7) +
  geom_text(aes(x = x_plot, label = n_loci), hjust = -0.5, size = 2.4, colour = COL_MARK) +
  scale_colour_manual(values = c("Present, JASPAR motif available" = COL_MARK,
                                 "Present, no JASPAR motif" = COL_B,
                                 "Not detected in the proteome" = "grey55"), name = NULL) +
  scale_x_log10(expand = expansion(mult = c(0.03, 0.12)),
                breaks = c(1, 3, 10, 30, 100, 300), labels = c("1", "3", "10", "30", "100", "300"),
                name = "D. laeve loci carrying the domain (log scale)") +
  labs(y = NULL, title = sprintf("TF DNA-binding-domain families in D. laeve: %d of %d present",
                                 n_present, nrow(dbd_comp)),
       subtitle = "One row per class-expected Pfam DBD family scanned across the proteome.\nLoci, not proteins: isoforms are collapsed. Open circles: family not detected.") +
  theme_pub() +
  theme(legend.position = "top", legend.justification = "left",
        axis.text.y = element_text(size = 6.8), panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.25),
        plot.subtitle = element_text(size = 7.6, colour = "grey30"))
save_supp(figD, "figS_tf_tf_family_complement", 7.0, 9.6)
cat(sprintf("  DBD family complement: %d of %d class-expected families present, %d TF loci total\n",
            n_present, nrow(dbd_comp), uniqueN(pd$locus)))

# ---- [10b] Figure: the method as a workflow diagram --------------------------
# What each stage asks, which tool answers it, how many motifs survive. Counts come
# from `fun`, so the diagram cannot drift from the numbers this run actually produced.
wf <- data.table(
  step  = 1:10,
  title = c("JASPAR CORE, four animal groups",
            "Motif to UniProt accession",
            "Accession to protein sequence",
            "Homology search, both directions",
            "Isoforms collapsed to locus; reciprocal best hits",
            "Pfam domains on both proteins",
            "DBD licensing: same class-expected domain",
            "DBD amino acid identity",
            "Motif-transfer threshold, per family",
            "The motif library"),
  tool  = c("JASPAR2024 sqlite, direct SQL", "MATRIX_PROTEIN table",
            "UniProt REST /uniprotkb/accessions", "DIAMOND blastp --ultra-sensitive",
            "best hit per locus, both legs required", "hmmscan and hmmsearch, Pfam-A --cut_ga",
            "JASPAR structural class to Pfam map", "MAFFT, identities / longer domain",
            sprintf("%s thresholds", if (ORTH_IDENTITY == "cisbp") "Cis-BP 3.10" else "Weirauch 2014"),
            sprintf("%s + %s", ORTH_HOMOLOGY, ORTH_IDENTITY)),
  n     = c(fun$n[1], fun$n[2], fun$n[3], NA, fun$n[5], fun$n[4], fun$n[6], NA, fun$n[7], fun$n[8]))
# leaves the background TRANSPARENT, which many viewers paint black: set white explicitly.
wf[, yc := -as.numeric(step)]
wf[, lab := ifelse(is.na(n), "", format(n, big.mark = ","))]
wf[, fill := fifelse(step == 1L, "grey85", fifelse(step == 10L, "#CDEBDF", "grey93"))]
figE <- ggplot(wf) +
  geom_rect(aes(xmin = 0, xmax = 6.6, ymin = yc - 0.34, ymax = yc + 0.34, fill = fill), colour = NA) +
  scale_fill_identity() +
  geom_segment(data = wf[step < 10], aes(x = 3.3, xend = 3.3, y = yc - 0.34, yend = yc - 0.66),
               colour = "grey40", linewidth = 0.5,
               arrow = arrow(length = unit(0.07, "in"), type = "closed")) +
  geom_text(aes(x = 0.28, y = yc + 0.115, label = paste0(step, ".  ", title)),
            hjust = 0, size = 2.9, colour = COL_MARK) +
  geom_text(aes(x = 0.28, y = yc - 0.145, label = tool),
            hjust = 0, size = 2.3, colour = "grey35", fontface = "italic") +
  geom_text(aes(x = 7.9, y = yc, label = lab), hjust = 1, size = 2.9, colour = COL_MARK) +
  annotate("text", x = 7.9, y = -0.42, label = "motifs in play", hjust = 1, size = 2.4, colour = "grey45") +
  scale_x_continuous(limits = c(-0.05, 8.0), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(add = c(0.3, 0.75))) +
  labs(title = "How the D. laeve motif library is built",
       subtitle = "01_genome_toolkit, TF annotation. Each stage is one filter; a JASPAR motif reaches the library only if every answer is yes.") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 10, hjust = 0),
        plot.subtitle = element_text(size = 7.8, colour = "grey30", hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 10, 10, 10))
save_supp(figE, "figS_tf_tf_annotation_workflow", 6.6, 5.6)

cat("[01_genome_toolkit TF-annotation section] done\n")

})

# ---- [9] Reproducibility: record the exact package versions this run used ----
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_01_genome_toolkit.txt"))
cat("[01_genome_toolkit] done\n")
