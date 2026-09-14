#!/usr/bin/env Rscript
set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges)
  library(Biostrings); library(bsseq); library(ggplot2); library(scales); library(patchwork)
})
# DESeq2, clusterProfiler, org.Hs.eg.db, DSS and svglite are attached later, where first used.

# ---- [0] Setup: pinned versions, paths, palette, figure helpers -------------
# Built and validated under R 4.4.1 / Bioconductor 3.20 (Fenix, /opt/apps/r/4.4.1-studio).
#   data.table 1.18.4  GenomicRanges 1.58.0  IRanges 2.40.1     Biostrings 2.74.1
#   bsseq 1.42.0         ggplot2 4.0.2         scales 1.4.0       patchwork 1.3.2
#   DESeq2 1.46.0        clusterProfiler 4.14.6  org.Hs.eg.db 3.20.0  DSS 2.54.0  svglite 2.2.2
# A full machine-checkable record is written to sessionInfo_03_promoters.txt at the end of the run.

PIPE <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
B01 <- file.path(PIPE, "01_genome_toolkit/objects"); B02 <- file.path(PIPE, "02_landscape/objects")
BATCH <- file.path(PIPE, "03_promoters")
DAT <- file.path(BATCH, "data"); FIGM <- file.path(BATCH, "figures/main")
FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(DAT, FIGM, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")

COL_WEBER <- c(HCP = "#117733", ICP = "#88CCEE", LCP = "#CC6677")
# Plot-only x-limit (Weber Fig 2 stops at 1.2; <0.3% of promoters exceed it).
# Classification and every TSV use the FULL data; drop_note() reports each drop.
OE_MAX <- 1.2
drop_note <- function(x, nm) {
  n <- sum(x > OE_MAX, na.rm = TRUE)
  cat(sprintf("  %s: %d/%d promoters (%.2f%%) above O/E %.1f not drawn\n",
              nm, n, length(x), 100*n/length(x), OE_MAX))
}
theme_pub <- function() theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"))
save_fig <- function(p, name, w, h) {
  # cairo_pdf so the β glyph in axis labels renders in the PDF
  ggsave(file.path(FIGM, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGM, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGM, paste0(name, ".svg")), p, width = w, height = h)   # vector (svg)
  cat(sprintf("  saved %s\n", name))
}
# same as save_fig but into figures/supplementary
save_supp <- function(p, name, w, h) {
  ggsave(file.path(FIGS, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGS, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGS, paste0(name, ".svg")), p, width = w, height = h)   # vector (svg)
  cat(sprintf("  saved %s\n", name))
}

# ---- [1] Load genome / GFF / bsseq (01_genome_toolkit + 02_landscape chrmt objects) --------
cat("[1] load genome / gff / bsseq (from lower batches)\n")
genome <- readRDS(file.path(B01, "genome_chrmt.rds"))
gff    <- readRDS(file.path(B01, "gff_chrmt.rds"))
bs     <- readRDS(file.path(B02, "bsseq_cov5_chrmt.rds"))
chr_len <- setNames(width(genome), names(genome))

# ---- [2] Weber promoter classification (sliding 500-bp, -1300/+200) ---------
# weber_sliding(): scan each promoter with 500-bp windows at 5-bp offset; return
# class + max window O/E + whole-promoter O/E. A single whole-window O/E averages
# a local CpG island away; the sliding max recovers it (the faithful Weber method).
weber_sliding <- function(gen, dt, w = 500L, off = 5L) {
  cls <- character(nrow(dt)); mxoe <- rep(NA_real_, nrow(dt)); whoe <- rep(NA_real_, nrow(dt))
  for (ac in intersect(unique(dt$seqid), names(gen))) {
    idx <- which(dt$seqid == ac); L <- length(gen[[ac]])
    v <- Views(gen[[ac]], start = pmax(1L, dt$ps[idx]), end = pmin(L, dt$pe[idx]))
    for (m in seq_along(idx)) {
      s <- v[[m]]; Ls <- length(s); if (Ls < w) next
      cw <- letterFrequencyInSlidingView(s, w, "C")[, 1]
      gw <- letterFrequencyInSlidingView(s, w, "G")[, 1]
      nwin <- Ls - w + 1L
      cg <- integer(Ls); st <- start(matchPattern("CG", s, fixed = TRUE)); if (length(st)) cg[st] <- 1L
      # Weber defines TWO O/E quantities: the sliding-window ratio (classification)
      # and the whole-promoter ratio -- the x-axis of Weber Fig 2 (their HCP histogram
      # extends BELOW 0.75, impossible for a sliding max). Classify max_oe; plot whole_oe.
      fr <- letterFrequency(s, c("C", "G"))
      whoe[idx[m]] <- (sum(cg) * Ls) / max(fr[["C"]] * fr[["G"]], 1)
      # window i spans [i, i+w-1] -> CpGs = ccg[i+w-1] - ccg[i-1]; use ccg[w:Ls], NOT
      # ccg[w:(Ls-1)], which recycles and corrupts the last window into a negative count
      ccg <- cumsum(cg); cgw <- ccg[w:Ls] - c(0, ccg[1:(nwin - 1)])
      oe <- (cgw * w) / pmax(cw * gw, 1); gc <- 100 * (cw + gw) / w
      sel <- seq(1L, nwin, by = off); oe <- oe[sel]; gc <- gc[sel]
      mxoe[idx[m]] <- max(oe, 0, na.rm = TRUE)
      cls[idx[m]]  <- if (any(oe > 0.75 & gc > 55)) "HCP" else if (!any(oe > 0.48)) "LCP" else "ICP"
    }
  }
  data.table(max_oe = mxoe, whole_oe = whoe,
             weber_class = factor(cls, levels = c("HCP", "ICP", "LCP")))
}

cat("[2] promoter Weber classification (sliding 500-bp, -1300/+200)\n")
gene_gr <- gff[gff$type == "gene"]
gid <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
tss <- ifelse(as.character(strand(gene_gr)) == "-", end(gene_gr), start(gene_gr))
neg <- as.character(strand(gene_gr)) == "-"
prom <- GRanges(seqnames(gene_gr),
                IRanges(ifelse(neg, tss - 200L, tss - 1300L), ifelse(neg, tss + 1300L, tss + 200L)),
                strand = strand(gene_gr), gene_id = gid,
                biotype = as.character(gene_gr$gene_biotype), tss = tss)
prom <- GenomeInfoDb::keepSeqlevels(prom, keep_chr, pruning.mode = "coarse")
seqlengths(prom) <- chr_len[keep_chr]
prom <- trim(prom); prom <- prom[width(prom) >= 500]           # need >=500 bp to slide
pdt <- as.data.table(prom)[, .(seqnames = as.character(seqnames), start, end,
                               strand = as.character(strand), gene_id, biotype, tss)]
pdt[, `:=`(seqid = seqnames, ps = start, pe = end)]
pdt <- cbind(pdt, weber_sliding(genome, pdt))                  # adds max_oe + whole_oe + weber_class
pdt <- pdt[is.finite(max_oe) & is.finite(whole_oe)]
# MAIN analysis (fig3a-g) is PROTEIN-CODING ONLY, matching Weber 2007's validated
# protein-coding set; lncRNA adds CpG-poor promoters that pile into LCP and
# flatten the human bimodality. pdt_all keeps all biotypes -- NOTE: no code below
# writes the by-biotype supplementary any more (the figS3_*_by_biotype files on
pdt_all <- copy(pdt)
pdt <- pdt[biotype == "protein_coding"]
class_n <- pdt[, .N, by = weber_class][order(weber_class)][, pct := 100*N/sum(N)][]
fwrite(pdt[, .(gene_id, max_oe, whole_oe, weber_class, biotype)],
       file.path(DAT, "promoter_weber_classification.tsv"), sep = "\t")
fwrite(class_n, file.path(DAT, "weber_class_composition.tsv"), sep = "\t")
cat(sprintf("  %d promoters: HCP %d / ICP %d / LCP %d\n", nrow(pdt),
            class_n[weber_class=="HCP", N], class_n[weber_class=="ICP", N],
            class_n[weber_class=="LCP", N]))

# -- fig3a: combined Weber histogram (whole_oe, protein coding) --
drop_note(pdt$whole_oe, "fig3a/b D. laeve")
pa <- ggplot(pdt[whole_oe <= OE_MAX], aes(whole_oe, fill = weber_class)) +
  geom_histogram(bins = 80, colour = "white", linewidth = 0.05) +
  scale_fill_manual(values = COL_WEBER, name = "Promoter class (Weber 2007)") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  # short x label: the full window string clips at this panel width; window in subtitle
  labs(x = "CpG observed/expected (entire promoter)", y = "Promoters",
       title = "Promoter CpG content of protein-coding genes",
       subtitle = sprintf("D. laeve, n = %s promoters, -1300/+200 window", comma(nrow(pdt)))) + theme_pub()
save_fig(pa, "fig3a_promoter_cpg_oe_histogram", 4.6, 3.0)

# -- fig3b: per-class faceted histogram --
nb3 <- pdt[, .N, by = weber_class]
pdt[, class_lab := factor(sprintf("%s (n = %s)", weber_class, comma(nb3$N[match(weber_class, nb3$weber_class)])),
                          levels = sprintf("%s (n = %s)", levels(weber_class),
                                           comma(nb3$N[match(levels(weber_class), nb3$weber_class)])))]
pb <- ggplot(pdt[whole_oe <= OE_MAX], aes(whole_oe, fill = weber_class)) +
  geom_histogram(binwidth = 0.05, boundary = 0, colour = "white", linewidth = 0.1) +
  facet_wrap(~ class_lab, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = COL_WEBER, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "CpG observed/expected (entire promoter)", y = "Count",
       title = "CpG O/E distribution by promoter class",
       subtitle = sprintf("D. laeve protein-coding genes, n = %s promoters", comma(nrow(pdt)))) +
  theme_pub() + theme(strip.background = element_blank())
save_fig(pb, "fig3b_promoter_oe_distribution_by_class", 4.0, 4.7)

# ---- [3] Pooled promoter methylation (all 4 WGBS samples) + fig3c -----------
# per-gene promoter beta, C1+C2+A1+A2 pooled; writes
# promoter_methylation_per_gene.tsv (reused by [7], [9] and 07_wgcna section 5).
chrs <- as.character(seqnames(bs)); bs <- bs[chrs %in% keep_chr, ]
gr <- granges(bs)
M_all <- rowSums(as.matrix(getCoverage(bs, type = "M")))
C_all <- rowSums(as.matrix(getCoverage(bs, type = "Cov")))
cpg_all <- data.table(chr = as.character(seqnames(gr)), start = start(gr), end = start(gr),
                      M = M_all, C = C_all)
setkey(cpg_all, chr, start, end)
prom_pool <- pdt[, .(chr = seqnames, start, end, gene_id)]
setkey(prom_pool, chr, start, end)
ovp <- foverlaps(cpg_all, prom_pool, nomatch = 0L)
prom_meth <- ovp[, .(mean_meth = sum(M)/pmax(sum(C),1), n_cpg = .N), by = gene_id]
fwrite(prom_meth, file.path(DAT, "promoter_methylation_per_gene.tsv"), sep = "\t")

# -- fig3c: HCP promoters ranked by pooled methylation --
g_id <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
g_note <- vapply(mcols(gene_gr)$Note, function(x) if (length(x)) as.character(x)[1] else NA_character_, character(1))
parse_name <- function(note, fb) {
  m <- regmatches(note, regexec("^Similar to ([^:]+):", note))
  vapply(seq_along(m), function(i) if (length(m[[i]]) >= 2 && nchar(m[[i]][2]) > 0) m[[i]][2] else fb[i], character(1))
}
hcp <- pdt[weber_class == "HCP"]
hcp[, gene_name := parse_name(g_note[match(gene_id, g_id)], gene_id)]
hcp <- merge(hcp, prom_meth, by = "gene_id", all.x = TRUE)
hcp <- hcp[!is.na(mean_meth)][order(-mean_meth)]
# ones with their LOC id, never merge. Also keeps `label` unique --
# factor(levels = rev(x)) dies on duplicates (the trap that broke fig3c once).
dup_sym <- duplicated(hcp$gene_name) | duplicated(hcp$gene_name, fromLast = TRUE)
hcp[dup_sym, gene_name := paste0(gene_name, " [", gene_id, "]")]
hcp[, label := paste0(gene_name, " (", round(100*mean_meth, 1), "%)")]
# Plot the top N only (one row per covered HCP gene is an unplottable axis); every
# gene stays in promoter_methylation_per_gene.tsv -- the figure shows the methylated extreme.
topn <- 40L
hcp_plot <- head(hcp, topn)
stopifnot(!any(duplicated(hcp_plot$label)))       # factor(levels=) needs unique levels
hcp_plot[, y := factor(label, levels = rev(label))]
cat(sprintf("  fig3c: %d HCP promoters with coverage; plotting top %d\n",
            nrow(hcp), nrow(hcp_plot)))
pc <- ggplot(hcp_plot, aes(100*mean_meth, y, fill = 100*mean_meth)) +
  geom_segment(aes(x = 0, xend = 100*mean_meth, yend = y), colour = "grey60", linewidth = 0.3) +
  geom_point(shape = 21, size = 3, colour = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#FDF2E9", high = "#C0392B", name = "β (%)", limits = c(0, 100)) +
  labs(x = "Mean promoter β (%)", y = NULL,
       title = sprintf("Top %d HCP promoters ranked by methylation", nrow(hcp_plot)),
       subtitle = sprintf("of %s HCP promoters with >= 1 CpG covered (cov >= 5 in all 4 samples); full list in promoter_methylation_per_gene.tsv",
                          comma(nrow(hcp)))) +
  theme_pub() + theme(axis.text.y = element_text(size = 6.5))
save_fig(pc, "fig3c_hcp_gene_list", 5.8, max(4.0, 0.18*nrow(hcp_plot) + 1))

# ---- [4] fig3d: TSS metagene by promoter class (all 4 samples, pooled) ------
M <- as.matrix(getCoverage(bs, type = "M")); Cv <- as.matrix(getCoverage(bs, type = "Cov"))
M_pool <- rowSums(M); C_pool <- rowSums(Cv)              # pool all 4 samples (C1,C2,A1,A2)
cpg <- data.table(chr = as.character(seqnames(gr)), pos = start(gr),
                  bg = M_pool/pmax(C_pool,1), cov = C_pool)
cpg[, `:=`(start = pos, end = pos)]
pm <- pdt[, .(chr = seqnames, tss = tss, strand, gene_id, weber_class)]
flank <- 2000L; nbin <- 50L; bw <- (2*flank)/nbin
pm[, `:=`(ws = tss - flank, we = tss + flank)]
setkey(pm, chr, ws, we); setkey(cpg, chr, start, end)
ov <- foverlaps(cpg, pm, by.x = c("chr","start","end"), by.y = c("chr","ws","we"), nomatch = 0L)
ov[, rel := ifelse(strand == "-", tss - pos, pos - tss)]
ov[, bin := pmin(nbin, pmax(1L, as.integer((rel + flank) %/% bw) + 1L))]
ov[, bc := -flank + (bin - 0.5)*bw]
agg <- ov[, .(beta = sum(bg*cov)/pmax(sum(cov),1)), by = .(weber_class, bc)]
fwrite(agg, file.path(DAT, "metagene_by_promoter_class.tsv"), sep = "\t")
pd <- ggplot(agg, aes(bc/1000, beta*100, colour = weber_class)) +
  geom_vline(xintercept = 0, colour = "#C0392B", linetype = "dashed", linewidth = 0.45) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 0.6, alpha = 0.55) +
  scale_colour_manual(values = COL_WEBER, name = "Promoter class") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Distance from TSS (kb)", y = "Mean CpG methylation β (%)",
       title = "TSS methylation by promoter class") +
  theme_pub()
save_fig(pd, "fig3d_metagene_by_promoter_class", 3.4, 2.2)

# ---- [5] Human GRCh38 comparison (fig3e): same classifier, same window ------
# RefSeq. The dataset/ cache is pre-staged on a login node; defq compute nodes
# have NO internet, so the wget fallback must never fire there.
# (ls -l before trusting: fna 972,898,531 B, gff 77,807,126 B), and the URL is the
# UNVERSIONED GRCh38_latest -- a re-fetch could silently change the human control.
cat("[5] human GRCh38 promoter CpG O/E vs D. laeve (fig3e)\n")
DS <- file.path(BATCH, "dataset"); dir.create(DS, showWarnings = FALSE, recursive = TRUE)
hs_fna <- file.path(DS, "GRCh38_latest_genomic.fna.gz")
hs_gff <- file.path(DS, "GRCh38_latest_genomic.gff.gz")
url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers"
if (!file.exists(hs_fna)) system2("wget", c("-q","-c", file.path(url,"GRCh38_latest_genomic.fna.gz"), "-O", hs_fna))
if (!file.exists(hs_gff)) system2("wget", c("-q","-c", file.path(url,"GRCh38_latest_genomic.gff.gz"), "-O", hs_gff))
# a failed wget -O leaves a ZERO-BYTE file that passes file.exists() forever; the good
# cache is 972,898,531 / 77,807,126 B, so any truncation fails these floors loudly
stopifnot(file.size(hs_fna) > 5e8, file.size(hs_gff) > 5e7)

# primary chromosomes only: chr1-22 (NC_000001..NC_000022), X (NC_000023), Y (NC_000024)
hs_primary <- c(sprintf("NC_0000%02d", 1:22), "NC_000023", "NC_000024")
# human genes -> TSS -> Weber -1300/+200 region; sliding-500bp classification
hg <- fread(cmd = sprintf("zcat '%s' | grep -v '^#'", hs_gff), sep = "\t", header = FALSE, quote = "",
            col.names = c("seqid","src","type","start","end","score","strand","phase","attr"))
# type=="gene" alone silently drops ~16.9k human pseudogenes, while D. laeve
# pseudogenes arrive as type=="gene". Keep both for cross-species parity (the
# main fig3a-e set is protein_coding-only, so this does not touch it).
hg <- hg[type %in% c("gene", "pseudogene") & sub("\\..*", "", seqid) %in% hs_primary]
hg[, biotype := sub(".*gene_biotype=([^;]+).*", "\\1", attr)]
hg[, tss := ifelse(strand == "-", end, start)]; hneg <- hg$strand == "-"
hg[, `:=`(ps = ifelse(hneg, tss - 200L, tss - 1300L), pe = ifelse(hneg, tss + 1300L, tss + 200L))]
cat(sprintf("  %d human loci on primary chromosomes (%d gene + %d pseudogene)\n",
            nrow(hg), sum(hg$type == "gene"), sum(hg$type == "pseudogene")))

hs_genome <- readDNAStringSet(hs_fna)
names(hs_genome) <- sub(" .*", "", names(hs_genome))          # header -> accession only
hs_genome <- hs_genome[sub("\\..*", "", names(hs_genome)) %in% hs_primary]  # drop scaffolds/alts
hoe <- cbind(hg[, .(biotype)], weber_sliding(hs_genome, hg))  # max_oe + whole_oe + weber_class
hoe <- hoe[is.finite(max_oe) & is.finite(whole_oe)]
hoe_all <- copy(hoe)                                    # all biotypes; currently unused (by-biotype supp not written)
hoe <- hoe[biotype == "protein_coding"]   # main: protein-coding only (Weber's set)
fwrite(hoe[, .(max_oe, whole_oe, weber_class)], file.path(DAT, "human_promoter_cpg_oe.tsv"), sep = "\t")
# Replication check vs Weber 2007 Fig 2 (n=15,609 hg17 validated promoters:
# HCP 64% / ICP 13% / LCP 23%); our GRCh38 RefSeq n differs, split should land close.
hpct <- function(cl) 100 * mean(hoe$weber_class == cl)
cat(sprintf("  %d human promoters -> HCP %.0f%% / ICP %.0f%% / LCP %.0f%%  (Weber 2007: 64/13/23)\n",
            nrow(hoe), hpct("HCP"), hpct("ICP"), hpct("LCP")))
cat(sprintf("  D. laeve            -> HCP %.0f%% / ICP %.0f%% / LCP %.0f%%\n",
            100*mean(pdt$weber_class == "HCP"), 100*mean(pdt$weber_class == "ICP"),
            100*mean(pdt$weber_class == "LCP")))

# -- fig3e: two-species histograms + class pies (Weber 2007 Fig 3c style) --
# protein-coding only; D. laeve left, human right
cmp <- rbind(data.table(species = "D. laeve",   whole_oe = pdt$whole_oe, weber_class = pdt$weber_class),
             data.table(species = "H. sapiens", whole_oe = hoe$whole_oe, weber_class = hoe$weber_class))
cmp[, species := factor(species, levels = c("D. laeve", "H. sapiens"))]     # slug left, human right
drop_note(cmp[species == "D. laeve", whole_oe],   "fig3e D. laeve")
drop_note(cmp[species == "H. sapiens", whole_oe], "fig3e H. sapiens")
# facet label carries n per species
n_sp <- cmp[, .N, by = species]
sp_lab <- setNames(sprintf("%s (n = %s)", levels(cmp$species),
                           comma(n_sp$N[match(levels(cmp$species), n_sp$species)])), levels(cmp$species))
cmp[, species_lab := factor(sp_lab[as.character(species)], levels = sp_lab)]
# they threshold the sliding-window ratio, not this axis -- a vline would imply a
# rule that does not apply (LOCKED convention, see header RULES).
pe_hist <- ggplot(cmp[whole_oe <= OE_MAX], aes(whole_oe, fill = weber_class)) +
  geom_histogram(bins = 80, colour = "white", linewidth = 0.05) +
  facet_wrap(~ species_lab, scales = "free_y") +
  scale_fill_manual(values = COL_WEBER, name = "Promoter class (Weber 2007)") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = seq(0, 1.2, 0.4)) +
  labs(x = "Promoter CpG observed/expected (entire promoter, -1300/+200)", y = "Promoters",
       title = "Promoter CpG content of protein-coding genes: D. laeve vs human") +
  theme_pub() + theme(strip.background = element_blank(),
                      plot.margin = margin(5.5, 5.5, 5.5, 10))
mk_pie <- function(sp) {
  # pies use ALL promoters of that species (not the OE_MAX-trimmed plotting subset),
  # so the percentages match the class split printed above. n is on the facet strip.
  tab <- cmp[species == sp, .N, by = weber_class][, frac := N/sum(N)][order(weber_class)]
  # a <8% wedge is too thin to hold its label: push those outside, and stagger the
  # radius so two ADJACENT small wedges (D. laeve HCP 4% + LCP 4%) cannot collide
  tab[, rad := 1]
  tab[frac >= 0.08 & frac < 0.30, rad := 1.35]
  tab[frac < 0.08, rad := ifelse(seq_len(.N) %% 2 == 1, 1.6, 2.0)]
  ggplot(tab, aes(x = 1, y = frac, fill = weber_class)) +
    geom_col(width = 1, colour = "white", linewidth = 0.3) + coord_polar(theta = "y") +
    geom_text(aes(x = rad, label = sprintf("%s %.0f%%", weber_class, 100*frac)),
              position = position_stack(vjust = 0.5), size = 2.8) +
    scale_x_continuous(limits = c(0.5, 2.3)) +
    scale_fill_manual(values = COL_WEBER, guide = "none") +
    labs(title = sp) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 9.5, face = "italic"))
}
# the wider canvas also un-clips the title, the facet strips and the x-axis label.
pe <- (pe_hist | (mk_pie("D. laeve") / mk_pie("H. sapiens"))) + plot_layout(widths = c(2, 1))
save_fig(pe, "fig3e_promoter_cpg_oe_human_vs_dlaeve", 7.2, 3.6)

# ---- [6] REMOVED: figS3_weber_tss1kb TSS+-1kb sensitivity ----

# ---- [7] fig3f: expression vs promoter methylation, per Weber class ---------
# Weber 2007 Fig 4e: promoter activity + hypermethylation incompatible for HCP
# and ICP but NOT LCP. NOT CIRCULAR: the class comes from SEQUENCE alone (O/E + GC).
# SCOPE: a CONSTITUTIVE cross-gene correlation -- methylation pooled over the 4
# WGBS samples, expression = mean of ALL 7 tail libraries (unpaired animals). It
# says NOTHING about the amputation response and does not contradict 07_wgcna's
# "DMR direction carries no expression information" (a differential claim).
cat("[7] fig3f expression vs promoter methylation, per Weber class\n")
suppressPackageStartupMessages(library(DESeq2))
HTSEQ <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"
# sample map rebuilt from filenames, no external metadata (matches 01_genome_toolkit §6):
# tail control C1S1..C4S4 ; tail amputated T2S6/T3S7/T4S8 (T1S5 excluded at gene level)
rna <- data.table(sample = c("C1S1","C2S2","C3S3","C4S4","T2S6","T3S7","T4S8"),
                  condition = c(rep("Control", 4), rep("Amputated", 3)))
rna[, file := file.path(HTSEQ, paste0(sample, "_htseq_gene_counts.txt"))]
stopifnot(all(file.exists(rna$file)))
cl <- lapply(rna$file, function(f) {
  x <- fread(f, header = FALSE, col.names = c("gene_id", "count")); x[!startsWith(gene_id, "__")]
})
stopifnot(sapply(cl, function(x) identical(x$gene_id, cl[[1]]$gene_id)))  # HTSeq row-order guard
cmat <- sapply(cl, function(x) x$count)
rownames(cmat) <- cl[[1]]$gene_id; colnames(cmat) <- rna$sample
cmat <- cmat[rownames(cmat) %in% pdt$gene_id, , drop = FALSE]   # protein-coding, keep_chr universe
cmat <- cmat[rowSums(cmat >= 5) >= 2, , drop = FALSE]           # expressed (same rule as 01_genome_toolkit)
dds  <- DESeqDataSetFromMatrix(cmat, data.frame(row.names = rna$sample), design = ~ 1)
nmat <- counts(estimateSizeFactors(dds), normalized = TRUE)     # median-of-ratios
expr <- data.table(gene_id = rownames(nmat), mean_expr = rowMeans(nmat))  # mean of ctrl AND amp

me <- merge(merge(pdt[, .(gene_id, weber_class)], prom_meth[, .(gene_id, mean_meth, n_cpg)],
                  by = "gene_id"), expr, by = "gene_id")
me[, `:=`(beta = 100 * mean_meth, log_expr = log2(mean_expr + 1))]
fwrite(me, file.path(DAT, "promoter_meth_vs_expression.tsv"), sep = "\t")

# Pearson (as requested) per class; Spearman reported alongside because the
# methylation-expression relationship is monotonic but not necessarily linear.
rstat <- me[, {
  ct <- cor.test(beta, log_expr, method = "pearson")
  sp <- suppressWarnings(cor.test(beta, log_expr, method = "spearman", exact = FALSE))
  .(n = .N, r = unname(ct$estimate), p = ct$p.value, rho = unname(sp$estimate), p_rho = sp$p.value)
}, by = weber_class][order(weber_class)]
# STAT TEST: Pearson and Spearman correlation tests (cor.test) of promoter beta vs log2 expression,
# per Weber class, plus the genome-wide Spearman the Results quote
rall <- me[, { sp <- suppressWarnings(cor.test(beta, log_expr, method = "spearman", exact = FALSE))
               ct <- cor.test(beta, log_expr, method = "pearson")
               .(weber_class = "All", n = .N, r = unname(ct$estimate), p = ct$p.value,
                 rho = unname(sp$estimate), p_rho = sp$p.value) }]
rstat <- rbind(rstat, rall)
fwrite(rstat, file.path(DAT, "promoter_meth_expr_correlation.tsv"), sep = "\t")
cat("  Pearson r (promoter beta vs log2 expression), per Weber class:\n")
for (i in seq_len(nrow(rstat)))
  cat(sprintf("    %-3s n=%6s  r=%+.3f  p=%-10.3g  (Spearman rho=%+.3f)\n",
              rstat$weber_class[i], comma(rstat$n[i]), rstat$r[i], rstat$p[i], rstat$rho[i]))

lab <- rstat[weber_class != "All", .(weber_class,
                 lab = sprintf("rho = %+.2f (P = %s)\nr = %+.2f\nn = %s", rho,
                               format.pval(p_rho, digits = 2, eps = 1e-16), r, comma(n)))]   # Spearman = the statistic the text reports
pf <- ggplot(me, aes(beta, log_expr)) +
  geom_point(aes(colour = weber_class), alpha = 0.10, size = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.5) +
  facet_wrap(~ weber_class) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = lab), hjust = 1.06, vjust = 1.12,
            size = 2.5, lineheight = 0.95) +      # plain text, no box (paper rule)
  scale_colour_manual(values = COL_WEBER, guide = "none") +
  labs(x = "Mean promoter methylation β (%)",
       y = expression(log[2]*"(mean normalized counts + 1)"),
       title = "Gene expression vs promoter methylation, by promoter class",
       subtitle = "D. laeve protein-coding genes; β pooled over 4 WGBS samples, expression = mean of 7 tail libraries") +
  theme_pub() + theme(strip.background = element_blank())
save_fig(pf, "fig3f_promoter_meth_vs_expression", 7.5, 3.2)

# ---- [8] GO/KEGG over-representation by Weber class, per species ------------
# Classes tested SEPARATELY (class genes vs all classified genes of the species).
# D. laeve = STRING v12 terms via enricher (no OrgDb exists); human = org.Hs.eg.db
cat("[8] GO over-representation by Weber promoter class (D. laeve + human)\n")
suppressPackageStartupMessages(library(clusterProfiler))
GO_CATS <- c(BP = "Biological Process (Gene Ontology)", MF = "Molecular Function (Gene Ontology)",
             CC = "Cellular Component (Gene Ontology)",
             KEGG = "KEGG (Kyoto Encyclopedia of Genes and Genomes)")
# Dotplot x = FOLD ENRICHMENT, not -log10(FDR): the FDR axis is dominated by set
# size and ranks huge low-fold housekeeping terms above specific ones. Colour =
# raw BH FDR, size = gene count; Count>=5 drops tiny noise terms.
CLS_GO <- c("HCP", "ICP", "LCP")
go_dot <- function(dt) {                         # shared aesthetic for every GO dotplot
  # colour = raw BH FDR (dark red = most significant); legend reads actual FDRs
  ggplot(dt, aes(FoldEnrichment, lab, size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#7B241C", high = "#F5CBA7", name = "BH FDR",
                          limits = c(0, 0.05), guide = guide_colourbar(reverse = TRUE)) +
    scale_size_continuous(name = "Genes", range = c(2, 6)) +
    labs(x = "Fold enrichment (observed / expected)", y = NULL) +
    theme_pub() + theme(axis.text.y = element_text(size = 6))
}
# go_class_figs(): ONE combined figure per species, plotted classes stacked on a
go_class_figs <- function(godt, tag, name, to_main = FALSE, classes = CLS_GO) {
  sig <- godt[p.adjust < 0.05 & Count >= 5 & is.finite(FoldEnrichment)]
  sig[, weber_class := factor(weber_class, levels = CLS_GO)]
  nper <- vapply(CLS_GO, function(cl) nrow(sig[weber_class == cl]), integer(1))
  cat(sprintf("  %s enriched terms (Count>=5): HCP=%d ICP=%d LCP=%d\n", name, nper[1], nper[2], nper[3]))
  sig <- sig[as.character(weber_class) %in% classes]
  CLS_GO <- classes                      # local: drives facets, placeholders and levels
  sig[, weber_class := factor(as.character(weber_class), levels = CLS_GO)]
  # A plain top-8 by fold plots ZERO KEGG on the human HCP panel (the best KEGG fold
  # sits just below the 8th GO term): take the best KEGG first (<= 1/4 of each class's
  # rows), then fill by fold. `sig` is already FDR<0.05 & Count>=5; ties broken by p.adjust.
  NPER <- 8L
  pick_terms <- function(d) {
    nk   <- min(sum(d$ontology == "KEGG"), NPER %/% 4L)
    keg  <- if (nk > 0) d[ontology == "KEGG"][order(-FoldEnrichment, p.adjust)][seq_len(nk)] else d[0]
    rest <- d[!(ID %in% keg$ID)][order(-FoldEnrichment, p.adjust)][seq_len(min(NPER - nrow(keg), .N))]
    rbind(keg, rest)[order(-FoldEnrichment, p.adjust)]
  }
  comb <- sig[, pick_terms(.SD), by = weber_class]
  if (!nrow(comb)) return(invisible())
  # a class with NO enriched term would leave an empty facet, which crashes
  # facet_grid(space="free_y") with a non-finite height. Give each empty class ONE
  # placeholder row (NA fold -> no point) so its facet reads "no term at FDR<0.05".
  miss <- setdiff(CLS_GO, as.character(unique(comb$weber_class)))
  if (length(miss)) comb <- rbind(comb, data.table(weber_class = factor(miss, levels = CLS_GO),
                        Description = paste0(miss, ": no term at FDR<0.05"),
                        FoldEnrichment = NA_real_, Count = NA_integer_, p.adjust = NA_real_), fill = TRUE)
  comb[, weber_class := factor(as.character(weber_class), levels = CLS_GO)]
  # Append the source ontology to each label ("Lysosome (KEGG)"), matching the 05_differential
  if (!"ontology" %in% names(comb)) comb[, ontology := NA_character_]
  comb[, Description := ifelse(is.na(ontology), Description,
                               sprintf("%s (%s)", Description, ontology))]
  # Wrap long term labels: the human MF/CC names run to ~70 characters, and unwrapped
  wrap_w <- if (to_main) 40 else 52
  comb[, Description := vapply(Description, function(s)
    paste(strwrap(s, width = wrap_w), collapse = "\n"), character(1))]
  comb[, lab := factor(make.unique(Description), levels = rev(make.unique(Description)))]
  p <- go_dot(comb) +
    facet_grid(weber_class ~ ., scales = "free_y", space = "free_y", drop = FALSE) +
    labs(title = sprintf("%s: GO and KEGG by Weber classes", name)) +
    theme(strip.background = element_blank(),
          plot.title = element_text(size = rel(1)),
          plot.title.position = "plot",  # "panel" anchoring pushed the title right of the wide y-label block, off the canvas
          legend.key.height = unit(10, "pt"))
  # edge and at 3.3 in the human figure's label block left the panel zero-width
  h <- max(2.4, 4.0 * length(CLS_GO) / 3)
  if (to_main) save_fig(p, sprintf("fig3g_%s_weber_class_go", tag), 4.6, h)
  else         save_supp(p, sprintf("figS3_%s_weber_class_go", tag), 5.4, max(2.6, 0.30 * nrow(comb) + 1.4))
}

# --- D. laeve: STRING v12 route ---
STRING_ENR <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/STRING.protein.enrichment.terms.v12.0.txt"
# KEGG tables cached once under tools/kegg/ and read BY PATH (compute nodes have
# no outbound internet); refresh via the curl calls recorded in tools/kegg/.
KEGG_LINK  <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/kegg/hsa_pathway_link.tsv"
KEGG_NAMES <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/kegg/hsa_pathway_names.tsv"
dl_cls <- pdt[biotype == "protein_coding" & !is.na(weber_class), .(gene_id, weber_class)]
sterms <- fread(STRING_ENR, col.names = c("string_id", "category", "term", "description"))
sterms[, gene_id := sub("^[^.]+\\.", "", string_id)]
run_enr_dl <- function(sig, uni_genes, cat_label) {              # one class x one category
  sub <- sterms[category == cat_label]
  uni <- intersect(unique(sub$gene_id), uni_genes)
  hit <- intersect(sig, uni)
  if (length(hit) < 5) return(data.table())
  # STAT TEST: hypergeometric ORA (clusterProfiler enricher), BH-adjusted.
  # never the whole STRING table (the inflated-universe bug that voided 05_differential's DMR terms).
  res <- suppressWarnings(enricher(gene = hit, universe = uni, TERM2GENE = sub[, .(term, gene_id)],
           TERM2NAME = unique(sub[, .(term, description)]), pAdjustMethod = "BH",
           pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
  if (is.null(res)) return(data.table())
  as.data.table(as.data.frame(res))
}
uni_dl <- unique(dl_cls$gene_id)
dl_go <- rbindlist(lapply(c("HCP", "ICP", "LCP"), function(cl) {
  sig <- dl_cls[weber_class == cl]$gene_id
  rbindlist(lapply(names(GO_CATS), function(o) {
    d <- run_enr_dl(sig, uni_dl, GO_CATS[[o]]); if (nrow(d)) d[, `:=`(ontology = o, weber_class = cl)]; d
  }), fill = TRUE)
}), fill = TRUE)
fwrite(dl_go, file.path(DAT, "dlaeve_weber_class_go_enrichment.tsv"), sep = "\t")
cat(sprintf("  D. laeve GO/KEGG at FDR<0.05: %s\n",
            paste(sprintf("%s=%d", c("HCP","ICP","LCP"),
              vapply(c("HCP","ICP","LCP"), function(cl) nrow(dl_go[weber_class==cl & p.adjust<0.05]), integer(1))),
              collapse = " ")))
go_class_figs(dl_go, "dlaeve", "D. laeve", to_main = TRUE)

# --- human: org.Hs.eg.db route (no human STRING file exists) ---
hg_go <- tryCatch({
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  hgc <- hg[type == "gene" & biotype == "protein_coding"]
  hgc[, symbol := sub(".*;gene=([^;]+).*", "\\1", attr)]
  hgc <- hgc[grepl(";gene=", attr)]
  hwc <- weber_sliding(hs_genome, hgc)                          # classify the protein-coding set
  hgc[, weber_class := as.character(hwc$weber_class)]
  hgc <- hgc[!is.na(weber_class) & !is.na(symbol) & symbol != ""]
  uni_h <- unique(hgc$symbol)
  cat(sprintf("  human classified protein-coding genes for GO: %d\n", length(uni_h)))
  # enrichGO ont="ALL" gives BP/MF/CC (ONTOLOGY renamed to `ontology` to match the
  # D. laeve table); KEGG comes from the CACHED tools/kegg tables via enricher, which
  # needs ENTREZ ids (mapped through org.Hs.eg.db). The KEGG arm is wrapped so any
  # failure degrades to GO-only instead of killing the batch.
  rbindlist(lapply(c("HCP", "ICP", "LCP"), function(cl) {
    cl_sym <- unique(hgc[weber_class == cl]$symbol)
    # STAT TEST: hypergeometric over-representation (clusterProfiler enrichGO), BH-adjusted
    e <- enrichGO(gene = cl_sym, OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL", ont = "ALL", universe = uni_h,
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
    d <- if (is.null(e)) data.table() else as.data.table(as.data.frame(e))
    if (nrow(d) && "ONTOLOGY" %in% names(d)) setnames(d, "ONTOLOGY", "ontology")
    k <- tryCatch({
      s2e <- function(v) unique(na.omit(AnnotationDbi::mapIds(
               org.Hs.eg.db, keys = v, keytype = "SYMBOL", column = "ENTREZID",
               multiVals = "first")))
      # with enricher() = the identical hypergeometric test, no network call.
      #   hsa_pathway_link.tsv : "path:hsa04973<TAB>hsa:57818"  (pathway -> Entrez gene)
      #   hsa_pathway_names.tsv: "hsa01522<TAB>Endocrine resistance - Homo sapiens (human)"
      kl <- fread(KEGG_LINK, header = FALSE, col.names = c("term", "gene"))
      kl[, `:=`(term = sub("^path:", "", term), gene = sub("^hsa:", "", gene))]
      kn <- fread(KEGG_NAMES, header = FALSE, col.names = c("term", "name"))
      kn[, name := sub(" - Homo sapiens \\(human\\)$", "", name)]
      ek <- enricher(gene = s2e(cl_sym), universe = s2e(uni_h),
                     TERM2GENE = kl[, .(term, gene)], TERM2NAME = kn[, .(term, name)],
                     pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
      if (is.null(ek)) data.table() else {
        kd <- as.data.table(as.data.frame(ek)); if (nrow(kd)) kd[, ontology := "KEGG"]; kd
      }
    }, error = function(e) { cat("  human KEGG skipped:", conditionMessage(e), "\n"); data.table() })
    d <- rbindlist(list(d, k), fill = TRUE)
    if (nrow(d)) d[, weber_class := cl]
    d
  }), fill = TRUE)
}, error = function(e) { cat("  human GO failed:", conditionMessage(e), "\n"); data.table() })
if (nrow(hg_go)) {
  fwrite(hg_go, file.path(DAT, "human_weber_class_go_enrichment.tsv"), sep = "\t")
  cat(sprintf("  human GO/KEGG at FDR<0.05: %s\n",
              paste(sprintf("%s=%d", c("HCP","ICP","LCP"),
                vapply(c("HCP","ICP","LCP"), function(cl) nrow(hg_go[weber_class==cl & p.adjust<0.05]), integer(1))),
                collapse = " ")))
  go_class_figs(hg_go, "human", "H. sapiens")
}

# ---- [9] Weber-class methylation & expression characterisation (supp only) --
# All three classes on: (a) methylation status, (b) basal expression level,
# (c) methylated vs silent. Reuses prom_meth/pdt + 01_genome_toolkit's DE table (an allowed
# upstream read). The HCP-with-DMR panels are drawn by 06_decoupling (they need 05_differential).
cat("[9] Weber-class promoter methylation & expression characterisation\n")
DE_TAIL3 <- file.path(PIPE, "01_genome_toolkit/data/gene_de_tail.tsv")
de3 <- fread(DE_TAIL3)[, .(gene_id, baseMean, log2FoldChange, padj)]
CLS <- c("HCP", "ICP", "LCP")
METH_CUT <- 0.2   # "methylated promoter" = mean beta > 0.2: a convention (stated in Methods), not a fitted value
cls_ids <- lapply(CLS, function(cl) pdt[biotype == "protein_coding" & weber_class == cl, gene_id])
names(cls_ids) <- CLS

# (a) methylation STATUS per class: no-coverage, covered-unmethylated, covered-methylated
status <- rbindlist(lapply(CLS, function(cl) {
  ids <- cls_ids[[cl]]; cm <- prom_meth[gene_id %in% ids]
  data.table(weber_class = cl, n_total = length(ids),
             `No coverage` = length(ids) - nrow(cm),
             `Covered, unmethylated` = sum(cm$mean_meth <= METH_CUT),
             `Covered, methylated`   = sum(cm$mean_meth > METH_CUT))
}))
fwrite(status, file.path(DAT, "weber_class_methylation_status.tsv"), sep = "\t")
print(status)
sl <- data.table::melt(status, id.vars = c("weber_class", "n_total"), variable.name = "state", value.name = "n")  # namespaced (masked-generic rule)
sl[, `:=`(weber_class = factor(weber_class, levels = CLS),
          state = factor(state, levels = c("No coverage", "Covered, unmethylated", "Covered, methylated")),
          pct = 100 * n / n_total)]
p_stat <- ggplot(sl, aes(weber_class, pct, fill = state)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5), size = 2.6, colour = "white") +
  scale_fill_manual(values = c("grey70", "#0072B2", "#D55E00"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Weber promoter class", y = "% of promoters",
       title = "Promoter methylation rises HCP < ICP < LCP",
       subtitle = "CpG-island (HCP) promoters stay unmethylated; CpG-poor (LCP) are mostly methylated") +
  theme_pub() + theme(legend.position = "bottom")
save_supp(p_stat, "figS3_weber_class_methylation_status", 6.5, 4.0)

# (b) BASAL expression level per class (answers "are HCP the basally-expressed ones?")
bas <- rbindlist(lapply(CLS, function(cl) de3[gene_id %in% cls_ids[[cl]], .(weber_class = cl, gene_id, baseMean)]))
bas[, weber_class := factor(weber_class, levels = CLS)]
# STAT TEST: Kruskal-Wallis rank-sum test (baseMean across the three Weber classes)
kw_bas <- kruskal.test(baseMean ~ weber_class, data = bas)
med_bas <- bas[, .(med = median(baseMean, na.rm = TRUE)), by = weber_class]
fwrite(med_bas, file.path(DAT, "weber_class_basal_expression.tsv"), sep = "\t")
p_bas <- ggplot(bas, aes(weber_class, log2(baseMean + 1), fill = weber_class)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  scale_fill_manual(values = COL_WEBER, guide = "none") +
  labs(x = "Weber promoter class", y = expression(log[2](baseMean + 1)),
       title = "Basal expression by promoter class",
       subtitle = sprintf("genes with a DESeq2 estimate; median baseMean HCP %.0f / ICP %.0f / LCP %.0f (Kruskal-Wallis p = %.2g)",
                          med_bas$med[1], med_bas$med[2], med_bas$med[3], kw_bas$p.value)) +
  theme_pub()
save_supp(p_bas, "figS3_weber_class_basal_expression", 5.5, 4.0)

# (c) are the METHYLATED promoters the silent ones? per class, methylated vs not
ms <- rbindlist(lapply(CLS, function(cl) {
  cm <- merge(prom_meth[gene_id %in% cls_ids[[cl]], .(gene_id, beta = mean_meth)], de3, by = "gene_id")
  if (!nrow(cm)) return(NULL)
  cm[, `:=`(weber_class = cl, meth_state = ifelse(beta > METH_CUT, "methylated", "unmethylated"))]; cm
}))
ms[, weber_class := factor(weber_class, levels = CLS)]
msum <- ms[, .(n = .N, median_baseMean = round(median(baseMean))), by = .(weber_class, meth_state)]
fwrite(msum, file.path(DAT, "weber_class_meth_state_expression.tsv"), sep = "\t")
print(msum)
# STAT TEST: per class, two-sided Mann-Whitney U (methylated vs unmethylated
# baseMean); effect size = rank-biserial r_rb = 1 - 2U/(n_meth*n_unmeth).
# the NEGATIVE of the textbook 2U/(n1 n2) - 1 form, matching the paper's positive values.
mw <- rbindlist(lapply(CLS, function(cl) {
  a <- ms[weber_class == cl & meth_state == "methylated",   baseMean]
  b <- ms[weber_class == cl & meth_state == "unmethylated", baseMean]
  if (length(a) < 2L || length(b) < 2L)
    return(data.table(weber_class = cl, n_meth = length(a), n_unmeth = length(b),
                      U = NA_real_, p = NA_real_, rank_biserial = NA_real_))
  wt  <- suppressWarnings(wilcox.test(a, b))          # Mann-Whitney U, two-sided
  U   <- unname(wt$statistic)
  data.table(weber_class = cl, n_meth = length(a), n_unmeth = length(b),
             U = U, p = wt$p.value, rank_biserial = round(1 - 2 * U / (length(a) * length(b)), 3))
}))
fwrite(mw, file.path(DAT, "weber_class_meth_state_mannwhitney.tsv"), sep = "\t")
cat("  Mann-Whitney meth vs unmeth baseMean, per class:\n"); print(mw)
p_ms <- ggplot(ms, aes(meth_state, log2(baseMean + 1), fill = meth_state)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~ weber_class, nrow = 1) +
  scale_fill_manual(values = c(methylated = "#D55E00", unmethylated = "#0072B2"), guide = "none") +
  labs(x = NULL, y = expression(log[2](baseMean + 1)),
       title = sprintf("Methylated promoters are lower expressed at %s", {
                 lo <- mw[!is.na(p) & p < 0.05 & rank_biserial > 0, weber_class]
                 if (length(lo)) paste(lo, collapse = " and ") else "no class" }),
       subtitle = paste0("beta > ", METH_CUT, " = methylated; ",
                         paste(mw[, sprintf("%s: n = %d/%d, P = %.2g, r = %+.2f", weber_class, n_meth, n_unmeth, p, rank_biserial)],
                               collapse = "; "))) +
  theme_pub() + theme(strip.background = element_blank(),
                      axis.text.x = element_text(angle = 20, hjust = 1, size = 7))
save_supp(p_ms, "figS3_weber_class_meth_state_expression", 7.0, 3.6)

# HCP-promoter-DMR panels need 05_differential's DMR calls, so 06_decoupling reads
# promoter_weber_classification.tsv and writes hcp_promoter_dmr_genes.tsv + the
# figS6_hcp_* panels. Nothing in this module depends on 05_differential.

# ---- [9a] AP-2 domain architecture ----
# The single D. laeve AP-2 gene (LOC_00012416; assembly symbol TFAP2B, JASPAR bait
# TFAP2A) dominates the HCP motif signal, so the supplement shows its protein's
# Pfam domain architecture, in the style of 01_genome_toolkit's figS_uhrf_domains. Domains
# come from a fresh hmmscan (project-local HMMER 3.4 + pressed Pfam-A, gathering
local({
  PROTE  <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.proteins.fasta"
  HMMSCAN <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/hmmer/bin/hmmscan"
  PFAM    <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/pfam/Pfam-A.hmm"
  AP2     <- "LOC_00012416"
  aa  <- readAAStringSet(PROTE)
  aa  <- aa[grepl(paste0("^", AP2, "-mRNA"), names(aa))]
  aa  <- aa[which.max(width(aa))]                      # longest isoform
  # intermediates live beside the other cached inputs in dataset/
  qfa <- file.path(DS, "ap2_protein.fa"); writeXStringSet(aa, qfa)
  dtb <- file.path(DS, "ap2_hmmscan.domtblout")
  stopifnot(system2(HMMSCAN, c("--cut_ga", "--domtblout", dtb, PFAM, qfa),
                    stdout = FALSE) == 0)
  dom <- fread(cmd = sprintf("/usr/bin/grep -v '^#' %s | awk '{print $1, $13, $20, $21}'", dtb),
               header = FALSE, col.names = c("domain", "i_evalue", "from", "to"))
  dom <- dom[i_evalue < 0.01][order(from)]
  fwrite(dom, file.path(DAT, "ap2_domains.tsv"), sep = "\t")
  plen <- width(aa)
  cat(sprintf("[9] AP-2 %s (%d aa): domains %s\n", names(aa), plen,
              paste(dom[, sprintf("%s %d-%d", domain, from, to)], collapse = ", ")))
  pd <- ggplot() +
    geom_segment(aes(x = 1, xend = plen, y = 0, yend = 0), linewidth = 2.2, colour = "grey70") +
    geom_rect(data = dom, aes(xmin = from, xmax = to, ymin = -0.28, ymax = 0.28, fill = domain),
              colour = "white", linewidth = 0.3) +
    geom_text(data = dom, aes(x = (from + to)/2, y = 0.52, label = domain), size = 2.8) +
    geom_text(data = dom, aes(x = (from + to)/2, y = -0.52,
                              label = sprintf("%d-%d", from, to)), size = 2.4, colour = "grey30") +
    scale_fill_manual(values = rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
                                   length.out = max(1, nrow(dom))), guide = "none") +
    scale_x_continuous(limits = c(1, plen), expand = expansion(mult = c(0.01, 0.01))) +
    ylim(-0.8, 0.8) +
    labs(x = sprintf("%s (AP-2), %d aa", AP2, plen), y = NULL,
         title = "Domain architecture of the D. laeve AP-2 protein") +
    theme_pub() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                        axis.line.y = element_blank(),
                        plot.title = element_text(size = rel(1)))
  save_supp(pd, "figS3_ap2_domain", 4.8, 1.7)
})

# ---- [10] Reproducibility: sessionInfo record -------------------------------
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_03_promoters.txt"))
cat("[03_promoters] done\n")
