#!/usr/bin/env Rscript
set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges)
  library(bsseq); library(DESeq2)
})
# DSS and svglite are attached later, where first used.

PIPE <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
HTSEQ <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"
B01 <- file.path(PIPE, "01_genome_toolkit/objects"); B01D <- file.path(PIPE, "01_genome_toolkit/data")
B02 <- file.path(PIPE, "02_landscape/objects"); B05D <- file.path(PIPE, "05_differential/data")
BATCH <- file.path(PIPE, "06_decoupling"); DAT <- file.path(BATCH, "data")
dir.create(DAT, showWarnings = FALSE, recursive = TRUE)
keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")

# ---- [1] Gene-body methylation (pooled + per condition) ----------------------
# Per-gene β = coverage-summed CpG counts over the gene body, pooled and per
# condition; genes covered at fewer than 5 CpGs are dropped (n_cpg >= 5 filter).
# dbeta = beta_amp - beta_ctrl (amputated minus control throughout).
cat("[1] gene-body methylation (pooled + per condition)\n")
bs <- readRDS(file.path(B02, "bsseq_cov5_chrmt.rds"))
chrs <- as.character(seqnames(bs)); bs <- bs[chrs %in% keep_chr, ]
gr <- granges(bs)
M <- as.matrix(getCoverage(bs, type = "M")); Cv <- as.matrix(getCoverage(bs, type = "Cov"))
cond <- ifelse(grepl("^A", sampleNames(bs)), "Amputated", "Control")
mc <- rowSums(M[, cond=="Control", drop=FALSE]); cc <- rowSums(Cv[, cond=="Control", drop=FALSE])
ma <- rowSums(M[, cond=="Amputated", drop=FALSE]); ca <- rowSums(Cv[, cond=="Amputated", drop=FALSE])
gff <- readRDS(file.path(B01, "gff_chrmt.rds"))
gene_gr <- gff[gff$type == "gene"]
gene_gr <- gene_gr[as.character(mcols(gene_gr)$gene_biotype) %in% c("protein_coding", "lncRNA")]
gid <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
gene_dt <- data.table(chr = as.character(seqnames(gene_gr)), start = start(gene_gr),
                      end = end(gene_gr), gene_id = gid)
cpg_dt <- data.table(chr = as.character(seqnames(gr)), start = start(gr), end = start(gr),
                     mc = mc, cc = cc, ma = ma, ca = ca)
setkey(gene_dt, chr, start, end); setkey(cpg_dt, chr, start, end)
ov <- foverlaps(cpg_dt, gene_dt, nomatch = 0L)
gm <- ov[, .(beta = (sum(mc)+sum(ma))/pmax(sum(cc)+sum(ca),1),
             beta_ctrl = sum(mc)/pmax(sum(cc),1), beta_amp = sum(ma)/pmax(sum(ca),1),
             n_cpg = .N), by = gene_id][n_cpg >= 5]
gm[, dbeta := beta_amp - beta_ctrl]

# ---- [2] Tail expression (HTSeq -> VST) + differential expression ------------
# 4 control (C1S1-C4S4) vs 3 amputated (T2S6-T4S8) tail libraries; genes kept if
# count >= 5 in >= 2 libraries; VST blind; expr_mean pools all 7 libraries.
# The DE table (01_genome_toolkit) and DMR list (05_differential) load behind file.exists fallbacks:
# an absent upstream file degrades the run to empty sets / NAs instead of stopping.
cat("[2] tail expression (HTSeq -> VST) + differential expression\n")
tail_s <- c("C1S1","C2S2","C3S3","C4S4","T2S6","T3S7","T4S8")
cl <- lapply(tail_s, function(s) {
  x <- fread(file.path(HTSEQ, paste0(s, "_htseq_gene_counts.txt")),
             header = FALSE, col.names = c("gene_id","count")); x[!startsWith(gene_id,"__")] })
ids <- cl[[1]]$gene_id; cm <- sapply(cl, function(x) x$count); rownames(cm) <- ids
cm <- cm[rownames(cm) %in% gene_dt$gene_id, , drop = FALSE]
cm <- cm[rowSums(cm >= 5) >= 2, , drop = FALSE]                  # expressed: >= 5 reads in >= 2 libraries
vsd <- vst(DESeqDataSetFromMatrix(cm, data.frame(s = tail_s), ~ 1), blind = TRUE)
expr <- data.table(gene_id = rownames(assay(vsd)), expr_mean = rowMeans(assay(vsd)))

de_path <- file.path(B01D, "gene_de_tail.tsv")
stopifnot(file.exists(de_path))                    # fail LOUD: a silent empty-DE fallback fakes a null
stopifnot(file.exists(de_path)); de <- fread(de_path)

dmr_path <- file.path(B05D, "dmrs_annotated.tsv")
stopifnot(file.exists(dmr_path))                   # same rule: sections 3b/4 read it unguarded anyway
stopifnot(file.exists(dmr_path))
# definition as the Venn and the Methods, not only the single display gene per DMR
dmr_asg_path <- file.path(B05D, "dmrs_gene_assignments.tsv")
dmr_genes <- if (file.exists(dmr_asg_path)) unique(fread(dmr_asg_path)$gene_id) else unique(fread(dmr_path)$gene_id)
dmr_genes <- dmr_genes[!is.na(dmr_genes) & dmr_genes != ""]

# ---- [3] Merge + decoupling stats (Fisher; length-adjusted CMH) --------------
# Per-gene table -> decoupling_per_gene.tsv; headline numbers -> decoupling_summary.tsv;
# per-stratum counts -> decoupling_length_strata.tsv.
cat("[3] merge + decoupling stats\n")
g <- merge(gm, expr, by = "gene_id")
g <- merge(g, de[, .(gene_id, log2FoldChange, padj)], by = "gene_id", all.x = TRUE)
g[, has_dmr := gene_id %in% dmr_genes]
g <- merge(g, gene_dt[, .(gene_id, gene_len = end - start + 1L)], by = "gene_id", all.x = TRUE)
fwrite(g, file.path(DAT, "decoupling_per_gene.tsv"), sep = "\t")

# baseline: gene-body methylation ~ expression
baseline_r2 <- summary(lm(expr_mean ~ beta, data = g))$r.squared
# differential: Δmethylation ~ Δexpression (only genes with a DE estimate)
gd <- g[is.finite(log2FoldChange)]
diff_r2 <- if (nrow(gd) > 10) summary(lm(log2FoldChange ~ dbeta, data = gd))$r.squared else NA_real_
gd[, is_de := !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1]
# STAT TEST: two-sided Fisher's exact test — is carrying a DMR associated with being DE?
fish <- tryCatch(fisher.test(table(gd$has_dmr, gd$is_de)), error = function(e) NULL)
# Length-adjusted version: DMR genes are ~3x longer, a confound for any DMR enrichment.
# STAT TEST: Cochran-Mantel-Haenszel (mantelhaen.test) on has_dmr x is_de stratified by gene-length tertile.
gd_l <- gd[!is.na(gene_len)]
gd_l[, ltert := cut(gene_len, quantile(gene_len, 0:3/3, na.rm = TRUE), include.lowest = TRUE)]
# exact = TRUE: only 5 genes are both DMR-carrying and DE and all fall in the longest
# length tertile, so two of three strata are uninformative; the asymptotic chi-square
# is not trustworthy that sparse. Per-stratum counts are written out to keep it visible.
cmh <- tryCatch(mantelhaen.test(table(gd_l$has_dmr, gd_l$is_de, gd_l$ltert), exact = TRUE),
                error = function(e) NULL)
strata_n <- gd_l[has_dmr == TRUE & is_de == TRUE, .N, by = ltert][order(ltert)]
fwrite(gd_l[, .(n_genes = .N, n_dmr = sum(has_dmr), n_de = sum(is_de),
                n_dmr_and_de = sum(has_dmr & is_de)), by = ltert][order(ltert)],
       file.path(DAT, "decoupling_length_strata.tsv"), sep = "\t")
cat("  DMR-and-DE genes per length tertile: ",
    paste(sprintf("%s=%d", strata_n$ltert, strata_n$N), collapse = ", "), "\n", sep = "")

summ <- data.table(
  n_genes = nrow(g),
  baseline_methylation_expr_R2 = baseline_r2,
  differential_dMeth_dExpr_R2 = diff_r2,
  de_definition = "padj<0.05 & |log2FC|>=1 (the paper's strict DE set)",
  dmr_de_OR = if (!is.null(fish)) unname(fish$estimate) else NA_real_,
  dmr_de_OR_lo = if (!is.null(fish)) fish$conf.int[1] else NA_real_,
  dmr_de_OR_hi = if (!is.null(fish)) fish$conf.int[2] else NA_real_,
  dmr_de_p  = if (!is.null(fish)) fish$p.value else NA_real_,
  dmr_de_OR_lengthadj = if (!is.null(cmh)) unname(cmh$estimate) else NA_real_,
  dmr_de_OR_lengthadj_lo = if (!is.null(cmh) && !is.null(cmh$conf.int)) cmh$conf.int[1] else NA_real_,
  dmr_de_OR_lengthadj_hi = if (!is.null(cmh) && !is.null(cmh$conf.int)) cmh$conf.int[2] else NA_real_,
  dmr_de_p_lengthadj  = if (!is.null(cmh)) cmh$p.value else NA_real_,
  headline = "asymmetry: baseline methylation tracks expression; Δmethylation barely predicts Δexpression")
fwrite(summ, file.path(DAT, "decoupling_summary.tsv"), sep = "\t")
cat(sprintf("  baseline R2=%.4f | differential R2=%s | DMR-DE OR=%s p=%s\n",
            baseline_r2, ifelse(is.na(diff_r2), "NA", sprintf("%.4f", diff_r2)),
            if (!is.null(fish)) sprintf("%.2f", fish$estimate) else "NA",
            if (!is.null(fish)) sprintf("%.3g", fish$p.value) else "NA"))

# ---- [3b] DMR ΔMeth vs gene log2FC correlation (Spearman) --------------------
# The manuscript's headline decoupling statistic, reported as (test, statistic):
# Spearman of per-DMR ΔMeth (amputated - control) against the host gene's DESeq2
# log2FC, overall and per genomic region -> dmr_de_correlation.tsv.
cat("[3b] DMR deltaMeth vs gene log2FC correlation (Spearman)\n")
dmr_de <- merge(fread(dmr_path)[, .(gene_id, dMeth = -diff.Methy, region)],   # diff.Methy = ctrl-amp in 05_differential; negate -> amputated - control (paper convention)
                de[, .(gene_id, log2FoldChange)], by = "gene_id")
dmr_de <- dmr_de[is.finite(dMeth) & is.finite(log2FoldChange)]
spear <- function(x, lab) {
  if (nrow(x) < 10) return(data.table(set = lab, n = nrow(x), rho = NA_real_, p = NA_real_))
  # STAT TEST: Spearman rank correlation test (DMR delta-methylation vs gene log2FC)
  ct <- suppressWarnings(cor.test(x$dMeth, x$log2FoldChange, method = "spearman"))
  data.table(set = lab, n = nrow(x), rho = unname(ct$estimate), p = ct$p.value)
}
corr <- rbindlist(list(
  spear(dmr_de,                                   "all DMRs"),
  spear(dmr_de[region == "Promoter"],             "promoter DMRs"),
  spear(dmr_de[region %in% c("Exon", "Intron")],  "gene-body DMRs"),
  spear(dmr_de[region == "Exon"],                 "exon DMRs"),
  spear(dmr_de[region == "Intron"],               "intron DMRs"),
  spear(dmr_de[, .SD[which.max(abs(dMeth))], by = gene_id], "per-gene (max|dMeth|)")))
fwrite(corr, file.path(DAT, "dmr_de_correlation.tsv"), sep = "\t")
print(corr[, .(set, n, rho = round(rho, 3), p = signif(p, 3))])
main_r <- corr[set == "all DMRs"]
cat(sprintf("  MANUSCRIPT VALUE: DMR ΔMeth vs gene log2FC — Spearman rho = %.3f, p = %.2g (n = %d DMR-gene pairs)\n",
            main_r$rho, main_r$p, main_r$n))

# ---- [3c] MAIN fig6a: the decoupling on its three axes --------------------------
# left, per-DMR methylation change against the host gene's log2 fold change (direction);
# centre, the DMR-bearing vs differentially-expressed odds ratios with 95% CIs, unadjusted
# Fisher and the length-stratified CMH (occurrence); right, the two R^2 values (magnitude).
# All numbers are the ones written to data/ above; nothing is recomputed here.
cat("[3c] fig6a: the decoupling on three axes\n")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
FIGM <- file.path(BATCH, "figures/main"); dir.create(FIGM, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(FIGM, "\\.(pdf|png|svg)$", full.names = TRUE))
theme_pub <- function() theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 9, face = "bold"),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGM, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)  # cairo: Δ glyphs
  ggsave(file.path(FIGM, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGM, paste0(name, ".svg")), p, width = w, height = h)                       # vector (svg)
  cat(sprintf("  saved %s\n", name))
}
p_dir <- ggplot(dmr_de, aes(dMeth, log2FoldChange)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(alpha = 0.35, size = 0.8, colour = "#0072B2") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.6,
           label = sprintf("Spearman rho = %.3f, P = %.2f\nn = %d DMR-gene pairs", main_r$rho, main_r$p, main_r$n)) +
  labs(x = "DMR methylation change (amputated - control)", y = "Gene log2 fold change",
       title = "Direction") + theme_pub()
ors <- data.table(test = c("Fisher, unadjusted", "CMH, gene-length tertiles"),
                  OR = c(summ$dmr_de_OR, summ$dmr_de_OR_lengthadj),
                  lo = c(summ$dmr_de_OR_lo, summ$dmr_de_OR_lengthadj_lo),
                  hi = c(summ$dmr_de_OR_hi, summ$dmr_de_OR_lengthadj_hi))
ors[, test := factor(test, levels = rev(test))]
p_occ <- ggplot(ors, aes(OR, test)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, yend = test), linewidth = 0.5, colour = "grey30") +
  geom_point(size = 2.4) +
  scale_x_log10() +
  labs(x = "Odds ratio, DMR-bearing vs differentially expressed (95% CI)", y = NULL,
       title = "Occurrence") + theme_pub()
r2 <- data.table(axis = c("Baseline methylation\nvs expression", "Methylation change\nvs expression change"),
                 R2 = c(summ$baseline_methylation_expr_R2, summ$differential_dMeth_dExpr_R2))
r2[, axis := factor(axis, levels = axis)]
p_mag <- ggplot(r2, aes(axis, R2)) +
  geom_col(width = 0.55, fill = "#0072B2") +
  geom_text(aes(label = sprintf("R² = %.3g", R2)), vjust = -0.4, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = NULL, y = "Variance explained (R²)", title = "Magnitude") +
  theme_pub() + theme(axis.text.x = element_text(size = 7))
save_fig(p_dir + p_occ + p_mag + plot_layout(widths = c(1.15, 1.15, 0.8)), "fig6a_decoupling_three_axes", 9.0, 2.8)

# ---- [4] Supplementary DSS tracks for the DMP∩DMR∩DE genes -------------------
# For genes simultaneously DMP-, DMR- and DE-associated (05_differential Venn, strict DE):
# (a) test whether any transcript TSS falls inside the DMR's host intron (a
# cryptic/alternative-promoter candidate) -> intersection_dmr_tss.tsv; (b) draw
# per-sample methylation across each DMR with DSS::showOneDMR (pdf + png + svg).
cat("[4] supplementary DMR methylation tracks (DSS::showOneDMR)\n")
suppressPackageStartupMessages({ library(DSS); library(svglite) })
FIGS <- file.path(BATCH, "figures/supplementary")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(FIGS, "\\.(pdf|png|svg)$", full.names = TRUE))   # no stale stems survive a rerun

isec_path <- file.path(B05D, "dmp_dmr_de_intersection.tsv")
dmr_all <- fread(dmr_path); dmr_all[, g := sub(";.*", "", gene_id)]
de[, g := sub(";.*", "", gene_id)]
stopifnot(file.exists(isec_path)); sel_genes <- unique(sub(";.*", "", fread(isec_path)$gene_id))
sel <- dmr_all[g %in% sel_genes]

# transcript TSS (strand-aware) + host-gene introns for the TSS-in-intron test
tx  <- gff[gff$type %in% c("mRNA", "lnc_RNA")]   # EviAnn types: mRNA + lnc_RNA (no "transcript")
tss <- ifelse(as.character(strand(tx)) == "-", end(tx), start(tx))
tss_gr  <- GRanges(seqnames(tx), IRanges(tss, tss), strand = strand(tx))
exon_gr <- gff[gff$type == "exon"]
gene_id_all <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
disp <- function(gg) { s <- de[g == gg, symbol][1]
  if (is.na(s) || s == "" || s == "Non annotated gene") gg else s }

annot <- rbindlist(lapply(seq_len(nrow(sel)), function(i) {
  r <- sel[i]; d <- GRanges(r$chr, IRanges(r$start, r$end))
  gg <- gene_gr[gene_id_all == r$g]
  intron <- if (length(gg)) { gaps <- setdiff(range(gg), reduce(exon_gr[overlapsAny(exon_gr, gg)]))
                              gaps[overlapsAny(gaps, d)] } else GRanges()
  nd <- distanceToNearest(d, tss_gr)
  data.table(gene_id = r$g, symbol = disp(r$g), chr = r$chr, start = r$start, end = r$end,
             region = r$region, direction = r$direction,
             log2FC = de[g == r$g, log2FoldChange][1],
             tss_in_dmr = sum(overlapsAny(tss_gr, d)),
             tss_in_host_intron = if (length(intron)) sum(overlapsAny(tss_gr, intron, ignore.strand = TRUE)) else 0L,
             nearest_tss_bp = if (length(nd)) mcols(nd)$distance else NA_integer_)
}))
fwrite(annot, file.path(DAT, "intersection_dmr_tss.tsv"), sep = "\t")

# keep only what is necessary — coordinates, the TSS-in-intron call and the
# nearest-TSS distance live in intersection_dmr_tss.tsv, not on the panel.
draw_one <- function(i) {
  r <- sel[i]; a <- annot[i]
  par(oma = c(0, 0, 2.4, 0))
  DSS::showOneDMR(data.frame(chr = r$chr, start = r$start, end = r$end), bs, ext = 2000)
  mtext(sprintf("%s  %s %s DMR  |  DE log2FC %+.2f",
                a$symbol, tolower(r$region), r$direction, a$log2FC),
        outer = TRUE, side = 3, line = 0.4, cex = 0.85, font = 2)
}
# (a) one combined multi-page PDF with all intersection DMRs
pdf(file.path(FIGS, "figS6_intersection_dmrs.pdf"), width = 6.5, height = 5)
for (i in seq_len(nrow(sel))) draw_one(i)
dev.off()
# (b) per-DMR png + svg (single-page formats)
for (i in seq_len(nrow(sel))) {
  nm <- sprintf("figS6_dmr_%s", gsub("[^A-Za-z0-9]+", "_", annot$symbol[i]))
  if (sum(annot$symbol == annot$symbol[i]) > 1) nm <- sprintf("%s_%d", nm, i)   # two DMRs in one gene keep separate files
  png(file.path(FIGS, paste0(nm, ".png")), width = 6.5, height = 5, units = "in", res = 150); draw_one(i); dev.off()
  svglite::svglite(file.path(FIGS, paste0(nm, ".svg")), width = 6.5, height = 5); draw_one(i); dev.off()
}
cat(sprintf("  %d intersection DMRs plotted (TSS-in-intron: %d) -> figS6_*\n",
            nrow(sel), sum(annot$tss_in_host_intron > 0)))

# ---- [5] Supplementary: HCP (CpG-island-like) promoters carrying a DMR -------
# Do any CpG-rich (HCP) protein-coding promoters change methylation after
# amputation? Each such DMR is drawn with DSS::showOneDMR; table + supp figures only.
# order; this needs 05_differential DMRs AND 03_promoters Weber classes, so 06_decoupling is its home).
cat("[5] HCP promoters carrying a DMR (DSS::showOneDMR)\n")
B03D <- file.path(PIPE, "03_promoters/data")
weber <- fread(file.path(B03D, "promoter_weber_classification.tsv"))      # 03_promoters (upstream)
hcp_ids <- weber[biotype == "protein_coding" & weber_class == "HCP", gene_id]
hdmr <- dmr_all[region == "Promoter" & g %in% hcp_ids]
hdmr <- merge(hdmr, de[, .(g, baseMean, log2FoldChange, padj)], by = "g", all.x = TRUE)
hdmr[, symbol := vapply(g, disp, character(1))]
hdmr[, dMeth := -diff.Methy]            # diff.Methy = ctrl - amp in 05_differential; negate -> amputated - control (paper convention, as in 3b)
hdmr[, ad := abs(dMeth)]; setorder(hdmr, -ad); hdmr[, ad := NULL]   # largest change first
fwrite(hdmr[, .(symbol, gene_id = g, chr, start, end, dMeth_amp_minus_ctrl = dMeth, direction,
                baseMean, log2FoldChange, padj)],
       file.path(DAT, "hcp_promoter_dmr_genes.tsv"), sep = "\t")
cat(sprintf("  %d of %d HCP promoters carry a DMR\n", nrow(hdmr), length(hcp_ids)))
if (nrow(hdmr)) {
  draw_hcp <- function(i) {                 # same legibility rules as draw_one above
    r <- hdmr[i]; par(oma = c(0, 0, 2.4, 0))
    DSS::showOneDMR(data.frame(chr = r$chr, start = r$start, end = r$end), bs, ext = 2000)
    mtext(sprintf("%s  %s promoter DMR  |  dMeth (amp - ctrl) %+.2f",
                  r$symbol, r$direction, r$dMeth),
          outer = TRUE, side = 3, line = 0.4, cex = 0.85, font = 2)
  }
  pdf(file.path(FIGS, "figS6_hcp_promoter_dmrs.pdf"), width = 6.5, height = 5)   # one page per DMR
  for (i in seq_len(nrow(hdmr))) draw_hcp(i)
  dev.off()
  for (i in seq_len(nrow(hdmr))) {                                             # single-page png + svg
    nm <- sprintf("figS6_hcp_dmr_%s", gsub("[^A-Za-z0-9]+", "_", hdmr$symbol[i]))
    png(file.path(FIGS, paste0(nm, ".png")), width = 6.5, height = 5, units = "in", res = 150); draw_hcp(i); dev.off()
    svglite::svglite(file.path(FIGS, paste0(nm, ".svg")), width = 6.5, height = 5); draw_hcp(i); dev.off()
  }
  cat(sprintf("  drew %d HCP promoter DMRs -> figS6_hcp_promoter_dmrs.pdf + figS6_hcp_dmr_*\n", nrow(hdmr)))
}

# ---- [6] Reproducibility: record the exact package versions this run used ----
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_06_decoupling.txt"))
cat("[06_decoupling] done\n")
