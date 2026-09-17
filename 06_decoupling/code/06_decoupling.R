#!/usr/bin/env Rscript
# 06_decoupling.R
# Methylation versus expression: baseline gene-body methylation against expression,
# differential methylation against differential expression, and DMR track figures.
# Inputs: 02_landscape/objects/bsseq_cov5_chrmt.rds; 01_genome_toolkit/objects/gff_chrmt.rds and
#   data/gene_de_tail.tsv; HTSeq tail gene counts; 05_differential/data/dmrs_annotated.tsv,
#   dmrs_gene_assignments.tsv and dmp_dmr_de_intersection.tsv; 03_promoters/data/
#   promoter_weber_classification.tsv
# Outputs: data/ decoupling_per_gene, decoupling_summary, decoupling_length_strata, dmr_de_correlation,
#   intersection_dmr_tss, hcp_promoter_dmr_genes (TSV); figures/main fig6a_decoupling_three_axes;
#   figures/supplementary figS6_* DSS DMR tracks; sessionInfo_06_decoupling.txt
# Run: sbatch 06_decoupling/code/06_decoupling.slurm  (from main/methylation_pipeline/)

# Step 1 - Seed, packages, paths and chromosome filter
set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges)
  library(bsseq); library(DESeq2)
})

PIPE <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
HTSEQ <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"
B01 <- file.path(PIPE, "01_genome_toolkit/objects"); B01D <- file.path(PIPE, "01_genome_toolkit/data")
B02 <- file.path(PIPE, "02_landscape/objects"); B05D <- file.path(PIPE, "05_differential/data")
BATCH <- file.path(PIPE, "06_decoupling"); DAT <- file.path(BATCH, "data")
dir.create(DAT, showWarnings = FALSE, recursive = TRUE)
keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")   # chr1-31 plus the mitochondrial scaffold

# Step 2 - Gene-body methylation per gene, pooled and per condition
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
# Gene universe: protein_coding and lncRNA genes only
gene_gr <- gene_gr[as.character(mcols(gene_gr)$gene_biotype) %in% c("protein_coding", "lncRNA")]
gid <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
gene_dt <- data.table(chr = as.character(seqnames(gene_gr)), start = start(gene_gr),
                      end = end(gene_gr), gene_id = gid)
cpg_dt <- data.table(chr = as.character(seqnames(gr)), start = start(gr), end = start(gr),
                     mc = mc, cc = cc, ma = ma, ca = ca)
setkey(gene_dt, chr, start, end); setkey(cpg_dt, chr, start, end)
ov <- foverlaps(cpg_dt, gene_dt, nomatch = 0L)
# Coverage-weighted beta per gene; genes with fewer than 5 covered CpGs are dropped
gm <- ov[, .(beta = (sum(mc)+sum(ma))/pmax(sum(cc)+sum(ca),1),
             beta_ctrl = sum(mc)/pmax(sum(cc),1), beta_amp = sum(ma)/pmax(sum(ca),1),
             n_cpg = .N), by = gene_id][n_cpg >= 5]
# dbeta = amputated minus control
gm[, dbeta := beta_amp - beta_ctrl]

# Step 3 - Tail expression (HTSeq counts -> VST) and the upstream DE and DMR tables
cat("[2] tail expression (HTSeq -> VST) + differential expression\n")
# 4 control and 3 amputated tail libraries
tail_s <- c("C1S1","C2S2","C3S3","C4S4","T2S6","T3S7","T4S8")
cl <- lapply(tail_s, function(s) {
  x <- fread(file.path(HTSEQ, paste0(s, "_htseq_gene_counts.txt")),
             header = FALSE, col.names = c("gene_id","count")); x[!startsWith(gene_id,"__")] })
ids <- cl[[1]]$gene_id; cm <- sapply(cl, function(x) x$count); rownames(cm) <- ids
cm <- cm[rownames(cm) %in% gene_dt$gene_id, , drop = FALSE]   # restrict to the filtered gene universe before the VST
cm <- cm[rowSums(cm >= 5) >= 2, , drop = FALSE]   # expressed: >= 5 reads in >= 2 libraries
vsd <- vst(DESeqDataSetFromMatrix(cm, data.frame(s = tail_s), ~ 1), blind = TRUE)
expr <- data.table(gene_id = rownames(assay(vsd)), expr_mean = rowMeans(assay(vsd)))

# Upstream tables must exist; the run stops rather than continuing with empty sets
de_path <- file.path(B01D, "gene_de_tail.tsv")
stopifnot(file.exists(de_path)); de <- fread(de_path)

dmr_path <- file.path(B05D, "dmrs_annotated.tsv")
stopifnot(file.exists(dmr_path))
# DMR-bearing genes = every gene a DMR is assigned to (long assignment table)
dmr_asg_path <- file.path(B05D, "dmrs_gene_assignments.tsv")
dmr_genes <- if (file.exists(dmr_asg_path)) unique(fread(dmr_asg_path)$gene_id) else unique(fread(dmr_path)$gene_id)
dmr_genes <- dmr_genes[!is.na(dmr_genes) & dmr_genes != ""]

# Step 4 - Merge and decoupling statistics (Fisher; length-adjusted CMH)
cat("[3] merge + decoupling stats\n")
g <- merge(gm, expr, by = "gene_id")
g <- merge(g, de[, .(gene_id, log2FoldChange, padj)], by = "gene_id", all.x = TRUE)
g[, has_dmr := gene_id %in% dmr_genes]
g <- merge(g, gene_dt[, .(gene_id, gene_len = end - start + 1L)], by = "gene_id", all.x = TRUE)
fwrite(g, file.path(DAT, "decoupling_per_gene.tsv"), sep = "\t")

# Baseline: gene-body methylation against expression (R2)
baseline_r2 <- summary(lm(expr_mean ~ beta, data = g))$r.squared
# Differential: methylation change against log2 fold change, genes with a DE estimate
gd <- g[is.finite(log2FoldChange)]
diff_r2 <- if (nrow(gd) > 10) summary(lm(log2FoldChange ~ dbeta, data = gd))$r.squared else NA_real_
# DE = padj < 0.05 and |log2FC| >= 1, the strict set used throughout the pipeline
gd[, is_de := !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1]
# STAT TEST: two-sided Fisher's exact test — is carrying a DMR associated with being DE?
fish <- tryCatch(fisher.test(table(gd$has_dmr, gd$is_de)), error = function(e) NULL)
# Length-adjusted version: DMR genes are longer, a confound for any DMR enrichment
# STAT TEST: Cochran-Mantel-Haenszel (mantelhaen.test) on has_dmr x is_de stratified by gene-length tertile.
gd_l <- gd[!is.na(gene_len)]
gd_l[, ltert := cut(gene_len, quantile(gene_len, 0:3/3, na.rm = TRUE), include.lowest = TRUE)]
# exact = TRUE: the DMR-and-DE cells are sparse, so the asymptotic test is not used
cmh <- tryCatch(mantelhaen.test(table(gd_l$has_dmr, gd_l$is_de, gd_l$ltert), exact = TRUE),
                error = function(e) NULL)
strata_n <- gd_l[has_dmr == TRUE & is_de == TRUE, .N, by = ltert][order(ltert)]
fwrite(gd_l[, .(n_genes = .N, n_dmr = sum(has_dmr), n_de = sum(is_de),
                n_dmr_and_de = sum(has_dmr & is_de)), by = ltert][order(ltert)],
       file.path(DAT, "decoupling_length_strata.tsv"), sep = "\t")
cat("  DMR-and-DE genes per length tertile: ",
    paste(sprintf("%s=%d", strata_n$ltert, strata_n$N), collapse = ", "), "\n", sep = "")

# Headline numbers in one row
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

# Step 5 - DMR methylation change vs host-gene log2FC (Spearman)
cat("[3b] DMR deltaMeth vs gene log2FC correlation (Spearman)\n")
dmr_de <- merge(fread(dmr_path)[, .(gene_id, dMeth = -diff.Methy, region)],   # diff.Methy is control minus amputated; negate to amputated minus control
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

# Step 6 - fig6a: the decoupling on three axes (direction, occurrence, magnitude)
cat("[3c] fig6a: the decoupling on three axes\n")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
FIGM <- file.path(BATCH, "figures/main"); dir.create(FIGM, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(FIGM, "\\.(pdf|png|svg)$", full.names = TRUE))
theme_pub <- function() theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 9, face = "bold"),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGM, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGM, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGM, paste0(name, ".svg")), p, width = w, height = h)
  cat(sprintf("  saved %s\n", name))
}
# Panels reuse the values written to data/ above; nothing is recomputed
p_dir <- ggplot(dmr_de, aes(dMeth, log2FoldChange)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(alpha = 0.35, size = 0.8, colour = "#0072B2") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.6,
           label = sprintf("Spearman rho = %.3f, P = %.2f\nn = %d DMR-gene pairs", main_r$rho, main_r$p, main_r$n)) +
  labs(x = "DMR methylation change (amputated - control)", y = "Gene log2 fold change",
       title = "Direction") + theme_pub()
# Occurrence: Fisher and CMH odds ratios with 95% CI
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
# Magnitude: the two R2 values
r2 <- data.table(axis = c("Baseline\nmethylation vs\nexpression", "Methylation\nchange vs\nexpression change"),
                 R2 = c(summ$baseline_methylation_expr_R2, summ$differential_dMeth_dExpr_R2))
r2[, axis := factor(axis, levels = axis)]
p_mag <- ggplot(r2, aes(axis, R2)) +
  geom_col(width = 0.55, fill = "#0072B2") +
  geom_text(aes(label = sprintf("R² = %.3g", R2)), vjust = -0.4, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = NULL, y = "Variance explained (R²)", title = "Magnitude") +
  theme_pub() + theme(axis.text.x = element_text(size = 7))
save_fig(p_dir + p_occ + p_mag + plot_layout(widths = c(1.1, 1.1, 1.0)), "fig6a_decoupling_three_axes", 9.5, 2.9)

# Step 7 - Supplementary DSS tracks for the DMP, DMR and DE intersection genes
cat("[4] supplementary DMR methylation tracks (DSS::showOneDMR)\n")
suppressPackageStartupMessages({ library(DSS); library(svglite) })
FIGS <- file.path(BATCH, "figures/supplementary")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(FIGS, "\\.(pdf|png|svg)$", full.names = TRUE))

# Intersection genes from 05_differential; their DMRs are drawn below
isec_path <- file.path(B05D, "dmp_dmr_de_intersection.tsv")
dmr_all <- fread(dmr_path); dmr_all[, g := sub(";.*", "", gene_id)]
de[, g := sub(";.*", "", gene_id)]
stopifnot(file.exists(isec_path)); sel_genes <- unique(sub(";.*", "", fread(isec_path)$gene_id))
sel <- dmr_all[g %in% sel_genes]

# Transcript TSS (strand-aware) and host-gene introns for the TSS-in-intron test
tx  <- gff[gff$type %in% c("mRNA", "lnc_RNA")]   # EviAnn transcript types
tss <- ifelse(as.character(strand(tx)) == "-", end(tx), start(tx))
tss_gr  <- GRanges(seqnames(tx), IRanges(tss, tss), strand = strand(tx))
exon_gr <- gff[gff$type == "exon"]
gene_id_all <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
disp <- function(gg) { s <- de[g == gg, symbol][1]
  if (is.na(s) || s == "" || s == "Non annotated gene") gg else s }

# Per DMR: TSS inside the DMR, TSS inside the host intron, distance to the nearest TSS
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

# showOneDMR draws text at a set size; a small canvas and a one-line header keep it legible
draw_one <- function(i) {
  r <- sel[i]; a <- annot[i]
  par(oma = c(0, 0, 2.4, 0))
  DSS::showOneDMR(data.frame(chr = r$chr, start = r$start, end = r$end), bs, ext = 2000)
  mtext(sprintf("%s  %s %s DMR  |  DE log2FC %+.2f",
                a$symbol, tolower(r$region), r$direction, a$log2FC),
        outer = TRUE, side = 3, line = 0.4, cex = 0.85, font = 2)
}
# (a) one multi-page PDF with all intersection DMRs
pdf(file.path(FIGS, "figS6_intersection_dmrs.pdf"), width = 6.5, height = 5)
for (i in seq_len(nrow(sel))) draw_one(i)
dev.off()
# (b) per-DMR png and svg
for (i in seq_len(nrow(sel))) {
  nm <- sprintf("figS6_dmr_%s", gsub("[^A-Za-z0-9]+", "_", annot$symbol[i]))
  if (sum(annot$symbol == annot$symbol[i]) > 1) nm <- sprintf("%s_%d", nm, i)   # two DMRs in one gene keep separate files
  png(file.path(FIGS, paste0(nm, ".png")), width = 6.5, height = 5, units = "in", res = 150); draw_one(i); dev.off()
  svglite::svglite(file.path(FIGS, paste0(nm, ".svg")), width = 6.5, height = 5); draw_one(i); dev.off()
}
cat(sprintf("  %d intersection DMRs plotted (TSS-in-intron: %d) -> figS6_*\n",
            nrow(sel), sum(annot$tss_in_host_intron > 0)))

# Step 8 - Supplementary: HCP (CpG-rich) promoters carrying a DMR
cat("[5] HCP promoters carrying a DMR (DSS::showOneDMR)\n")
# HCP = CpG-rich promoter class from 03_promoters (Weber classification)
B03D <- file.path(PIPE, "03_promoters/data")
weber <- fread(file.path(B03D, "promoter_weber_classification.tsv"))
hcp_ids <- weber[biotype == "protein_coding" & weber_class == "HCP", gene_id]
hdmr <- dmr_all[region == "Promoter" & g %in% hcp_ids]
hdmr <- merge(hdmr, de[, .(g, baseMean, log2FoldChange, padj)], by = "g", all.x = TRUE)
hdmr[, symbol := vapply(g, disp, character(1))]
hdmr[, dMeth := -diff.Methy]   # amputated minus control, as in Step 5
hdmr[, ad := abs(dMeth)]; setorder(hdmr, -ad); hdmr[, ad := NULL]   # largest change first
fwrite(hdmr[, .(symbol, gene_id = g, chr, start, end, dMeth_amp_minus_ctrl = dMeth, direction,
                baseMean, log2FoldChange, padj)],
       file.path(DAT, "hcp_promoter_dmr_genes.tsv"), sep = "\t")
cat(sprintf("  %d of %d HCP promoters carry a DMR\n", nrow(hdmr), length(hcp_ids)))
if (nrow(hdmr)) {
  draw_hcp <- function(i) {   # same layout as draw_one
    r <- hdmr[i]; par(oma = c(0, 0, 2.4, 0))
    DSS::showOneDMR(data.frame(chr = r$chr, start = r$start, end = r$end), bs, ext = 2000)
    mtext(sprintf("%s  %s promoter DMR  |  dMeth (amp - ctrl) %+.2f",
                  r$symbol, r$direction, r$dMeth),
          outer = TRUE, side = 3, line = 0.4, cex = 0.85, font = 2)
  }
  pdf(file.path(FIGS, "figS6_hcp_promoter_dmrs.pdf"), width = 6.5, height = 5)   # one page per DMR
  for (i in seq_len(nrow(hdmr))) draw_hcp(i)
  dev.off()
  for (i in seq_len(nrow(hdmr))) {   # single-page png and svg
    nm <- sprintf("figS6_hcp_dmr_%s", gsub("[^A-Za-z0-9]+", "_", hdmr$symbol[i]))
    if (sum(hdmr$symbol == hdmr$symbol[i]) > 1) nm <- sprintf("%s_%d", nm, i)   # two DMRs in one gene keep separate files
    png(file.path(FIGS, paste0(nm, ".png")), width = 6.5, height = 5, units = "in", res = 150); draw_hcp(i); dev.off()
    svglite::svglite(file.path(FIGS, paste0(nm, ".svg")), width = 6.5, height = 5); draw_hcp(i); dev.off()
  }
  cat(sprintf("  drew %d HCP promoter DMRs -> figS6_hcp_promoter_dmrs.pdf + figS6_hcp_dmr_*\n", nrow(hdmr)))
}

# Step 9 - Record package versions (sessionInfo)
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_06_decoupling.txt"))
cat("[06_decoupling] done\n")
