#!/usr/bin/env Rscript
# 07_wgcna.R
# Multi-tissue WGCNA co-expression network with a methylation overlay: module-level DMP/DMR
# enrichment (Fisher, with a length-stratified sensitivity test), GO/KEGG of the enriched modules, tail eigengene scores,
# and hub-gene analyses.
# Inputs: HTSeq gene counts (all tissues); 01_genome_toolkit (gff_chrmt.rds, gene_de_tail.tsv);
#   02_landscape (gene-body methylation table, bsseq object); 03_promoters (Weber classes);
#   05_differential (DMP and DMR tables); STRING v12 terms; JASPAR2024 sqlite; eggNOG annotations.
# Outputs: data/*.tsv (sample map, module assignments, enrichment and hub tables);
#   objects/wgcna.rds (cached network); figures/main/fig7{a,b,c,d}_*; figures/supplementary/fig7_s*, figS7_*
# Run: sbatch 07_wgcna/code/07_wgcna.slurm from main/methylation_pipeline/

# Step 1 - Setup: seed, packages, WGCNA thread count
set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(WGCNA)
  library(GenomicRanges); library(ggplot2)
})
options(stringsAsFactors = FALSE)
# Thread count follows the SLURM allocation; falls back to 4 outside SLURM
enableWGCNAThreads(nThreads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

# Step 2 - Input paths and output directories
PIPE   <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
HTSEQ  <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"
STRING <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/STRING.protein.enrichment.terms.v12.0.txt"
GFF_RDS <- file.path(PIPE, "01_genome_toolkit/objects/gff_chrmt.rds")
DMP_TSV <- file.path(PIPE, "05_differential/data/dmps_annotated.tsv")
DMR_TSV <- file.path(PIPE, "05_differential/data/dmrs_annotated.tsv")
# JASPAR is queried for TF names only (Step 25); all motif scanning is done in 08_motifs
JASPAR_SQLITE <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/jaspar/JASPAR2024.sqlite"
EMAPPER <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/eggnog_mapper/dlasi_proteome.emapper.annotations"
BATCH <- file.path(PIPE, "07_wgcna")
OBJ  <- file.path(BATCH, "objects")
DAT  <- file.path(BATCH, "data")
FIGM <- file.path(BATCH, "figures/main")
FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(OBJ, DAT, FIGM, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
# Figure directories are emptied so no stale panel survives a rerun
unlink(list.files(c(FIGM, FIGS), full.names = TRUE))

# Step 3 - Shared plot theme and figure savers (pdf, png, svg)
COL_COND <- c(Control = "#0072B2", Amputated = "#D55E00")  # Okabe-Ito, paper-wide
theme_pub <- function() theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"))
save_gg <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)  # cairo renders the Greek glyphs
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(dir, paste0(name, ".svg")), p, width = w, height = h)
  cat(sprintf("  saved %s\n", name))
}
# Base-graphics saver for the WGCNA dendrograms and heatmaps
save_base <- function(dir, name, w, h, FUN) {
  pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h); FUN(); dev.off()
  png(file.path(dir, paste0(name, ".png")), width = w, height = h,
      units = "in", res = 150); FUN(); dev.off()
  svglite::svglite(file.path(dir, paste0(name, ".svg")), width = w, height = h); FUN(); dev.off()
  cat(sprintf("  saved %s\n", name))
}

# Step 4 - Sample map reconstructed from the HTSeq count file names
cat("[1] sample map\n")
files <- list.files(HTSEQ, pattern = "_htseq_gene_counts.txt$")
base  <- sub("\\.Aligned\\.out\\.bam_htseq_gene_counts\\.txt$", "", files)
base  <- sub("_htseq_gene_counts\\.txt$", "", base)
# Library prefix -> tissue, condition, experiment group. The group (third element) drives the
# metadata and colours: 'condition' alone would merge controls of different experiments
classify <- function(b) {  # order matters (C#S# before C#)
  if (grepl("^C[0-9]+S[0-9]+", b))  return(c("tail",      "control",    "TailControl"))
  if (grepl("^T[0-9]+S[0-9]+", b))  return(c("tail",      "amputated",  "TailAmputated"))
  if (grepl("^R[0-9]+$", b))        return(c("eye",       "amputated",  "EyeAmputated"))
  if (grepl("^C[0-9]+$", b))        return(c("eye",       "control",    "EyeControl"))
  if (grepl("^dcrep", b))           return(c("bodywall",  "control",    "IrradiatedControl"))
  if (grepl("^fungicide_l0", b))    return(c("bodywall",  "control",    "FungicideControl"))
  if (grepl("^fungicide_l30", b))   return(c("bodywall",  "fungicide",  "FungicideTreated"))
  if (grepl("^irrep", b))           return(c("bodywall",  "irradiated", "IrradiatedTreated"))
  if (grepl("Head", b))             return(c("head",      "control",    "Head"))
  if (grepl("Juv", b))              return(c("juvenile",  "control",    "Juvenile"))
  if (grepl("Ovo", b))              return(c("ovotestis", "control",    "Ovotestis"))
  c("other", "unknown", "Other")
}
cls  <- t(vapply(base, classify, character(3)))
meta <- data.table(file = files, orig_sample = base,
                   tissue = cls[, 1], condition = cls[, 2], group = cls[, 3])
meta <- meta[tissue != "other"]  # keep only recognised libraries
# Readable sample names <Group><n>, numbered within each group in library-name order
setorder(meta, group, orig_sample)
meta[, sample := paste0(group, seq_len(.N)), by = group]
fwrite(meta, file.path(DAT, "sample_map.tsv"), sep = "\t")
cat(sprintf("  %d libraries across %d tissues, %d experiment groups\n",
            nrow(meta), uniqueN(meta$tissue), uniqueN(meta$group)))
print(meta[, .N, by = group][order(group)])
# Plain data.frame with row names so samples can be indexed by name below
meta <- as.data.frame(meta); rownames(meta) <- meta$sample

# Step 5 - Count matrix, gene universe, low-count filter and VST
cat("[2] counts + VST\n")
cl  <- lapply(meta$file, function(f) {
  x <- fread(file.path(HTSEQ, f), header = FALSE, col.names = c("gene_id", "count"))
  x[!startsWith(gene_id, "__")]  # drop HTSeq __no_feature etc.
})
ids <- cl[[1]]$gene_id
cm  <- sapply(cl, function(x) x$count[match(ids, x$gene_id)])
rownames(cm) <- ids; colnames(cm) <- meta$sample

# Gene universe from the 01_genome_toolkit GFF (chr1-31 + mito), so the chromosome filter is inherited
gff <- readRDS(GFF_RDS)
chr_genes <- unique(gff$ID[gff$type == "gene"])
cm <- cm[rownames(cm) %in% chr_genes, , drop = FALSE]
cat(sprintf("  %d genes on chr1-31\n", nrow(cm)))

cm  <- cm[rowSums(cm >= 10) >= 3, , drop = FALSE]  # >=10 reads in >=3 libraries
cat(sprintf("  %d genes pass the low-count filter\n", nrow(cm)))
vsd <- vst(DESeqDataSetFromMatrix(cm, meta[colnames(cm), ], ~ 1), blind = TRUE)
datExpr <- t(assay(vsd))  # samples x genes

# Step 6 - Sample clustering and outlier removal
cat("[3] sample clustering / outliers\n")
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# Colour by experiment group, not condition (see Step 4)
trait_colors <- function(m) data.frame(
  tissue = labels2colors(as.numeric(factor(m$tissue))),
  group  = labels2colors(as.numeric(factor(m$group))))

tree1 <- hclust(dist(datExpr), method = "average")  # Euclidean distance, average linkage
save_base(FIGS, "fig7_s1_sample_clustering", 12, 6, function()
  plotDendroAndColors(tree1, trait_colors(meta[rownames(datExpr), ]),
    groupLabels = c("tissue", "group"),
    main = "Sample clustering (outliers marked)", cex.dendroLabels = 0.7))

# Outliers are given by original library name and translated to the new sample names
outlier_orig <- c("T1S5", "dcrep4", "R6", "irrep7")
outliers <- meta$sample[match(outlier_orig, meta$orig_sample)]
stopifnot(!anyNA(outliers))
cat(sprintf("  removing outliers: %s\n",
            paste(sprintf("%s (was %s)", outliers, outlier_orig), collapse = ", ")))
datExpr <- datExpr[!rownames(datExpr) %in% outliers, ]
meta    <- meta[rownames(datExpr), ]

tree2 <- hclust(dist(datExpr), method = "average")
save_base(FIGS, "fig7_s2_sample_clustering_no_outliers", 12, 6, function()
  plotDendroAndColors(tree2, trait_colors(meta[rownames(datExpr), ]),
    groupLabels = c("tissue", "group"),
    main = "Sample clustering (outliers removed)", cex.dendroLabels = 0.7))

# Step 7 - Variance pre-filter sweep: largest cutoff that keeps scale-free topology
# Per-gene variance after outlier removal; for each candidate quantile the soft-threshold fit is
# recomputed and the largest quantile with R^2 >= 0.85 at some power <= 20 is chosen
cat("[3b] variance pre-filter sweep (largest cutoff keeping scale-free topology)\n")
SFT_CUT  <- 0.85
powers   <- c(1:10, seq(12, 20, 2))
gene_var <- apply(datExpr, 2, var)
q_grid   <- seq(0.20, 0.60, by = 0.05)
sweep <- rbindlist(lapply(q_grid, function(q) {
  keep <- gene_var > quantile(gene_var, q)
  s  <- pickSoftThreshold(datExpr[, keep, drop = FALSE], powerVector = powers,
                          networkType = "signed", verbose = 0)
  fi <- s$fitIndices
  ok <- which(fi$SFT.R.sq >= SFT_CUT)
  data.table(q = q, n_genes = sum(keep),
             best_r2 = max(fi$SFT.R.sq, na.rm = TRUE),
             power   = if (length(ok)) fi$Power[ok[1]] else NA_integer_)
}))
fwrite(sweep, file.path(DAT, "variance_filter_sweep.tsv"), sep = "\t")
print(sweep)
good <- sweep[!is.na(power)]
if (!nrow(good)) stop("no variance cutoff reaches the scale-free criterion (R^2 >= 0.85)")
sel     <- good[which.max(q)]  # most aggressive filter that still works
var_q   <- sel$q
var_cut <- quantile(gene_var, var_q)
cat(sprintf("  chosen variance quantile = %.2f -> %d genes (best R2 = %.3f at power %d)\n",
            var_q, sel$n_genes, sel$best_r2, sel$power))
save_base(FIGS, "fig7_s2b_variance_cutoff", 10, 4.5, function() {
  par(mfrow = c(1, 2))
  hist(log10(gene_var), breaks = 50, col = "lightblue", border = "white",
       main = "Gene-variance distribution (log10)", xlab = "log10(variance)")
  abline(v = log10(var_cut), col = "red", lty = 2, lwd = 2)
  legend("topright", legend = sprintf("chosen cutoff (q = %.2f)", var_q),
         col = "red", lty = 2, lwd = 2, bty = "n")
  plot(sweep$q, sweep$best_r2, type = "b", pch = 19, col = "grey30",
       xlab = "variance quantile dropped", ylab = "best scale-free R^2",
       main = "Scale-free fit vs variance filter",
       ylim = range(c(sweep$best_r2, SFT_CUT), na.rm = TRUE))
  abline(h = SFT_CUT, col = "red", lty = 2)
  points(var_q, sel$best_r2, col = "red", pch = 19, cex = 1.6)
  legend("bottomleft", legend = c("scale-free criterion", "chosen"),
         col = "red", lty = c(2, NA), pch = c(NA, 19), bty = "n")
})
keep_var <- gene_var > var_cut  # strictly greater
cat(sprintf("  %d -> %d genes after variance filter\n", ncol(datExpr), sum(keep_var)))
expr_all <- datExpr  # pre-variance-filter matrix, used in Step 14
datExpr <- datExpr[, keep_var, drop = FALSE]

# Step 8 - Soft-threshold power
cat("[4] soft threshold\n")
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
soft_power <- sft$powerEstimate
if (is.na(soft_power)) soft_power <- 14  # fallback when no power reaches R^2 0.85
cat(sprintf("  chosen soft power = %d\n", soft_power))
fi <- sft$fitIndices
save_base(FIGS, "fig7_s3_soft_threshold", 10, 5, function() {
  par(mfrow = c(1, 2))
  plot(fi[, 1], -sign(fi[, 3]) * fi[, 2], type = "n",
       xlab = "Soft threshold (power)", ylab = expression("Scale-free topology " * R^2),
       main = "Scale independence")
  text(fi[, 1], -sign(fi[, 3]) * fi[, 2], labels = powers, col = "red", cex = 0.9)
  abline(h = 0.85, col = "red", lty = 2)
  plot(fi[, 1], fi[, 5], type = "n",
       xlab = "Soft threshold (power)", ylab = "Mean connectivity",
       main = "Mean connectivity")
  text(fi[, 1], fi[, 5], labels = powers, col = "red", cex = 0.9)
})

# Step 9 - Module detection (blockwiseModules), cached in objects/wgcna.rds
# The cache is reused only if gene identity, matrix dimensions, matrix sum and soft power match
cat("[5] modules (reuse cached wgcna.rds if present)\n")
wgcna_rds <- file.path(OBJ, "wgcna.rds")
if (file.exists(wgcna_rds)) {
  W <- readRDS(wgcna_rds); net <- W$net; modColors <- W$modColors; MEs <- W$MEs
  stopifnot(identical(W$genes, colnames(datExpr)), identical(W$dim, dim(datExpr)),
            isTRUE(all.equal(W$expr_sum, sum(datExpr))), identical(W$soft_power, soft_power))
  cat("  reused cached network — blockwiseModules skipped (fast)\n")
} else {
  # blockwiseModules calls cor() unqualified; WGCNA::cor is masked in temporarily
  cor <- WGCNA::cor
  net <- blockwiseModules(datExpr,
    power = soft_power, networkType = "signed", TOMType = "signed",
    minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    maxBlockSize = 25000, saveTOMs = FALSE, verbose = 0)  # single block
  cor <- stats::cor
  stopifnot(length(net$dendrograms) == 1)  # if this fails, raise maxBlockSize
  modColors <- labels2colors(net$colors)
  MEs <- orderMEs(moduleEigengenes(datExpr, modColors)$eigengenes)
  MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
  saveRDS(list(net = net, modColors = modColors, MEs = MEs, genes = colnames(datExpr), meta = meta,
               soft_power = soft_power, dim = dim(datExpr), expr_sum = sum(datExpr)), wgcna_rds)
}
mod_dt <- data.table(gene_id = colnames(datExpr), module = modColors)
fwrite(mod_dt, file.path(DAT, "module_assignments.tsv"), sep = "\t")
mods <- setdiff(sort(unique(modColors)), "grey")  # real modules (grey = unassigned)
cat(sprintf("  %d modules (+ grey); %d genes unassigned\n",
            length(mods), sum(modColors == "grey")))

# Gene dendrogram with module colours
save_base(FIGS, "fig7_s4_module_dendrogram", 10, 6, function()
  plotDendroAndColors(net$dendrograms[[1]], modColors[net$blockGenes[[1]]],
    "Module", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
    guideHang = 0.05, main = "Gene dendrogram and module colours"))

# Step 10 - Module DMP-burden enrichment (Fisher, unadjusted)
cat("[6] module DMP-burden Fisher test\n")
dmp_genes <- unique(fread(DMP_TSV)$gene_id)
dmp_genes <- dmp_genes[!is.na(dmp_genes) & dmp_genes != ""]
univ <- mod_dt$gene_id  # universe = ALL network genes (grey/unassigned included)
enr <- rbindlist(lapply(mods, function(m) {
  inmod <- mod_dt$gene_id[mod_dt$module == m]
  a <- sum(inmod %in% dmp_genes)  # in-module & DMP
  b <- length(inmod) - a  # in-module & not-DMP
  c_ <- sum(univ %in% dmp_genes) - a  # out-module & DMP
  d <- length(univ) - length(inmod) - c_  # out-module & not-DMP
  # STAT TEST: two-sided Fisher's exact test on the 2x2 (module members vs rest of network)
  ft <- fisher.test(matrix(c(a, b, c_, d), 2))
  data.table(module = m, n_genes = length(inmod), n_dmp = a,
             OR = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2],
             p = ft$p.value)
}))
enr[, fdr := p.adjust(p, "BH")]
enr <- enr[order(-OR)]
fwrite(enr, file.path(DAT, "module_dmp_enrichment.tsv"), sep = "\t")

# Step 11 - Length-stratified CMH test per module (sensitivity analysis) and main panel fig7a
# Longer genes carry more CpGs and collect DMPs by size alone, so the 2x2 test is repeated as a
# Cochran-Mantel-Haenszel common odds ratio across gene-length quintiles
cat("[6b] module DMP enrichment, stratified by gene-length quintile\n")
gn8b <- gff[gff$type == "gene"]
glen <- data.table(gene_id = sub(";.*", "", as.character(gn8b$ID)), len = width(gn8b))
mod_len <- merge(mod_dt[, .(gene_id, module)], glen, by = "gene_id")
mod_len[, lenq := cut(len, quantile(len, seq(0, 1, 0.2)), include.lowest = TRUE, labels = FALSE)]
mod_len[, is_dmp := gene_id %in% dmp_genes]
enr_adj <- rbindlist(lapply(mods, function(m) {
  x <- copy(mod_len)[, inmod := module == m]
  tb <- table(factor(x$inmod, c(TRUE, FALSE)), factor(x$is_dmp, c(TRUE, FALSE)), x$lenq)
  # Strata with an empty margin cannot enter the CMH test
  keep <- apply(tb, 3, function(s) all(rowSums(s) > 0) && all(colSums(s) > 0))
  if (sum(keep) < 2) return(data.table(module = m, OR_adj = NA_real_, lo_adj = NA_real_,
                                       hi_adj = NA_real_, p_adj_test = NA_real_))
  ct <- mantelhaen.test(tb[, , keep, drop = FALSE])
  data.table(module = m, OR_adj = unname(ct$estimate), lo_adj = ct$conf.int[1],
             hi_adj = ct$conf.int[2], p_adj_test = ct$p.value)
}))
enr_adj[, fdr_adj := p.adjust(p_adj_test, "BH")]
enr2 <- merge(enr, enr_adj, by = "module")[order(-OR)]
fwrite(enr2, file.path(DAT, "module_dmp_enrichment_length_adjusted.tsv"), sep = "\t")
cat("  module            raw OR   raw FDR | length-adj OR  adj FDR   survives?\n")
for (i in seq_len(nrow(enr2))) with(enr2[i], cat(sprintf(
  "  %-16s %6.2f %9.2g | %11.2f %9.2g   %s\n", module, OR, fdr, OR_adj, fdr_adj,
  if (!is.na(fdr_adj) && fdr_adj < 0.05) "YES" else "no")))

# fig7a plots the unadjusted Fisher test (odds ratio, 95% CI, FDR); the length-stratified
# CMH in enr2 is kept only as a sensitivity table
enr_fig <- copy(enr)
setorder(enr_fig, OR)
enr_fig[, module := factor(module, levels = module)]
p_fisher <- ggplot(enr_fig, aes(module, OR)) +
  geom_segment(aes(xend = module, y = 1, yend = OR), colour = "grey85", linewidth = 0.3) +
  geom_segment(aes(xend = module, y = lo, yend = hi), colour = "grey45", linewidth = 0.45) +
  geom_point(aes(size = n_dmp, colour = fdr < 0.05), alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  scale_colour_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey65"), name = "FDR < 0.05") +
  scale_size_continuous(name = "Genes with ≥1 DMP", range = c(2, 9)) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_flip() +
  labs(x = "WGCNA module",
       y = "Odds ratio and 95% CI, Fisher's exact test\n(module vs rest of network)",
       title = "DMP enrichment per WGCNA module") +
  theme_pub() + theme(axis.text.y = element_text(size = 7),
                      axis.title.x = element_text(size = 8),
                      plot.title = element_text(size = 9, face = "bold"))
save_gg(p_fisher, FIGM, "fig7a_module_dmp_fisher", 4.2, 2.8)

# Step 12 - CMH stratified by gene length x baseline methylation, and by movable CpGs
# A CpG can only change by >= 0.10 if it is not at 0 or 1; 'movable' CpGs have pooled beta in
# [0.10, 0.90] within the gene body + 2 kb upstream (the DMP assignment window)
cat("[6c] module DMP enrichment stratified by baseline methylation and by movable CpGs\n")
GB_TSV <- file.path(PIPE, "02_landscape/data/genebody_methylation_per_gene.tsv")
BSSEQ  <- file.path(PIPE, "02_landscape/objects/bsseq_cov5_chrmt.rds")
gb <- fread(GB_TSV)[, .(gene_id, beta)]
suppressPackageStartupMessages(library(bsseq))
bs7 <- readRDS(BSSEQ)
beta_pool <- rowSums(as.matrix(getCoverage(bs7, type = "M"))) / pmax(rowSums(as.matrix(getCoverage(bs7, type = "Cov"))), 1)
cpg7 <- granges(bs7); rm(bs7); invisible(gc())
mov <- cpg7[beta_pool >= 0.10 & beta_pool <= 0.90]  # CpGs that can move by >= 0.10 either way
win7 <- gn8b; s7 <- as.character(strand(gn8b))  # gene body + 2 kb upstream
start(win7)[s7 != "-"] <- pmax(1L, start(win7)[s7 != "-"] - 2000L)
end(win7)[s7 == "-"] <- end(win7)[s7 == "-"] + 2000L
win7 <- suppressWarnings(trim(win7))  # clip 2 kb overhangs at chromosome ends
movdt <- data.table(gene_id = sub(";.*", "", as.character(gn8b$ID)),
                    n_covered = countOverlaps(win7, cpg7), n_movable = countOverlaps(win7, mov))
fwrite(movdt, file.path(DAT, "movable_cpgs_per_gene.tsv"), sep = "\t")
ml <- merge(mod_len, gb, by = "gene_id", all.x = TRUE)
ml <- merge(ml, movdt, by = "gene_id", all.x = TRUE)
ml[!is.na(beta), betaq := cut(beta, unique(quantile(beta, seq(0, 1, 0.2))), include.lowest = TRUE, labels = FALSE)]
ml[!is.na(betaq), len_beta := (lenq - 1L) * 5L + betaq]  # 25 joint strata
ml[!is.na(n_movable) & n_covered > 0, movq := cut(n_movable, unique(quantile(n_movable, seq(0, 1, 0.2))), include.lowest = TRUE, labels = FALSE)]
cmh_by <- function(col, label) {
  out <- rbindlist(lapply(mods, function(m) {
    x <- ml[!is.na(get(col))]; x[, inmod := module == m]
    tb <- table(factor(x$inmod, c(TRUE, FALSE)), factor(x$is_dmp, c(TRUE, FALSE)), x[[col]])
    keep <- apply(tb, 3, function(s) all(rowSums(s) > 0) && all(colSums(s) > 0))
    if (sum(keep) < 2) return(data.table(module = m, stratification = label, OR = NA_real_, lo = NA_real_,
                                         hi = NA_real_, p = NA_real_, n_genes = nrow(x), n_strata = sum(keep)))
    # STAT TEST: Cochran-Mantel-Haenszel common odds ratio (mantelhaen.test) across the strata
    ct <- mantelhaen.test(tb[, , keep, drop = FALSE])
    data.table(module = m, stratification = label, OR = unname(ct$estimate), lo = ct$conf.int[1],
               hi = ct$conf.int[2], p = ct$p.value, n_genes = nrow(x), n_strata = sum(keep))
  }))
  out[, fdr := p.adjust(p, "BH")]; out
}
strat <- rbind(
  enr_adj[, .(module, stratification = "Gene length (quintiles)", OR = OR_adj, lo = lo_adj, hi = hi_adj,
              p = p_adj_test, n_genes = nrow(mod_len), n_strata = NA_integer_, fdr = fdr_adj)],
  cmh_by("len_beta", "Gene length x gene-body methylation (5 x 5)"),
  cmh_by("movq", "Movable CpGs per gene (quintiles)"))
fwrite(strat, file.path(DAT, "module_dmp_enrichment_length_meth_adjusted.tsv"), sep = "\t")
cat("  module            length-adj OR  FDR   | length x beta OR  FDR   | movable-CpG OR  FDR\n")
for (m in enr2$module) {
  r1 <- strat[module == m & stratification == "Gene length (quintiles)"]
  r2 <- strat[module == m & stratification == "Gene length x gene-body methylation (5 x 5)"]
  r3 <- strat[module == m & stratification == "Movable CpGs per gene (quintiles)"]
  cat(sprintf("  %-16s %8.2f %8.2g | %8.2f %8.2g | %8.2f %8.2g\n", m, r1$OR, r1$fdr, r2$OR, r2$fdr, r3$OR, r3$fdr))
}
# Forest plot of all three stratifications, modules ordered by the length-adjusted OR
sf <- strat[!is.na(OR)]
sf[, module := factor(module, levels = levels(enr_fig$module))]
sf[, stratification := factor(stratification, levels = unique(strat$stratification))]
p_forest <- ggplot(sf, aes(OR, module, colour = stratification, group = stratification)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), linewidth = 0.45) +
  geom_point(aes(shape = fdr < 0.05), position = position_dodge(width = 0.6), size = 1.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "FDR < 0.05") +
  scale_colour_manual(values = c("#2C7FB8", "#D55E00", "#009E73"), name = "Stratified by") +
  scale_x_log10() +
  labs(x = "CMH common odds ratio, DMP-bearing genes (95% CI)", y = "WGCNA module",
       title = "Module DMP enrichment under three stratifications") +
  theme_pub() + theme(axis.text.y = element_text(size = 7), legend.position = "right",
                      legend.text = element_text(size = 7), plot.title = element_text(size = 9, face = "bold"))
save_gg(p_forest, FIGS, "figS7_module_dmp_enrichment_forest", 7.0, 4.6)

# Step 13 - Enrichment split by DMP direction (hyper / hypo), length-stratified CMH
cat("[6d] module DMP enrichment split by direction (length-stratified CMH)\n")
dmp_dir <- fread(DMP_TSV)[!is.na(gene_id) & gene_id != "", .(gene_id, direction)]
dir_sets <- list(Hyper = unique(dmp_dir[direction == "Hyper", gene_id]),
                 Hypo  = unique(dmp_dir[direction == "Hypo",  gene_id]))
enr_dir <- rbindlist(lapply(names(dir_sets), function(dn) {
  x0 <- copy(mod_len)[, is_dir := gene_id %in% dir_sets[[dn]]]
  rbindlist(lapply(mods, function(m) {
    x <- copy(x0)[, inmod := module == m]
    tb <- table(factor(x$inmod, c(TRUE, FALSE)), factor(x$is_dir, c(TRUE, FALSE)), x$lenq)
    keep <- apply(tb, 3, function(s) all(rowSums(s) > 0) && all(colSums(s) > 0))
    if (sum(keep) < 2) return(data.table(module = m, direction = dn, n_genes_dir = sum(x$is_dir),
                                         n_in_module = sum(x$inmod & x$is_dir), OR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
    # STAT TEST: Cochran-Mantel-Haenszel common odds ratio across gene-length quintiles
    ct <- mantelhaen.test(tb[, , keep, drop = FALSE])
    data.table(module = m, direction = dn, n_genes_dir = sum(x$is_dir), n_in_module = sum(x$inmod & x$is_dir),
               OR = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2], p = ct$p.value)
  }))
}))
enr_dir[, fdr := p.adjust(p, "BH"), by = direction]
fwrite(enr_dir[order(direction, -OR)], file.path(DAT, "module_dmp_enrichment_by_direction.tsv"), sep = "\t")
cat("  module            hyper OR  FDR     | hypo OR  FDR\n")
for (m in enr2$module) {
  h <- enr_dir[module == m & direction == "Hyper"]; l <- enr_dir[module == m & direction == "Hypo"]
  cat(sprintf("  %-16s %7.2f %8.2g | %7.2f %8.2g\n", m, h$OR, h$fdr, l$OR, l$fdr))
}
ed <- enr_dir[!is.na(OR)]
ed[, module := factor(module, levels = levels(enr_fig$module))]
p_dir <- ggplot(ed, aes(OR, module, colour = direction, group = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55), linewidth = 0.45) +
  geom_point(aes(shape = fdr < 0.05), position = position_dodge(width = 0.55), size = 1.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "FDR < 0.05") +
  scale_colour_manual(values = c(Hyper = "#C0392B", Hypo = "#0072B2"), name = "DMP direction") +
  scale_x_log10() +
  labs(x = "CMH common odds ratio, stratified by gene length (95% CI)", y = "WGCNA module",
       title = "Module enrichment of hyper- and hypomethylated DMP genes") +
  theme_pub() + theme(axis.text.y = element_text(size = 7), legend.text = element_text(size = 7),
                      plot.title = element_text(size = 9, face = "bold"))
save_gg(p_dir, FIGS, "figS7_module_dmp_direction", 6.4, 4.4)

# Step 14 - Expression breadth (tau) across the atlas vs gene-body methylation
# tau = sum(1 - x_i / max x) / (n - 1) over mean VST per experiment group, on values shifted to a
# zero floor (0 = uniform, 1 = one group only); computed on the pre-variance-filter matrix
cat("[6e] expression breadth (tau) vs gene-body methylation\n")
grp <- meta[rownames(expr_all), "group"]
gm  <- t(apply(expr_all, 2, function(v) tapply(v, grp, mean)))  # genes x groups
gm0 <- gm - min(gm, na.rm = TRUE)  # zero floor for tau
tau <- apply(gm0, 1, function(v) { mx <- max(v); if (!is.finite(mx) || mx <= 0) return(NA_real_); sum(1 - v / mx) / (length(v) - 1) })
cv  <- apply(expr_all, 2, function(v) sd(v) / mean(v))
br <- data.table(gene_id = colnames(expr_all), tau = tau, cv = cv)
br <- merge(br, gb, by = "gene_id")
br[, is_dmp := gene_id %in% dmp_genes]
br[, module := mod_dt$module[match(gene_id, mod_dt$gene_id)]]
# STAT TEST: Spearman rank correlation of tau (and CV) with gene-body beta; Mann-Whitney U
# (rank-biserial r) of tau in DMP-bearing vs other genes.
ct_tau <- cor.test(br$tau, br$beta, method = "spearman", exact = FALSE)
ct_cv  <- cor.test(br$cv,  br$beta, method = "spearman", exact = FALSE)
mw <- wilcox.test(tau ~ is_dmp, data = br)
rb <- 1 - 2 * unname(mw$statistic) / (sum(!br$is_dmp) * sum(br$is_dmp))
br[, tau_decile := cut(tau, quantile(tau, seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
breadth_summ <- data.table(
  n_genes = nrow(br), spearman_tau_beta = unname(ct_tau$estimate), p_tau = ct_tau$p.value,
  spearman_cv_beta = unname(ct_cv$estimate), p_cv = ct_cv$p.value,
  median_tau_dmp = median(br[is_dmp == TRUE, tau], na.rm = TRUE), median_tau_other = median(br[is_dmp == FALSE, tau], na.rm = TRUE),
  mannwhitney_p = mw$p.value, rank_biserial_r = rb)
fwrite(breadth_summ, file.path(DAT, "expression_breadth_vs_genebody_methylation.tsv"), sep = "\t")
fwrite(br[, .(gene_id, module, tau, cv, beta, is_dmp)], file.path(DAT, "expression_breadth_per_gene.tsv"), sep = "\t")
fwrite(br[!is.na(tau_decile), .(n = .N, median_beta = median(beta), mean_beta = mean(beta),
                                 pct_dmp = 100 * mean(is_dmp)), by = tau_decile][order(tau_decile)],
       file.path(DAT, "expression_breadth_deciles.tsv"), sep = "\t")
print(breadth_summ)
p_br1 <- ggplot(br[!is.na(tau_decile)], aes(factor(tau_decile), beta)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, fill = "grey92") +
  labs(x = "Expression breadth decile (tau; 1 = broadest, 10 = most group specific)", y = "Gene-body methylation (β)",
       title = sprintf("Spearman rho = %.2f", ct_tau$estimate)) + theme_pub()
p_br2 <- ggplot(br, aes(is_dmp, tau, fill = is_dmp)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, width = 0.55) +
  scale_x_discrete(labels = c(`FALSE` = "Other genes", `TRUE` = "DMP-bearing")) +
  scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = "#D55E00"), guide = "none") +
  labs(x = NULL, y = "Expression breadth (tau)", title = sprintf("rank-biserial r = %.2f, P = %.2g", rb, mw$p.value)) +
  theme_pub()
suppressPackageStartupMessages(library(patchwork))  # the `+` between two ggplots needs patchwork attached
save_gg(p_br1 + p_br2 + plot_layout(widths = c(1.6, 1)), FIGS, "figS7_expression_breadth_vs_genebody_meth", 7.2, 3.0)

# Step 15 - Raw Fisher vs length-adjusted CMH side by side, and why gene length matters
cat("[6f] raw Fisher vs length-adjusted CMH per module + why\n")
cmp <- rbind(enr2[, .(module, test = "Fisher, unadjusted", OR, lo, hi, fdr)],
             enr2[, .(module, test = "CMH, gene-length quintiles", OR = OR_adj, lo = lo_adj, hi = hi_adj, fdr = fdr_adj)])
cmp <- cmp[!is.na(OR)]
med_len <- mod_len[, .(median_len = as.numeric(median(len)), n_genes = .N, pct_dmp = 100 * mean(is_dmp)), by = module]
eff <- merge(enr2[, .(module, OR_raw = OR, fdr_raw = fdr, OR_adj, fdr_adj)], med_len, by = "module")
eff[, verdict_raw := fifelse(fdr_raw < 0.05, fifelse(OR_raw > 1, "enriched", "depleted"), "n.s.")]
eff[, verdict_adj := fifelse(!is.na(fdr_adj) & fdr_adj < 0.05, fifelse(OR_adj > 1, "enriched", "depleted"), "n.s.")]
eff[, verdict_changes := verdict_raw != verdict_adj]
fwrite(eff[order(-OR_adj)], file.path(DAT, "length_adjustment_effect.tsv"), sep = "\t")
cat(sprintf("  verdict changes with the length adjustment for %d of %d modules: %s\n", sum(eff$verdict_changes), nrow(eff),
            paste(eff[verdict_changes == TRUE, sprintf("%s (%s -> %s)", module, verdict_raw, verdict_adj)], collapse = "; ")))
share_q <- mod_len[, .(pct_with_dmp = 100 * mean(is_dmp), n = .N, median_len = as.numeric(median(len))), by = lenq][order(lenq)]
fwrite(share_q, file.path(DAT, "dmp_share_by_length_quintile_network.tsv"), sep = "\t")
cmp[, module := factor(module, levels = levels(enr_fig$module))]
cmp[, test := factor(test, levels = c("Fisher, unadjusted", "CMH, gene-length quintiles"))]
p_cmp <- ggplot(cmp, aes(OR, module, colour = test, group = test)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), linewidth = 0.45) +
  geom_point(aes(shape = fdr < 0.05), position = position_dodge(width = 0.6), size = 2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "FDR < 0.05") +
  scale_colour_manual(values = c(`Fisher, unadjusted` = "grey45", `CMH, gene-length quintiles` = "#2C7FB8"), name = NULL) +
  scale_x_log10() +
  labs(x = "Odds ratio, DMP-bearing genes (95% CI)", y = "WGCNA module", title = "Same test, without and with gene length matching") +
  theme_pub() + theme(axis.text.y = element_text(size = 7), legend.position = "bottom", legend.text = element_text(size = 7),
                      plot.title = element_text(size = 9, face = "bold"))
ml_plot <- copy(mod_len)[, module := factor(module, levels = levels(enr_fig$module))][!is.na(module)]
p_len <- ggplot(ml_plot, aes(module, len / 1000, fill = module)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  geom_hline(yintercept = median(mod_len$len) / 1000, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = setNames(levels(ml_plot$module), levels(ml_plot$module)), guide = "none") +  # module names are valid R colours
  scale_y_log10() + coord_flip() +
  labs(x = NULL, y = "Gene length (kb, log scale)", title = "Modules differ in gene length") +
  theme_pub() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(size = 9, face = "bold"))
p_share <- ggplot(share_q, aes(factor(lenq), pct_with_dmp)) +
  geom_col(width = 0.65, fill = "#D55E00") +
  geom_text(aes(label = sprintf("%.1f%%", pct_with_dmp)), vjust = -0.4, size = 2.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Gene-length quintile (1 = shortest)", y = "Genes with a DMP (%)", title = "DMP carriage rises with gene length") +
  theme_pub() + theme(plot.title = element_text(size = 9, face = "bold"))
save_gg(p_cmp + p_len + p_share + plot_layout(widths = c(1.4, 1, 0.9)), FIGS, "figS7_length_adjustment_effect", 10.5, 4.4)

# Step 16 - Module DMR-burden enrichment (Fisher); main panel fig7d
cat("[6b] module DMR-burden Fisher test\n")
dmr_genes <- unique(fread(DMR_TSV)$gene_id); dmr_genes <- dmr_genes[!is.na(dmr_genes) & dmr_genes != ""]
enr_dmr <- rbindlist(lapply(mods, function(m) {
  inmod <- mod_dt$gene_id[mod_dt$module == m]
  a <- sum(inmod %in% dmr_genes); b <- length(inmod) - a
  c_ <- sum(univ %in% dmr_genes) - a; d <- length(univ) - length(inmod) - c_
  # STAT TEST: two-sided Fisher's exact test on the 2x2 (module members vs rest of network)
  ft <- fisher.test(matrix(c(a, b, c_, d), 2))
  data.table(module = m, n_genes = length(inmod), n_dmr = a,
             OR = unname(ft$estimate), p = ft$p.value)
}))
enr_dmr[, fdr := p.adjust(p, "BH")]
fwrite(enr_dmr[order(fdr)], file.path(DAT, "module_dmr_enrichment.tsv"), sep = "\t")
setorder(enr_dmr, OR)
enr_dmr[, module := factor(module, levels = module)]
p_dmr <- ggplot(enr_dmr, aes(module, OR)) +
  geom_segment(aes(xend = module, y = 1, yend = OR), colour = "grey85", linewidth = 0.3) +
  geom_point(aes(size = n_dmr, colour = fdr < 0.05), alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  scale_colour_manual(values = c(`TRUE` = "#1B9E9E", `FALSE` = "grey65"), name = "FDR < 0.05") +
  scale_size_continuous(name = "Genes with ≥1 DMR", range = c(2, 9)) +
  coord_flip() +
  labs(x = "WGCNA module", y = "Odds ratio (module vs rest of network)",
       title = "DMR enrichment per WGCNA module") +
  theme_pub() + theme(axis.text.y = element_text(size = 7))
save_gg(p_dmr, FIGM, "fig7d_module_dmr_fisher", 7.0, 5.0)

# Step 17 - HCP (CpG-island promoter) gene enrichment per module; selection of enriched modules
cat("[6c] HCP-gene module enrichment (Fisher)\n")
WEBER_TSV <- file.path(PIPE, "03_promoters/data/promoter_weber_classification.tsv")
hcp_genes <- unique(fread(WEBER_TSV)[weber_class == "HCP", gene_id])
hcp_genes <- intersect(hcp_genes, univ)  # HCP genes that are in the network
enr_hcp <- rbindlist(lapply(mods, function(m) {
  inmod <- mod_dt$gene_id[mod_dt$module == m]
  a <- sum(inmod %in% hcp_genes); b <- length(inmod) - a
  c_ <- sum(univ %in% hcp_genes) - a; d <- length(univ) - length(inmod) - c_
  # STAT TEST: two-sided Fisher's exact test on the 2x2 (module members vs rest of network)
  ft <- fisher.test(matrix(c(a, b, c_, d), 2))
  data.table(module = m, n_genes = length(inmod), n_hcp = a,
             OR = unname(ft$estimate), p = ft$p.value)
}))
enr_hcp[, fdr := p.adjust(p, "BH")]
fwrite(enr_hcp[order(fdr)], file.path(DAT, "module_hcp_enrichment.tsv"), sep = "\t")
cat(sprintf("  %d HCP genes in the network; modules enriched (FDR<0.05, OR>1): %s\n",
            length(hcp_genes), paste(enr_hcp[fdr < 0.05 & OR > 1]$module, collapse = ", ")))
setorder(enr_hcp, OR)
enr_hcp[, module := factor(module, levels = module)]
p_hcp <- ggplot(enr_hcp, aes(module, OR)) +
  geom_segment(aes(xend = module, y = 1, yend = OR), colour = "grey85", linewidth = 0.3) +
  geom_point(aes(size = n_hcp, colour = fdr < 0.05), alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  scale_colour_manual(values = c(`TRUE` = "#117733", `FALSE` = "grey65"), name = "FDR < 0.05") +
  scale_size_continuous(name = "HCP genes", range = c(2, 9)) +
  coord_flip() +
  labs(x = "WGCNA module", y = "Odds ratio (module vs rest of network)",
       title = "HCP (CpG-island) gene enrichment per WGCNA module") +
  theme_pub() + theme(axis.text.y = element_text(size = 7))
save_gg(p_hcp, FIGS, "fig7_s13_module_hcp_enrichment", 7.0, 5.0)

# Enriched modules for Steps 18-19: FDR < 0.05 and OR > 1 on the unadjusted Fisher test;
# falls back to the single most significant module so the set is never empty
enriched <- as.character(enr$module[which(enr$fdr < 0.05 & enr$OR > 1)])
if (!length(enriched)) enriched <- as.character(enr$module[which.min(enr$fdr)])
cat(sprintf("  enriched modules: %s\n", paste(enriched, collapse = ", ")))

# Step 18 - GO/KEGG enrichment of the enriched modules; main panel fig7b
# No OrgDb exists for D. laeve; terms come from the STRING v12 file as TERM2GENE/TERM2NAME
cat("[7] GO enrichment of enriched modules\n")
go_raw <- fread(STRING, header = TRUE, sep = "\t", quote = "")
setnames(go_raw, 1:4, c("protein", "category", "term", "description"))
go_raw[, gene := sub("^[^.]+\\.", "", protein)]  # STRG0A31YWK.LOC_x -> LOC_x
onto_map <- c("Biological Process (Gene Ontology)" = "BP",
              "Molecular Function (Gene Ontology)" = "MF",
              "Cellular Component (Gene Ontology)" = "CC",
              "KEGG (Kyoto Encyclopedia of Genes and Genomes)" = "KEGG")
go <- unique(go_raw[category %in% names(onto_map),
                    .(gene, term, description, ontology = onto_map[category])])
universe  <- intersect(colnames(datExpr), unique(go$gene))  # network genes with a GO
go_u      <- go[gene %in% universe]
N         <- length(universe)
term_size <- table(go_u$term)
term_desc <- unique(go_u[, .(term, description, ontology)])
cat(sprintf("  %d/%d network genes carry a GO term\n", N, ncol(datExpr)))

# enricher output is reshaped to k (hits), n (module genes with a term), K (term size), N (universe)
TERM2GENE <- go_u[, .(term, gene)]
TERM2NAME <- unique(go_u[, .(term, description)])
enrich_module <- function(mod_genes, min_term = 5, max_term = 2000, min_hits = 3) {
  mg <- intersect(mod_genes, universe)
  if (length(mg) < 10) return(NULL)
  # STAT TEST: hypergeometric over-representation (clusterProfiler enricher), BH-adjusted
  e <- tryCatch(clusterProfiler::enricher(
         gene = mg, universe = universe,
         TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
         minGSSize = min_term, maxGSSize = max_term,
         pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1),
       error = function(err) NULL)
  df <- if (is.null(e)) NULL else as.data.frame(e)
  if (is.null(df) || !nrow(df)) return(NULL)
  r <- as.data.table(df)
  r[, `:=`(k = as.integer(sub("/.*", "", GeneRatio)),  # hits in module
           n = as.integer(sub(".*/", "", GeneRatio)),  # module genes with any term
           K = as.integer(sub("/.*", "", BgRatio)),  # term size in universe
           N = as.integer(sub(".*/", "", BgRatio)))]  # universe size
  r <- r[k >= min_hits]
  if (!nrow(r)) return(NULL)
  r[, fold := (k / n) / (K / N)]
  out <- r[, .(term = ID, k, K, n, N, fold, p = pvalue, padj = p.adjust,
               description = Description)]
  merge(out, unique(term_desc[, .(term, ontology)]), by = "term")[order(p)]
}

go_all <- rbindlist(lapply(enriched, function(m) {
  r <- enrich_module(mod_dt$gene_id[mod_dt$module == m])
  if (is.null(r)) return(NULL)
  cbind(module = m, r)
}), fill = TRUE)

if (nrow(go_all)) {
  fwrite(go_all, file.path(DAT, "module_go_enrichment.tsv"), sep = "\t")
  # Facets: modules (rows) x ontology (columns); x = fold enrichment on a log axis
  go_all[, ratio := k / n]
  go_all[, log_q := pmin(-log10(pmax(padj, 1e-300)), 50)]
  # Only terms with padj < 0.05, at most three per module x ontology; KEGG must be a factor level
  # or its rows drop from the facet; make.unique guards against duplicated labels
  top <- go_all[padj < 0.05][order(padj)][, head(.SD, 3), by = .(module, ontology)]
  top[, ont_lab := factor(ontology, levels = c("BP", "MF", "CC", "KEGG"))]
  setorder(top, ont_lab, module, padj)
  top[, lab_txt := sprintf("%s (%s)", description, ontology)]
  top[, label_key := paste(module, lab_txt, sep = "||")]
  top[, label := factor(label_key, levels = rev(unique(label_key)))]
  p_go <- ggplot(top, aes(fold, label, colour = padj, size = k)) +
    geom_point() +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    facet_grid(module ~ ont_lab, scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#7B241C", high = "#F5CBA7", name = "BH FDR",
                          guide = guide_colourbar(reverse = TRUE)) +
    scale_size_continuous(range = c(2, 6), name = "Module genes") +
    scale_x_log10(breaks = c(3, 30), labels = function(x) paste0(x, "x")) +
    labs(x = "Fold enrichment (observed / expected, log scale)", y = NULL,
         title = "GO and KEGG enrichment per module") +
    theme_pub() +
    theme(axis.text.y = element_text(size = 7), strip.text.x = element_text(size = 8, face = "bold"),
          strip.text.y = element_text(size = 7, face = "bold", angle = 0),
          strip.background = element_blank(),
          plot.margin = margin(5.5, 12, 5.5, 5.5),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))
  save_gg(p_go, FIGM, "fig7b_module_go_enrichment", 7.8, max(3.0, 1.5 * length(enriched)))
} else {
  cat("  no GO terms evaluable for the enriched modules\n")
}

# Step 19 - Module eigengene scores, tail Control vs Amputated; main panel fig7c
# Tail samples are selected by experiment group so no control of another experiment is picked up
cat("[8] eigengene scores (tail only)\n")
tail_s <- rownames(meta)[meta$group %in% c("TailControl", "TailAmputated")]
me_cols <- paste0("ME", enriched)  # eigengenes of the enriched modules
me_cols <- intersect(me_cols, colnames(MEs))
eig <- data.table(sample = tail_s, condition = meta[tail_s, "group"],
                  MEs[tail_s, me_cols, drop = FALSE])
eig_long <- melt(eig, id.vars = c("sample", "condition"),
                 variable.name = "module", value.name = "eigengene")
eig_long[, module := sub("^ME", "", module)]
eig_long[, condition := factor(ifelse(condition == "TailControl", "Control", "Amputated"),
                               levels = c("Control", "Amputated"))]
fwrite(eig_long, file.path(DAT, "module_eigengene_scores_tail.tsv"), sep = "\t")

p_eig <- ggplot(eig_long, aes(condition, eigengene, colour = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey40") +
  geom_jitter(position = position_jitter(width = 0.12, seed = 20260426), size = 2) +  # seeded: identical pdf/png/svg
  facet_wrap(~ module, scales = "free_y") +
  scale_colour_manual(values = COL_COND, guide = "none") +
  labs(x = NULL, y = "Module eigengene (tail samples)",
       title = "Eigengene scores of DMP-enriched modules") +
  theme_pub() + theme(strip.background = element_blank(),
                      plot.title = element_text(size = 9, face = "bold"))
save_gg(p_eig, FIGM, "fig7c_eigengene_scores_tail",
        3.8, 1.36 + 1.05 * ceiling(length(me_cols) / 2))

# Step 20 - Module palette; module-trait heatmap over one-hot experiment groups
mod_pal <- setNames(sort(unique(modColors)), sort(unique(modColors)))

cat("[9] module-trait heatmap (tissue x condition groups)\n")
s   <- rownames(meta)
grp <- meta$group
stopifnot(!anyNA(grp))
ord <- intersect(c("Ovotestis","Juvenile","Head","TailControl","TailAmputated",
                   "EyeControl","EyeAmputated","FungicideControl","FungicideTreated",
                   "IrradiatedControl","IrradiatedTreated"), unique(grp))
traits <- model.matrix(~ 0 + factor(grp, levels = ord)); colnames(traits) <- ord
rownames(traits) <- s; traits <- traits[rownames(MEs), , drop = FALSE]
mtc <- cor(MEs, traits, use = "p")
# STAT TEST: module-trait association = Pearson correlation (cor, above) with WGCNA's
# Student asymptotic p-value (corPvalueStudent, n = #samples); shown on the heatmap.
mtp <- corPvalueStudent(mtc, nrow(MEs))
# BH adjustment across the whole module x trait grid (one family of tests)
mtp <- matrix(p.adjust(mtp, "BH"), nrow(mtp), dimnames = dimnames(mtp))
data.table::fwrite(data.table::data.table(module = rownames(mtc), as.data.frame(mtc)),
                   file.path(DAT, "module_trait_cor.tsv"), sep = "\t")
data.table::fwrite(data.table::data.table(module = rownames(mtp), as.data.frame(mtp)),
                   file.path(DAT, "module_trait_fdr.tsv"), sep = "\t")
txt <- paste0(signif(mtc, 2), "\n(", signif(mtp, 1), ")"); dim(txt) <- dim(mtc)
save_base(FIGS, "fig7_s5_module_trait_heatmap", 8, 7.5, function() {
  par(mar = c(9, 8, 3, 1))
  labeledHeatmap(Matrix = mtc, xLabels = colnames(mtc), yLabels = rownames(mtc),
    ySymbols = rownames(mtc), colorLabels = FALSE, colors = blueWhiteRed(50),
    textMatrix = txt, setStdMargins = FALSE, cex.text = 0.45, cex.lab = 0.7,
    zlim = c(-1, 1), main = "Module-trait relationships (tissue + condition)")
})
fwrite(data.table(module = rownames(mtc), mtc), file.path(DAT, "module_trait_cor.tsv"), sep = "\t")

# Step 21 - Module sizes
cat("[10] module sizes\n")
sz <- as.data.table(sort(table(modColors), decreasing = TRUE)); setnames(sz, c("module", "n_genes"))
fwrite(sz, file.path(DAT, "module_sizes.tsv"), sep = "\t")
sz[, module := factor(module, levels = module)]
p_sz <- ggplot(sz, aes(module, n_genes, fill = module)) + geom_col() +
  geom_text(aes(label = n_genes), vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = mod_pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of genes", title = "Module sizes (all genes per module)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_gg(p_sz, FIGS, "fig7_s6_module_sizes", 8, 3.5)

# Step 22 - Module membership (kME), intramodular connectivity and hub genes
cat("[11] module membership (kME) + hub genes\n")
nGenes <- ncol(datExpr)
kME <- as.matrix(signedKME(datExpr, MEs, outputColumnName = "kME"))
mm_idx <- match(paste0("kME", modColors), colnames(kME))
MM <- rep(NA_real_, nGenes); ok <- !is.na(mm_idx)
MM[ok] <- kME[cbind(which(ok), mm_idx[ok])]  # kME of a gene to its OWN module
kWithin <- setNames(rep(NA_real_, nGenes), colnames(datExpr)); kWs <- kWithin
for (m in mods) {  # intramodular connectivity per module
  idx <- which(modColors == m); g <- colnames(datExpr)[idx]
  if (length(g) < 3) next
  adj <- adjacency(datExpr[, g, drop = FALSE], power = soft_power, type = "signed")
  kin <- rowSums(adj) - 1; rm(adj); gc(verbose = FALSE)
  kWithin[idx] <- kin; kWs[idx] <- kin / max(kin)
}
# Hub = top 10 percent of scaled intramodular connectivity within its module and |kME| >= 0.80
kin_thr <- tapply(kWs, modColors, quantile, probs = 0.90, na.rm = TRUE)
is_hub <- (modColors != "grey") & !is.na(MM) & (kWs >= kin_thr[modColors]) & (abs(MM) >= 0.80)
is_hub[is.na(is_hub)] <- FALSE
mm_tbl <- data.table(gene_id = colnames(datExpr), module = modColors,
                     ModuleMembership = MM, kWithin = kWithin, kWithin_scaled = kWs, is_hub = is_hub)
fwrite(mm_tbl, file.path(DAT, "module_membership_hubs.tsv"), sep = "\t")
cat(sprintf("  %d hub genes (top 10%% kWithin & |kME|>0.8)\n", sum(is_hub)))

hc <- as.data.table(table(module = factor(mm_tbl$module[mm_tbl$is_hub], levels = mods)))
setnames(hc, c("module", "n_hub")); hc[, module := factor(module, levels = mods)]
p_hub <- ggplot(hc, aes(module, n_hub, fill = module)) + geom_col() +
  geom_text(aes(label = n_hub), vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = mod_pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Hub genes", title = "Hub genes per module",
       subtitle = "hub = top 10% intramodular connectivity & |kME| > 0.8") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_gg(p_hub, FIGS, "fig7_s7_hub_counts", 8, 3.5)

mmp <- mm_tbl[module != "grey" & is.finite(kWithin)]
p_mm <- ggplot(mmp, aes(ModuleMembership, kWithin, colour = is_hub)) +
  geom_point(size = 0.4, alpha = 0.5) + facet_wrap(~ module, scales = "free_y") +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70"), name = "Hub") +
  labs(x = "Module membership (kME)", y = "Intramodular connectivity (kWithin)",
       title = "Module membership vs connectivity (hubs in orange)") +
  theme_pub() + theme(axis.text = element_text(size = 5), strip.background = element_blank())
save_gg(p_mm, FIGS, "fig7_s8_module_membership", 11, 8)

# Step 23 - Do hub genes carry the DMPs / DMRs? Fisher plus length-tertile CMH
cat("[11b] hub genes vs DMP / DMR carriage\n")
gene_gr8 <- gff[gff$type == "gene"]
gid8  <- sub(";.*", "", as.character(mcols(gene_gr8)$ID))
glen8 <- setNames(width(gene_gr8), gid8)
note8 <- vapply(mcols(gene_gr8)$Note, function(x) if (length(x)) as.character(x)[1] else NA_character_, character(1))
sym8  <- setNames(sub("^Similar to ([^:]+):.*$", "\\1", note8), gid8)
sym8[!grepl("^Similar to [^:]+:", note8)] <- NA
disp8 <- function(g) ifelse(is.na(sym8[g]), g, sym8[g])  # symbol, else the LOC id

hubs <- copy(mm_tbl)
hubs[, `:=`(has_dmp = gene_id %in% dmp_genes, has_dmr = gene_id %in% dmr_genes,
            gene_len = glen8[gene_id], symbol = disp8(gene_id))]
fwrite(hubs[is_hub == TRUE][order(module, -abs(ModuleMembership))][
         , .(gene_id, symbol, module, ModuleMembership, kWithin_scaled, gene_len, has_dmp, has_dmr)],
       file.path(DAT, "hub_genes_dmp_dmr.tsv"), sep = "\t")

hub_test <- function(col) {  # hub vs non-hub, network-wide
  x <- hubs[module != "grey" & !is.na(gene_len)]
  # STAT TEST: two-sided Fisher's exact (hub vs feature); plus a Cochran-Mantel-Haenszel
  # test (mantelhaen.test) stratified by gene-length tertile to control the length confound.
  ft <- fisher.test(table(factor(x$is_hub, c(TRUE, FALSE)), factor(x[[col]], c(TRUE, FALSE))))
  x[, ltert := cut(gene_len, quantile(gene_len, 0:3/3, na.rm = TRUE), include.lowest = TRUE)]
  cmh <- tryCatch(mantelhaen.test(table(factor(x$is_hub, c(TRUE, FALSE)),
                                        factor(x[[col]], c(TRUE, FALSE)), x$ltert)),
                  error = function(e) NULL)
  data.table(mark = toupper(sub("has_", "", col)),
             n_hub = sum(x$is_hub), n_hub_mark = sum(x$is_hub & x[[col]]),
             pct_hub = 100 * mean(x[[col]][x$is_hub]), pct_nonhub = 100 * mean(x[[col]][!x$is_hub]),
             OR = unname(ft$estimate), p = ft$p.value,
             OR_len_adj = if (is.null(cmh)) NA_real_ else unname(cmh$estimate),
             p_len_adj  = if (is.null(cmh)) NA_real_ else cmh$p.value,
             med_len_hub = median(x$gene_len[x$is_hub]), med_len_nonhub = median(x$gene_len[!x$is_hub]))
}
hub_overall <- rbindlist(lapply(c("has_dmp", "has_dmr"), hub_test))
fwrite(hub_overall, file.path(DAT, "hub_dmp_dmr_enrichment.tsv"), sep = "\t")
print(hub_overall)

# Per module: fraction of hubs carrying a DMP / a DMR
hub_mod <- hubs[is_hub == TRUE, .(n_hub = .N, n_dmp = sum(has_dmp), n_dmr = sum(has_dmr),
                                  pct_dmp = 100 * mean(has_dmp), pct_dmr = 100 * mean(has_dmr)),
                by = module][order(-n_hub)]
fwrite(hub_mod, file.path(DAT, "hub_dmp_dmr_by_module.tsv"), sep = "\t")
cat("  hub genes carrying a DMR (all modules):\n")
print(hubs[is_hub == TRUE & has_dmr == TRUE][order(module)][
        , .(symbol, gene_id, module, kME = round(ModuleMembership, 2))])

hm_long <- melt(hub_mod[, .(module, DMP = pct_dmp, DMR = pct_dmr)], id.vars = "module",
                variable.name = "mark", value.name = "pct")
p_hubmeth <- ggplot(hm_long, aes(reorder(module, -pct), pct, fill = mark)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c(DMP = "#2C7FB8", DMR = "#1B9E9E"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "WGCNA module", y = "% of hub genes carrying the mark",
       title = "Methylation changes in hub genes",
       subtitle = sprintf("DMP: hubs %.1f%% vs other module genes %.1f%% (OR %.2f, p = %.2g; length-adjusted OR %.2f)",
                          hub_overall$pct_hub[1], hub_overall$pct_nonhub[1], hub_overall$OR[1],
                          hub_overall$p[1], hub_overall$OR_len_adj[1])) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_gg(p_hubmeth, FIGS, "fig7_s9_hub_dmp_dmr", 8, 4)

# Step 24 - Hub DMRs: methylation change vs expression change
# dMeth = control - amputated (positive = loss of methylation on amputation); two DE definitions
# are reported: padj < 0.05 with |log2FC| >= 1, and padj < 0.05 alone
cat("[11c] hub DMRs: dMeth vs log2FC\n")
DE_TAIL8 <- file.path(PIPE, "01_genome_toolkit/data/gene_de_tail.tsv")
dmr_full <- fread(DMR_TSV)[, .(gene_id, region, nCG, dmr_len = length,
                               dmr_chr = chr, dmr_start = start, dmr_end = end,
                               meth_ctrl = meanMethy1, meth_amp = meanMethy2,
                               dMeth = diff.Methy, direction)]
de8 <- fread(DE_TAIL8)[, .(gene_id, baseMean, log2FoldChange, padj, symbol)]
hub_dmr <- merge(dmr_full[gene_id %in% hubs[is_hub == TRUE]$gene_id],
                 hubs[is_hub == TRUE, .(gene_id, module, kME = ModuleMembership, gene_len)],
                 by = "gene_id")
hub_dmr <- merge(hub_dmr, de8, by = "gene_id", all.x = TRUE)
hub_dmr[, `:=`(dmp_enriched_module = module %in% enriched,
               de_paper = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1,
               de_padj_only = !is.na(padj) & padj < 0.05)]
setorder(hub_dmr, module, -dMeth)
fwrite(hub_dmr, file.path(DAT, "hub_dmr_meth_vs_expression.tsv"), sep = "\t")
cat(sprintf("  %d DMRs in %d hub genes; DE at the paper rule: %d; at padj<0.05 only: %d\n",
            nrow(hub_dmr), uniqueN(hub_dmr$gene_id), sum(hub_dmr$de_paper), sum(hub_dmr$de_padj_only)))
print(hub_dmr[, .(symbol, module, region, nCG, meth_ctrl = round(meth_ctrl, 3),
                  meth_amp = round(meth_amp, 3), dMeth = round(dMeth, 3),
                  LFC = round(log2FoldChange, 2), padj = signif(padj, 2))])

# Gene-body DMRs of all genes as background, hubs highlighted
gb <- merge(dmr_full[region %in% c("Exon", "Intron")], de8, by = "gene_id")
gb <- merge(gb, hubs[, .(gene_id, is_hub, module)], by = "gene_id", all.x = TRUE)
gb <- gb[is.finite(log2FoldChange) & is.finite(dMeth)]
gb[, is_hub := !is.na(is_hub) & is_hub]
rho <- suppressWarnings(cor(gb$dMeth, gb$log2FoldChange, method = "spearman", use = "complete.obs"))
rho_h <- suppressWarnings(cor(gb[is_hub == TRUE]$dMeth, gb[is_hub == TRUE]$log2FoldChange,
                              method = "spearman", use = "complete.obs"))
lab_gb <- gb[is_hub == TRUE][order(-abs(log2FoldChange))][seq_len(min(12, .N))]
p_gb <- ggplot(gb, aes(dMeth, log2FoldChange)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(data = gb[is_hub == FALSE], colour = "grey78", size = 0.7, alpha = 0.6) +
  geom_point(data = gb[is_hub == TRUE], colour = "#D55E00", size = 1.6) +
  ggrepel::geom_text_repel(data = lab_gb, aes(label = ifelse(is.na(symbol) | symbol == "", gene_id, symbol)),
                           size = 2.3, max.overlaps = 20, segment.colour = "grey70") +
  labs(x = expression(Delta*"Meth (control - amputated)"), y = expression(log[2]~"fold change"),
       title = "Gene-body DMRs: methylation change vs expression change",
       subtitle = sprintf("hub genes in orange (n = %d DMRs); Spearman rho = %.3f all, %.3f hubs",
                          nrow(gb[is_hub == TRUE]), rho, rho_h)) +
  theme_pub()
save_gg(p_gb, FIGS, "fig7_s10_hub_dmr_meth_vs_expression", 7, 5)

# Step 25 - Transcription-factor content of the hubs
# TF = JASPAR2024 CORE animal TF name matched to the locus symbol (EviAnn 'Similar to', else
# eggNOG Preferred_name); a name-level match restricted to TFs that have a JASPAR motif
cat("[11e] TF content of the WGCNA hubs\n")
suppressPackageStartupMessages(library(RSQLite))
con <- dbConnect(SQLite(), JASPAR_SQLITE)
jn  <- dbGetQuery(con, paste0("SELECT DISTINCT m.NAME FROM MATRIX m ",
        "JOIN MATRIX_ANNOTATION a ON m.ID=a.ID WHERE m.COLLECTION='CORE' AND ",
        "a.TAG='tax_group' AND a.VAL IN ('vertebrates','insects','nematodes','urochordates')"))
dbDisconnect(con)
tf_syms <- unique(toupper(unlist(strsplit(jn$NAME, "::|/"))))  # split heterodimer names
cat(sprintf("  %d JASPAR2024 CORE animal TF symbols\n", length(tf_syms)))
# Symbol per locus: EviAnn 'Similar to SYM:' first, else eggNOG Preferred_name
note8b <- vapply(mcols(gene_gr8)$Note, function(x) if (length(x)) as.character(x)[1] else NA_character_, character(1))
evi8   <- toupper(sub("^Similar to ([^:]+):.*$", "\\1", note8b)); evi8[!grepl("^Similar to [^:]+:", note8b)] <- NA
egg8   <- fread(EMAPPER, sep = "\t", quote = "", header = TRUE, skip = "#query", na.strings = c("-","","NA"), fill = TRUE)
setnames(egg8, 1, "query"); egg8 <- egg8[!startsWith(query, "##")]; egg8[, locus := sub("-mRNA-.*$", "", query)]
egg8m  <- unique(egg8[locus %in% gid8 & !is.na(Preferred_name), .(gene_id = locus, egg = toupper(Preferred_name))], by = "gene_id")
tfmap  <- merge(data.table(gene_id = gid8, evi = evi8), egg8m, by = "gene_id", all.x = TRUE)
tfmap[, sym_use := fifelse(!is.na(evi), evi, egg)]
tfmap[, is_tf := !is.na(sym_use) & sym_use %in% tf_syms]
cat(sprintf("  %d D. laeve loci flagged as TFs\n", sum(tfmap$is_tf)))

# net is reassigned here; the blockwiseModules object is no longer needed
net <- merge(mm_tbl, tfmap[, .(gene_id, sym_use, is_tf)], by = "gene_id", all.x = TRUE)
net[, is_tf := !is.na(is_tf) & is_tf]  # genes with no ortholog symbol -> not a TF
net <- net[module != "grey"]
# STAT TEST: two-sided Fisher's exact test — are hub genes enriched for transcription factors?
ft_tf <- fisher.test(table(factor(net$is_hub, c(TRUE, FALSE)), factor(net$is_tf, c(TRUE, FALSE))))
cat(sprintf("  network: %d genes (%d TFs, %.1f%%);  hubs: %d (%d TFs, %.1f%%);  Fisher OR %.2f p %.3g\n",
            nrow(net), sum(net$is_tf), 100*mean(net$is_tf), sum(net$is_hub),
            sum(net$is_hub & net$is_tf), 100*mean(net$is_tf[net$is_hub]),
            unname(ft_tf$estimate), ft_tf$p.value))
tf_by_mod <- net[, .(n_genes = .N, n_tf = sum(is_tf), n_hub = sum(is_hub),
                     n_hub_tf = sum(is_hub & is_tf)), by = module][order(-n_hub_tf)]
fwrite(tf_by_mod, file.path(DAT, "module_tf_content.tsv"), sep = "\t")
fwrite(net[is_hub == TRUE & is_tf == TRUE][order(module, -ModuleMembership),
           .(gene_id, symbol = sym_use, module, ModuleMembership, kWithin_scaled)],
       file.path(DAT, "hub_transcription_factors.tsv"), sep = "\t")
cat("  hub TFs:\n"); print(net[is_hub == TRUE & is_tf == TRUE][order(module), .(sym_use, gene_id, module)])
p_tf <- {
  mm <- rbind(`all module genes` = 100 * tf_by_mod$n_tf / tf_by_mod$n_genes,
              hubs = 100 * tf_by_mod$n_hub_tf / pmax(tf_by_mod$n_hub, 1))
  colnames(mm) <- tf_by_mod$module
  function() { par(mar = c(6, 4.5, 3, 1))
    barplot(mm, beside = TRUE, col = c("grey70", "#D55E00"), border = NA, las = 2, cex.names = 0.7,
            ylab = "% transcription factors", ylim = c(0, max(mm, na.rm = TRUE) * 1.15),
            main = "TF content of modules and their hubs")
    legend("topright", legend = rownames(mm), fill = c("grey70", "#D55E00"), bty = "n", cex = 0.8) }
}
save_base(FIGS, "fig7_s12_hub_tf_content", 8, 4.5, p_tf)

# Step 26 - GO enrichment for all modules (one dot-plot page per module)
cat("[12] GO for all modules\n")
go_all_mods <- rbindlist(lapply(mods, function(m) {
  r <- enrich_module(mod_dt$gene_id[mod_dt$module == m]); if (is.null(r)) return(NULL); cbind(module = m, r)
}), fill = TRUE)
if (nrow(go_all_mods)) {
  fwrite(go_all_mods, file.path(DAT, "module_go_all.tsv"), sep = "\t")
  pdf(file.path(FIGS, "fig7_s9_go_all_modules.pdf"), width = 8, height = 6)
  for (m in mods) {
    r <- go_all_mods[module == m]
    if (!nrow(r)) { plot.new(); title(main = paste0("module ", m, ": no GO terms")); next }
    top <- head(r[order(p)], 15)
    top[, desc_w := vapply(description, function(s) paste(strwrap(s, 34), collapse = "\n"), character(1))]
    top[, label := factor(make.unique(desc_w), levels = rev(make.unique(desc_w)))]
    print(ggplot(top, aes(fold, label, colour = padj, size = k)) +
      geom_point() + geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
      scale_colour_gradient(low = "#7B241C", high = "#F5CBA7", name = "BH FDR",
                            guide = guide_colourbar(reverse = TRUE)) +
      scale_size_continuous(range = c(2, 7), name = "Genes") +
      labs(x = "Fold enrichment", y = NULL,
           title = sprintf("GO enrichment — module %s", m)) + theme_pub())
  }
  dev.off()
  cat("  saved fig8_s9_go_all_modules.pdf (one dot-plot page per module)\n")
} else cat("  no GO terms evaluable across modules\n")

# Step 27 - Session info
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_07_wgcna.txt"))
cat("[07_wgcna] done\n")
