#!/usr/bin/env Rscript

set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(WGCNA)
  library(GenomicRanges); library(ggplot2)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

PIPE   <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
HTSEQ  <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"
STRING <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/STRING.protein.enrichment.terms.v12.0.txt"
GFF_RDS <- file.path(PIPE, "batch01/objects/gff_chrmt.rds")
DMP_TSV <- file.path(PIPE, "batch05/data/dmps_annotated.tsv")
DMR_TSV <- file.path(PIPE, "batch05/data/dmrs_annotated.tsv")
JASPAR_SQLITE <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/jaspar/JASPAR2024.sqlite"
EMAPPER <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/eggnog_mapper/dlasi_proteome.emapper.annotations"
BATCH <- file.path(PIPE, "batch07")
OBJ  <- file.path(BATCH, "objects")
DAT  <- file.path(BATCH, "data")
FIGM <- file.path(BATCH, "figures/main")
FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(OBJ, DAT, FIGM, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(c(FIGM, FIGS), full.names = TRUE))

COL_COND <- c(Control = "#0072B2", Amputated = "#D55E00")
theme_pub <- function() theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"))
save_gg <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(dir, paste0(name, ".svg")), p, width = w, height = h)
  cat(sprintf("  saved %s\n", name))
}
save_base <- function(dir, name, w, h, FUN) {
  pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h); FUN(); dev.off()
  png(file.path(dir, paste0(name, ".png")), width = w, height = h,
      units = "in", res = 150); FUN(); dev.off()
  svglite::svglite(file.path(dir, paste0(name, ".svg")), width = w, height = h); FUN(); dev.off()
  cat(sprintf("  saved %s\n", name))
}

cat("[1] sample map\n")
files <- list.files(HTSEQ, pattern = "_htseq_gene_counts.txt$")
base  <- sub("\\.Aligned\\.out\\.bam_htseq_gene_counts\\.txt$", "", files)
base  <- sub("_htseq_gene_counts\\.txt$", "", base)
classify <- function(b) {
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
meta <- meta[tissue != "other"]
setorder(meta, group, orig_sample)
meta[, sample := paste0(group, seq_len(.N)), by = group]
fwrite(meta, file.path(DAT, "sample_map.tsv"), sep = "\t")
cat(sprintf("  %d libraries across %d tissues, %d experiment groups\n",
            nrow(meta), uniqueN(meta$tissue), uniqueN(meta$group)))
print(meta[, .N, by = group][order(group)])
meta <- as.data.frame(meta); rownames(meta) <- meta$sample

cat("[2] counts + VST\n")
cl  <- lapply(meta$file, function(f) {
  x <- fread(file.path(HTSEQ, f), header = FALSE, col.names = c("gene_id", "count"))
  x[!startsWith(gene_id, "__")]
})
ids <- cl[[1]]$gene_id
cm  <- sapply(cl, function(x) x$count[match(ids, x$gene_id)])
rownames(cm) <- ids; colnames(cm) <- meta$sample

gff <- readRDS(GFF_RDS)
chr_genes <- unique(gff$ID[gff$type == "gene"])
cm <- cm[rownames(cm) %in% chr_genes, , drop = FALSE]
cat(sprintf("  %d genes on chr1-31\n", nrow(cm)))

cm  <- cm[rowSums(cm >= 10) >= 3, , drop = FALSE]
cat(sprintf("  %d genes pass the low-count filter\n", nrow(cm)))
vsd <- vst(DESeqDataSetFromMatrix(cm, meta[colnames(cm), ], ~ 1), blind = TRUE)
datExpr <- t(assay(vsd))

cat("[3] sample clustering / outliers\n")
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

trait_colors <- function(m) data.frame(
  tissue = labels2colors(as.numeric(factor(m$tissue))),
  group  = labels2colors(as.numeric(factor(m$group))))

tree1 <- hclust(dist(datExpr), method = "average")
save_base(FIGS, "fig7_s1_sample_clustering", 12, 6, function()
  plotDendroAndColors(tree1, trait_colors(meta[rownames(datExpr), ]),
                      groupLabels = c("tissue", "group"),
                      main = "Sample clustering (outliers marked)", cex.dendroLabels = 0.7))

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
sel     <- good[which.max(q)]
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
keep_var <- gene_var > var_cut
cat(sprintf("  %d -> %d genes after variance filter\n", ncol(datExpr), sum(keep_var)))
datExpr <- datExpr[, keep_var, drop = FALSE]

cat("[4] soft threshold\n")
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
soft_power <- sft$powerEstimate
if (is.na(soft_power)) soft_power <- 14
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

cat("[5] modules (reuse cached wgcna.rds if present)\n")
wgcna_rds <- file.path(OBJ, "wgcna.rds")
if (file.exists(wgcna_rds)) {
  W <- readRDS(wgcna_rds); net <- W$net; modColors <- W$modColors; MEs <- W$MEs
  soft_power <- W$soft_power
  stopifnot(identical(W$genes, colnames(datExpr)))
  cat("  reused cached network — blockwiseModules skipped (fast)\n")
} else {
  cor <- WGCNA::cor
  net <- blockwiseModules(datExpr,
                          power = soft_power, networkType = "signed", TOMType = "signed",
                          minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          maxBlockSize = 25000, saveTOMs = FALSE, verbose = 0)
  cor <- stats::cor
  stopifnot(length(net$dendrograms) == 1)
  modColors <- labels2colors(net$colors)
  MEs <- orderMEs(moduleEigengenes(datExpr, modColors)$eigengenes)
  MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
  saveRDS(list(net = net, modColors = modColors, MEs = MEs,
               genes = colnames(datExpr), meta = meta, soft_power = soft_power), wgcna_rds)
}
mod_dt <- data.table(gene_id = colnames(datExpr), module = modColors)
fwrite(mod_dt, file.path(DAT, "module_assignments.tsv"), sep = "\t")
mods <- setdiff(sort(unique(modColors)), "grey")
cat(sprintf("  %d modules (+ grey); %d genes unassigned\n",
            length(mods), sum(modColors == "grey")))

save_base(FIGS, "fig7_s4_module_dendrogram", 10, 6, function()
  plotDendroAndColors(net$dendrograms[[1]], modColors[net$blockGenes[[1]]],
                      "Module", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                      guideHang = 0.05, main = "Gene dendrogram and module colours"))

cat("[6] module DMP-burden Fisher test\n")
dmp_genes <- unique(fread(DMP_TSV)$gene_id)
dmp_genes <- dmp_genes[!is.na(dmp_genes) & dmp_genes != ""]
univ <- mod_dt$gene_id
enr <- rbindlist(lapply(mods, function(m) {
  inmod <- mod_dt$gene_id[mod_dt$module == m]
  a <- sum(inmod %in% dmp_genes)
  b <- length(inmod) - a
  c_ <- sum(univ %in% dmp_genes) - a
  d <- length(univ) - length(inmod) - c_
  ft <- fisher.test(matrix(c(a, b, c_, d), 2))
  data.table(module = m, n_genes = length(inmod), n_dmp = a,
             OR = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2],
             p = ft$p.value)
}))
enr[, fdr := p.adjust(p, "BH")]
enr <- enr[order(-OR)]
fwrite(enr, file.path(DAT, "module_dmp_enrichment.tsv"), sep = "\t")

cat("[6b] module DMP enrichment, stratified by gene-length quintile\n")
gn8b <- gff[gff$type == "gene"]
glen <- data.table(gene_id = sub(";.*", "", as.character(gn8b$ID)), len = width(gn8b))
mod_len <- merge(mod_dt[, .(gene_id, module)], glen, by = "gene_id")
mod_len[, lenq := cut(len, quantile(len, seq(0, 1, 0.2)), include.lowest = TRUE, labels = FALSE)]
mod_len[, is_dmp := gene_id %in% dmp_genes]
enr_adj <- rbindlist(lapply(mods, function(m) {
  x <- copy(mod_len)[, inmod := module == m]
  tb <- table(factor(x$inmod, c(TRUE, FALSE)), factor(x$is_dmp, c(TRUE, FALSE)), x$lenq)
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

enr_fig <- enr2[!is.na(OR_adj) & !is.na(fdr_adj)]
setorder(enr_fig, OR_adj)
enr_fig[, module := factor(module, levels = module)]
p_fisher <- ggplot(enr_fig, aes(module, OR_adj)) +
  geom_segment(aes(xend = module, y = 1, yend = OR_adj), colour = "grey85", linewidth = 0.3) +
  geom_point(aes(size = n_dmp, colour = fdr_adj < 0.05), alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  scale_colour_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey65"), name = "FDR < 0.05") +
  scale_size_continuous(name = "Genes with ≥1 DMP", range = c(2, 9)) +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5)) +
  coord_flip() +
  labs(x = "WGCNA module",
       y = "CMH common odds ratio, stratified by gene length\n(module vs rest of network)",
       title = "DMP enrichment per WGCNA module,\nmatched for gene length") +
  theme_pub() + theme(axis.text.y = element_text(size = 7),
                      axis.title.x = element_text(size = 8),
                      plot.title = element_text(size = 9, face = "bold"))
save_gg(p_fisher, FIGM, "fig7a_module_dmp_fisher", 4.2, 2.8)

cat("[6b] module DMR-burden Fisher test\n")
dmr_genes <- unique(fread(DMR_TSV)$gene_id); dmr_genes <- dmr_genes[!is.na(dmr_genes) & dmr_genes != ""]
enr_dmr <- rbindlist(lapply(mods, function(m) {
  inmod <- mod_dt$gene_id[mod_dt$module == m]
  a <- sum(inmod %in% dmr_genes); b <- length(inmod) - a
  c_ <- sum(univ %in% dmr_genes) - a; d <- length(univ) - length(inmod) - c_
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

cat("[6c] HCP-gene module enrichment (Fisher)\n")
WEBER_TSV <- file.path(PIPE, "batch03/data/promoter_weber_classification.tsv")
hcp_genes <- unique(fread(WEBER_TSV)[weber_class == "HCP", gene_id])
hcp_genes <- intersect(hcp_genes, univ)
enr_hcp <- rbindlist(lapply(mods, function(m) {
  inmod <- mod_dt$gene_id[mod_dt$module == m]
  a <- sum(inmod %in% hcp_genes); b <- length(inmod) - a
  c_ <- sum(univ %in% hcp_genes) - a; d <- length(univ) - length(inmod) - c_
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

enriched <- as.character(enr2$module[which(enr2$fdr_adj < 0.05 & enr2$OR_adj > 1)])
if (!length(enriched)) enriched <- as.character(enr2$module[which.min(enr2$fdr_adj)])
cat(sprintf("  enriched modules: %s\n", paste(enriched, collapse = ", ")))

cat("[7] GO enrichment of enriched modules\n")
go_raw <- fread(STRING, header = TRUE, sep = "\t", quote = "")
setnames(go_raw, 1:4, c("protein", "category", "term", "description"))
go_raw[, gene := sub("^[^.]+\\.", "", protein)]
onto_map <- c("Biological Process (Gene Ontology)" = "BP",
              "Molecular Function (Gene Ontology)" = "MF",
              "Cellular Component (Gene Ontology)" = "CC",
              "KEGG (Kyoto Encyclopedia of Genes and Genomes)" = "KEGG")
go <- unique(go_raw[category %in% names(onto_map),
                    .(gene, term, description, ontology = onto_map[category])])
universe  <- intersect(colnames(datExpr), unique(go$gene))
go_u      <- go[gene %in% universe]
N         <- length(universe)
term_size <- table(go_u$term)
term_desc <- unique(go_u[, .(term, description, ontology)])
cat(sprintf("  %d/%d network genes carry a GO term\n", N, ncol(datExpr)))

TERM2GENE <- go_u[, .(term, gene)]
TERM2NAME <- unique(go_u[, .(term, description)])
enrich_module <- function(mod_genes, min_term = 5, max_term = 2000, min_hits = 3) {
  mg <- intersect(mod_genes, universe)
  if (length(mg) < 10) return(NULL)
  e <- tryCatch(clusterProfiler::enricher(
    gene = mg, universe = universe,
    TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
    minGSSize = min_term, maxGSSize = max_term,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1),
    error = function(err) NULL)
  df <- if (is.null(e)) NULL else as.data.frame(e)
  if (is.null(df) || !nrow(df)) return(NULL)
  r <- as.data.table(df)
  r[, `:=`(k = as.integer(sub("/.*", "", GeneRatio)),
           n = as.integer(sub(".*/", "", GeneRatio)),
           K = as.integer(sub("/.*", "", BgRatio)),
           N = as.integer(sub(".*/", "", BgRatio)))]
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
  go_all[, ratio := k / n]
  go_all[, log_q := pmin(-log10(pmax(padj, 1e-300)), 50)]
  top <- go_all[padj < 0.05][order(padj)][, head(.SD, 3), by = .(module, ontology)]
  top[, ont_lab := factor(ontology, levels = c("BP", "MF", "CC", "KEGG"))]
  setorder(top, ont_lab, module, padj)
  top[, lab_txt := sprintf("%s (%s)", description, ontology)]
  top[, label := factor(make.unique(lab_txt), levels = rev(make.unique(lab_txt)))]
  p_go <- ggplot(top, aes(fold, label, colour = padj, size = k)) +
    geom_point() +
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

cat("[8] eigengene scores (tail only)\n")
tail_s <- rownames(meta)[meta$group %in% c("TailControl", "TailAmputated")]
me_cols <- paste0("ME", enriched)
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
  geom_jitter(position = position_jitter(width = 0.12, seed = 20260426), size = 2) +
  facet_wrap(~ module, scales = "free_y") +
  scale_colour_manual(values = COL_COND, guide = "none") +
  labs(x = NULL, y = "Module eigengene (tail samples)",
       title = "Eigengene scores of DMP-enriched modules") +
  theme_pub() + theme(strip.background = element_blank(),
                      plot.title = element_text(size = 9, face = "bold"))
save_gg(p_eig, FIGM, "fig7c_eigengene_scores_tail",
        3.8, 1.36 + 1.05 * ceiling(length(me_cols) / 2))

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
mtp <- corPvalueStudent(mtc, nrow(MEs))
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

cat("[11] module membership (kME) + hub genes\n")
nGenes <- ncol(datExpr)
kME <- as.matrix(signedKME(datExpr, MEs, outputColumnName = "kME"))
mm_idx <- match(paste0("kME", modColors), colnames(kME))
MM <- rep(NA_real_, nGenes); ok <- !is.na(mm_idx)
MM[ok] <- kME[cbind(which(ok), mm_idx[ok])]
kWithin <- setNames(rep(NA_real_, nGenes), colnames(datExpr)); kWs <- kWithin
for (m in mods) {
  idx <- which(modColors == m); g <- colnames(datExpr)[idx]
  if (length(g) < 3) next
  adj <- adjacency(datExpr[, g, drop = FALSE], power = soft_power, type = "signed")
  kin <- rowSums(adj) - 1; rm(adj); gc(verbose = FALSE)
  kWithin[idx] <- kin; kWs[idx] <- kin / max(kin)
}
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

cat("[11b] hub genes vs DMP / DMR carriage\n")
gene_gr8 <- gff[gff$type == "gene"]
gid8  <- sub(";.*", "", as.character(mcols(gene_gr8)$ID))
glen8 <- setNames(width(gene_gr8), gid8)
note8 <- vapply(mcols(gene_gr8)$Note, function(x) if (length(x)) as.character(x)[1] else NA_character_, character(1))
sym8  <- setNames(sub("^Similar to ([^:]+):.*$", "\\1", note8), gid8)
sym8[!grepl("^Similar to [^:]+:", note8)] <- NA
disp8 <- function(g) ifelse(is.na(sym8[g]), g, sym8[g])

hubs <- copy(mm_tbl)
hubs[, `:=`(has_dmp = gene_id %in% dmp_genes, has_dmr = gene_id %in% dmr_genes,
            gene_len = glen8[gene_id], symbol = disp8(gene_id))]
fwrite(hubs[is_hub == TRUE][order(module, -abs(ModuleMembership))][
  , .(gene_id, symbol, module, ModuleMembership, kWithin_scaled, gene_len, has_dmp, has_dmr)],
  file.path(DAT, "hub_genes_dmp_dmr.tsv"), sep = "\t")

hub_test <- function(col) {
  x <- hubs[module != "grey" & !is.na(gene_len)]
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

cat("[11c] hub DMRs: dMeth vs log2FC\n")
DE_TAIL8 <- file.path(PIPE, "batch01/data/gene_de_tail.tsv")
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


cat("[11e] TF content of the WGCNA hubs\n")
suppressPackageStartupMessages(library(RSQLite))
con <- dbConnect(SQLite(), JASPAR_SQLITE)
jn  <- dbGetQuery(con, paste0("SELECT DISTINCT m.NAME FROM MATRIX m ",
                              "JOIN MATRIX_ANNOTATION a ON m.ID=a.ID WHERE m.COLLECTION='CORE' AND ",
                              "a.TAG='tax_group' AND a.VAL IN ('vertebrates','insects','nematodes','urochordates')"))
dbDisconnect(con)
tf_syms <- unique(toupper(unlist(strsplit(jn$NAME, "::|/"))))
cat(sprintf("  %d JASPAR2024 CORE animal TF symbols\n", length(tf_syms)))
note8b <- vapply(mcols(gene_gr8)$Note, function(x) if (length(x)) as.character(x)[1] else NA_character_, character(1))
evi8   <- toupper(sub("^Similar to ([^:]+):.*$", "\\1", note8b)); evi8[!grepl("^Similar to [^:]+:", note8b)] <- NA
egg8   <- fread(EMAPPER, sep = "\t", quote = "", header = TRUE, skip = "#query", na.strings = c("-","","NA"), fill = TRUE)
setnames(egg8, 1, "query"); egg8 <- egg8[!startsWith(query, "##")]; egg8[, locus := sub("-mRNA-.*$", "", query)]
egg8m  <- unique(egg8[locus %in% gid8 & !is.na(Preferred_name), .(gene_id = locus, egg = toupper(Preferred_name))], by = "gene_id")
tfmap  <- merge(data.table(gene_id = gid8, evi = evi8), egg8m, by = "gene_id", all.x = TRUE)
tfmap[, sym_use := fifelse(!is.na(evi), evi, egg)]
tfmap[, is_tf := !is.na(sym_use) & sym_use %in% tf_syms]
cat(sprintf("  %d D. laeve loci flagged as TFs\n", sum(tfmap$is_tf)))

net <- merge(mm_tbl, tfmap[, .(gene_id, sym_use, is_tf)], by = "gene_id", all.x = TRUE)
net[, is_tf := !is.na(is_tf) & is_tf]
net <- net[module != "grey"]
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

writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_batch07.txt"))
cat("[batch07] done\n")