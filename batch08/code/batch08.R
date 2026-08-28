#!/usr/bin/env Rscript

set.seed(20260426)

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges); library(Biostrings)
})
PIPE   <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
B01    <- file.path(PIPE, "batch01/objects"); B02 <- file.path(PIPE, "batch02/objects")
B03D   <- file.path(PIPE, "batch03/data")             
BATCH  <- file.path(PIPE, "batch08")
OBJ <- file.path(BATCH, "objects"); DAT <- file.path(BATCH, "data")
FIGM <- file.path(BATCH, "figures/main"); FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(OBJ, DAT, FIGM, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(FIGM, full.names = TRUE)); unlink(list.files(FIGS, full.names = TRUE))
keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")  

MINCOV <- 10L
cov10_keep <- function(bs) rowSums(as.matrix(getCoverage(bs, type = "Cov")) >= MINCOV) == ncol(bs)

cat("[0] MethylSeekR LMR/UMR segmentation (before enrichment stack)\n")
suppressPackageStartupMessages({ library(MethylSeekR); library(bsseq) })

suppressMessages(setMethod("getSeq", "DNAStringSet", function(x, names, as.character = FALSE, ...) {
  gr <- names; chrs <- as.character(seqnames(gr))
  sl <- setNames(as.integer(width(x)), names(x))              
  a <- pmax(1L, as.integer(start(gr)))                        
  b <- pmin(as.integer(sl[chrs]), as.integer(end(gr)))        
  out <- DNAStringSet(mapply(function(ch, s, e)
           if (is.null(x[[ch]]) || e < s) DNAString("") else subseq(x[[ch]], s, e),
           chrs, a, b, SIMPLIFY = FALSE))
  st <- as.character(strand(gr)); if (any(st == "-")) out[st == "-"] <- reverseComplement(out[st == "-"])
  out
}))
genome  <- readRDS(file.path(B01, "genome_chrmt.rds"))
chr_len <- setNames(width(genome), names(genome))

gff <- readRDS(file.path(B01, "gff_chrmt.rds")); gene_gr <- gff[gff$type == "gene"]

bs <- readRDS(file.path(B02, "bsseq_cov5_chrmt.rds"))
chrs <- as.character(GenomeInfoDb::seqnames(bs)); bs <- bs[chrs %in% keep_chr, ]
n0 <- nrow(bs); bs <- bs[cov10_keep(bs), ]                      
cat(sprintf("  CpGs: %s -> %s after >=%dx-in-all-4 (%.1f%% kept)\n",
            format(n0, big.mark=","), format(nrow(bs), big.mark=","), MINCOV, 100*nrow(bs)/n0))
gr <- granges(bs)
M  <- as.matrix(getCoverage(bs, type = "M")); Cv <- as.matrix(getCoverage(bs, type = "Cov"))
meth <- rowSums(M); cov <- rowSums(Cv)                         
sl <- chr_len[keep_chr]

nuc  <- as.character(seqnames(gr)) != "HiC_scaffold_1563"
m <- GRanges(as.character(seqnames(gr))[nuc],
             IRanges(start(gr)[nuc], start(gr)[nuc]), T = cov[nuc], M = meth[nuc])
GenomeInfoDb::seqlevels(m) <- keep_chr; GenomeInfoDb::seqlengths(m) <- sl


alpha_qc <- file.path(OBJ, "qc_alpha_distribution_pmd_check.pdf")
pdf(alpha_qc, width = 6, height = 5)
try(MethylSeekR::plotAlphaDistributionOneChr(m = m, chr.sel = "chr1", num.cores = 1)); dev.off()
cat("  PMD alpha check written to objects/qc_alpha_distribution_pmd_check.pdf (QC, not a figure)\n")


seg_cache <- file.path(OBJ, "methylseekr_segments_cov10_chrmt_gr.rds")   
if (file.exists(seg_cache)) {
  seg <- readRDS(seg_cache); cat("  reuse cached MethylSeekR segments\n")
} else {

  seg_pdf <- file.path(OBJ, "methylseekr_segmentation_qc.pdf")
  seg <- segmentUMRsLMRs(m = m, meth.cutoff = 0.5, nCpG.cutoff = 3L,
                         myGenomeSeq = genome, seqLengths = sl, num.cores = 1, pdfFilename = seg_pdf)
  saveRDS(seg, seg_cache)

  fwrite(data.table(chr = as.character(seqnames(seg)), start = start(seg), end = end(seg),
                    width = width(seg), type = as.character(seg$type),
                    nCG = seg$nCG, nCG_segmentation = seg$nCG.segmentation,
                    T = seg$T, M = seg$M, pmeth = seg$pmeth, median_meth = seg$median.meth),
         file.path(DAT, "methylseekr_segments.tsv"), sep = "\t")
}
lmr_gr0 <- seg[seg$type == "LMR"]
cat(sprintf("  segments: %d (LMR %d, UMR %d)\n",
            length(seg), sum(seg$type == "LMR"), sum(seg$type == "UMR")))

cat("[0a2] Takai-Jones CGI annotation + calculateFDRs\n")
TE_TJ <- "/mnt/data/alfredvar/30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv"
W_CGI <- 200L; GC_MIN <- 0.55; OE_MIN <- 0.65; LEN_MIN <- 500L; GAP_MAX <- 100L
cgi_rds <- file.path(OBJ, "cgi_takai_jones_chrmt.rds")
if (file.exists(cgi_rds)) {
  cgi_dt <- readRDS(cgi_rds)
} else {
  scan_chr_tj <- function(seqv, ch) {
    L <- length(seqv)
    if (L < W_CGI) return(data.table())
    s <- chartr("acgtn", "ACGTN", as.character(seqv))
    v <- charToRaw(s); rm(s)
    isC <- v == charToRaw("C"); isG <- v == charToRaw("G")
    cg  <- isC[-L] & isG[-1]
    cC  <- c(0, cumsum(isC)); cG <- c(0, cumsum(isG)); cCG <- c(0, cumsum(cg))
    rm(v, isC, isG, cg)
    n   <- L - W_CGI + 1L
    nC  <- cC[(W_CGI + 1):(L + 1)] - cC[1:n]
    nG  <- cG[(W_CGI + 1):(L + 1)] - cG[1:n]
    nCG <- cCG[W_CGI:L] - cCG[1:n]
    gc_ok <- (nC + nG) / W_CGI >= GC_MIN
    oe    <- ifelse(nC * nG > 0, nCG * W_CGI / (nC * nG), 0)
    pass  <- gc_ok & oe >= OE_MIN
    rm(nC, nG, nCG, gc_ok, oe)
    if (!any(pass)) return(data.table())
    r <- rle(pass)
    ends_i <- cumsum(r$lengths); starts_i <- ends_i - r$lengths + 1L
    isl <- data.table(start = starts_i[r$values], end = ends_i[r$values] + W_CGI - 1L)
    # ⚠️ data.table:: qualified — IRanges is attached after data.table, so a bare
    # shift() dispatches to IRanges::shift (no fill=) and dies (cost job 259208)
    isl[, gap := start - data.table::shift(end, fill = -1e9L)]
    isl[, grp := cumsum(gap >= GAP_MAX)]
    isl <- isl[, .(start = min(start), end = max(end)), by = grp][, grp := NULL][]
    isl[, len := end - start + 1L]
    isl[, nC  := cC[end + 1L] - cC[start]]
    isl[, nG  := cG[end + 1L] - cG[start]]
    isl[, nCG := cCG[end] - cCG[start]]          
    isl[, gc  := (nC + nG) / len]
    isl[, oe  := ifelse(nC * nG > 0, nCG * len / (nC * nG), 0)]
    isl <- isl[len >= LEN_MIN & gc >= GC_MIN & oe >= OE_MIN]
    if (nrow(isl)) isl[, chr := ch]
    isl
  }
  cgi_dt <- rbindlist(lapply(names(genome), function(ch) scan_chr_tj(genome[[ch]], ch)),
                      fill = TRUE)
  setcolorder(cgi_dt, c("chr", "start", "end", "len", "gc", "oe", "nCG"))
  saveRDS(cgi_dt, cgi_rds)
}
cgi_dt[, .(chr, start = start - 1L, end,                    
           name = sprintf("CGI_%05d", .I), score = 0L, strand = ".")] |>
  fwrite(file.path(DAT, "cgi_takai_jones.bed"), sep = "\t", col.names = FALSE)
fwrite(cgi_dt, file.path(DAT, "cgi_takai_jones.tsv"), sep = "\t")

te_tj  <- fread(TE_TJ)[chrom %in% unique(cgi_dt$chr)]
te_gr  <- reduce(GRanges(te_tj$chrom, IRanges(te_tj$start, te_tj$end)))
cgi_gr <- GRanges(cgi_dt$chr, IRanges(cgi_dt$start, cgi_dt$end))
ov_bp  <- sum(width(GenomicRanges::intersect(cgi_gr, te_gr, ignore.strand = TRUE)))
data.table(n_cgi = nrow(cgi_dt), total_bp = sum(cgi_dt$len),
           pct_genome = 100 * sum(cgi_dt$len) / sum(as.numeric(width(genome))),
           median_len = median(cgi_dt$len),
           n_mito = sum(cgi_dt$chr == "HiC_scaffold_1563"),
           pct_cgi_bp_in_te = 100 * ov_bp / sum(cgi_dt$len)) |>
  fwrite(file.path(DAT, "cgi_summary.tsv"), sep = "\t")
cat(sprintf("  %s CGIs, %.1f Mb (%.2f%% of genome), median %d bp, %.1f%% bp in TEs\n",
            format(nrow(cgi_dt), big.mark = ","), sum(cgi_dt$len) / 1e6,
            100 * sum(cgi_dt$len) / sum(as.numeric(width(genome))),
            as.integer(median(cgi_dt$len)), 100 * ov_bp / sum(cgi_dt$len)))


GenomeInfoDb::seqlevels(cgi_gr) <- keep_chr
GenomeInfoDb::seqlengths(cgi_gr) <- sl
ncpu_tj <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
fdr_tj <- calculateFDRs(m = m, CGIs = cgi_gr, num.cores = ncpu_tj,
                        pdfFilename = file.path(OBJ, "qc_methylseekr_fdr_grid.pdf"))
fdr_dt <- as.data.table(as.table(fdr_tj$FDRs))
setnames(fdr_dt, c("meth_cutoff", "n_cpg", "fdr_pct"))
fwrite(fdr_dt, file.path(DAT, "methylseekr_fdr_table.tsv"), sep = "\t")

cat(sprintf("  calculateFDRs: min FDR %.0f%% (>100%% everywhere) -> nCpG.cutoff = 3 retained\n",
            min(fdr_dt$fdr_pct, na.rm = TRUE)))
rm(te_tj, te_gr, ov_bp, fdr_tj)

suppressPackageStartupMessages({
  library(TFBSTools); library(monaLisa)               
  library(SummarizedExperiment); library(BiocParallel); library(BSgenome)
})
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

BPP  <- if (ncpu > 1) MulticoreParam(workers = min(ncpu, 8), RNGseed = 20260426L) else SerialParam(RNGseed = 20260426L)

cat("[1] JASPAR2024, 4 animal taxa + SEQUENCE-level ortholog filter\n")
JASPAR_SQLITE <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/jaspar/JASPAR2024.sqlite"
stopifnot(file.exists(JASPAR_SQLITE))
pwms <- getMatrixSet(JASPAR_SQLITE, opts = list(matrixtype = "PWM",
          tax_group = c("vertebrates", "insects", "nematodes", "urochordates")))
cat(sprintf("  %d JASPAR2024 animal PWMs (before ortholog filter)\n", length(pwms)))

SEQ_BRIDGE <- file.path(PIPE, "batch01", "data", "jaspar_ortholog_bridge.tsv")
stopifnot(file.exists(SEQ_BRIDGE))
seqb <- fread(SEQ_BRIDGE, sep = "\t")
pwms_all     <- pwms                                    
tf_names_all <- vapply(pwms_all, name, character(1))    
mids_all     <- vapply(pwms_all, ID,   character(1))

stopifnot(all(mids_all %in% seqb$motif_id))
setkey(seqb, motif_id)
keep <- seqb[.(mids_all)]$has_ortholog %in% c(TRUE, "TRUE")

fwrite(seqb[.(mids_all)], file.path(DAT, "jaspar_ortholog_bridge.tsv"), sep = "\t")
pwms     <- pwms_all[keep]
tf_names <- tf_names_all[keep]
cat(sprintf("  %d PWMs kept (TF has a SEQUENCE-verified D. laeve ortholog)\n", length(pwms)))
stopifnot(length(pwms) >= 20)


suppressPackageStartupMessages(library(universalmotif))
homer_motifs <- file.path(OBJ, "jaspar_ortholog_homer.motif")
homer_mots <- convert_motifs(pwms)


homer_thr <- vapply(homer_mots, function(m) {
  ppm <- convert_type(m, "PPM")["motif"]              
  0.8 * sum(log(pmax(apply(ppm, 2, max), 1e-3) / 0.25))   
}, numeric(1))
write_homer(homer_mots, homer_motifs, overwrite = TRUE,
            threshold = homer_thr, threshold.type = "logodds.abs")
cat(sprintf("  wrote %d ortholog-filtered motifs to HOMER format (%s); threshold %.1f-%.1f (natural log)\n",
            length(pwms), basename(homer_motifs), min(homer_thr), max(homer_thr)))


prep_gr <- function(gr, w = NULL) {
  gr <- gr[as.character(seqnames(gr)) %in% keep_chr]
  seqlevels(gr) <- keep_chr; seqlengths(gr) <- chr_len[keep_chr]
  if (is.null(w)) w <- as.integer(median(width(gr)))
  gr <- trim(resize(gr, width = max(w, 200L), fix = "center"))
  gr[width(gr) >= 100]
}
get_seqs <- function(gr, prefix = "LMR") {
  out <- DNAStringSet(rep("", length(gr)))
  for (ch in intersect(keep_chr, as.character(unique(seqnames(gr))))) {
    idx <- which(as.character(seqnames(gr)) == ch)
    out[idx] <- DNAStringSet(Views(genome[[ch]], start = start(gr)[idx], end = end(gr)[idx]))
  }
  names(out) <- sprintf("%s_%05d", prefix, seq_along(out)); out
}


cat("[3] per-LMR methylation change (Amp - Ctrl)\n")
bs <- readRDS(file.path(B02, "bsseq_cov5_chrmt.rds"))
bchr <- as.character(GenomeInfoDb::seqnames(bs)); bs <- bs[bchr %in% keep_chr, ]
bs <- bs[cov10_keep(bs), ]                             
cpg <- granges(bs)
Mm <- as.matrix(getCoverage(bs, type = "M")); Cc <- as.matrix(getCoverage(bs, type = "Cov"))
isA <- grepl("^A", sampleNames(bs))                    
mc <- rowSums(Mm[, !isA, drop = FALSE]); cc <- rowSums(Cc[, !isA, drop = FALSE])
ma <- rowSums(Mm[,  isA, drop = FALSE]); ca <- rowSums(Cc[,  isA, drop = FALSE])
ov <- findOverlaps(cpg, lmr_gr0); q <- queryHits(ov); s <- subjectHits(ov)
pl <- data.table(lmr = s, mc = mc[q], cc = cc[q], ma = ma[q], ca = ca[q])
pl <- pl[, .(beta_ctrl = sum(mc)/pmax(sum(cc), 1), beta_amp = sum(ma)/pmax(sum(ca), 1),
             ncpg = .N), by = lmr][ncpg >= 3]           
pl[, delta := beta_amp - beta_ctrl]
lmr_use <- lmr_gr0[pl$lmr]; mcols(lmr_use) <- NULL
lmr_use$beta_ctrl <- pl$beta_ctrl; lmr_use$beta_amp <- pl$beta_amp
lmr_use$delta <- pl$delta; lmr_use$ncpg <- pl$ncpg
cat(sprintf("  %d LMRs (>=3 CpGs); deltaMeth %.3f .. %.3f (median %.3f)\n",
            length(lmr_use), min(pl$delta), max(pl$delta), median(pl$delta)))
fwrite(as.data.frame(lmr_use), file.path(DAT, "lmr_methylation_change.tsv"), sep = "\t")


cat("[4] bin LMRs by deltaMeth + calcBinnedMotifEnrR\n")
lmr_p   <- prep_gr(lmr_use)                             
lmrseqs <- get_seqs(lmr_p, "LMR")
nEl  <- max(100L, as.integer(length(lmr_p) / 7L))       
bins <- bin(x = lmr_p$delta, binmode = "equalN", nElement = nEl)
cat("  bins along the methylation-change gradient:\n"); print(table(bins))

se_rds <- file.path(OBJ, "lmr_deltameth_cov10_chrmt_se.rds")   
if (file.exists(se_rds)) { se <- readRDS(se_rds); cat("  reuse cached SE\n") } else {
  se <- calcBinnedMotifEnrR(seqs = lmrseqs, bins = bins, pwmL = pwms,
                            background = "otherBins", BPPARAM = BPP, verbose = FALSE)
  saveRDS(se, se_rds)
}
enr_dt <- as.data.table(assay(se, "log2enr")); setnames(enr_dt, paste0("log2enr.", colnames(se)))
enr_dt[, `:=`(motif = rownames(se), tf = rowData(se)$motif.name)]
fwrite(enr_dt, file.path(DAT, "lmr_motif_enrichment.tsv"), sep = "\t")


cat("[5] figures\n")
suppressPackageStartupMessages({ library(ComplexHeatmap); library(grid) })
save_plot <- function(dir, name, draw, w = 7, h = 9) {

    grDevices::cairo_pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h); try(draw()); dev.off()
  png(file.path(dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 150); try(draw()); dev.off()
  svglite::svglite(file.path(dir, paste0(name, ".svg")), width = w, height = h); try(draw()); dev.off()
  cat("  saved", name, "\n")
}
show_obj <- function(p) { if (inherits(p, c("Heatmap", "HeatmapList"))) ComplexHeatmap::draw(p)
  else if (!is.null(p) && inherits(p, "gg")) print(p) else invisible(NULL) }

theme_pub <- function() {
  ggplot2::theme_classic(base_size = 9, base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, colour = "grey90"))
}
pick_sig <- function(se_) {
  pj <- assay(se_, "negLog10Padj")
  sg <- which(apply(abs(pj), 1, function(x) max(x, 0, na.rm = TRUE)) > -log10(0.05))
  if (length(sg) < 2) sg <- order(apply(pj, 1, max, na.rm = TRUE), decreasing = TRUE)[seq_len(min(25, nrow(se_)))]
  if (length(sg) > 80) sg <- sg[order(apply(pj[sg, , drop = FALSE], 1, max, na.rm = TRUE), decreasing = TRUE)][1:80]
  se_[sg, ]
}

TOP_PROMOTER <- 10L
top_n_sig <- function(se_, n = TOP_PROMOTER) {
  if (nrow(se_) <= n) return(se_)
  best <- apply(assay(se_, "negLog10Padj"), 1, max, na.rm = TRUE)
  se_[order(best, decreasing = TRUE)[seq_len(n)], ]
}

wrap_title <- function(s, width = 46) paste(strwrap(s, width = width), collapse = "\n")


draw_motif_hm <- function(se_t, title, cluster = TRUE, dendro = FALSE) {
  ms <- max(2, ceiling(max(assay(se_t, "negLog10Padj"), na.rm = TRUE)))
  hl <- plotMotifHeatmaps(x = se_t, which.plots = c("log2enr", "negLog10Padj"), width = 1.6,
                          cluster = cluster, show_dendrogram = dendro, show_seqlogo = TRUE,
                          show_motif_GC = TRUE,
                          maxEnr = 2, maxSig = ms, width.seqlogo = 1.2, doPlot = FALSE)

    ComplexHeatmap::draw(Reduce(`+`, hl), column_title = wrap_title(title),
                       column_title_gp = grid::gpar(fontface = "bold", fontsize = 11),
                       padding = grid::unit(c(2, 10, 4, 10), "mm"))
}


label_ids <- function(se_) {
  nm <- paste0(rowData(se_)$motif.name, " (", rownames(se_), ")")
  if (!is.null(rowData(se_)$n_similar)) nm <- ifelse(rowData(se_)$n_similar > 1,
    sprintf("%s x%d", nm, rowData(se_)$n_similar), nm)
  rowData(se_)$motif.name <- nm; se_
}
collapse_similar <- function(se_, cutoff = 0.95) {
  if (nrow(se_) < 3) return(se_)
  SMx <- tryCatch(monaLisa::motifSimilarity(rowData(se_)$motif.pfm, BPPARAM = BPP),
                  error = function(e) NULL)
  if (is.null(SMx)) return(se_)
  cl   <- cutree(hclust(as.dist(1 - SMx), method = "average"), h = 1 - cutoff)
  best <- apply(assay(se_, "negLog10Padj"), 1, max, na.rm = TRUE)   # cluster representative
  keep <- sort(vapply(split(seq_len(nrow(se_)), cl), function(i) i[which.max(best[i])], integer(1)))
  out  <- se_[keep, ]
  rowData(out)$n_similar <- as.integer(table(cl)[as.character(cl[keep])])
  cat(sprintf("  motif redundancy: %d rows -> %d after collapsing at similarity %.2f\n",
              nrow(se_), nrow(out), cutoff))
  out
}

dedupe_by_tf <- function(se_) {
  if (nrow(se_) < 2) return(se_)
  nm   <- rowData(se_)$motif.name
  best <- apply(assay(se_, "negLog10Padj"), 1, max, na.rm = TRUE)
  keep <- sort(vapply(split(seq_len(nrow(se_)), nm), function(i) i[which.max(best[i])], integer(1)))
  cat(sprintf("  TF dedupe: %d matrices -> %d TFs (one matrix per TF name)\n", nrow(se_), length(keep)))
  se_[keep, ]
}
motif_comp <- function(se_) {                     
  pfmL <- rowData(se_)$motif.pfm
  rbindlist(lapply(seq_along(pfmL), function(i) {
    m <- TFBSTools::Matrix(pfmL[[i]]); p <- sweep(m, 2, pmax(colSums(m), 1), "/")
    L <- ncol(p)
    data.table(motif = rownames(se_)[i], tf = rowData(se_)$motif.name[i],
               exp_cpg = if (L > 1) sum(p["C", 1:(L-1)] * p["G", 2:L]) else 0,
               gc = mean(p["C", ] + p["G", ]), len = L)
  }))
}

seSel <- pick_sig(se); hh <- max(6, 0.16 * nrow(seSel) + 2)
SM  <- tryCatch(monaLisa::motifSimilarity(rowData(seSel)$motif.pfm, BPPARAM = BPP), error = function(e) NULL)
hcl <- if (!is.null(SM) && nrow(seSel) >= 3) hclust(as.dist(1 - SM), method = "average") else TRUE
save_plot(FIGM, "fig8_lmr_motif_enrichment",
          function() draw_motif_hm(label_ids(seSel), "LMR sequences — TF motifs across the methylation-change gradient",
                                   cluster = hcl, dendro = !isTRUE(hcl)), w = 8.6, h = hh * 0.8)   


save_plot(FIGS, "figS8_lmr_bin_density",
          function() plotBinDensity(lmr_p$delta, bins,
                                    xlab = "Methylation change (amputated - control)",
                                    main = "LMRs per ΔMeth bin"),
          w = 7, h = 5)

if ("plotBinDiagnostics" %in% getNamespaceExports("monaLisa")) {
  save_plot(FIGS, "figS8_lmr_bindiag_GC",    function() show_obj(plotBinDiagnostics(lmrseqs, bins, "GCfrac")),    w = 6, h = 4.5)
  save_plot(FIGS, "figS8_lmr_bindiag_dinuc", function() show_obj(plotBinDiagnostics(lmrseqs, bins, "dinucfreq")), w = 7, h = 6)
}

save_plot(FIGS, "figS8_lmr_full_enrichment",
          function() draw_motif_hm(label_ids(seSel), "LMR sequences — all significant motifs (FDR<0.05)"), w = 7.5, h = hh)

save_plot(FIGS, "figS8_lmr_overview", function() {
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1)); on.exit(par(op))
  barplot(table(factor(seg$type, levels = c("UMR", "LMR"))), col = c(UMR = "#0072B2", LMR = "#009E73"),
          ylab = "n segments", main = sprintf("MethylSeekR segments (%d)", length(seg)))
  hist(lmr_use$delta, breaks = 40, col = "#009E73", border = "white",
       main = "LMR methylation change", xlab = "deltaMeth (Amp - Ctrl)")
  hist(log10(pmax(width(lmr_use), 1)), breaks = 40, col = "#009E73", border = "white",
       main = "LMR width", xlab = "log10(width, bp)")
}, w = 10, h = 4)


cat("[5b] LMR / UMR genomic location (bar + pie)\n")
exon_r <- reduce(granges(gff[gff$type == "exon"]))
prom_r <- reduce(trim(suppressWarnings(promoters(gene_gr, 2000, 0))))     
intr_r <- GenomicRanges::setdiff(reduce(granges(gene_gr)), exon_r)
FEAT_LV <- c("Promoter", "Exon", "Intron", "Intergenic")
feat_of <- function(x) {                       
  f <- rep("Intergenic", length(x))
  f[overlapsAny(x, intr_r)] <- "Intron"
  f[overlapsAny(x, exon_r)] <- "Exon"
  f[overlapsAny(x, prom_r)] <- "Promoter"
  factor(f, levels = FEAT_LV)
}
umr_use <- seg[seg$type == "UMR"]              
lfeat <- feat_of(lmr_use); ufeat <- feat_of(umr_use)
fwrite(rbind(data.table(region = "LMR", feature = FEAT_LV, n = as.integer(table(lfeat))),
             data.table(region = "UMR", feature = FEAT_LV, n = as.integer(table(ufeat)))),
       file.path(DAT, "lmr_umr_genomic_location.tsv"), sep = "\t")
COL_LFEAT <- c(Promoter = "#2C7FB8", Exon = "#1B9E9E", Intron = "#6A51A3", Intergenic = "#7FB3D5")
save_plot(FIGS, "figS8_lmr_genomic_location", function() {
  cnt <- table(lfeat); par(mar = c(4, 5, 3, 1))
  bp <- barplot(cnt, col = COL_LFEAT[names(cnt)], border = NA, ylab = "Number of LMRs",
                ylim = c(0, max(cnt) * 1.12), main = "LMR genomic location")
  text(bp, as.integer(cnt), labels = as.integer(cnt), pos = 3, cex = 0.9, xpd = NA)
  mtext(sprintf("n = %s", format(length(lmr_use), big.mark = ",")), side = 3, line = 0.2, cex = 0.85)
}, w = 5, h = 4)

pie_feat <- function(f, title) {
  cnt <- table(f); pct <- 100 * as.numeric(cnt) / sum(cnt)
  pie(as.numeric(cnt), labels = sprintf("%s\n%.1f%% (%s)", names(cnt), pct,
                                        format(as.integer(cnt), big.mark = ",")),

            col = COL_LFEAT[names(cnt)], border = "white", radius = 0.70, cex = 0.75,
      main = sprintf("%s  (n = %s)", title, format(length(f), big.mark = ",")))
}
save_plot(FIGS, "figS8_lmr_umr_region_pie", function() {

  op <- par(mfrow = c(1, 2), mar = c(1, 3, 3, 3)); on.exit(par(op))
  pie_feat(lfeat, "LMRs"); pie_feat(ufeat, "UMRs")
}, w = 10, h = 5)

weber_sliding <- function(gen, dt, w = 500L, off = 5L) {
  cls <- character(nrow(dt)); mxoe <- rep(NA_real_, nrow(dt)); whoe <- rep(NA_real_, nrow(dt))
  for (ac in intersect(unique(dt$seqid), names(gen))) {
    idx <- which(dt$seqid == ac); L <- length(gen[[ac]])
    v <- Views(gen[[ac]], start = pmax(1L, dt$ps[idx]), end = pmin(L, dt$pe[idx]))
    for (mm in seq_along(idx)) {
      s <- v[[mm]]; Ls <- length(s); if (Ls < w) next
      cw <- letterFrequencyInSlidingView(s, w, "C")[, 1]
      gw <- letterFrequencyInSlidingView(s, w, "G")[, 1]
      nwin <- Ls - w + 1L
      cg <- integer(Ls); st <- start(matchPattern("CG", s, fixed = TRUE)); if (length(st)) cg[st] <- 1L
      fr <- letterFrequency(s, c("C", "G"))
      whoe[idx[mm]] <- (sum(cg) * Ls) / max(fr[["C"]] * fr[["G"]], 1)   
      ccg <- cumsum(cg); cgw <- ccg[w:Ls] - c(0, ccg[1:(nwin - 1)])
      oe <- (cgw * w) / pmax(cw * gw, 1); gc <- 100 * (cw + gw) / w
      sel <- seq(1L, nwin, by = off); oe <- oe[sel]; gc <- gc[sel]
      mxoe[idx[mm]] <- max(oe, 0, na.rm = TRUE)                          
      cls[idx[mm]]  <- if (any(oe > 0.75 & gc > 55)) "HCP" else if (!any(oe > 0.48)) "LCP" else "ICP"
    }
  }
  data.table(max_oe = mxoe, whole_oe = whoe,
             weber_class = factor(cls, levels = c("HCP", "ICP", "LCP")))
}


cat("[6b] Weber-class promoter motif enrichment (categorical bins)\n")
COL_WEBER <- c(HCP = "#117733", ICP = "#88CCEE", LCP = "#CC6677")   
KEEP_BIOTYPE <- "protein_coding"          
gid9  <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
neg9  <- as.character(strand(gene_gr)) == "-"
tss9  <- ifelse(neg9, end(gene_gr), start(gene_gr))                 
bio9  <- as.character(mcols(gene_gr)$gene_biotype)
prom9 <- GRanges(seqnames(gene_gr),                                
                 IRanges(ifelse(neg9, tss9 - 200L, tss9 - 1300L),
                         ifelse(neg9, tss9 + 1300L, tss9 + 200L)),
                 strand = strand(gene_gr), gene_id = gid9, biotype = bio9)
prom9 <- prom9[bio9 %in% KEEP_BIOTYPE & as.character(seqnames(prom9)) %in% keep_chr]
seqlevels(prom9) <- keep_chr; seqlengths(prom9) <- chr_len[keep_chr]
prom9 <- trim(prom9); prom9 <- prom9[width(prom9) >= 100]          
cat(sprintf("  promoter universe: %d genes (%s; pseudogenes dropped)\n",
            length(prom9), paste(sprintf("%s %d", names(table(prom9$biotype)),
                                         as.integer(table(prom9$biotype))), collapse = " + ")))

wcls <- weber_sliding(genome, data.table(seqid = as.character(seqnames(prom9)),
                                         ps = start(prom9), pe = end(prom9)))
prom9$weber_class <- wcls$weber_class
prom9 <- prom9[!is.na(prom9$weber_class)]
wcls <- wcls[!is.na(weber_class)]
weber_dt <- data.table(gene_id = prom9$gene_id, biotype = prom9$biotype,
                       max_oe = wcls$max_oe, whole_oe = wcls$whole_oe,
                       weber_class = as.character(prom9$weber_class))
fwrite(weber_dt, file.path(DAT, "promoter_weber_classification_pc.tsv"), sep = "\t")
wbins <- prom9$weber_class
cat("  Weber promoter bins:\n"); print(table(wbins))
cat("  by biotype:\n"); print(table(prom9$biotype, wbins))
promseqs <- get_seqs(prom9, "PROM")
se_w_rds <- file.path(OBJ, "promoter_weber_pc_chrmt_se.rds")          
if (file.exists(se_w_rds)) { se_w <- readRDS(se_w_rds); cat("  reuse cached SE\n") } else {
  se_w <- calcBinnedMotifEnrR(seqs = promseqs, bins = wbins, pwmL = pwms,
                              background = "otherBins", BPPARAM = BPP, verbose = FALSE)
  saveRDS(se_w, se_w_rds)
}
enr_w <- as.data.table(assay(se_w, "log2enr")); setnames(enr_w, paste0("log2enr.", colnames(se_w)))
enr_w[, `:=`(motif = rownames(se_w), tf = rowData(se_w)$motif.name)]
fwrite(enr_w, file.path(DAT, "promoter_weber_motif_enrichment.tsv"), sep = "\t")
seW  <- pick_sig(se_w)
nW   <- table(wbins)
ttlW <- sprintf("Promoter TF motifs across Weber CpG classes (HCP %s / ICP %s / LCP %s)",
                format(nW[["HCP"]], big.mark=","), format(nW[["ICP"]], big.mark=","),
                format(nW[["LCP"]], big.mark=","))
seWc <- top_n_sig(dedupe_by_tf(seW))                         

save_plot(FIGS, "figS8_promoter_weber_motif_enrichment",
          function() draw_motif_hm(label_ids(seWc), ttlW),
          w = 10, h = max(5, 0.30 * nrow(seWc) + 2.5))

gcf <- rowSums(letterFrequency(promseqs, c("G", "C"), as.prob = TRUE))   
save_plot(FIGS, "figS8_promoter_weber_gc", function() {
  par(mar = c(4, 4, 3, 1))
  boxplot(gcf ~ wbins, col = COL_WEBER[levels(wbins)], border = "grey30",
          ylab = "Promoter GC fraction", xlab = "Weber class",
          main = "Promoter GC by Weber class (HCP is high-GC by definition)")
}, w = 5, h = 4)


cat("[6c] LMR / UMR overlap with Weber promoter classes\n")
p_lmr <- overlapsAny(prom9, lmr_use); p_umr <- overlapsAny(prom9, umr_use)
ov_dt <- rbindlist(lapply(c("LMR", "UMR"), function(k) {
  hit <- if (k == "LMR") p_lmr else p_umr
  rbindlist(lapply(levels(wbins), function(cl) {
    inc <- wbins == cl
    ft  <- fisher.test(table(factor(inc, c(TRUE, FALSE)), factor(hit, c(TRUE, FALSE))))
    data.table(region = k, weber_class = cl, n_promoters = sum(inc), n_hit = sum(hit & inc),
               pct_hit = 100 * mean(hit[inc]), or = unname(ft$estimate), p = ft$p.value)
  }))
}))
ov_dt[, fdr := p.adjust(p, "BH")]
fwrite(ov_dt, file.path(DAT, "weber_class_lmr_umr_overlap.tsv"), sep = "\t")
print(ov_dt[, .(region, weber_class, n_promoters, n_hit, pct_hit = round(pct_hit, 2),
                or = round(or, 2), fdr = signif(fdr, 2))])
save_plot(FIGS, "figS8_weber_class_lmr_umr", function() {
  mm <- matrix(ov_dt$pct_hit, nrow = 2, byrow = TRUE,
               dimnames = list(c("LMR", "UMR"), levels(wbins)))
  par(mar = c(4, 4.5, 3, 1))
  bp <- barplot(mm, beside = TRUE, col = c(LMR = "#009E73", UMR = "#0072B2"), border = NA,
                ylab = "% of promoters overlapping", xlab = "Weber promoter class",
                ylim = c(0, max(mm) * 1.2), main = "LMRs and UMRs by promoter CpG class")
  text(bp, mm, sprintf("%.1f%%", mm), pos = 3, cex = 0.7, xpd = NA)
  legend("topright", legend = c("LMR", "UMR"), fill = c("#009E73", "#0072B2"), bty = "n")
}, w = 6, h = 4.2)


cat("[6e] HOMER motif enrichment (UMR promoters, LMRs, DMRs): known motifs\n")
HOMER <- "/mnt/data/alfredvar/rlopezt/meth_paper/tools/homer/bin/findMotifsGenome.pl"
if (!file.exists(HOMER)) cat("  SKIPPED — HOMER not installed at tools/homer\n") else {
  Sys.setenv(PATH = paste(dirname(HOMER), Sys.getenv("PATH"), sep = ":"))  
  ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
  genome_fa <- file.path(OBJ, "genome_chrmt.fa")
  if (!file.exists(genome_fa)) { Biostrings::writeXStringSet(genome, genome_fa); cat("  wrote genome FASTA\n") }


  HOMER_PITCH <- 25
  svg_inner <- function(f) {                      
    s <- paste(readLines(f, warn = FALSE), collapse = "\n")
    list(w = as.numeric(sub('.*<svg[^>]*width="([0-9.]+)".*',  "\\1", s)),
         h = as.numeric(sub('.*<svg[^>]*height="([0-9.]+)".*', "\\1", s)),
         body = sub("</svg>\\s*$", "", sub("^.*?<svg[^>]*>", "", s)))
  }

  xml_esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", x))

  RSVG <- "/usr/bin/rsvg-convert"

  save_homer_svg <- function(rows, nm, title, ncol = 2, outdir_fig = FIGS) {

    COLW <- 300; LOGO_H <- 46; LOGO_W <- 258
    ROWH <- LOGO_H + 34; TOP <- 62          
    nr <- ceiling(length(rows) / ncol)
    body <- vapply(seq_along(rows), function(i) {
      r  <- rows[[i]]
      x0 <- 30 + ((i - 1) %% ncol) * COLW
      y0 <- TOP + floor((i - 1) / ncol) * ROWH
      s  <- min(LOGO_H / r$lg$h, LOGO_W / r$lg$w)   
      pitch <- HOMER_PITCH * s; lh <- r$lg$h * s
      rect <- if (length(r$cg)) paste(sprintf(
        '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#9a9a9a" opacity="0.32"/>',
        x0 + (r$cg - 1) * pitch, y0, 2 * pitch, lh), collapse = "\n") else ""
      sprintf('%s\n<g transform="translate(%.2f,%.2f) scale(%.4f)">%s</g>
<text x="%.1f" y="%.1f" font-family="Arial" font-size="14" font-weight="bold">%s</text>
<text x="%.1f" y="%.1f" font-family="Arial" font-size="10" fill="#4d4d4d">%s</text>',
        rect, x0, y0, s, r$lg$body, x0, y0 - 15, xml_esc(r$tf), x0, y0 - 5, xml_esc(r$st))
    }, character(1))
    svg <- sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d">
<rect width="100%%" height="100%%" fill="white"/>
<text x="22" y="28" font-family="Arial" font-size="16" font-weight="bold">%s</text>
%s</svg>', 30 + ncol * COLW, as.integer(TOP + nr * ROWH),
      xml_esc(title), paste(body, collapse = "\n"))
    f_svg <- file.path(outdir_fig, paste0(nm, ".svg")); writeLines(svg, f_svg)
    if (file.exists(RSVG)) {
      system2(RSVG, c("-f","pdf","-o", shQuote(file.path(outdir_fig, paste0(nm,".pdf"))), shQuote(f_svg)))
      system2(RSVG, c("-f","png","-d","300","-p","300",
                      "-o", shQuote(file.path(outdir_fig, paste0(nm,".png"))), shQuote(f_svg)))
    } else cat("  rsvg-convert missing - wrote .svg only\n")
    cat(sprintf("  saved %s (%d HOMER logos)\n", nm, length(rows)))
  }

  
  FOLD_MIN <- 1.5; UBIQ <- 95; TOPN <- 10; CORE <- 6
  pick_homer <- function(kr) {
    ok <- kr[!is.na(fold) & q < 0.05 & fold >= FOLD_MIN & pct_target_num <= UBIQ][order(-fold)]

    keep <- integer(0)
    for (i in seq_len(nrow(ok))) {
      cons <- ok$consensus[i]
      subs <- if (nchar(cons) >= CORE)
        unique(substring(cons, 1:(nchar(cons)-CORE+1), CORE:nchar(cons))) else cons
      dup <- length(keep) > 0 && any(vapply(ok$consensus[keep], function(k)
               any(vapply(subs, grepl, logical(1), x = k, fixed = TRUE)), logical(1)))
      if (!dup) keep <- c(keep, i)
      if (length(keep) >= TOPN) break
    }
    ok[keep]
  }

  run_homer <- function(gr, tag, title, cpg = FALSE, main = FALSE) {
    bed <- file.path(OBJ, sprintf("homer_%s_chrmt.bed", tag))
    fwrite(data.table(chr = as.character(seqnames(gr)), start = start(gr) - 1L, end = end(gr),
                      id = sprintf("%s_%05d", tag, seq_along(gr)), score = 0L, strand = "+"),
           bed, sep = "\t", col.names = FALSE)
    outdir <- file.path(OBJ, sprintf("homer_%s_chrmt", tag)); log <- file.path(OBJ, sprintf("homer_%s_chrmt.log", tag))
    args <- c(bed, genome_fa, outdir, "-mknown", homer_motifs, "-p", ncpu,   # ortholog-filtered JASPAR motifs, NOT -mset vertebrates
              "-size", "given", "-len", "6,8,10,12", "-nomotif",   # KNOWN motifs only (de novo dropped 2026-07-24)
              "-preparsedDir", file.path(OBJ, "homer_preparsed"))
    if (cpg) args <- c(args, "-cpg")
    kr_file <- file.path(outdir, "knownResults.txt")

    if (file.exists(kr_file)) {
      cat(sprintf("  [%s] reuse cached HOMER scan (%s)\n", tag, basename(outdir)))
    } else {
      cat(sprintf("  [%s] findMotifsGenome.pl on %s regions (known motifs only, %s norm)...\n",
                  tag, format(length(gr), big.mark=","), if (cpg) "CpG%" else "GC%"))
      st <- system2(HOMER, args, stdout = log, stderr = log)
      if (!file.exists(kr_file)) { cat(sprintf("  [%s] no knownResults.txt (exit %s); see %s\n", tag, st, basename(log))); return(invisible()) }
    }
    kr <- fread(kr_file)
    setnames(kr, c("motif","consensus","p","logp","q","n_target","pct_target","n_bg","pct_bg")[seq_len(ncol(kr))])

    kr[, `:=`(rank = .I, tf = sub("/.*", "", motif),
              pct_target_num = as.numeric(sub("%", "", pct_target)),
              pct_bg_num = as.numeric(sub("%", "", pct_bg)))]
    kr[, `:=`(P = exp(logp),
              fold = fifelse(pct_bg_num > 0, pct_target_num / pct_bg_num, NA_real_),
              cpg_motif = grepl("CG", consensus, fixed = TRUE))]
    fwrite(kr, file.path(DAT, sprintf("%s_homer_known_motifs.tsv", tag)), sep = "\t")
    cat(sprintf("  [%s] HOMER known motifs q<0.05: %d; top: %s\n", tag,
                sum(kr$q < 0.05, na.rm = TRUE), paste(head(kr[order(logp)]$tf, 6), collapse = ", ")))

    if (sum(kr$q < 0.05, na.rm = TRUE) == 0)
      cat(sprintf("  [%s] NULL: min q = %.3f; %d/%d motifs have target hits, %d/%d background hits%s\n",
                  tag, min(kr$q, na.rm = TRUE), sum(kr$n_target > 0), nrow(kr),
                  sum(kr$n_bg > 0), nrow(kr),
                  if (sum(kr$n_target > 0) == 0) "  <-- ZERO hits: thresholds are broken, NOT a real null" else "  (real null)"))

    sel <- pick_homer(kr)
    n_sig  <- sum(kr$q < 0.05, na.rm = TRUE)
    n_ubiq <- sum(kr$q < 0.05 & kr$pct_target_num > UBIQ, na.rm = TRUE)
    n_weak <- sum(kr$q < 0.05 & kr$pct_target_num <= UBIQ &
                  (is.na(kr$fold) | kr$fold < FOLD_MIN), na.rm = TRUE)
    if (!nrow(sel)) {

      cat(sprintf("  [%s] no motif passes q<0.05 & fold>=%.1f -- NO figure written (null reported in %s_homer_known_motifs.tsv)\n",
                  tag, FOLD_MIN, tag)); return(invisible())
    }
    miss <- 0L
    rows <- lapply(seq_len(nrow(sel)), function(i) {
      f <- file.path(outdir, "knownResults", sprintf("known%d.logo.svg", sel$rank[i]))
      if (!file.exists(f)) { miss <<- miss + 1L; return(NULL) }   
      cons <- sel$consensus[i]; m <- gregexpr("CG", cons, fixed = TRUE)[[1]]
      list(lg = svg_inner(f), cg = if (m[1] > 0) as.integer(m) else integer(0), tf = sel$tf[i],
           st = sprintf("fold %.2fx | P = %.1e, q %s | %s vs %s bg", sel$fold[i], sel$P[i],
                        if (sel$q[i] < 1e-4) "< 1e-4" else sprintf("= %.4f", sel$q[i]),
                        sel$pct_target[i], sel$pct_bg[i]))
    })
    rows <- Filter(Negate(is.null), rows)
    if (miss) cat(sprintf("  [%s] %d selected motifs had no HOMER logo file\n", tag, miss))
    if (!length(rows)) { cat(sprintf("  [%s] no HOMER logos available -- no figure\n", tag)); return(invisible()) }
    cat(sprintf("  [%s] %d motifs q<0.05; showing top %d by fold (dropped %d ubiquitous >%d%% of regions, %d with fold<%.1f; families collapsed)\n",
                tag, n_sig, length(rows), n_ubiq, UBIQ, n_weak, FOLD_MIN))
    save_homer_svg(rows,
      sprintf(if (main) "fig8_%s_homer_motifs" else "figS8_%s_homer_motifs", tag), title,
      outdir_fig = if (main) FIGM else FIGS)
  }

  ph <- trim(suppressWarnings(promoters(gene_gr, upstream = 2000, downstream = 200)))
  strand(ph) <- "*"; prom_hp <- reduce(ph)
  umr5 <- umr_use[width(umr_use) < 5000L]
  ov   <- findOverlaps(umr5, prom_hp)
  frac <- width(pintersect(umr5[queryHits(ov)], prom_hp[subjectHits(ov)])) / width(umr5[queryHits(ov)])
  umr_prom <- umr5[unique(queryHits(ov)[frac >= 0.8])]
  cat(sprintf("  UMRs <5kb: %s -> %s with >=80%% promoter overlap\n",
              format(length(umr5), big.mark=","), format(length(umr_prom), big.mark=",")))
  run_homer(umr_prom, "umr_promoter", "TF motifs at unmethylated promoters",
            cpg = TRUE, main = TRUE)     

  run_homer(lmr_use, "lmr", "TF motifs at low-methylated regions", cpg = FALSE)

  dmr_hg <- with(fread(file.path(PIPE, "batch05/data/dmrs_annotated.tsv"))[chr %in% keep_chr],
                 GRanges(chr, IRanges(start, end)))
  run_homer(dmr_hg, "dmr", "TF motifs at differentially methylated regions", cpg = FALSE)
}


cat("[6f] HUMAN GRCh38 Weber-class promoter motif enrichment (JASPAR human motifs)\n")
pwms_h <- getMatrixSet(JASPAR_SQLITE, opts = list(matrixtype = "PWM", species = "9606"))
cat(sprintf("  %d JASPAR2024 CORE human (tax_id 9606) PWMs for the human run\n", length(pwms_h)))
B03DS  <- file.path(PIPE, "batch03/dataset")
hs_fna <- file.path(B03DS, "GRCh38_latest_genomic.fna.gz")
hs_gff <- file.path(B03DS, "GRCh38_latest_genomic.gff.gz")
human_ok <- file.exists(hs_fna) && file.exists(hs_gff)
if (!human_ok) cat("  SKIPPED — human GRCh38 reference not staged in batch03/dataset\n")
if (human_ok) {
 
  hs_primary <- c(sprintf("NC_0000%02d", 1:22), "NC_000023", "NC_000024")
  hg <- fread(cmd = sprintf("zcat '%s' | grep -v '^#'", hs_gff), sep = "\t", header = FALSE, quote = "",
              col.names = c("seqid","src","type","start","end","score","strand","phase","attr"))
  hg <- hg[type == "gene" & sub("\\..*", "", seqid) %in% hs_primary]
  hg[, biotype := sub(".*gene_biotype=([^;]+).*", "\\1", attr)]
  hg <- hg[biotype %in% KEEP_BIOTYPE]   # SAME universe as the slug: protein_coding only (Weber's set)
  hg[, symbol := sub(".*;gene=([^;]+).*", "\\1", attr)]
  hg[, tss := ifelse(strand == "-", end, start)]; hneg <- hg$strand == "-"
  hg[, `:=`(ps = ifelse(hneg, tss - 200L, tss - 1300L), pe = ifelse(hneg, tss + 1300L, tss + 200L))]
  hs_genome <- readDNAStringSet(hs_fna)
  names(hs_genome) <- sub(" .*", "", names(hs_genome))
  hs_genome <- hs_genome[sub("\\..*", "", names(hs_genome)) %in% hs_primary]
  hcls <- weber_sliding(hs_genome, hg)
  hg[, weber_class := hcls$weber_class]
  hg <- hg[!is.na(weber_class)]
  cat(sprintf("  %d human protein-coding promoters -> HCP %.0f%% / ICP %.0f%% / LCP %.0f%%\n",
              nrow(hg), 100*mean(hg$weber_class=="HCP"), 100*mean(hg$weber_class=="ICP"),
              100*mean(hg$weber_class=="LCP")))

    hg[, `:=`(ps = pmax(1L, ps), pe = pmin(as.integer(width(hs_genome))[match(seqid, names(hs_genome))], pe))]
  hg <- hg[pe - ps + 1L >= 100]
  hseqs <- DNAStringSet(unlist(lapply(split(seq_len(nrow(hg)), hg$seqid), function(i) {
    ac <- hg$seqid[i[1]]
    as.list(DNAStringSet(Views(hs_genome[[ac]], start = hg$ps[i], end = hg$pe[i])))
  }), use.names = FALSE))
  hord  <- unlist(split(seq_len(nrow(hg)), hg$seqid), use.names = FALSE)  
  hbins <- hg$weber_class[hord]
  names(hseqs) <- sprintf("HSPROM_%05d", seq_along(hseqs))
  cat("  human Weber promoter bins:\n"); print(table(hbins))
  rm(hs_genome); invisible(gc())
  se_h_rds <- file.path(OBJ, "human_promoter_weber_hs9606_se.rds")   
  if (file.exists(se_h_rds)) { se_h <- readRDS(se_h_rds); cat("  reuse cached SE\n") } else {
    se_h <- calcBinnedMotifEnrR(seqs = hseqs, bins = hbins, pwmL = pwms_h,  
                                background = "otherBins", BPPARAM = BPP, verbose = FALSE)
    saveRDS(se_h, se_h_rds)
  }
  enr_h <- as.data.table(assay(se_h, "log2enr")); setnames(enr_h, paste0("log2enr.", colnames(se_h)))
  enr_h[, `:=`(motif = rownames(se_h), tf = rowData(se_h)$motif.name)]
  fwrite(enr_h, file.path(DAT, "human_promoter_weber_motif_enrichment.tsv"), sep = "\t")
  nH <- table(hbins)
  seH <- top_n_sig(dedupe_by_tf(pick_sig(se_h)))
  save_plot(FIGS, "figS8_human_promoter_weber_motif_enrichment",
            function() draw_motif_hm(label_ids(seH),
              sprintf("HUMAN promoter TF motifs across Weber classes (HCP %s / ICP %s / LCP %s)",
                      format(nH[["HCP"]], big.mark=","), format(nH[["ICP"]], big.mark=","),
                      format(nH[["LCP"]], big.mark=","))),
            w = 10, h = max(5, 0.30 * nrow(seH) + 2.5))   


  cmp9 <- merge(
    data.table(motif = rownames(se_w), tf = rowData(se_w)$motif.name,
               dlaeve_HCP = assay(se_w, "log2enr")[, "HCP"],
               dlaeve_padj = assay(se_w, "negLog10Padj")[, "HCP"]),
    data.table(motif = rownames(se_h),
               human_HCP = assay(se_h, "log2enr")[, "HCP"],
               human_padj = assay(se_h, "negLog10Padj")[, "HCP"]), by = "motif")
  cmp9 <- cmp9[is.finite(dlaeve_HCP) & is.finite(human_HCP)]        # intersection, both scored
  cmp9[, delta := dlaeve_HCP - human_HCP]
  cmp9[, sig_both := dlaeve_padj > -log10(0.05) & human_padj > -log10(0.05)]
  setorder(cmp9, -delta)
  fwrite(cmp9, file.path(DAT, "hcp_motif_enrichment_human_vs_dlaeve.tsv"), sep = "\t")
  cat(sprintf("  HCP enrichment correlates across species: Pearson r = %.2f (n = %d motifs)\n",
              cor(cmp9$dlaeve_HCP, cmp9$human_HCP, use = "complete.obs"), nrow(cmp9)))
  cat("  most SLUG-specific HCP motifs:\n"); print(head(cmp9[, .(tf, motif, dlaeve_HCP, human_HCP, delta)], 8))
  cat("  most HUMAN-specific HCP motifs:\n"); print(tail(cmp9[, .(tf, motif, dlaeve_HCP, human_HCP, delta)], 8))
  lab9 <- cmp9[order(-abs(delta))][1:15]
  save_plot(FIGS, "figS8_hcp_motif_human_vs_dlaeve", function() {
    print(ggplot2::ggplot(cmp9, ggplot2::aes(human_HCP, dlaeve_HCP)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
      ggplot2::geom_point(ggplot2::aes(colour = sig_both), size = 1.1, alpha = 0.8) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70"),
                                   name = "FDR<0.05 in both") +
      ggrepel::geom_text_repel(data = lab9, ggplot2::aes(label = sprintf("%s (%s)", tf, motif)),
                               size = 2.2, max.overlaps = 20, segment.colour = "grey70") +
      ggplot2::labs(x = "log2 enrichment in HUMAN HCP promoters",
                    y = "log2 enrichment in D. laeve HCP promoters",
                    title = "Which TF motifs mark CpG-island promoters in each species?",
                    subtitle = sprintf("motifs scored in both (D. laeve-ortholog set vs JASPAR human set); n = %d, Pearson r = %.2f",
                                       nrow(cmp9), cor(cmp9$dlaeve_HCP, cmp9$human_HCP, use = "complete.obs"))) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey30")))
  }, w = 7, h = 6)

}

cat("[6g] REMOVED — GC-matched genome-background check dropped (author, 2026-08-26)\n")

cat("[7] DMR motif enrichment: REMOVED — DMRs are tested with HOMER in §6e (see comment)\n")

writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_batch08.txt"))
cat("[batch08] done\n")
