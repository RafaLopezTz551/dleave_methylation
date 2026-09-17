#!/usr/bin/env Rscript
# 02_landscape.R
# Global CpG methylation landscape of the D. laeve tail: per sample, per chromosome and per genomic
# region, its relation to gene expression (gene bodies, TSS, discrete gene segments), and the
# agreement between PacBio HiFi and WGBS methylation calls.
# Inputs: four Bismark CpG reports (C1, C2 control; A1, A2 amputated), EviAnn GFF, HTSeq counts (tail
#   and bodywall), HiFi 5mC beds (00_data_pacbio_hifi), 01_genome_toolkit objects/gff_chrmt.rds and
#   data/gene_de_tail.tsv.
# Outputs: data/    per-sample, window, region, decile, metagene and agreement tables (.tsv)
#          objects/ bsseq_cov5_chrmt.rds, hifi_bodywall_cpg_persample_chrmt.rds
#          figures/ main fig2a, b, d, e, f, g, j; supplementary figS_* and fig2i (.pdf, .png, .svg)
# Run: sbatch 02_landscape/code/02_landscape.slurm from main/methylation_pipeline/

# Step 1 - Seed, packages, paths, chromosome set, palettes, theme and figure saver
set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges)
  library(Biostrings); library(bsseq); library(ggplot2); library(scales)
  library(patchwork)
})

PIPE   <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
GFF    <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff"
HTSEQ  <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"
MCALL  <- "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls"
CPG <- c(C1 = file.path(MCALL, "C1.CpG_report.txt.gz"),
         C2 = file.path(MCALL, "C2.CpG_report.txt.gz"),
         A1 = file.path(MCALL, "A1.CpG_report.txt.gz"),
         A2 = file.path(MCALL, "A2.CpG_report.txt.gz"))
B01    <- file.path(PIPE, "01_genome_toolkit/objects")   # 01_genome_toolkit objects
BATCH  <- file.path(PIPE, "02_landscape")
OBJ <- file.path(BATCH, "objects"); DAT <- file.path(BATCH, "data")
FIGM <- file.path(BATCH, "figures/main"); FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(OBJ, DAT, FIGM, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")   # chr1-31 + mito scaffold
`%||%` <- function(a, b) if (is.null(a)) b else a

COL_COND <- c(Control = "#2166AC", Amputated = "#B2182B")
COL_REGION <- c(Promoter = "#2C7FB8", Exon = "#1B9E9E", Intron = "#6A51A3", Intergenic = "#7FB3D5")
DECILE_PAL <- grDevices::colorRampPalette(c("#E07B6A","#6FAE6F","#5FB3C4","#7A5BAA","#3D2A66"))(10)

theme_pub <- function(base_size = 9) theme_classic(base_size = base_size, base_family = "sans") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"))
# cairo_pdf so Unicode glyphs in axis labels render in the PDF
save_fig <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(dir, paste0(name, ".svg")), p, width = w, height = h)   # vector (svg)
  cat(sprintf("  saved %s\n", name))
}

# Step 2 - Strand-collapsed bsseq from the four CpG reports, cov >= 5 in all (cached)
cat("[1] bsseq (strand-collapsed CpG, cov>=5 all samples)\n")
bs_rds <- file.path(OBJ, "bsseq_cov5_chrmt.rds")
stale_cache <- function(cache, inputs) !file.exists(cache) || any(file.mtime(inputs) > file.mtime(cache))   # rebuild when any input is newer than the cache
if (!stale_cache(bs_rds, CPG)) {
  bs <- readRDS(bs_rds)
} else {
  read_one <- function(f) {
    d <- fread(f, header = FALSE,
               col.names = c("chr","pos","strand","meth","unmeth","context","tri"))
    d <- d[chr %in% keep_chr & context == "CG"]
    d[strand == "-", pos := pos - 1L]                 # collapse to + strand CpG
    d[, .(meth = sum(meth), unmeth = sum(unmeth)), by = .(chr, pos)]
  }
  dl <- lapply(names(CPG), function(s) { x <- read_one(CPG[[s]]); x[, sample := s]; x })
  dt <- rbindlist(dl)
  wide <- dcast(dt, chr + pos ~ sample, value.var = c("meth","unmeth"))
  for (s in names(CPG)) wide[, paste0("cov_", s) := get(paste0("meth_", s)) + get(paste0("unmeth_", s))]
  keep <- wide[, Reduce(`&`, lapply(names(CPG), function(s) get(paste0("cov_", s)) >= 5))]
  wide <- wide[keep]
  M  <- as.matrix(wide[, lapply(names(CPG), function(s) get(paste0("meth_", s)))])
  Cv <- as.matrix(wide[, lapply(names(CPG), function(s) get(paste0("cov_",  s)))])
  colnames(M) <- colnames(Cv) <- names(CPG)
  bs <- BSseq(chr = wide$chr, pos = wide$pos, M = M, Cov = Cv, sampleNames = names(CPG))
  pData(bs)$condition <- c("Control","Control","Amputated","Amputated")
  saveRDS(bs, bs_rds)
}
chrs <- as.character(seqnames(bs)); bs <- bs[chrs %in% keep_chr, ]   # defensive re-filter of the cached object
gr <- granges(bs)
M  <- as.matrix(getCoverage(bs, type = "M")); Cv <- as.matrix(getCoverage(bs, type = "Cov"))
samp <- sampleNames(bs); cond <- as.character(pData(bs)$condition)   # design taken from the object
stopifnot(identical(sort(unique(cond)), c("Amputated", "Control")))
cat(sprintf("  %s CpGs x %d samples\n", format(nrow(M), big.mark = ","), ncol(M)))

# Pooled per-condition counts and per-sample mean beta, used by several panels
M_ctrl <- rowSums(M[, cond == "Control", drop = FALSE]); C_ctrl <- rowSums(Cv[, cond == "Control", drop = FALSE])
M_amp  <- rowSums(M[, cond == "Amputated", drop = FALSE]); C_amp  <- rowSums(Cv[, cond == "Amputated", drop = FALSE])
beta_cpg_per <- colMeans(M / pmax(Cv, 1L), na.rm = TRUE)
gmean_global <- mean(beta_cpg_per)

# Step 3 - GFF gene models (01_genome_toolkit object when present, otherwise imported)
cat("[2] gff\n")
gff_rds <- file.path(B01, "gff_chrmt.rds")
if (file.exists(gff_rds)) {
  gff <- readRDS(gff_rds)
} else {
  suppressPackageStartupMessages(library(rtracklayer))
  gff <- import(GFF); gff <- gff[as.character(seqnames(gff)) %in% keep_chr]
  gff <- GenomeInfoDb::keepSeqlevels(gff, intersect(keep_chr, GenomeInfoDb::seqlevels(gff)),
                                     pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(gff) <- keep_chr
}
gene_gr <- gff[gff$type == "gene"]

# Step 4 - Tail expression deciles (HTSeq -> DESeq2 VST) and the per-CpG pooled table
cat("[3] tail expression deciles\n")
suppressPackageStartupMessages(library(DESeq2))
tail_s <- c("C1S1","C2S2","C3S3","C4S4","T2S6","T3S7","T4S8")
cl <- lapply(tail_s, function(s) {
  x <- fread(file.path(HTSEQ, paste0(s, "_htseq_gene_counts.txt")),
             header = FALSE, col.names = c("gene_id","count"))
  x[!startsWith(gene_id, "__")]
})
gids <- cl[[1]]$gene_id
stopifnot(sapply(cl, function(x) identical(x$gene_id, gids)))  # HTSeq files must share row order or counts misassign silently
cm <- sapply(cl, function(x) x$count); rownames(cm) <- gids; colnames(cm) <- tail_s
gid_chr <- sub(";.*", "", as.character(mcols(gff[gff$type == "gene"])$ID))   # keep_chr gene universe (gff is already filtered)
cm <- cm[rownames(cm) %in% gid_chr, , drop = FALSE]   # keep_chr universe before the VST
cm <- cm[rowSums(cm >= 5) >= 2, , drop = FALSE]        # expressed: >= 5 reads in >= 2 libraries (same prefilter as 01's DE)
vsd <- vst(DESeqDataSetFromMatrix(cm, data.frame(s = tail_s), ~ 1), blind = TRUE)
expr <- rowMeans(assay(vsd))
# Per-condition expression means, so methylation can be compared to the matching transcriptome
ctrl_rna <- c("C1S1","C2S2","C3S3","C4S4"); amp_rna <- c("T2S6","T3S7","T4S8")
ex_dt <- data.table(gene_id = names(expr), expr_mean = as.numeric(expr),
                    expr_ctrl = as.numeric(rowMeans(assay(vsd)[, ctrl_rna, drop = FALSE])),
                    expr_amp  = as.numeric(rowMeans(assay(vsd)[, amp_rna,  drop = FALSE])))

gid_all <- sub(";.*", "", as.character(mcols(gene_gr)$ID))
gene_dt <- data.table(gene_id = gid_all, chr = as.character(seqnames(gene_gr)),
                      start = start(gene_gr), end = end(gene_gr),
                      strand = as.character(strand(gene_gr)))
gene_dt[, tss := ifelse(strand == "-", end, start)]
gene_dt <- merge(gene_dt, ex_dt, by = "gene_id")
gene_dt <- gene_dt[is.finite(expr_mean)]
gene_dt[, decile := cut(expr_mean, quantile(expr_mean, seq(0, 1, 0.1), na.rm = TRUE),
                        labels = 1:10, include.lowest = TRUE)]

# Per-CpG pooled coverage and beta for control and amputated, keyed for overlaps
cpg_r <- data.table(chr = as.character(seqnames(gr)), pos = start(gr),
                    cv_c = C_ctrl, cv_a = C_amp,
                    bg_c = M_ctrl / pmax(C_ctrl, 1), bg_a = M_amp / pmax(C_amp, 1))
cpg_r[, `:=`(start = pos, end = pos)]

# Step 5 - figS_per_chromosome_methylation: per-sample beta in 1 Mb windows
cat("[4] per-chromosome methylation (supplementary)\n")
chr_vec <- as.character(seqnames(gr)); win_idx <- (start(gr) %/% 1e6) * 1e6
per <- rbindlist(lapply(seq_len(ncol(M)), function(i)
  data.table(chr = chr_vec, win = win_idx, M = M[, i], Cv = Cv[, i]
  )[, .(beta = sum(M)/pmax(sum(Cv),1), n_cpg = sum(Cv > 0)), by = .(chr, win)
  ][, `:=`(sample = samp[i], condition = cond[i])]))
per <- per[n_cpg >= 50]   # window floor: >=50 covered CpGs
per[, chr := factor(chr, levels = keep_chr)]; per[, win_mb := win/1e6 + 0.5]
cond_mean <- per[, .(beta = sum(beta*n_cpg)/pmax(sum(n_cpg),1)), by = .(chr, win_mb, condition)]
pa <- ggplot(per, aes(win_mb, beta*100)) +
  geom_hline(yintercept = gmean_global*100, colour = "#C0392B", linetype = "dashed", linewidth = 0.35) +
  geom_point(aes(colour = condition), size = 0.35, alpha = 0.45) +
  geom_line(data = cond_mean, aes(colour = condition, group = condition), linewidth = 0.45) +
  facet_wrap(~ chr, ncol = 4, scales = "free_x") +
  scale_colour_manual(values = COL_COND, name = NULL) +
  labs(x = "Position (Mb)", y = "Mean CpG methylation β (%)",
       title = "Genome-wide methylation, per chromosome") +
  theme_pub() + theme(strip.background = element_blank(),
                      strip.text = element_text(size = 6, face = "bold"),
                      axis.text = element_text(size = 5), legend.position = "bottom")
save_fig(pa, FIGS, "figS_per_chromosome_methylation", 8.5, 10.5)
fwrite(data.table(sample = samp, condition = cond, beta_cpg_mean = beta_cpg_per),
       file.path(DAT, "per_sample_methylation.tsv"), sep = "\t")

# Step 6 - fig2a: global CpG methylation per sample
cat("[4b] fig2a global methylation per sample\n")
gm <- data.table(sample = factor(samp, levels = samp),
                 condition = factor(cond, levels = c("Control","Amputated")),
                 beta = beta_cpg_per * 100)
pa0 <- ggplot(gm, aes(sample, beta, fill = condition)) +
  geom_col(width = 0.58, alpha = 0.92) +
  geom_text(aes(label = sprintf("%.1f", beta)), vjust = -0.55, size = 3,
            fontface = "bold", colour = "grey15") +
  scale_fill_manual(values = COL_COND, name = NULL) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Mean CpG methylation β (%)",
       title = "Global CpG methylation per sample") +
  theme_pub() +
  theme(legend.position = "top",
        legend.key.size = unit(9, "pt"),
        plot.title = element_text(size = 9.5, face = "bold"),
        plot.margin = margin(5.5, 9, 5.5, 5.5),
        axis.line.x = element_line(linewidth = 0.3, colour = "grey30"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
save_fig(pa0, FIGM, "fig2a_global_methylation_per_sample", 3.1, 2.3)

# Step 7 - fig2b: pooled genome-wide methylation in 1 Mb windows
cat("[5] fig2b genome-wide 1 Mb\n")
cpg <- data.table(chr = chr_vec, pos = start(gr), Mt = M_ctrl + M_amp, Ct = C_ctrl + C_amp)
cpg[, win := (pos %/% 1e6) * 1e6]
agg <- cpg[, .(beta = sum(Mt)/pmax(sum(Ct),1)), by = .(chr, win)]
agg[, chr := factor(chr, levels = keep_chr)]; agg[, end := win + 1e6 - 1]
clen <- agg[, .(maxpos = max(end)), by = chr][order(chr)]
clen[, offset := cumsum(c(0, head(maxpos, -1)))]
agg <- merge(agg, clen[, .(chr, offset)], by = "chr")
agg[, x := (win + end)/2 + offset]; agg[, band := as.integer(chr) %% 2]
chr_mid <- agg[, .(mid = mean(x)), by = chr][order(chr)]
# Mito scaffold labelled 'Mt'; a few odd right-end labels blanked to avoid overlap
chr_lab <- ifelse(chr_mid$chr == "HiC_scaffold_1563", "Mt", sub("chr", "", chr_mid$chr))
chr_lab[chr_lab %in% c("27", "29", "31")] <- ""   # thin the crowded right-end tick labels
pb <- ggplot(agg, aes(x, beta*100, colour = factor(band))) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_hline(yintercept = gmean_global*100, colour = "#C0392B", linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = c("0" = "#2C3E50", "1" = "#7F8C8D"), guide = "none") +
  scale_x_continuous(breaks = chr_mid$mid, labels = chr_lab) +
  labs(x = "Chromosome", y = "Mean CpG methylation β (%)",
       title = "Pooled genome-wide methylation (1 Mb)") +
  theme_pub() + theme(axis.text.x = element_text(size = 6))
save_fig(pb, FIGM, "fig2b_genomewide_1mb", 6.8, 1.96)
fwrite(agg[, .(chr, win_start = as.integer(win), beta_pct = beta*100)],   # integer, never scientific notation
       file.path(DAT, "genomewide_methylation_1mb.tsv"), sep = "\t")

# Step 8 - Region sets: promoter (2 kb upstream of the TSS), exon, intron, genic
exons_gr <- gff[gff$type == "exon"]
seqlengths(gene_gr) <- NA
prom <- trim(suppressWarnings(promoters(gene_gr, 2000, 0)))
# Strand-neutralise before reduce()/setdiff(): GRanges set operations are strand-aware
gu <- gene_gr; strand(gu) <- "*"; eu <- exons_gr; strand(eu) <- "*"; pu <- prom; strand(pu) <- "*"
body <- reduce(gu); exon_r <- reduce(eu)
intron <- GenomicRanges::setdiff(body, exon_r)
genic <- reduce(c(granges(body), granges(pu))); strand(genic) <- "*"; genic <- reduce(genic)
region_gr <- list(Promoter = reduce(granges(pu)), Exon = exon_r, Intron = intron)

# Step 9 - fig2e: mean methylation per region (regions may overlap)
cat("[6] fig2e region mean methylation\n")
setkey(cpg_r, chr, start, end)
region_beta <- function(rgr) {
  dt <- data.table(chr = as.character(seqnames(rgr)), start = start(rgr), end = end(rgr))
  setkey(dt, chr, start, end)
  ov <- foverlaps(cpg_r, dt, nomatch = 0L)
  meth <- sum(ov$bg_c*ov$cv_c + ov$bg_a*ov$cv_a)          # total methylated reads in region
  c(beta = meth / pmax(sum(ov$cv_c + ov$cv_a),1), n = nrow(ov), meth = meth)
}
reg_rows <- rbindlist(lapply(names(region_gr), function(nm) {
  v <- region_beta(region_gr[[nm]])
  data.table(region = nm, mean_beta = v["beta"], n = v["n"], meth = v["meth"]) }))
# Intergenic = CpGs not overlapping any gene body or promoter
genic_dt <- data.table(chr = as.character(seqnames(genic)), start = start(genic), end = end(genic))
setkey(genic_dt, chr, start, end)
ov_g <- foverlaps(cpg_r, genic_dt, nomatch = NA, which = TRUE)
inter_mask <- is.na(ov_g$yid)
ig <- cpg_r[inter_mask]
reg_rows <- rbind(reg_rows, data.table(region = "Intergenic",
  mean_beta = sum(ig$bg_c*ig$cv_c + ig$bg_a*ig$cv_a)/pmax(sum(ig$cv_c+ig$cv_a),1), n = nrow(ig),
  meth = sum(ig$bg_c*ig$cv_c + ig$bg_a*ig$cv_a)))
reg_rows[, region := factor(region, levels = names(COL_REGION))]
fwrite(reg_rows, file.path(DAT, "region_methylation.tsv"), sep = "\t")
pe <- ggplot(reg_rows, aes(region, mean_beta*100, fill = region)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", mean_beta*100)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = COL_REGION, guide = "none") +
  labs(x = NULL, y = "Mean CpG methylation β (%)", title = "Methylation by genomic region") +
  theme_pub()
save_fig(pe, FIGM, "fig2e_region_methylation", 3.7, 2.63)

# Step 10 - fig2j: share of methylated reads per region, exclusive partition
cat("[6b] fig2j region methylation-signal pie\n")
cpg_gr <- GRanges(cpg_r$chr, IRanges(cpg_r$start, width = 1))
# Priority Promoter > Exon > Intron > Intergenic so each CpG is counted once
ra <- rep("Intergenic", length(cpg_gr))
ra[overlapsAny(cpg_gr, body)]                <- "Intron"
ra[overlapsAny(cpg_gr, exon_r)]              <- "Exon"
ra[overlapsAny(cpg_gr, reduce(granges(pu)))] <- "Promoter"
reg_pie <- data.table(region = ra, m = cpg_r$bg_c*cpg_r$cv_c + cpg_r$bg_a*cpg_r$cv_a)[
             , .(meth = sum(m)), by = region][, frac := meth / sum(meth)][]
reg_pie[, region := factor(region, levels = names(COL_REGION))]
pj <- ggplot(reg_pie, aes("", frac, fill = region)) +
  geom_col(width = 1, colour = "white") + coord_polar(theta = "y") +
  geom_text(aes(label = sprintf("%.0f%%", 100*frac)),        # only the number on the slice
            position = position_stack(vjust = 0.5), size = 3.0) +
  geom_point(aes(colour = region), x = 1, y = 0, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = COL_REGION, guide = "none") +
  scale_colour_manual(values = COL_REGION, name = NULL,
                      guide = guide_legend(override.aes = list(alpha = 1, size = 3.5))) +
  labs(title = "Share of methylation signal by region",
       subtitle = "Fraction of total methylated-CpG reads", x = NULL, y = NULL) +
  theme_pub() + theme(axis.line = element_blank(), axis.text = element_blank(),
                      axis.ticks = element_blank(), panel.grid = element_blank(),
                      legend.position = "right")
save_fig(pj, FIGM, "fig2j_region_methylation_pie", 4.2, 4.0)
fwrite(reg_pie, file.path(DAT, "region_methylation_signal.tsv"), sep = "\t")

# Step 11 - Gene-body methylation per gene by tail expression decile
cat("[7] fig2f gene-body methylation by decile\n")
gene_full <- gene_dt[, .(gene_id, chr, start, end, strand, decile)]
setkey(gene_full, chr, start, end)
ov2 <- foverlaps(cpg_r, gene_full, nomatch = 0L)
gb <- ov2[, .(beta = (sum(bg_c*cv_c)+sum(bg_a*cv_a))/pmax(sum(cv_c)+sum(cv_a),1), n_cpg = .N),
          by = .(gene_id, decile)][n_cpg >= 5 & !is.na(decile)]   # per-gene floor: >=5 covered CpGs; expressed genes only
fwrite(gb, file.path(DAT, "genebody_methylation_per_gene.tsv"), sep = "\t")

# Step 12 - Methylation vs expression statistics and the decile-10 split
cat("[7b] methylation vs expression statistics + the decile-10 dip\n")
gbx <- merge(gb, gene_dt[, .(gene_id, expr_mean, length_bp = end - start + 1L)], by = "gene_id")
ct <- cor.test(gbx$beta, gbx$expr_mean, method = "spearman", exact = FALSE)   # STAT TEST: Spearman rank correlation
# The Spearman P underflows at this n; the machine floor is written instead of 0
fwrite(data.table(statistic = c("spearman_rho", "spearman_P", "n_genes"),
                  value = c(unname(ct$estimate), max(ct$p.value, .Machine$double.xmin), nrow(gbx))),   # floor, not a literal 0
       file.path(DAT, "genebody_expression_correlation.tsv"), sep = "\t")
dec_sum <- gbx[, .(n = .N, median_beta = median(beta), q1_beta = quantile(beta, 0.25),
                   q3_beta = quantile(beta, 0.75), pct_beta_below_0.10 = 100 * mean(beta < 0.10)),
               by = decile][order(as.integer(as.character(decile)))]
fwrite(dec_sum, file.path(DAT, "genebody_decile_summary.tsv"), sep = "\t")
cat(sprintf("  Spearman rho(beta, expression) = %.3f (n = %s)\n", ct$estimate, format(nrow(gbx), big.mark = ",")))
print(dec_sum)
# Decile-10 split at beta < 0.10; compare length, CpG count and density, symbols, ribosomal genes, DE
d10 <- gbx[decile == "10"]
d10[, state := ifelse(beta < 0.10, "unmethylated", "methylated")]
note_all <- vapply(mcols(gene_gr)$Note, function(x) if (length(x)) as.character(x)[1] else NA_character_, character(1))
sym_all  <- sub("^Similar to ([^:]+):.*$", "\\1", note_all); sym_all[!grepl("^Similar to [^:]+:", note_all)] <- NA
d10[, symbol := setNames(sym_all, gid_all)[gene_id]]
# Ribosomal protein genes by symbol (^Rp[ls], excluding Rps6k*)
d10[, ribosomal := !is.na(symbol) & grepl("^Rp[ls]", symbol, ignore.case = TRUE) & !grepl("^Rps6k", symbol, ignore.case = TRUE)]
d10[, cpg_per_kb := 1000 * n_cpg / length_bp]
DE_TAIL2 <- file.path(PIPE, "01_genome_toolkit/data/gene_de_tail.tsv")                 # 01_genome_toolkit (upstream)
de2 <- fread(DE_TAIL2)[, .(gene_id, log2FoldChange, padj)]
de2[, de_strict := !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1]   # the paper's DE definition
d10 <- merge(d10, de2, by = "gene_id", all.x = TRUE)
fwrite(d10[order(state, -beta), .(gene_id, symbol, state, beta, n_cpg, length_bp, cpg_per_kb, expr_mean,
                                   ribosomal, log2FoldChange, padj, de_strict)],
       file.path(DAT, "decile10_genes.tsv"), sep = "\t")
u <- d10[state == "unmethylated"]; m <- d10[state == "methylated"]
mw <- function(x, y) { w <- wilcox.test(x, y); list(P = w$p.value, rb = 2 * as.numeric(w$statistic) / (as.numeric(length(x)) * length(y)) - 1) }   # STAT TEST: Mann-Whitney
t_len <- mw(u$length_bp, m$length_bp); t_cpg <- mw(u$n_cpg, m$n_cpg); t_den <- mw(u$cpg_per_kb, m$cpg_per_kb)
named_u <- sum(!is.na(u$symbol)); named_m <- sum(!is.na(m$symbol))
f_rib <- fisher.test(matrix(c(sum(u$ribosomal), named_u - sum(u$ribosomal),
                              sum(m$ribosomal), named_m - sum(m$ribosomal)), 2))   # STAT TEST: Fisher, among named genes
# Strict DE among unmethylated decile-10 genes vs all other genes in the DE table
de_bg <- de2[!gene_id %in% u$gene_id]
f_de <- fisher.test(matrix(c(sum(u$de_strict, na.rm = TRUE), sum(!u$de_strict, na.rm = TRUE),
                             sum(de_bg$de_strict), sum(!de_bg$de_strict)), 2))   # STAT TEST: Fisher
row7 <- function(statistic, scope, value) data.table(statistic = statistic, scope = scope, value = as.numeric(value))
dip <- rbindlist(list(
  row7("n_genes", c("unmethylated", "methylated"), c(nrow(u), nrow(m))),
  row7("median_length_bp", c("unmethylated", "methylated"), c(median(u$length_bp), median(m$length_bp))),
  row7("median_n_cpg", c("unmethylated", "methylated"), c(median(u$n_cpg), median(m$n_cpg))),
  row7("median_cpg_per_kb", c("unmethylated", "methylated"), c(median(u$cpg_per_kb), median(m$cpg_per_kb))),
  row7("mann_whitney_P_length", "unmethylated vs methylated", t_len$P), row7("rank_biserial_length", "unmethylated vs methylated", t_len$rb),
  row7("mann_whitney_P_n_cpg", "unmethylated vs methylated", t_cpg$P), row7("rank_biserial_n_cpg", "unmethylated vs methylated", t_cpg$rb),
  row7("mann_whitney_P_cpg_per_kb", "unmethylated vs methylated", t_den$P), row7("rank_biserial_cpg_per_kb", "unmethylated vs methylated", t_den$rb),
  row7("n_with_symbol", c("unmethylated", "methylated"), c(named_u, named_m)),
  row7("n_ribosomal_protein_genes", c("unmethylated", "methylated"), c(sum(u$ribosomal), sum(m$ribosomal))),
  row7("fisher_P_ribosomal_among_named", "unmethylated vs methylated", f_rib$p.value),
  row7("n_de_strict", c("unmethylated", "methylated"), c(sum(u$de_strict, na.rm = TRUE), sum(m$de_strict, na.rm = TRUE))),
  row7("n_de_strict_up", "unmethylated", sum(u$de_strict & u$log2FoldChange > 0, na.rm = TRUE)),
  row7("pct_de_strict", c("unmethylated", "all other genes in the DE table"),
       100 * c(mean(u$de_strict, na.rm = TRUE), mean(de_bg$de_strict))),
  row7("fisher_OR_de_strict", "unmethylated vs all other genes", unname(f_de$estimate)),
  row7("fisher_P_de_strict", "unmethylated vs all other genes", f_de$p.value)))
fwrite(dip, file.path(DAT, "decile10_unmethylated_vs_methylated.tsv"), sep = "\t")
cat("  decile-10 dip:\n"); print(dip)

# Step 13 - fig2f: per-gene gene-body methylation by expression decile
pf <- ggplot(gb, aes(decile, beta*100, fill = decile)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, outlier.alpha = 0.3, linewidth = 0.25) +
  scale_fill_manual(values = setNames(DECILE_PAL, 1:10), guide = "none") +
  labs(x = "Tail absolute expression decile (1=low, 10=high)", y = "Gene-body β (%)",
       title = "Gene-body methylation by expression decile") + theme_pub()
save_fig(pf, FIGM, "fig2f_genebody_methylation_decile", 4.6, 3.0)

# Step 14 - figS_genebody_metagene_decile: 100-bin gene-body metagene by decile
cat("[8] gene-body metagene by decile (supplementary)\n")
nb <- 100L   # 100 bins along the gene body
setkey(gene_full, chr, start, end)
ovc <- foverlaps(cpg_r, gene_full, nomatch = 0L)
ovc[, rel := (pos - start)/pmax(end - start, 1)]
ovc[strand == "-", rel := 1 - rel]
ovc[, bin := pmin(nb, pmax(1L, as.integer(rel*nb) + 1L))]
ovc[, beta := (bg_c*cv_c + bg_a*cv_a)/pmax(cv_c+cv_a,1)]
mgc <- ovc[!is.na(decile), .(mean_beta = mean(beta)), by = .(decile, bin)]
fwrite(mgc, file.path(DAT, "metagene_genebody_decile.tsv"), sep = "\t")
pc <- ggplot(mgc, aes(bin, mean_beta*100, colour = decile, group = decile)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.5, alpha = 0.55) +
  scale_colour_manual(values = setNames(DECILE_PAL, 1:10), name = "Expr decile",
                      guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Position along gene body (bin, TSS → TES)",
       y = "Mean CpG methylation β (%)",
       title = "Gene-body metagene by expression decile") + theme_pub() +
  theme(legend.position = "right")
save_fig(pc, FIGS, "figS_genebody_metagene_decile", 6.5, 3.2)

# Step 15 - fig2d: TSS +/- 5 kb metagene by decile, tail control, tail blastema, bodywall
cat("[9] fig2d TSS metagene (Tail control / Tail blastema / Bodywall)\n")

# Step 15.1 - HiFi bodywall CpG methylome, cov >= 5 in both samples (cached)
# methbat beds are already strand-combined (one row per CpG dyad); they are not collapsed again
HIFI <- file.path(PIPE, "00_data_pacbio_hifi")
hifi_rds <- file.path(OBJ, "hifi_bodywall_cpg_persample_chrmt.rds")
if (file.exists(hifi_rds)) {
  hifi <- readRDS(hifi_rds)
} else {
  read_hifi <- function(f) {              # bedMethyl-style: col2 is 0-based start
    d <- fread(cmd = sprintf("zcat '%s' | grep -v '^#'", f), header = FALSE,
               select = c(1, 2, 8, 9, 10),
               col.names = c("chr", "start0", "type", "cov", "mod"))
    d <- d[chr %in% keep_chr & type == "Total"]
    d[, .(chr, pos = start0 + 1L, cov, mod)]
  }
  h1 <- read_hifi(file.path(HIFI, "sample1.5mC.bed.gz"))
  h2 <- read_hifi(file.path(HIFI, "sample2.5mC.bed.gz"))
  hifi <- merge(h1, h2, by = c("chr", "pos"), suffixes = c("_1", "_2"))
  hifi <- hifi[cov_1 >= 5 & cov_2 >= 5]   # cov>=5 in both samples
  saveRDS(hifi, hifi_rds)                        # per-sample cov/mod kept in the cache
}
hifi[, `:=`(bg_bw = (mod_1 + mod_2) / (cov_1 + cov_2), cv_bw = cov_1 + cov_2)]
hifi[, `:=`(start = pos, end = pos)]
cat(sprintf("  HiFi bodywall CpGs (cov>=5 in both samples): %s\n",
            format(nrow(hifi), big.mark = ",")))

# Step 15.2 - Bodywall expression deciles from the four control libraries (tail recipe)
bws <- c("dcrep1", "dcrep2", "dcrep3", "dcrep6")
bwl <- lapply(bws, function(s) {
  x <- fread(file.path(HTSEQ, paste0(s, "_htseq_gene_counts.txt")),
             header = FALSE, col.names = c("gene_id", "count"))
  x[!startsWith(gene_id, "__")]
})
stopifnot(sapply(bwl, function(x) identical(x$gene_id, bwl[[1]]$gene_id)))  # same row-order guard as the tail matrix
bwm <- sapply(bwl, function(x) x$count)
rownames(bwm) <- bwl[[1]]$gene_id; colnames(bwm) <- bws
bwm <- bwm[rownames(bwm) %in% gid_chr, , drop = FALSE]   # keep_chr universe before the VST
bwm <- bwm[rowSums(bwm >= 5) >= 2, , drop = FALSE]        # expressed: >= 5 reads in >= 2 libraries
bw_vsd <- vst(DESeqDataSetFromMatrix(bwm, data.frame(s = bws), ~ 1), blind = TRUE)
bw_ex <- data.table(gene_id = rownames(bw_vsd), expr_bw = rowMeans(assay(bw_vsd)))
gene_bw <- merge(gene_dt[, .(gene_id, chr, tss, strand)], bw_ex, by = "gene_id")
gene_bw <- gene_bw[is.finite(expr_bw)]
gene_bw[, decile := cut(expr_bw, quantile(expr_bw, seq(0, 1, 0.1), na.rm = TRUE),
                        labels = 1:10, include.lowest = TRUE)]

# Step 15.3 - Bin the three methylomes on the same TSS +/- 5 kb grid and draw fig2d
flank <- 5000L; nbt <- 100L; bw <- (2*flank)/nbt
tss_bins <- function(cpgs, genes) {       # cpgs: chr/start/end + signal columns
  gwx <- copy(genes); gwx[, `:=`(ws = tss - flank, we = tss + flank)]
  setkey(gwx, chr, ws, we)
  ov <- foverlaps(cpgs, gwx, by.x = c("chr","start","end"),
                  by.y = c("chr","ws","we"), nomatch = 0L)
  ov[, rel := ifelse(strand == "-", tss - pos, pos - tss)]
  ov[, bin := pmin(nbt, pmax(1L, as.integer((rel + flank) %/% bw) + 1L))]
  ov[, bc := -flank + (bin - 0.5)*bw]
  ov[!is.na(decile)]
}
gw  <- gene_dt[, .(gene_id, chr, tss, strand, decile)]
ovt <- tss_bins(cpg_r, gw)
mgt <- ovt[, .(`Tail control`  = sum(bg_c*cv_c)/pmax(sum(cv_c),1),
               `Tail blastema` = sum(bg_a*cv_a)/pmax(sum(cv_a),1)), by = .(decile, bc)]
ovb <- tss_bins(hifi, gene_bw[, .(gene_id, chr, tss, strand, decile)])
mgb <- ovb[, .(condition = "Bodywall",
               beta = sum(bg_bw*cv_bw)/pmax(sum(cv_bw),1)), by = .(decile, bc)]
mgl <- rbind(data.table::melt(mgt, id.vars = c("decile","bc"), variable.name = "condition",   # data.table::melt, namespaced
                  value.name = "beta"),
             mgb[, .(decile, bc, condition, beta)])
mgl[, condition := factor(condition,
                          levels = c("Tail control", "Tail blastema", "Bodywall"))]
fwrite(mgl, file.path(DAT, "metagene_tss5kb_decile.tsv"), sep = "\t")
pd <- ggplot(mgl, aes(bc/1000, beta*100, colour = decile, group = decile)) +
  geom_vline(xintercept = 0, colour = "#C0392B", linetype = "dashed", linewidth = 0.45) +
  geom_line(linewidth = 0.6) + facet_wrap(~ condition, ncol = 3) +
  scale_colour_manual(values = setNames(DECILE_PAL, 1:10), name = "Expr decile",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "Distance from TSS (kb)", y = "Mean β (%)",
       title = "TSS ± 5 kb metagene by expression decile") + theme_pub() +
  theme(legend.position = "right", strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))
save_fig(pd, FIGM, "fig2d_tss5kb_metagene_decile", 7.0, 2.6)

# Step 15.4 - figS_tss5kb_delta_metagene: blastema minus control on the same grid
mgd <- mgt[, .(decile, bc, dbeta = `Tail blastema` - `Tail control`)]
fwrite(mgd, file.path(DAT, "metagene_tss5kb_delta.tsv"), sep = "\t")
p_dlt <- ggplot(mgd, aes(bc/1000, dbeta*100, colour = decile, group = decile)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#C0392B", linetype = "dashed", linewidth = 0.45) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = setNames(DECILE_PAL, 1:10), name = "Expr decile",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "Distance from TSS (kb)", y = "Δβ, blastema − control (%)",
       title = "TSS ± 5 kb methylation change by expression decile") +
  theme_pub() +
  theme(legend.position = "right",
        plot.title = element_text(size = 9, face = "bold"),
        plot.title.position = "plot",
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))
save_fig(p_dlt, FIGS, "figS_tss5kb_delta_metagene", 5.2, 2.8)

# Step 16 - figS mitochondrial and global methylation for the three groups
COL_GRP <- c(`Tail control` = "#2166AC", `Tail blastema` = "#B2182B",
             Bodywall = "#009E73")
# Same estimator everywhere: unweighted mean of per-CpG beta over each platform's cov >= 5 set
grp_of <- ifelse(cond == "Control", "Tail control", "Tail blastema")

# Step 16.1 - figS_mt_methylation_3groups: mitochondrial CpG beta per sample
cat("[9d] figS mitochondrial methylation, 3 groups\n")
MT <- "HiC_scaffold_1563"
mt_idx <- which(chr_vec == MT)
mt_wgbs <- rbindlist(lapply(seq_along(samp), function(i)
  data.table(grp = grp_of[i], sample = samp[i],
             beta = M[mt_idx, i] / Cv[mt_idx, i])))
mt_hifi <- rbind(
  hifi[chr == MT, .(grp = "Bodywall", sample = "BW1", beta = mod_1 / cov_1)],
  hifi[chr == MT, .(grp = "Bodywall", sample = "BW2", beta = mod_2 / cov_2)])
mt3 <- rbind(mt_wgbs, mt_hifi)
mt3[, grp := factor(grp, levels = names(COL_GRP))]
mt3[, .(n_cpg = .N, mean_beta_pct = 100 * mean(beta)), by = .(grp, sample)] |>
  fwrite(file.path(DAT, "mt_methylation_3groups.tsv"), sep = "\t")
p_mt <- ggplot(mt3, aes(sample, beta * 100, colour = grp)) +
  geom_jitter(position = position_jitter(width = 0.18, seed = 20260426), size = 0.8, alpha = 0.55) +   # seeded: identical in pdf/png/svg
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.4,
               colour = "black") +
  facet_grid(~ grp, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = COL_GRP, guide = "none") +
  labs(x = NULL, y = "Mitochondrial CpG methylation β (%)",
       title = "Mitochondrial methylation across the three groups",
       subtitle = paste(sprintf("Mitochondrial scaffold; points = CpG sites (WGBS n=%d, cov>=5 in all 4;", length(mt_idx)),
                        paste("HiFi cov≥5 both slugs);",
                              "black bar = unweighted mean"), sep = "\n")) +
  theme_pub() + theme(strip.background = element_blank(),
                      strip.text = element_text(face = "bold"))
save_fig(p_mt, FIGS, "figS_mt_methylation_3groups", 6.5, 3.2)

# Step 16.2 - figS_global_methylation_3groups: global mean beta per sample
cat("[9e] figS global mean methylation, 3 groups\n")
g3 <- rbind(
  data.table(grp = grp_of, sample = samp, beta = beta_cpg_per),
  data.table(grp = "Bodywall", sample = c("BW1", "BW2"),
             beta = c(hifi[, mean(mod_1 / cov_1)], hifi[, mean(mod_2 / cov_2)])))
g3[, grp := factor(grp, levels = names(COL_GRP))]
fwrite(g3[, .(grp, sample, beta_pct = 100 * beta)],
       file.path(DAT, "global_methylation_3groups.tsv"), sep = "\t")
g3[, lab_hjust := rep(c(1.35, -0.35), length.out = .N), by = grp]
p_g3 <- ggplot(g3, aes(grp, beta * 100, fill = grp)) +
  stat_summary(fun = mean, geom = "col", width = 0.65, alpha = 0.9) +
  geom_point(shape = 21, fill = "white", size = 2.2,
             position = position_jitter(width = 0.06, seed = 1)) +
  geom_text(aes(label = sample, hjust = lab_hjust), size = 2.4,
            position = position_jitter(width = 0.06, seed = 1)) +
  scale_fill_manual(values = COL_GRP, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Mean CpG methylation β (%)",
       title = "Global methylation across the three groups",
       subtitle = paste("Unweighted mean of per CpG β; WGBS cov≥5 in all 4,",
                        "HiFi cov≥5 in both; points = samples", sep = "\n")) +
  theme_pub()
save_fig(p_g3, FIGS, "figS_global_methylation_3groups", 4.6, 3.2)

# Step 17 - HiFi vs WGBS agreement on shared CpGs and 1 Mb windows
cat("[9f] HiFi vs WGBS agreement\n")
wgbs_cpg <- data.table(chr = chr_vec, pos = start(gr),
                       bg = (M_ctrl + M_amp) / pmax(C_ctrl + C_amp, 1))
jn <- merge(wgbs_cpg, hifi[, .(chr, pos, bg_bw)], by = c("chr", "pos"))
jn[, win := (pos %/% 1e6) * 1e6]
w1 <- jn[, .(n_cpg = .N, beta_wgbs = mean(bg), beta_hifi = mean(bg_bw)),
         by = .(chr, win)][n_cpg >= 50]
agree <- data.table(
  n_shared_cpg     = nrow(jn),
  pearson_r_cpg    = cor(jn$bg, jn$bg_bw),
  spearman_rho_cpg = cor(jn$bg, jn$bg_bw, method = "spearman"),
  pearson_r_1mb    = cor(w1$beta_wgbs, w1$beta_hifi),
  n_windows_1mb    = nrow(w1),
  median_offset_1mb = median(w1$beta_hifi - w1$beta_wgbs))
fwrite(agree, file.path(DAT, "hifi_wgbs_agreement.tsv"), sep = "\t")
cat(sprintf("  shared CpGs %s | r(CpG) %.3f | r(1Mb) %.3f | median offset %+.3f\n",
            format(agree$n_shared_cpg, big.mark = ","), agree$pearson_r_cpg,
            agree$pearson_r_1mb, agree$median_offset_1mb))
rm(wgbs_cpg, jn, w1); invisible(gc(FALSE))

# Step 18 - figS_region_decile_metagene_wgbs_tail: discrete-region metagene, WGBS tail
cat("[8b] figS discrete-region metagene by decile, WGBS tail\n")

# Step 18.1 - Longest mRNA per gene
mrna <- as.data.table(gff[gff$type == "mRNA"])
mrna[, mrna_id := sub(";.*", "", as.character(ID))]
mrna[, g_id := vapply(Parent, function(x) if (length(x)) x[1] else NA_character_, character(1))]
mrna[, mlen := end - start + 1]; setorder(mrna, g_id, -mlen)
mk_keep <- mrna[!duplicated(g_id), .(mrna_id, g_id, chr = as.character(seqnames),
                                     mstart = start, mend = end, mstrand = as.character(strand))]

# Step 18.2 - Exons of those mRNAs in TSS-to-TES order, classified first/internal/last
ex <- as.data.table(gff[gff$type == "exon"])
ex[, parent := vapply(Parent, function(x) if (length(x)) x[1] else NA_character_, character(1))]
ex <- ex[parent %in% mk_keep$mrna_id, .(parent, chr = as.character(seqnames), start, end)]
ex <- merge(ex, mk_keep[, .(mrna_id, g_id, mstrand)], by.x = "parent", by.y = "mrna_id")
setorder(ex, g_id, start)
ex[, idx := seq_len(.N), by = g_id]; ex[, n := .N, by = g_id]
ex[mstrand == "-", idx := n - idx + 1L]
# Genes need >= 3 exons so that first, internal and last classes all exist
ex <- ex[n >= 3]                                          # need first/internal/last
ex[, cls := fifelse(idx == 1L, "first_exon", fifelse(idx == n, "last_exon", "internal_exon"))]

# Step 18.3 - Introns as gaps between consecutive exons, classified the same way
intr <- ex[, { o <- order(start); s <- start[o]; e <- end[o]
               if (length(s) >= 2L) list(istart = e[-length(e)] + 1L, iend = s[-1] - 1L,
                                         chr = chr[1], mstrand = mstrand[1]) else NULL },
           by = g_id][iend >= istart]
setorder(intr, g_id, istart)
intr[, idx := seq_len(.N), by = g_id]; intr[, n := .N, by = g_id]
intr[mstrand == "-", idx := n - idx + 1L]
intr[, cls := fifelse(idx == 1L, "first_intron", fifelse(idx == n, "last_intron", "internal_intron"))]

# Step 18.4 - Promoter (2 kb), distal upstream (5 kb) and downstream (5 kb), strand-aware
mu <- mk_keep[g_id %in% unique(ex$g_id)]
mu[, tss := fifelse(mstrand == "-", mend, mstart)]; mu[, tes := fifelse(mstrand == "-", mstart, mend)]
mu[, `:=`(prom_s = fifelse(mstrand=="-", tss+1L, tss-2000L), prom_e = fifelse(mstrand=="-", tss+2000L, tss-1L),
          dist_s = fifelse(mstrand=="-", tss+2001L, tss-7000L), dist_e = fifelse(mstrand=="-", tss+7000L, tss-2001L),
          down_s = fifelse(mstrand=="-", tes-5000L, tes+1L), down_e = fifelse(mstrand=="-", tes-1L, tes+5000L))]
mk <- function(rg, s, e) mu[, .(region = rg, chr, start = pmax(pmin(get(s), get(e)), 1L),
                                end = pmax(get(s), get(e)), strand = mstrand, gene_id = g_id)]

# Nine segments in TSS-to-TES order
regs <- rbindlist(list(
  mk("distal_upstream","dist_s","dist_e"), mk("promoter","prom_s","prom_e"),
  ex[cls=="first_exon",      .(region="first_exon",      chr, start, end, strand=mstrand, gene_id=g_id)],
  intr[cls=="first_intron",  .(region="first_intron",    chr, start=istart, end=iend, strand=mstrand, gene_id=g_id)],
  intr[cls=="internal_intron", .(region="internal_intron", chr, start=istart, end=iend, strand=mstrand, gene_id=g_id)],
  ex[cls=="internal_exon",   .(region="internal_exon",   chr, start, end, strand=mstrand, gene_id=g_id)],
  intr[cls=="last_intron",   .(region="last_intron",     chr, start=istart, end=iend, strand=mstrand, gene_id=g_id)],
  ex[cls=="last_exon",       .(region="last_exon",       chr, start, end, strand=mstrand, gene_id=g_id)],
  mk("downstream","down_s","down_e")))[end >= start]

gene_dec <- unique(gene_full[, .(gene_id, decile)])[!is.na(decile)]
regs <- regs[gene_id %in% gene_dec$gene_id]

# Step 18.5 - Overlap pooled CpGs with each segment, 10 bins per segment, mean per decile
setkey(regs, chr, start, end)
ovr <- foverlaps(cpg_r, regs, by.x = c("chr","start","end"), nomatch = 0L)
nbg <- 10L
ovr[, rel := (pos - start)/pmax(end - start, 1)]; ovr[strand == "-", rel := 1 - rel]
ovr[, bin := pmin(nbg, pmax(1L, as.integer(rel*nbg) + 1L))]
ovr[, beta := (bg_c*cv_c + bg_a*cv_a)/pmax(cv_c+cv_a, 1)]; ovr[, cvt := cv_c + cv_a]
ovr <- merge(ovr, gene_dec, by = "gene_id")
gbin <- ovr[cvt > 0, .(beta = sum(beta*cvt)/sum(cvt)), by = .(gene_id, decile, region, bin)]
aggg <- gbin[, .(beta = mean(beta)), by = .(decile, region, bin)]
fwrite(aggg, file.path(DAT, "metagene_region_decile.tsv"), sep = "\t")

# Step 18.6 - Lay the nine segments on one x axis and draw the figure
rlev <- c("distal_upstream","promoter","first_exon","first_intron","internal_intron",
          "internal_exon","last_intron","last_exon","downstream")
rlab <- c("Distal\nupstream","Promoter\n(2 kb)","First\nexon","First\nintron",
          "Internal\nintron","Internal\nexon","Last\nintron","Last\nexon","Down-\nstream")
rw <- setNames(c(1.3,0.7,0.4,0.9,1.0,0.4,0.9,0.4,1.3), rlev)
rstart <- setNames(head(c(0, cumsum(rw)), -1), rlev); rend <- setNames(cumsum(rw), rlev)
rmid <- (rstart + rend)/2
aggg[, region := factor(region, levels = rlev)]; aggg[, decile := factor(decile, levels = 1:10)]
aggg[, x := rstart[as.character(region)] + (bin - 0.5)/nbg * rw[as.character(region)]]

pg <- ggplot(aggg, aes(x, beta*100, colour = decile, group = interaction(decile, region))) +
  geom_vline(xintercept = rend[-length(rend)], colour = "grey60", linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = rend["promoter"], colour = "#C0392B", linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.35, se = FALSE, linewidth = 0.7) +
  scale_colour_manual(values = setNames(DECILE_PAL, 1:10), name = "Expr decile",
                      guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(breaks = rmid, labels = rlab, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Mean CpG methylation β (%)",
       title = "Methylation across discrete gene regions, Tail (WGBS)",
       subtitle = sprintf("n = %s protein-coding genes (mRNA models with >= 3 exons); tail expression deciles", comma(uniqueN(gbin$gene_id)))) +
  theme_pub(base_size = 11) + theme(panel.grid = element_blank(),
                      axis.text.x = element_text(size = 9), legend.position = "right")
save_fig(pg, FIGS, "figS_region_decile_metagene_wgbs_tail", 9.5, 4.0)

# Step 19 - fig2g: the same discrete-region metagene on HiFi bodywall methylation
cat("[8c] fig2g discrete-region metagene by decile, HiFi bodywall (MAIN)\n")
# Deciles from bodywall expression (Step 15.2); segments and binning reused from Step 18
gene_dec_bw <- unique(gene_bw[!is.na(decile), .(gene_id, decile)])
regs_bw <- regs[gene_id %in% gene_dec_bw$gene_id]
setkey(regs_bw, chr, start, end)
ovr_bw <- foverlaps(hifi, regs_bw, by.x = c("chr", "start", "end"), nomatch = 0L)
ovr_bw[, rel := (pos - start) / pmax(end - start, 1)]
ovr_bw[strand == "-", rel := 1 - rel]
ovr_bw[, bin := pmin(nbg, pmax(1L, as.integer(rel * nbg) + 1L))]
ovr_bw <- merge(ovr_bw, gene_dec_bw, by = "gene_id")
gbin_bw <- ovr_bw[cv_bw > 0, .(beta = sum(bg_bw * cv_bw) / sum(cv_bw)),
                  by = .(gene_id, decile, region, bin)]
agg_bw <- gbin_bw[, .(beta = mean(beta)), by = .(decile, region, bin)]
fwrite(agg_bw, file.path(DAT, "metagene_region_decile_bodywall.tsv"), sep = "\t")
agg_bw[, region := factor(region, levels = rlev)]
agg_bw[, decile := factor(decile, levels = 1:10)]
agg_bw[, x := rstart[as.character(region)] + (bin - 0.5) / nbg * rw[as.character(region)]]
pg_bw <- ggplot(agg_bw, aes(x, beta * 100, colour = decile,
                            group = interaction(decile, region))) +
  geom_vline(xintercept = rend[-length(rend)], colour = "grey60",
             linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = rend["promoter"], colour = "#C0392B",
             linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.35, se = FALSE, linewidth = 0.7) +
  scale_colour_manual(values = setNames(DECILE_PAL, 1:10), name = "Expr decile",
                      guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(breaks = rmid, labels = rlab,
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Mean CpG methylation β (%)",
       title = "Methylation across discrete gene regions, bodywall",
       subtitle = sprintf(
         "n = %s protein-coding genes (mRNA models with >= 3 exons); bodywall expression deciles",
         comma(uniqueN(gbin_bw$gene_id)))) +
  theme_pub() + theme(panel.grid = element_blank(),
                      axis.text.x = element_text(size = 8),
                      legend.position = "right")
save_fig(pg_bw, FIGM, "fig2g_region_decile_metagene_bodywall", 7.0, 2.7)

# Step 20 - fig2i: sample PCA and correlation on 1 Mb-window beta
cat("[10] fig2i sample PCA (supp)\n")
wmat <- dcast(per, chr + win ~ sample, value.var = "beta")
wmat <- wmat[complete.cases(wmat)]
pmat <- as.matrix(wmat[, ..samp]); pmat <- pmat[apply(pmat, 1, sd) > 0, , drop = FALSE]
pcd <- prcomp(t(pmat), scale. = TRUE); ve <- round(100*pcd$sdev^2/sum(pcd$sdev^2), 1)
pca_dt <- data.table(sample = samp, condition = cond, PC1 = pcd$x[,1], PC2 = pcd$x[,2])
p_pca <- ggplot(pca_dt, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 4) + geom_text(aes(label = sample), vjust = -1.1, size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = COL_COND) +
  labs(title = "Sample PCA (1 Mb-window β)", x = sprintf("PC1 (%.1f%%)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_pub()
cmat <- cor(pmat); cm_dt <- as.data.table(as.table(cmat)); setnames(cm_dt, c("s1","s2","r"))
p_cor <- ggplot(cm_dt, aes(s1, s2, fill = r)) + geom_tile() +
  geom_text(aes(label = sprintf("%.3f", r)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = mean(cmat[lower.tri(cmat)])) +
  labs(title = "Sample β correlation", x = NULL, y = NULL) + theme_pub()
save_fig(p_pca | p_cor, FIGS, "fig2i_sample_pca_correlation", 10, 4.2)
fwrite(pca_dt, file.path(DAT, "sample_pca_coords.tsv"), sep = "\t")

# Step 21 - Record the package versions used by this run
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_02_landscape.txt"))
cat("[02_landscape] done\n")
