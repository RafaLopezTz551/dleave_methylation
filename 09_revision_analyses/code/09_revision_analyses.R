#!/usr/bin/env Rscript
# ==============================================================================
# 09_revision_analyses — analyses added during the co-author review (October 2026)
# ------------------------------------------------------------------------------
# QUESTION  Figure and statistics requests from the internal review of the
#           manuscript: regenerative capacity next to CpG O/E (Fig. 1B), E-values
#           in the toolkit panel (Fig. 1E), per-feature violins with tests
#           (Fig. 2B), whether high methylation calls exceed the bisulfite
#           non-conversion rate, whether the lower HiFi upper mode is a platform
#           or a tissue effect, a rolling-mean genome-wide profile, and the
#           number of TE copies plotted per class.
#
# READS     02_landscape/objects/bsseq_cov5_chrmt.rds, hifi_bodywall_cpg_persample_chrmt.rds
#           02_landscape/data/genomewide_methylation_1mb{,_bodywall}.tsv
#           01_genome_toolkit/data/mollusc_dinucleotide_oe.tsv, toolkit_presence.tsv
#           01_genome_toolkit/objects/gff_chrmt.rds; eggNOG-mapper annotations
#           04_TEs/data/non_conversion_rate.tsv, te_methylation_per_copy{,_hifi}.tsv
#           06_decoupling/data/decoupling_per_gene.tsv; HTSeq bodywall counts
# WRITES    data/*.tsv, figures/{main,supplementary}/*.{svg,png}; sessionInfo_09_revision_analyses.txt
#
# RULES     fixed seed 20260426; chr1-31 only; reads only lower-numbered modules;
#           figures carry no title (titles are in the LaTeX captions).
# ==============================================================================

# ---- [0] Setup ----------------------------------------------------------------
set.seed(20260426)
suppressPackageStartupMessages({ library(data.table); library(bsseq); library(GenomicRanges); library(DESeq2); library(ggplot2); library(svglite) })
PIPE  <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
HTSEQ <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"
BATCH <- Sys.getenv("BATCH_DIR", file.path(PIPE, "09_revision_analyses"))
OUT   <- file.path(BATCH, "data"); FIGM <- file.path(BATCH, "figures/main"); FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(OUT, FIGM, FIGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
FONT <- "Liberation Sans"; chrs <- paste0("chr", 1:31)
theme_t <- function(bs = 8) theme_classic(base_size = bs, base_family = FONT) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(), axis.text = element_text(colour = "black"))
sv <- function(p, name, w, h, dir = FIGM) { ggsave(file.path(dir, paste0(name, ".svg")), p, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200); cat("saved", name, "\n") }
summ <- list(); add <- function(k, v) { summ[[length(summ) + 1]] <<- data.table(item = k, value = as.character(v)); cat(k, "=", v, "\n") }

# ---- [1] Fig. 1B: genome-wide CpG O/E across five mollusks, with regenerative capacity ------------------------------
oe <- fread(file.path(PIPE, "01_genome_toolkit/data/mollusc_dinucleotide_oe.tsv"))[dinucleotide == "CG"]
regen <- data.table(species = c("Deroceras laeve","Pomacea canaliculata","Elysia atroviridis","Aplysia californica","Octopus bimaculoides"),
  regen = c("tail","eye","whole body\nfrom head","neurons","arm"))
oe <- merge(oe, regen, by = "species"); setorder(oe, ratio_obs_exp)
oe[, lab := sub("^([A-Z])[a-z]+ ", "\\1. ", species)]; oe[, lab := factor(lab, levels = lab)]
for (i in seq_len(nrow(oe))) add(paste0("CpG O/E ", oe$species[i]), sprintf("%.3f (GC %.1f%%)", oe$ratio_obs_exp[i], 100 * oe$gc[i]))
p1 <- ggplot(oe, aes(lab, ratio_obs_exp, fill = group)) + geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", ratio_obs_exp)), vjust = -0.4, size = 2.6, family = FONT) +
  geom_text(aes(y = 0.04, label = sprintf("%.0f%%", 100 * gc)), size = 2.0, colour = "white", family = FONT) +
  geom_text(aes(y = -0.10, label = regen), size = 2.1, lineheight = 0.85, colour = "grey25", family = FONT, vjust = 1) +
  scale_fill_manual(values = c(Gastropoda = "#4C72B0", Cephalopoda = "#DD8452"), name = NULL) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) + coord_cartesian(ylim = c(-0.24, 1.08), clip = "off") +
  labs(x = NULL, y = "Genome-wide CpG O/E") + theme_t() +
  theme(axis.text.x = element_text(face = "italic", angle = 20, hjust = 1, vjust = 1.1, margin = margin(t = 22)),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.key.size = unit(3, "mm"))
sv(p1, "fig1e_mollusc_cpg_oe", 3.4, 2.5)

# ---- [2] Fig. 1E: methylation toolkit, absent genes in grey, eggNOG E-value beside the bit score ------------
tk <- fread(file.path(PIPE, "01_genome_toolkit/data/toolkit_presence.tsv"))
em <- fread(cmd = "grep -v '^##' /mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/eggnog_mapper/dlasi_proteome.emapper.annotations", select = 1:4, sep = "\t", quote = "")
setnames(em, c("query","seed","evalue","score")); em[, gene_id := sub("-mRNA-.*|\\.t[0-9]+$|_t[0-9]+$", "", query)]
emg <- em[, .SD[which.max(score)], by = gene_id][, .(gene_id, evalue, score_em = score)]
tk <- merge(tk, emg, by = "gene_id", all.x = TRUE, sort = FALSE)
add("toolkit rows with E-value matched", sprintf("%d of %d present", tk[present == TRUE & !is.na(evalue), .N], tk[present == TRUE, .N]))
fwrite(tk, file.path(OUT, "toolkit_presence_with_evalue.tsv"), sep = "\t")
bl <- c("Writer","Eraser","Reader","Support","Cofactor","Demethylation pathway")
tk[, cat_broad := factor(sub(" *\\(.*", "", category), levels = bl)]
tk[, label := ifelse(present, ortholog, gene_symbol)]; setorder(tk, -present, -bit, gene_symbol)
tk[, label := factor(make.unique(label), levels = rev(make.unique(label)))]
fmtE <- function(e) ifelse(is.na(e), "", ifelse(e == 0, "E < 1e-300", sprintf("E = %s", formatC(e, format = "e", digits = 0))))
tk[, txt := ifelse(present, paste0(format(round(bit), big.mark = ",", trim = TRUE), "  (", fmtE(evalue), ")"), "absent")]
COL <- c(Writer = "#1B9E77", Eraser = "#D95F02", Reader = "#7570B3", Support = "#E7298A", Cofactor = "#66A61E", `Demethylation pathway` = "#E6AB02")
p2 <- ggplot(tk, aes(bit, label, fill = cat_broad)) + geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_text(aes(x = -60, label = as.character(label), colour = present), hjust = 1, size = 2.3, family = FONT) +
  geom_text(aes(x = ifelse(present, bit, 0) + 60, label = txt, colour = present), hjust = 0, size = 2.3, family = FONT) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"), guide = "none") +
  facet_grid(cat_broad ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_fill_manual(values = COL, guide = "none") +
  scale_x_continuous(limits = c(-1150, 4300), breaks = c(0, 1000, 2000, 3000), expand = c(0, 0)) +
  labs(x = "eggNOG ortholog bit score", y = NULL) + theme_t() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.5, colour = "grey25"))
sv(p2, "fig2a_methylation_toolkit_presence", 4.4, 5.1)

# ---- [3] Fig. 2B: per-feature methylation by genomic region (violins, Kruskal-Wallis, pairwise Wilcoxon) ----------------------------
bs <- readRDS(file.path(PIPE, "02_landscape/objects/bsseq_cov5_chrmt.rds"))
chrs <- paste0("chr", 1:31); bs <- bs[as.character(seqnames(bs)) %in% chrs]
M <- as.matrix(getCoverage(bs, type = "M")); Cv <- as.matrix(getCoverage(bs, type = "Cov")); sn <- colnames(bs)
ic <- grep("^C", sn); ia <- grep("^A", sn)
cg <- data.table(chr = as.character(seqnames(bs)), pos = start(bs), m = rowSums(M), cov = rowSums(Cv),
                 mc = rowSums(M[, ic]), cc = rowSums(Cv[, ic]), ma = rowSums(M[, ia]), ca = rowSums(Cv[, ia]))
Mi <- M; Ci <- Cv; rm(bs); gc()
gff <- readRDS(file.path(PIPE, "01_genome_toolkit/objects/gff_chrmt.rds")); gff <- gff[as.character(seqnames(gff)) %in% chrs]
gene <- gff[gff$type == "gene"]; ex <- gff[gff$type == "exon"]; strand(ex) <- "*"
prom <- trim(suppressWarnings(promoters(gene, 2000, 0))); strand(prom) <- "*"
gu <- gene; strand(gu) <- "*"; body <- reduce(gu); exr <- reduce(ex); intr <- GenomicRanges::setdiff(body, exr)
genic <- reduce(c(granges(body), granges(prom))); sl <- cg[, .(len = max(pos) + 1L), by = chr]
gen <- GRanges(sl$chr, IRanges(1, sl$len)); inter <- GenomicRanges::setdiff(gen, genic)
inter <- unlist(tile(inter, width = 2000)); inter <- inter[width(inter) >= 500]
feat <- list(Promoter = reduce(prom), Exon = exr, Intron = intr, Intergenic = inter)
gcg <- GRanges(cg$chr, IRanges(cg$pos, width = 1))
fb <- rbindlist(lapply(names(feat), function(nm) { f <- feat[[nm]]; ov <- findOverlaps(gcg, f)
  d <- data.table(fid = subjectHits(ov), m = cg$m[queryHits(ov)], cov = cg$cov[queryHits(ov)])
  d <- d[, .(n_cpg = .N, beta = sum(m) / sum(cov)), by = fid][n_cpg >= 3]; d[, region := nm]; d }))
fb[, region := factor(region, levels = names(feat))]
kw <- kruskal.test(beta ~ region, data = fb); pw <- pairwise.wilcox.test(fb$beta, fb$region, p.adjust.method = "BH")
rs <- fb[, .(n_features = .N, median_beta = median(beta), iqr_lo = quantile(beta, .25), iqr_hi = quantile(beta, .75), mean_beta = mean(beta),
             pct_gt_0.5 = 100 * mean(beta > 0.5)), by = region][order(region)]
print(rs); fwrite(rs, file.path(OUT, "region_feature_beta_summary.tsv"), sep = "\t")
pwd <- as.data.table(as.table(pw$p.value))[!is.na(N)]; setnames(pwd, c("a","b","padj")); fwrite(pwd, file.path(OUT, "region_pairwise_wilcoxon.tsv"), sep = "\t")
add("region Kruskal-Wallis", sprintf("H = %.0f, df = %d, P = %.2g", kw$statistic, kw$parameter, kw$p.value))
COLR <- c(Promoter = "#2E86C1", Exon = "#17A589", Intron = "#7D3C98", Intergenic = "#85C1E9")
p3 <- ggplot(fb, aes(region, 100 * beta, fill = region)) + geom_violin(scale = "width", linewidth = 0.2, colour = "grey30", adjust = 1.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "white") +
  geom_text(data = rs, aes(region, 104, label = sprintf("n = %s", format(n_features, big.mark = ","))), size = 2.1, family = FONT, inherit.aes = FALSE) +
  scale_fill_manual(values = COLR, guide = "none") + scale_y_continuous(breaks = seq(0, 100, 25)) + coord_cartesian(ylim = c(0, 106)) +
  labs(x = NULL, y = "CpG methylation per feature (%)") + theme_t()
sv(p3, "fig2e_region_methylation", 3.7, 2.63)

# ---- [4] High-methylation calls vs the bisulfite non-conversion rate (one-sided binomial, BH) ---
nc <- fread(file.path(PIPE, "04_TEs/data/non_conversion_rate.tsv")); p0 <- max(nc$chh_pct) / 100
add("non-conversion rate used (max CHH across samples)", sprintf("%.3f%%", 100 * p0))
cg[, beta := m / cov]
bt <- function(x, n, p) pbinom(x - 1, n, p, lower.tail = FALSE)
cg[, p_nc := bt(m, cov, p0)]; cg[, p_nc1 := bt(m, cov, 0.01)]
cg[, fdr := p.adjust(p_nc, "BH")]; cg[, fdr1 := p.adjust(p_nc1, "BH")]
# per-sample: methylated in every sample?
pers <- sapply(seq_len(ncol(Mi)), function(j) p.adjust(bt(Mi[, j], Ci[, j], p0), "BH") < 0.05)
cg[, sig_all4 := rowSums(pers) == ncol(Mi)]
hb <- cg[beta >= 0.8]
fp <- data.table(set = c("beta >= 0.8","beta 0.5-0.8","beta 0.2-0.5","all"),
  n = c(nrow(hb), cg[beta >= 0.5 & beta < 0.8, .N], cg[beta >= 0.2 & beta < 0.5, .N], nrow(cg)),
  pct_FDR05_pooled = c(100 * mean(hb$fdr < 0.05), 100 * cg[beta >= 0.5 & beta < 0.8, mean(fdr < 0.05)], 100 * cg[beta >= 0.2 & beta < 0.5, mean(fdr < 0.05)], 100 * mean(cg$fdr < 0.05)),
  pct_FDR05_pooled_p1pct = c(100 * mean(hb$fdr1 < 0.05), 100 * cg[beta >= 0.5 & beta < 0.8, mean(fdr1 < 0.05)], 100 * cg[beta >= 0.2 & beta < 0.5, mean(fdr1 < 0.05)], 100 * mean(cg$fdr1 < 0.05)),
  pct_FDR05_each_sample = c(100 * mean(hb$sig_all4), 100 * cg[beta >= 0.5 & beta < 0.8, mean(sig_all4)], 100 * cg[beta >= 0.2 & beta < 0.5, mean(sig_all4)], 100 * mean(cg$sig_all4)),
  median_cov_pooled = c(median(hb$cov), cg[beta >= 0.5 & beta < 0.8, median(cov)], cg[beta >= 0.2 & beta < 0.5, median(cov)], median(cg$cov)))
q99 <- quantile(cg$cov, 0.99); fp[, pct_cov_above_q99 := c(100 * mean(hb$cov > q99), 100 * cg[beta >= 0.5 & beta < 0.8, mean(cov > q99)], 100 * cg[beta >= 0.2 & beta < 0.5, mean(cov > q99)], 1)]
print(fp); fwrite(fp, file.path(OUT, "high_methylation_false_positive_check.tsv"), sep = "\t")
# in-sample consistency of high CpGs: SD of beta across the four libraries
cg[, sd4 := apply(Mi / pmax(Ci, 1), 1, sd)]
add("high CpGs (beta>=0.8): median SD of beta across 4 libraries", sprintf("%.3f", median(cg[beta >= 0.8, sd4])))

# ---- [5] Upper mode of HiFi vs WGBS at CpGs shared by both platforms -------------------------------
h <- readRDS(file.path(PIPE, "02_landscape/objects/hifi_bodywall_cpg_persample_chrmt.rds"))[chr %in% chrs & cov_1 >= 5 & cov_2 >= 5]
h[, beta_h := (mod_1 + mod_2) / (cov_1 + cov_2)]; h[, cov_h := cov_1 + cov_2]
best <- 0L; bn <- -1; for (o in c(-1L, 0L, 1L)) { n <- nrow(merge(cg[, .(chr, pos)], h[, .(chr, pos = pos + o)], by = c("chr","pos"))); if (n > bn) { bn <- n; best <- o } }
h[, pos := pos + best]; sh <- merge(cg[, .(chr, pos, bc = mc / cc, ba = ma / ca, cov)], h[, .(chr, pos, beta_h, cov_h)], by = c("chr","pos"))
add("CpGs shared WGBS and HiFi", nrow(sh))
umode <- function(x) { x <- x[x >= 0.5]; d <- density(x, from = 0.5, to = 1, n = 512, adjust = 1.5); d$x[which.max(d$y)] }
um <- data.table(set = c("WGBS tail control","WGBS amputated tail","HiFi body wall","HiFi body wall, cov >= 30"),
  upper_mode = c(umode(sh$bc), umode(sh$ba), umode(sh$beta_h), umode(sh[cov_h >= 30, beta_h])),
  pct_ge_0.8 = c(100 * mean(sh$bc >= .8), 100 * mean(sh$ba >= .8), 100 * mean(sh$beta_h >= .8), 100 * sh[cov_h >= 30, mean(beta_h >= .8)]))
hi <- sh[bc >= 0.8]
um2 <- data.table(set = "CpGs with WGBS control beta >= 0.8", n = nrow(hi), median_wgbs = median(hi$bc), median_hifi = median(hi$beta_h),
  median_diff = median(hi$bc - hi$beta_h), median_hifi_cov30 = hi[cov_h >= 30, median(beta_h)])
# genes expressed in both tissues: restrict to gene bodies with high methylation in both platforms is circular; instead compare at mito-free,
# platform-agnostic CpGs: HiFi upper mode among CpGs fully methylated by WGBS in BOTH conditions (beta >= 0.95)
um3 <- sh[bc >= 0.95 & ba >= 0.95, .(n = .N, median_hifi = median(beta_h), pct_hifi_ge_0.9 = 100 * mean(beta_h >= 0.9))]
print(um); print(um2); print(um3)
fwrite(um, file.path(OUT, "upper_mode_by_platform.tsv"), sep = "\t"); fwrite(um2, file.path(OUT, "hifi_at_wgbs_high_cpgs.tsv"), sep = "\t"); fwrite(um3, file.path(OUT, "hifi_at_wgbs_full_cpgs.tsv"), sep = "\t")

# ---- [6] Platform vs tissue: HiFi values by coverage, and the WGBS - HiFi gap by gene
# [6a] coverage dependence of HiFi values at CpGs fully methylated in WGBS
fu <- sh[bc >= 0.95 & ba >= 0.95]
fu[, covbin := cut(cov_h, c(9, 19, 29, 39, 49, Inf), labels = c("10-19","20-29","30-39","40-49",">=50"))]
cvt <- fu[!is.na(covbin), .(n = .N, median_hifi = median(beta_h), mean_hifi = mean(beta_h), pct_ge_0.9 = 100 * mean(beta_h >= 0.9)), by = covbin][order(covbin)]
print(cvt); fwrite(cvt, file.path(OUT, "hifi_at_wgbs_full_by_coverage.tsv"), sep = "\t")

# [6b] gene level: does the WGBS-HiFi gap track tissue-specific expression?
gff <- readRDS(file.path(PIPE, "01_genome_toolkit/objects/gff_chrmt.rds")); gene <- gff[gff$type == "gene" & as.character(seqnames(gff)) %in% chrs]
gene <- gene[as.character(mcols(gene)$gene_biotype) %in% c("protein_coding", "lncRNA")]
gdt <- data.table(chr = as.character(seqnames(gene)), start = start(gene), end = end(gene), gene_id = sub(";.*", "", as.character(mcols(gene)$ID)))
sh[, `:=`(start = pos, end = pos)]; setkey(gdt, chr, start, end)
ov <- foverlaps(sh, gdt, by.x = c("chr","start","end"), nomatch = 0L)
pg <- ov[, .(n = .N, w = mean(bc), hb = mean(beta_h), w_hi = .SD[bc >= 0.8, mean(bc)], h_hi = .SD[bc >= 0.8, mean(beta_h)], n_hi = sum(bc >= 0.8)), by = gene_id][n >= 10]
tail <- fread(file.path(PIPE, "06_decoupling/data/decoupling_per_gene.tsv"), select = c("gene_id","expr_mean"))
bws <- c("dcrep1","dcrep2","dcrep3","dcrep6")
bwl <- lapply(bws, function(s) { x <- fread(file.path(HTSEQ, paste0(s, "_htseq_gene_counts.txt")), header = FALSE, col.names = c("gene_id","count")); x[!startsWith(gene_id, "__")] })
bwm <- sapply(bwl, function(x) x$count); rownames(bwm) <- bwl[[1]]$gene_id; colnames(bwm) <- bws
bwm <- bwm[rownames(bwm) %in% gdt$gene_id, , drop = FALSE]; bwm <- bwm[rowSums(bwm >= 5) >= 2, , drop = FALSE]
bwv <- vst(DESeqDataSetFromMatrix(bwm, data.frame(s = bws), ~ 1), blind = TRUE)
bwe <- data.table(gene_id = rownames(bwv), expr_bw = rowMeans(assay(bwv)))
gx <- merge(merge(pg, tail, by = "gene_id"), bwe, by = "gene_id")
gx[, `:=`(r_tail = frank(expr_mean) / .N, r_bw = frank(expr_bw) / .N)]; gx[, disc := abs(r_tail - r_bw)]
gx[, gap := w - hb]; gx[, gap_hi := w_hi - h_hi]
gm <- gx[w >= 0.5]
st <- cor.test(gm$gap, gm$disc, method = "spearman", exact = FALSE)
add("methylated genes (WGBS >= 50%), n", nrow(gm))
add("median gene gap WGBS - HiFi (pp)", sprintf("%.1f", 100 * median(gm$gap)))
add("Spearman gap vs expression discordance", sprintf("rho = %.3f, P = %.2g", st$estimate, st$p.value))
gm[, grp := fifelse(disc <= 0.1 & r_tail >= 0.5 & r_bw >= 0.5, "concordant, expressed in both", fifelse(disc >= 0.4, "discordant (tissue-specific)", "other"))]
gt <- gm[, .(n = .N, median_gap_pp = 100 * median(gap), median_gap_hiCpG_pp = 100 * median(gap_hi, na.rm = TRUE), median_wgbs = 100 * median(w), median_hifi = 100 * median(hb)), by = grp]
wt <- wilcox.test(gm[grp == "concordant, expressed in both", gap], gm[grp == "discordant (tissue-specific)", gap])
print(gt); fwrite(gt, file.path(OUT, "gap_by_expression_concordance.tsv"), sep = "\t")
add("Wilcoxon gap concordant vs discordant", sprintf("P = %.2g", wt$p.value))
# share of gene-to-gene variance of the gap explained by tissue expression discordance
fit <- lm(gap ~ disc + I(r_tail - r_bw), data = gm); add("R2 gap ~ expression discordance", sprintf("%.4f", summary(fit)$r.squared))
fwrite(gx, file.path(OUT, "gene_platform_gap.tsv"), sep = "\t")

# ---- [7] Fig. S (genome-wide profile), panel B: rolling mean over five 1 Mb bins
wg <- fread(file.path(PIPE, "02_landscape/data/genomewide_methylation_1mb.tsv"))[chr %in% chrs]
bw <- fread(file.path(PIPE, "02_landscape/data/genomewide_methylation_1mb_bodywall.tsv"))[chr %in% chrs]
d <- merge(wg[, .(chr, win_start, control = beta_pct_control, amputated = beta_pct_blastema)], bw[, .(chr, win_start, bodywall = beta_pct)], by = c("chr","win_start"), all = TRUE)
d[, chr := factor(chr, levels = chrs)]; setorder(d, chr, win_start)
for (v in c("control","amputated","bodywall")) d[, (paste0(v, "_rm")) := frollmean(get(v), 5, align = "center", na.rm = TRUE), by = chr]
off <- d[, .(len = max(win_start) + 1e6), by = chr][, off := data.table::shift(cumsum(len), fill = 0)]
d <- merge(d, off[, .(chr, off)], by = "chr"); d[, x := (win_start + off) / 1e6]
ml <- melt(d, id.vars = c("chr","x"), measure.vars = c("control_rm","amputated_rm","bodywall_rm"), variable.name = "set", value.name = "b")
ml[, set := factor(set, levels = c("control_rm","amputated_rm","bodywall_rm"), labels = c("Tail control (WGBS)","Amputated tail (WGBS)","Body wall (PacBio HiFi)"))]
mids <- d[, .(m = mean(range(x))), by = chr]; shade <- d[, .(x0 = min(x), x1 = max(x) + 1), by = chr][as.integer(chr) %% 2 == 0]
p <- ggplot(ml, aes(x, b, colour = set, group = interaction(chr, set))) + geom_rect(data = shade, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey94") +
  geom_line(linewidth = 0.3) + geom_hline(yintercept = 18.1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  scale_colour_manual(values = c("#2166AC","#B2182B","#1B7837"), name = NULL) + scale_x_continuous(breaks = mids$m, labels = sub("chr", "", mids$chr), expand = c(0.005, 0)) +
  labs(x = "Chromosome", y = "Mean CpG methylation (%)") + theme_t(7) + theme(legend.position = "top", axis.text.x = element_text(size = 5.5))
sv(p, "figS_genomewide_1mb_rolling", 7.0, 2.6, dir = FIGS)
cr <- d[, .(r_ctrl_bw = cor(control, bodywall, use = "complete"), r_ctrl_amp = cor(control, amputated, use = "complete"))]
add("1 Mb Pearson control vs body wall", sprintf("%.3f", cr$r_ctrl_bw)); add("1 Mb Pearson control vs amputated", sprintf("%.3f", cr$r_ctrl_amp))

# ---- [8] TE copies per class and location (Fig. 6 legend) ------------------------------------------------------
th <- fread(file.path(PIPE, "04_TEs/data/te_methylation_per_copy_hifi.tsv"), select = c("class","loc"))
tw <- fread(file.path(PIPE, "04_TEs/data/te_methylation_per_copy.tsv"), select = c("class","loc"))
tc <- merge(th[, .(hifi_bodywall = .N), by = .(class, loc)], tw[, .(wgbs_tail = .N), by = .(class, loc)], by = c("class","loc"), all = TRUE)
tc <- rbind(tc, tc[, .(loc = "All", hifi_bodywall = sum(hifi_bodywall, na.rm = TRUE), wgbs_tail = sum(wgbs_tail, na.rm = TRUE)), by = class])
print(tc[order(class, loc)]); fwrite(tc[order(class, loc)], file.path(OUT, "te_copies_per_class_plotted.tsv"), sep = "\t")

fwrite(rbindlist(summ), file.path(OUT, "revision_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_09_revision_analyses.txt"))
cat("[09_revision_analyses] done\n")
