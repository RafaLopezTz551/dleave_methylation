#!/usr/bin/env Rscript

set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges)
  library(bsseq); library(ggplot2); library(ggridges)
})

PIPE <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
TE   <- "/mnt/data/alfredvar/30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv"
B01 <- file.path(PIPE, "01_genome_toolkit/objects"); B02 <- file.path(PIPE, "02_landscape/objects")
BATCH <- file.path(PIPE, "04_TEs")
DAT <- file.path(BATCH, "data"); FIGM <- file.path(BATCH, "figures/main")
FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(DAT, FIGM, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")  # chr1-31 + mito scaffold (2026-08-17)

COL_CLASS <- c(LTR = "#2471A3", RC = "#F39C12", LINE = "#8E44AD",
               DNA = "#1ABC9C", SINE = "#D81B60")
COL_LOC   <- c(`Gene-overlapping` = "#2CA25F", Intergenic = "#E78A8A")
theme_pub <- function() theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
        strip.background = element_blank())   # no boxes around facet labels (rule)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGM, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGM, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGM, paste0(name, ".svg")), p, width = w, height = h)   # vector (svg)
  cat(sprintf("  saved %s\n", name))
}
save_supp <- function(p, name, w, h) {   # mirrors save_fig(), writes to supplementary
  ggsave(file.path(FIGS, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGS, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGS, paste0(name, ".svg")), p, width = w, height = h)
  cat(sprintf("  saved %s (supplementary)\n", name))
}

cat("[1] TE table + bsseq + gff\n")
te <- fread(TE)
te <- te[chrom %in% keep_chr]

gen_len <- sum(width(readRDS(file.path(B01, "genome_chrmt.rds"))))
rep_bp  <- sum(width(reduce(GRanges(te$chrom, IRanges(te$start, te$end)))))   # same 1-based convention as te_gr below
fwrite(data.table(genome_bp = gen_len, merged_repeat_bp = rep_bp,
                  repeat_fraction = rep_bp / gen_len),
       file.path(DAT, "te_genome_fraction.tsv"), sep = "\t")
cat(sprintf("  repeats cover %.1f%% of the chr1-31+mt assembly (%d of %.0f bp, merged)\n",
            100 * rep_bp / gen_len, rep_bp, as.numeric(gen_len)))
cf <- as.character(te$class_family)
te[, class := ifelse(grepl("^LINE", cf), "LINE", ifelse(grepl("^SINE", cf), "SINE",
              ifelse(grepl("^LTR", cf), "LTR", ifelse(grepl("^DNA", cf), "DNA",
              ifelse(grepl("^RC", cf), "RC", "Unknown")))))]
te <- te[class %in% names(COL_CLASS)]
te[, te_uid := .I]
te_gr <- GRanges(te$chrom, IRanges(te$start, te$end), te_uid = te$te_uid)

bs <- readRDS(file.path(B02, "bsseq_cov5_chrmt.rds"))
chrs <- as.character(seqnames(bs)); bs <- bs[chrs %in% keep_chr, ]
gr <- granges(bs)
M <- as.matrix(getCoverage(bs, type = "M")); Cv <- as.matrix(getCoverage(bs, type = "Cov"))
samp <- sampleNames(bs); cond <- ifelse(grepl("^A", samp), "Amputated", "Control")
mc <- rowSums(M[, cond == "Control", drop = FALSE]); cc <- rowSums(Cv[, cond == "Control", drop = FALSE])
ma <- rowSums(M[, cond == "Amputated", drop = FALSE]); ca <- rowSums(Cv[, cond == "Amputated", drop = FALSE])
gff <- readRDS(file.path(B01, "gff_chrmt.rds"))
gene_gr <- reduce(gff[gff$type == "gene"])

cat("[2] per-TE methylation\n")
cpg_dt <- data.table(chr = as.character(seqnames(gr)), start = start(gr), end = start(gr),
                     mc = mc, cc = cc, ma = ma, ca = ca)
te_dt <- data.table(chr = te$chrom, start = te$start, end = te$end, te_uid = te$te_uid)
setkey(te_dt, chr, start, end); setkey(cpg_dt, chr, start, end)
ov <- foverlaps(cpg_dt, te_dt, nomatch = 0L)
per_te <- ov[, .(beta_ctrl = sum(mc)/pmax(sum(cc),1), beta_amp = sum(ma)/pmax(sum(ca),1),
                 n_cpg = .N), by = te_uid][n_cpg >= 3]
te_loc <- te[, .(te_uid, class, kimura = kimura_div_CpG, kimura_unadj = kimura_div)]
te_loc[, loc := ifelse(overlapsAny(te_gr, gene_gr)[match(te_uid, te$te_uid)],
                       "Gene-overlapping", "Intergenic")]
per_te <- merge(per_te, te_loc, by = "te_uid")
per_te <- merge(te[, .(te_uid, chrom, start, end, te_name, class_family)], per_te, by = "te_uid")   # coordinates, so a reviewer can locate any copy
per_te[, class := factor(class, levels = names(COL_CLASS))]
per_te[, loc := factor(loc, levels = names(COL_LOC))]
per_te[, beta_pooled := (beta_ctrl*1 + beta_amp*1)/2]   # simple average of the two condition betas (equal weight)
fwrite(per_te, file.path(DAT, "te_methylation_per_copy.tsv"), sep = "\t")
cat(sprintf("  %s TE copies with >=3 CpGs\n", format(nrow(per_te), big.mark = ",")))

ridge_classes <- c("LTR", "RC", "LINE", "DNA", "SINE")

long <- melt(per_te, id.vars = c("te_uid","class","loc"),
             measure.vars = c("beta_ctrl","beta_amp"),
             variable.name = "condition", value.name = "beta")
long[, condition := factor(ifelse(condition == "beta_ctrl","Control","Amputated"),
                           levels = c("Control","Amputated"))]
long[, class := factor(class, levels = rev(ridge_classes))]
SPLIT_DIR <- "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls"
chh <- vapply(c("C1","C2","A1","A2"), function(s) {
  f <- file.path(SPLIT_DIR, sprintf("%s_paired_bismark_bt2_pe.deduplicated_splitting_report.txt", s))
  if (!file.exists(f)) return(NA_real_)
  l <- grep("C methylated in CHH context", readLines(f), value = TRUE)[1]
  as.numeric(sub("%.*", "", sub(".*:\\s*", "", l)))
}, numeric(1))
nonconv <- mean(chh, na.rm = TRUE)                       # percent
fwrite(data.table(sample = names(chh), chh_pct = chh), file.path(DAT, "non_conversion_rate.tsv"), sep = "\t")
cat(sprintf("  non-conversion proxy (mean CHH): %.2f%%  [%s]\n", nonconv,
            paste(sprintf("%s %.1f", names(chh), chh), collapse = ", ")))

pa <- ggplot(long, aes(beta*100, class, fill = class)) +
  geom_density_ridges(alpha = 0.85, scale = 1.4, colour = "white", linewidth = 0.2) +
  facet_wrap(~ condition, ncol = 2) +
  scale_fill_manual(values = COL_CLASS, guide = "none") +
  scale_x_continuous(limits = c(0, 100)) +
  labs(x = "Per-copy mean CpG methylation β (%)", y = "TE class",
       title = "TE methylation by class") +
  theme_pub()
save_supp(pa, "figS4_te_methylation_by_class_wgbs_tail", 6.0, 3.4)

per_te[, class := factor(class, levels = rev(ridge_classes))]
pb_dt <- copy(per_te)
pb_dt[, loc := factor(loc, levels = c("Intergenic", "Gene-overlapping"))]  # gene-overlapping on top
pb <- ggplot(pb_dt, aes(beta_pooled*100, class, fill = loc)) +
  geom_density_ridges(alpha = 0.7, scale = 1.3, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = COL_LOC, name = "Genomic location",
                    breaks = c("Gene-overlapping", "Intergenic")) +  # legend: green first
  scale_x_continuous(limits = c(-3, 103), breaks = c(0, 25, 50, 75, 100)) +
  labs(x = "Per-copy mean CpG methylation β (%)", y = "TE class",
       title = "TE methylation: gene-overlapping vs intergenic") +
  theme_pub() + theme(legend.position = "right")
save_supp(pb, "figS4_te_methylation_by_class_location_wgbs_tail", 6.5, 3.6)

age <- copy(per_te)[!is.na(kimura)]
qlo <- quantile(age$kimura, probs = seq(0.0, 0.8, 0.2))
qhi <- quantile(age$kimura, probs = seq(0.2, 1.0, 0.2))
qbrk <- c(-Inf, qhi[1:4], Inf)                       # 6 breaks -> 5 quintile bins
qlab <- sprintf("Q%d", 1:5)                          # x-axis: quintile label only (Q1=youngest)
qrng <- sprintf("Q%d %.1f-%.1f%%", 1:5, qlo, qhi)    # quintile -> actual Kimura % range
age[, age_bin := cut(kimura, breaks = qbrk, labels = qlab, include.lowest = TRUE)]
cat(sprintf("  Kimura quintile ranges: %s\n", paste(qrng, collapse = "; ")))
kim_cor <- rbind(
  data.table(class = "All", pearson_r = cor(age$kimura, age$beta_pooled, use = "complete.obs"),
             spearman_rho = cor(age$kimura, age$beta_pooled, method = "spearman", use = "complete.obs"), n = nrow(age)),
  age[, .(pearson_r = cor(kimura, beta_pooled, use = "complete.obs"),
          spearman_rho = cor(kimura, beta_pooled, method = "spearman", use = "complete.obs"), n = .N), by = class])
fwrite(kim_cor, file.path(DAT, "te_kimura_methylation_correlation.tsv"), sep = "\t")
cat("  Kimura vs methylation correlation (Pearson r / Spearman rho):\n"); print(kim_cor)
agg <- age[!is.na(age_bin), .(beta = mean(beta_pooled, na.rm = TRUE),
                              se = sd(beta_pooled, na.rm = TRUE)/sqrt(.N), n = .N),
           by = .(class, loc, age_bin)]
fwrite(agg, file.path(DAT, "te_age_by_class_location.tsv"), sep = "\t")
agg[, class := factor(class, levels = ridge_classes)]
pc <- ggplot(agg, aes(age_bin, beta*100, fill = loc)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = (beta-se)*100, ymax = (beta+se)*100),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.25) +
  facet_wrap(~ class, ncol = 5, scales = "free_x") +
  scale_fill_manual(values = COL_LOC, name = "Genomic location") +
  labs(x = "TE age (Kimura-divergence quintile)", y = "Per-copy mean CpG methylation β (%)",
       title = "TE methylation and Kimura Divergence",
       caption = NULL) +
  theme_pub() + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),   # 45 deg tilt (co-author request)
                      plot.caption = element_text(hjust = 1, size = 6.5, colour = "grey30"))
save_supp(pc, "figS4_te_age_by_class_location_wgbs_tail", 9.0, 3.0)

cat("[4] PacBio HiFi bodywall TE arm\n")
hifi <- readRDS(file.path(B02, "hifi_bodywall_cpg_persample_chrmt.rds"))
hifi <- hifi[chr %in% keep_chr]
hifi_dt <- hifi[, .(chr, start = pos, end = pos, mod = mod_1 + mod_2, cov = cov_1 + cov_2)]
setkey(hifi_dt, chr, start, end)
ovh <- foverlaps(hifi_dt, te_dt, nomatch = 0L)
per_te_h <- ovh[, .(beta_bw = sum(mod) / pmax(sum(cov), 1), n_cpg_hifi = .N),
                by = te_uid][n_cpg_hifi >= 3]                  # same >=3-CpG rule as the WGBS arm
per_te_h <- merge(per_te_h, te_loc, by = "te_uid")
per_te_h <- merge(te[, .(te_uid, chrom, start, end, te_name, class_family)], per_te_h, by = "te_uid")
per_te_h[, loc := factor(loc, levels = names(COL_LOC))]
fwrite(per_te_h, file.path(DAT, "te_methylation_per_copy_hifi.tsv"), sep = "\t")
cat(sprintf("  %s TE copies with >=3 HiFi CpGs vs %s on WGBS cov5 (%.2fx)\n",
            format(nrow(per_te_h), big.mark = ","), format(nrow(per_te), big.mark = ","),
            nrow(per_te_h) / nrow(per_te)))

kb   <- c(0, 5, 10, 15, 20, 25, Inf)
klab <- c("0-5", "5-10", "10-15", "15-20", "20-25", ">25")
cov_kim <- rbind(
  te[!is.na(kimura_div), .(platform = "annotated_total", n = .N),
     by = .(kim_bin = cut(kimura_div, kb, labels = klab, include.lowest = TRUE))],
  per_te[!is.na(kimura), .(platform = "WGBS_tail", n = .N),
         by = .(kim_bin = cut(kimura, kb, labels = klab, include.lowest = TRUE))],
  per_te_h[!is.na(kimura), .(platform = "HiFi_bodywall", n = .N),
           by = .(kim_bin = cut(kimura, kb, labels = klab, include.lowest = TRUE))])
cov_kim <- dcast(cov_kim, kim_bin ~ platform, value.var = "n", fill = 0L)
setcolorder(cov_kim, c("kim_bin", "annotated_total", "WGBS_tail", "HiFi_bodywall"))
cov_kim[, `:=`(pct_wgbs = 100 * WGBS_tail / annotated_total,
               pct_hifi = 100 * HiFi_bodywall / annotated_total)]
fwrite(cov_kim, file.path(DAT, "te_platform_coverage_by_kimura.tsv"), sep = "\t")
cat("  quantifiable copies per Kimura bin (annotated / WGBS / HiFi):\n"); print(cov_kim)

joint <- merge(per_te[, .(te_uid, beta_pooled, n_cpg)],
               per_te_h[, .(te_uid, beta_bw, n_cpg_hifi, class, kimura, loc)], by = "te_uid")
cmp <- rbind(
  data.table(class = "All", n_wgbs = nrow(per_te), n_hifi = nrow(per_te_h), n_both = nrow(joint)),
  Reduce(function(a, b) merge(a, b, by = "class", all = TRUE),
         list(per_te[, .(n_wgbs = .N), by = .(class = as.character(class))],
              per_te_h[, .(n_hifi = .N), by = .(class = as.character(class))],
              joint[, .(n_both = .N), by = .(class = as.character(class))])))
fwrite(cmp, file.path(DAT, "te_platform_coverage_comparison.tsv"), sep = "\t")
agree <- rbind(
  data.table(class = "All",
             pearson_r = joint[, cor(beta_pooled, beta_bw, use = "complete.obs")],
             spearman_rho = joint[, cor(beta_pooled, beta_bw, method = "spearman", use = "complete.obs")],
             n = nrow(joint)),
  joint[, .(pearson_r = cor(beta_pooled, beta_bw, use = "complete.obs"),
            spearman_rho = cor(beta_pooled, beta_bw, method = "spearman", use = "complete.obs"),
            n = .N), by = .(class = as.character(class))])
fwrite(agree, file.path(DAT, "te_platform_percopy_agreement.tsv"), sep = "\t")
cat("  cross-platform per-copy agreement (WGBS tail vs HiFi bodywall; cross-tissue):\n")
print(agree)

ph_dt <- copy(per_te_h)
ph_dt[, class := factor(class, levels = rev(ridge_classes))]
pha <- ggplot(ph_dt, aes(beta_bw * 100, class, fill = class)) +
  geom_density_ridges(alpha = 0.85, scale = 1.4, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = COL_CLASS, guide = "none") +
  scale_x_continuous(limits = c(0, 100)) +
  labs(x = "Per-copy mean CpG methylation β (%)", y = "TE class",
       title = "TE methylation by class") +
  theme_pub()
save_fig(pha, "fig4a_te_methylation_by_class_bodywall", 3.9, 3.16)

phb_dt <- copy(per_te_h)
phb_dt[, class := factor(class, levels = rev(ridge_classes))]
phb_dt[, loc := factor(loc, levels = c("Intergenic", "Gene-overlapping"))]
phb <- ggplot(phb_dt, aes(beta_bw * 100, class, fill = loc)) +
  geom_density_ridges(alpha = 0.7, scale = 1.3, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = COL_LOC, name = "Genomic location",
                    breaks = c("Gene-overlapping", "Intergenic")) +
  scale_x_continuous(limits = c(-3, 103), breaks = c(0, 25, 50, 75, 100)) +
  labs(x = "Per-copy mean CpG methylation β (%)", y = "TE class",
       title = "TE methylation: gene-overlapping vs intergenic") +
  theme_pub() +
  theme(legend.position = "right", plot.title = element_text(size = 9, face = "bold"))
save_fig(phb, "fig4b_te_methylation_by_class_location_bodywall", 3.6, 1.66)

age_h <- copy(per_te_h)[!is.na(kimura)]
qlo_h <- quantile(age_h$kimura, probs = seq(0.0, 0.8, 0.2))
qhi_h <- quantile(age_h$kimura, probs = seq(0.2, 1.0, 0.2))
qbrk_h <- c(-Inf, qhi_h[1:4], Inf)
age_h[, age_bin := cut(kimura, breaks = qbrk_h, labels = qlab, include.lowest = TRUE)]
cat(sprintf("  HiFi Kimura quintile ranges: %s\n",
            paste(sprintf("Q%d %.1f-%.1f%%", 1:5, qlo_h, qhi_h), collapse = "; ")))
kim_cor_h <- rbind(
  data.table(class = "All",
             pearson_r = cor(age_h$kimura, age_h$beta_bw, use = "complete.obs"),
             spearman_rho = cor(age_h$kimura, age_h$beta_bw, method = "spearman", use = "complete.obs"),
             n = nrow(age_h)),
  age_h[, .(pearson_r = cor(kimura, beta_bw, use = "complete.obs"),
            spearman_rho = cor(kimura, beta_bw, method = "spearman", use = "complete.obs"),
            n = .N), by = .(class = as.character(class))])
fwrite(kim_cor_h, file.path(DAT, "te_kimura_methylation_correlation_hifi.tsv"), sep = "\t")
cat("  HiFi Kimura vs methylation correlation (Pearson r / Spearman rho):\n"); print(kim_cor_h)
agg_h <- age_h[!is.na(age_bin), .(beta = mean(beta_bw, na.rm = TRUE),
                                  se = sd(beta_bw, na.rm = TRUE) / sqrt(.N), n = .N),
               by = .(class, loc, age_bin)]
fwrite(agg_h, file.path(DAT, "te_age_by_class_location_hifi.tsv"), sep = "\t")
agg_h[, class := factor(as.character(class), levels = ridge_classes)]
phc <- ggplot(agg_h, aes(age_bin, beta * 100, fill = loc)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = (beta - se) * 100, ymax = (beta + se) * 100),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.25) +
  facet_wrap(~ class, ncol = 5, scales = "free_x") +
  scale_fill_manual(values = COL_LOC, name = "Genomic location") +
  labs(x = "TE age (Kimura-divergence quintile)", y = "Per-copy mean CpG methylation β (%)",
       title = "TE methylation and Kimura divergence") +
  theme_pub() + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
save_fig(phc, "fig4c_te_age_by_class_location_bodywall", 7.0, 2.4)

cat("[5] TE statistics for the text\n")
te_stats <- function(dt, beta_col, agg_tab, binned) {
  d <- copy(dt)[!is.na(get(beta_col))]
  d[, beta := get(beta_col)]
  d[, high := beta > 0.5]
  g <- d[loc == "Gene-overlapping"]; i <- d[loc == "Intergenic"]
  w  <- wilcox.test(g$beta, i$beta)                       # Mann-Whitney U (W = U of group 1)
  rb <- 2 * as.numeric(w$statistic) / (as.numeric(nrow(g)) * nrow(i)) - 1   # double: n1*n2 overflows a 32-bit integer
  ft <- fisher.test(table(d$loc == "Gene-overlapping", d$high))
  row <- function(statistic, scope, value) data.table(statistic = statistic, scope = scope, value = as.numeric(value))
  loc2 <- c("Gene-overlapping", "Intergenic")
  out <- rbindlist(list(
    if (all(c("chrom", "start", "end") %in% names(d)))
      row("pct_copies_overlapping_another", "All",
          100 * mean(countOverlaps(GRanges(d$chrom, IRanges(d$start, d$end))) > 1L)),
    row("n_copies", "All", nrow(d)),
    row("n_copies", loc2, c(nrow(g), nrow(i))),
    row("mean_beta", "All", mean(d$beta)),
    row("mean_beta", loc2, c(mean(g$beta), mean(i$beta))),
    row("median_beta", "All", median(d$beta)),
    row("median_beta", loc2, c(median(g$beta), median(i$beta))),
    row("pct_beta_gt_0.5", "All", 100 * mean(d$high)),
    row("pct_beta_gt_0.5", loc2, 100 * c(mean(g$high), mean(i$high))),
    row("n_copies_beta_gt_0.5", "All", sum(d$high)),
    row("pct_of_high_copies_gene_overlapping", "All", 100 * mean(d[high == TRUE, loc == "Gene-overlapping"])),
    row("mann_whitney_P", "Gene-overlapping vs Intergenic", w$p.value),
    row("rank_biserial_r", "Gene-overlapping vs Intergenic", rb),
    row("fisher_OR_beta_gt_0.5", "Gene-overlapping vs Intergenic", unname(ft$estimate)),
    row("fisher_OR_ci_low", "Gene-overlapping vs Intergenic", ft$conf.int[1]),
    row("fisher_OR_ci_high", "Gene-overlapping vs Intergenic", ft$conf.int[2]),
    row("fisher_P", "Gene-overlapping vs Intergenic", ft$p.value),
    d[, row("n_copies", paste(.BY$class, .BY$loc, sep = " | "), .N), by = .(class, loc)][, .(statistic, scope, value)],
    d[, row("mean_beta", paste(.BY$class, .BY$loc, sep = " | "), mean(beta)), by = .(class, loc)][, .(statistic, scope, value)],
    d[, row("median_beta", paste(.BY$class, .BY$loc, sep = " | "), median(beta)), by = .(class, loc)][, .(statistic, scope, value)],
    d[, row("pct_beta_gt_0.5", paste(.BY$class, .BY$loc, sep = " | "), 100 * mean(high)), by = .(class, loc)][, .(statistic, scope, value)],
    d[, row("n_copies", as.character(.BY$class), .N), by = class][, .(statistic, scope, value)],
    d[, row("mean_beta", as.character(.BY$class), mean(beta)), by = class][, .(statistic, scope, value)],
    d[, row("median_beta", as.character(.BY$class), median(beta)), by = class][, .(statistic, scope, value)],
    d[, row("pct_beta_gt_0.5", as.character(.BY$class), 100 * mean(high)), by = class][, .(statistic, scope, value)]
  ))
  out <- rbind(out,
    row("point_biserial_r", "Gene-overlapping vs Intergenic", cor(as.numeric(d$loc == "Gene-overlapping"), d$beta)),
    row("median_beta_of_copies_gt_0.5", "All", median(d[high == TRUE, beta])),
    row("pct_kimura_gt_50", "All", 100 * mean(d$kimura > 50, na.rm = TRUE)))   # saturated K2P estimates (short fragments)
  k <- d[!is.na(kimura)]
  sp <- function(x) { ct <- cor.test(x$kimura, x$beta, method = "spearman", exact = FALSE)
                      list(rho = unname(ct$estimate), P = ct$p.value, n = nrow(x)) }
  s_all <- sp(k)
  out <- rbind(out,
    row("kimura_spearman_rho", "All", s_all$rho), row("kimura_spearman_P", "All", s_all$P), row("kimura_n", "All", s_all$n),
    k[, { r <- sp(.SD); sc <- as.character(.BY$class)
          rbind(row("kimura_spearman_rho", sc, r$rho), row("kimura_spearman_P", sc, r$P), row("kimura_n", sc, r$n)) }, by = class][, .(statistic, scope, value)],
    k[, { r <- sp(.SD); sc <- paste(.BY$class, .BY$loc, sep = " | ")
          rbind(row("kimura_spearman_rho", sc, r$rho), row("kimura_spearman_P", sc, r$P), row("kimura_n", sc, r$n)) },
      by = .(class, loc)][, .(statistic, scope, value)],
    row("kimura_spearman_rho_range_by_class", "min", min(k[, cor(kimura, beta, method = "spearman"), by = class]$V1)),
    row("kimura_spearman_rho_range_by_class", "max", max(k[, cor(kimura, beta, method = "spearman"), by = class]$V1)),
    if ("kimura_unadj" %in% names(k)) rbindlist(list(
      row("kimura_spearman_rho_UNADJUSTED", "All", cor(k$kimura_unadj, k$beta, method = "spearman", use = "complete.obs")),
      k[, row("kimura_spearman_rho_UNADJUSTED", as.character(.BY$class),
              cor(kimura_unadj, beta, method = "spearman", use = "complete.obs")), by = class][, .(statistic, scope, value)]))
    else NULL)
  q <- agg_tab[, .(value = beta[age_bin == "Q5"] - beta[age_bin == "Q1"],
                   monotone = as.numeric(all(diff(beta[order(as.character(age_bin))]) <= 0)),   # 1 = falls from Q1 to Q5 without a rise
                   n_q1 = n[age_bin == "Q1"], n_q5 = n[age_bin == "Q5"]), by = .(class, loc)]
  bb <- copy(binned)[!is.na(age_bin)]
  bb[, b := get(beta_col)]
  bb <- bb[!is.na(b) & age_bin %in% c("Q1", "Q5")]
  qtest <- bb[, { b5 <- b[age_bin == "Q5"]; b1 <- b[age_bin == "Q1"]
                  if (length(b5) < 20L || length(b1) < 20L) {
                    list(rb = NA_real_, P = NA_real_, dmed = NA_real_)
                  } else {
                    wt <- wilcox.test(b5, b1)
                    list(rb = 2 * as.numeric(wt$statistic) / (as.numeric(length(b5)) * length(b1)) - 1,
                         P = wt$p.value, dmed = median(b5) - median(b1))
                  } }, by = .(class, loc)]
  out <- rbind(out,
    qtest[, row("quintile_Q5_vs_Q1_rank_biserial_r", paste(.BY$class, .BY$loc, sep = " | "), rb), by = .(class, loc)][, .(statistic, scope, value)],
    qtest[, row("quintile_Q5_vs_Q1_mannwhitney_P",   paste(.BY$class, .BY$loc, sep = " | "), P),  by = .(class, loc)][, .(statistic, scope, value)],
    qtest[, row("quintile_Q5_minus_Q1_median_beta",  paste(.BY$class, .BY$loc, sep = " | "), dmed), by = .(class, loc)][, .(statistic, scope, value)])
  out <- rbind(out,
    q[, row("quintile_Q5_minus_Q1_beta", paste(.BY$class, .BY$loc, sep = " | "), value), by = .(class, loc)][, .(statistic, scope, value)],
    q[, row("quintile_monotone_decline", paste(.BY$class, .BY$loc, sep = " | "), monotone), by = .(class, loc)][, .(statistic, scope, value)],
    q[, row("quintile_Q1_n", paste(.BY$class, .BY$loc, sep = " | "), n_q1), by = .(class, loc)][, .(statistic, scope, value)],
    q[, row("quintile_Q5_n", paste(.BY$class, .BY$loc, sep = " | "), n_q5), by = .(class, loc)][, .(statistic, scope, value)])
  out[]
}
stats_h <- te_stats(per_te_h, "beta_bw", agg_h, age_h)
fwrite(stats_h, file.path(DAT, "te_statistics_bodywall.tsv"), sep = "\t")
stats_w <- te_stats(per_te, "beta_pooled", agg, age)
fwrite(stats_w, file.path(DAT, "te_statistics_wgbs_tail.tsv"), sep = "\t")
counts <- data.table(statistic = c("n_classified_copies", "pct_classified_copies_quantified_wgbs_tail",
                                   "pct_classified_copies_quantified_hifi_bodywall", "n_copies_quantified_on_both_platforms"),
                     scope = "All",
                     value = c(nrow(te), 100 * nrow(per_te) / nrow(te), 100 * nrow(per_te_h) / nrow(te), nrow(joint)))
fwrite(counts, file.path(DAT, "te_copy_counts.tsv"), sep = "\t")
fl <- joint[beta_pooled == 0, beta_bw]
floor_dt <- data.table(
  statistic = c("n_copies_wgbs_beta_0", "hifi_median_beta_at_wgbs_beta_0", "hifi_q25_beta_at_wgbs_beta_0",
                "hifi_q75_beta_at_wgbs_beta_0", "n_copies_wgbs_beta_gt_0.9", "hifi_median_beta_at_wgbs_beta_gt_0.9",
                "pct_concordant_at_beta_0.5", "n_gt_0.5_both", "n_gt_0.5_wgbs_only", "n_gt_0.5_hifi_only",
                "pct_wgbs_copies_below_conversion_floor"),
  scope = c(rep("Jointly quantified copies", 10), "All WGBS-quantified copies"),
  value = c(length(fl), median(fl), quantile(fl, 0.25), quantile(fl, 0.75),
            joint[, sum(beta_pooled > 0.9)], joint[beta_pooled > 0.9, median(beta_bw)],
            100 * joint[, mean((beta_pooled > 0.5) == (beta_bw > 0.5))],
            joint[, sum(beta_pooled > 0.5 & beta_bw > 0.5)], joint[, sum(beta_pooled > 0.5 & beta_bw <= 0.5)],
            joint[, sum(beta_pooled <= 0.5 & beta_bw > 0.5)],
            100 * mean(per_te$beta_pooled < mean(chh, na.rm = TRUE) / 100)))
fwrite(floor_dt, file.path(DAT, "te_platform_floor_concordance.tsv"), sep = "\t")
cat("  cross-platform floor / concordance:\n"); print(floor_dt)
cat("  bodywall (HiFi) headline statistics:\n")
print(stats_h[scope %in% c("All", "Gene-overlapping", "Intergenic", "Gene-overlapping vs Intergenic")])
cat("  WGBS tail headline statistics:\n")
print(stats_w[scope %in% c("All", "Gene-overlapping", "Intergenic", "Gene-overlapping vs Intergenic")])

cc_dt <- per_te[!is.na(beta_ctrl) & !is.na(beta_amp)]
cc_dt[, delta := beta_amp - beta_ctrl]
cond_all <- cc_dt[, .(class = "All", n = .N, mean_beta_ctrl = mean(beta_ctrl), mean_beta_amp = mean(beta_amp),
                      mean_delta = mean(delta), median_delta = median(delta),
                      pct_abs_delta_gt_0.1 = 100 * mean(abs(delta) > 0.1),
                      pct_gain_gt_0.1 = 100 * mean(delta > 0.1), pct_loss_gt_0.1 = 100 * mean(delta < -0.1),
                      pct_abs_delta_gt_0.2 = 100 * mean(abs(delta) > 0.2),
                      pct_gain_gt_0.2 = 100 * mean(delta > 0.2), pct_loss_gt_0.2 = 100 * mean(delta < -0.2),
                      wilcoxon_paired_P = wilcox.test(beta_amp, beta_ctrl, paired = TRUE)$p.value)]
cond_cls <- cc_dt[, .(n = .N, mean_beta_ctrl = mean(beta_ctrl), mean_beta_amp = mean(beta_amp),
                      mean_delta = mean(delta), median_delta = median(delta),
                      pct_abs_delta_gt_0.1 = 100 * mean(abs(delta) > 0.1),
                      pct_gain_gt_0.1 = 100 * mean(delta > 0.1), pct_loss_gt_0.1 = 100 * mean(delta < -0.1),
                      pct_abs_delta_gt_0.2 = 100 * mean(abs(delta) > 0.2),
                      pct_gain_gt_0.2 = 100 * mean(delta > 0.2), pct_loss_gt_0.2 = 100 * mean(delta < -0.2),
                      wilcoxon_paired_P = wilcox.test(beta_amp, beta_ctrl, paired = TRUE)$p.value),
                  by = .(class = as.character(class))]
cond_tab <- rbind(cond_all, cond_cls)
fwrite(cond_tab, file.path(DAT, "te_condition_contrast_wgbs_tail.tsv"), sep = "\t")
cat("  WGBS tail control vs amputated per copy:\n"); print(cond_tab)

writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_04_TEs.txt"))
cat("[04_TEs] done\n")
