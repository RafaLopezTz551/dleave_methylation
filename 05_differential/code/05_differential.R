#!/usr/bin/env Rscript
# 05_differential.R
# Differential methylation after tail amputation: DSS DMP/DMR calling, gene assignment,
# region and GO/KEGG enrichment, DMR-in-TE permutation test and control analyses.
# Inputs: 02_landscape/objects/bsseq_cov5_chrmt.rds; 01_genome_toolkit/objects/gff_chrmt.rds;
#   01_genome_toolkit/data/gene_de_tail.tsv and argonaute_census.tsv; STRING v12 per-protein
#   GO/KEGG terms; TE copy table with Kimura ages; OMARK conserved-unknown protein fasta.
# Outputs: data/ DMP/DMR annotation, gene-assignment, enrichment and control tables (TSV);
#   objects/ per-chromosome DMLtest caches; figures/main fig5a-fig5h and
#   figures/supplementary fig5h/fig5i hyper-hypo panels, fig5j-fig5n; sessionInfo_05_differential.txt
# Run: sbatch 05_differential/code/05_differential.slurm  (from main/methylation_pipeline/)

# Step 1 - Seed, packages, paths, palettes, plot theme and figure savers
set.seed(20260426)
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges)
  library(bsseq); library(DSS); library(ggplot2); library(scales); library(patchwork)
})

# External annotation: STRING v12 GO/KEGG terms and the TE copy table used in Step 19
STRING_ENR <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/STRING.protein.enrichment.terms.v12.0.txt"
TE <- "/mnt/data/alfredvar/30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv"
PIPE <- "/mnt/data/alfredvar/rlopezt/meth_paper/main/methylation_pipeline"
B01 <- file.path(PIPE, "01_genome_toolkit/objects")
B02 <- file.path(PIPE, "02_landscape/objects")
BATCH <- file.path(PIPE, "05_differential")
OBJ <- file.path(BATCH, "objects"); DAT <- file.path(BATCH, "data")
FIGM <- file.path(BATCH, "figures/main"); FIGS <- file.path(BATCH, "figures/supplementary")
for (d in c(OBJ, DAT, FIGM, FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
# Remove figures from earlier runs; cached objects and data tables are kept
unlink(c(list.files(FIGM, "\\.(pdf|png|svg)$", full.names = TRUE),
         list.files(FIGS, "\\.(pdf|png|svg)$", full.names = TRUE)))
keep_chr <- c(paste0("chr", 1:31), "HiC_scaffold_1563")   # chr1-31 plus the mitochondrial scaffold

# Colours: direction (Hyper/Hypo in amputated), genomic region, gene biotype
COL_DIR <- c(Hyper = "#C0392B", Hypo = "#0072B2", NS = "#9E9E9E")
COL_REGION <- c(Promoter = "#2C7FB8", Exon = "#1B9E9E", Intron = "#6A51A3", Intergenic = "#7FB3D5")
# EviAnn assigns every gene one of three biotypes; all three are kept explicit
BIOTYPES <- c("protein_coding", "lncRNA", "processed_pseudogene")
COL_BT   <- c(protein_coding = "#4C72B0", lncRNA = "#DD8452", processed_pseudogene = "#55A868")
LAB_BT   <- c(protein_coding = "Protein coding", lncRNA = "lncRNA",
              processed_pseudogene = "Pseudogenes")
# Plot theme and savers (pdf via cairo, png, svg) for main and supplementary figures
theme_pub <- function() theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGM, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGM, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGM, paste0(name, ".svg")), p, width = w, height = h)
  cat(sprintf("  saved %s\n", name))
}
save_supp <- function(p, name, w, h) {
  ggsave(file.path(FIGS, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGS, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(FIGS, paste0(name, ".svg")), p, width = w, height = h)
  cat(sprintf("  saved supp/%s\n", name))
}

# Step 2 - Load the CpG object and annotation; derive gene symbols and biotypes
cat("[1] load bsseq + gff\n")
bs <- readRDS(file.path(B02, "bsseq_cov5_chrmt.rds"))
chrs <- as.character(seqnames(bs)); bs <- bs[chrs %in% keep_chr, ]
gff <- readRDS(file.path(B01, "gff_chrmt.rds"))
gene_gr <- gff[gff$type == "gene"]; exon_gr <- gff[gff$type == "exon"]
prom_gr <- suppressWarnings(trim(promoters(gene_gr, 2000, 0)))
gid_of_gene <- sub(";.*", "", as.character(mcols(gene_gr)$ID))

# Gene symbol parsed from the GFF Note field ('Similar to SYM: ...'), NA when absent
g_note <- vapply(mcols(gene_gr)$Note,
                 function(x) if (length(x)) as.character(x)[1] else NA_character_, character(1))
sym_all <- sub("^Similar to ([^:]+):.*$", "\\1", g_note)
sym_all[!grepl("^Similar to [^:]+:", g_note)] <- NA_character_
names(sym_all) <- gid_of_gene
disp_name <- function(gid) { s <- unname(sym_all[gid]); ifelse(is.na(s), gid, s) }
# EviAnn gene_biotype per gene (summary tables only)
biotype_all <- as.character(mcols(gene_gr)$gene_biotype); names(biotype_all) <- gid_of_gene

# Step 3 - DSS DMLtest per chromosome (cached), then callDML (DMPs) and callDMR (DMRs)
cat("[2] DSS DMLtest (per chromosome — memory-safe) + callDML/callDMR\n")
# DMLtest runs per chromosome (smoothing does not cross chromosomes) and is cached
# per chromosome; a cache older than the bsseq object it derives from stops the run
dml_rds <- file.path(OBJ, "dmltest_chrmt.rds")
BS_MTIME <- file.mtime(file.path(B02, "bsseq_cov5_chrmt.rds"))
if (file.exists(dml_rds) && file.mtime(dml_rds) < BS_MTIME)
  stop("objects/dmltest_chrmt.rds is OLDER than the bsseq object it was computed from: delete the dmltest_* caches and rerun", call. = FALSE)
if (file.exists(dml_rds)) dml_test <- readRDS(dml_rds) else {
  chrs <- intersect(keep_chr, as.character(unique(seqnames(bs))))
  dml_list <- lapply(chrs, function(ch) {
    cf <- file.path(OBJ, sprintf("dmltest_chrmt_%s.rds", ch))
    if (file.exists(cf)) { cat(sprintf("  reuse DMLtest %s\n", ch)); return(readRDS(cf)) }
    cat(sprintf("  DMLtest %s\n", ch))
    r <- DMLtest(bs[as.character(seqnames(bs)) == ch, ],
                 group1 = c("C1","C2"), group2 = c("A1","A2"), smoothing = TRUE)
    saveRDS(r, cf); r
  })
  dml_test <- do.call(rbind, dml_list)
  saveRDS(dml_test, dml_rds)
}
# STAT TEST: DSS callDML = Wald test on the smoothed, dispersion-shrunk beta-binomial
# methylation model, per CpG (Wu et al., DSS). A DMP requires p < 0.05 AND |delta-beta| >= 0.10.
dmps <- as.data.table(callDML(dml_test, p.threshold = 0.05, delta = 0.10))
dmrs_raw <- callDMR(dml_test, p.threshold = 0.05, delta = 0.10, minlen = 50, minCG = 3)
dmrs <- if (is.null(dmrs_raw)) data.table() else as.data.table(dmrs_raw)
# diff = mu_control - mu_amputated, so diff < 0 means hypermethylated in amputated
dmps[, direction := ifelse(diff < 0, "Hyper", "Hypo")]
if (nrow(dmrs)) dmrs[, direction := ifelse(diff.Methy < 0, "Hyper", "Hypo")]
cat(sprintf("  DMPs %s | DMRs %s\n", format(nrow(dmps), big.mark=","), format(nrow(dmrs), big.mark=",")))

# Step 4 - Region annotation and multi-gene assignment (2 kb promoter + gene body)
# Region label precedence: Promoter > Exon > Intron > Intergenic
annotate_region <- function(chr, pos) {
  q <- GRanges(chr, IRanges(pos, width = 1))
  r <- rep("Intergenic", length(q))
  r[overlapsAny(q, gene_gr)] <- "Intron"
  r[overlapsAny(q, exon_gr)] <- "Exon"
  r[overlapsAny(q, prom_gr)] <- "Promoter"
  factor(r, levels = names(COL_REGION))
}
# Same precedence applied to a DMR's whole interval
annotate_region_range <- function(chr, start, end) {
  q <- GRanges(chr, IRanges(start, end))
  r <- rep("Intergenic", length(q))
  r[overlapsAny(q, gene_gr)] <- "Intron"
  r[overlapsAny(q, exon_gr)] <- "Exon"
  r[overlapsAny(q, prom_gr)] <- "Promoter"
  factor(r, levels = names(COL_REGION))
}
# Feature universe: 2 kb strand-aware promoter and gene body, each tagged with its gene_id
gene_body <- gene_gr; mcols(gene_body) <- DataFrame(gene_id = gid_of_gene, feature = "Body")
gene_prom <- prom_gr; mcols(gene_prom) <- DataFrame(gene_id = gid_of_gene, feature = "Promoter")
feat_gr   <- c(gene_body, gene_prom)

# Long-format overlap: one row per (feature, gene) pair, so a feature can map to several genes
assign_long <- function(gr, subj) {
  h <- findOverlaps(gr, subj)
  data.table(row = queryHits(h),
             gene_id = mcols(subj)$gene_id[subjectHits(h)],
             feature = as.character(mcols(subj)$feature[subjectHits(h)]))
}
dmp_gr <- GRanges(dmps$chr, IRanges(dmps$pos, width = 1))
dmr_gr <- if (nrow(dmrs)) GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end)) else GRanges()

# Single region label per feature (for the pies): DMP by position, DMR by interval
dmps[, region := annotate_region(chr, pos)]; dmps[, row := .I]
if (nrow(dmrs)) { dmrs[, region := annotate_region_range(chr, start, end)]; dmrs[, row := .I] }

# Primary assignment in long format, with direction and region carried onto each row
dmp_genes <- merge(assign_long(dmp_gr, feat_gr),
                   dmps[, .(row, chr, pos, diff, fdr, direction, region)], by = "row")
dmp_genes[, symbol := disp_name(gene_id)]
dmr_genes <- if (nrow(dmrs))
  merge(assign_long(dmr_gr, feat_gr),
        dmrs[, .(row, chr, start, end, diff.Methy, direction, region)], by = "row")[, symbol := disp_name(gene_id)][] else data.table()

# One display gene per feature: Body before Promoter, then a named gene, else the first
pick_display <- function(asg) {
  a <- copy(asg); a[, feat_rank := fifelse(feature == "Body", 1L, 2L)]
  a[, named := as.integer(!is.na(sym_all[gene_id]))]
  setorder(a, row, feat_rank, -named); a[, .(gene_id = gene_id[1L]), by = row]
}
dmps[, gene_id := NA_character_]; dp <- pick_display(dmp_genes); dmps[dp$row, gene_id := dp$gene_id]
if (nrow(dmrs)) { dmrs[, gene_id := NA_character_]
  if (nrow(dmr_genes)) { dr <- pick_display(dmr_genes); dmrs[dr$row, gene_id := dr$gene_id] } }

# Per-feature tables and the long multi-gene assignment tables
fwrite(dmps, file.path(DAT, "dmps_annotated.tsv"), sep = "\t")
fwrite(dmrs, file.path(DAT, "dmrs_annotated.tsv"), sep = "\t")
fwrite(dmp_genes, file.path(DAT, "dmps_gene_assignments.tsv"), sep = "\t")
if (nrow(dmr_genes)) fwrite(dmr_genes, file.path(DAT, "dmrs_gene_assignments.tsv"), sep = "\t")

# Step 5 - Orphan features overlapping no promoter or gene body
orphan_dmp <- setdiff(seq_len(nrow(dmps)), unique(dmp_genes$row))
fwrite(dmps[orphan_dmp], file.path(DAT, "dmps_orphan_intergenic.tsv"), sep = "\t")
n_dmr_orphan <- 0L
if (nrow(dmrs)) {
  orphan_dmr <- setdiff(seq_len(nrow(dmrs)), unique(dmr_genes$row))
  fwrite(dmrs[orphan_dmr], file.path(DAT, "dmrs_orphan_intergenic.tsv"), sep = "\t")
  n_dmr_orphan <- length(orphan_dmr)
}
cat(sprintf("  primary assign: DMP %d genes / %d rows | DMR %d genes / %d rows | orphans DMP %d DMR %d\n",
            uniqueN(dmp_genes$gene_id), nrow(dmp_genes),
            if (nrow(dmr_genes)) uniqueN(dmr_genes$gene_id) else 0L, nrow(dmr_genes),
            length(orphan_dmp), n_dmr_orphan))

# Step 6 - Gene-annotation summary: known vs unknown symbol, by biotype
gene_anno_summary <- function(gids, set_name) {
  gids <- unique(gids[!is.na(gids)])
  sym  <- sym_all[gids]; bt <- biotype_all[gids]
  # The biotype set is closed; an unexpected biotype stops the run
  stopifnot(all(bt[!is.na(bt)] %in% BIOTYPES))
  tb <- table(factor(bt, levels = BIOTYPES))
  data.table(set = set_name, n_genes = length(gids),
             known   = sum(!is.na(sym)), unknown = sum(is.na(sym)),
             protein_coding       = as.integer(tb[["protein_coding"]]),
             lncRNA               = as.integer(tb[["lncRNA"]]),
             processed_pseudogene = as.integer(tb[["processed_pseudogene"]]))
}
gene_summary <- rbind(gene_anno_summary(dmp_genes$gene_id, "DMP genes"),
                      if (nrow(dmr_genes)) gene_anno_summary(dmr_genes$gene_id, "DMR genes"))
fwrite(gene_summary, file.path(DAT, "dmp_dmr_gene_annotation_summary.tsv"), sep = "\t")
cat("  gene-annotation summary (known/unknown, by EviAnn gene biotype):\n")
print(gene_summary)

# Step 7 - Feature counts by biotype and on conserved genes without annotation
# Conserved-unknown genes: OMARK proteins with gastropod orthologs but no functional annotation
OMARK_FA <- "/mnt/data/alfredvar/rlopezt/OMARK/gasteropods_unknown_ortos_dlaeve.fasta"
stopifnot(file.exists(OMARK_FA))
cons_unknown <- unique(sub("-mRNA-.*$", "", sub("^.*\\|", "",
                       grep("^>", readLines(OMARK_FA), value = TRUE))))
feat_counts <- function(gtab, set_name) {
  gtab <- gtab[!is.na(gene_id)]
  bt  <- biotype_all[gtab$gene_id]; sym <- sym_all[gtab$gene_id]
  data.table(set = set_name,
             n_features          = nrow(gtab),
             n_genes             = uniqueN(gtab$gene_id),
             lnc_genes           = uniqueN(gtab$gene_id[bt == "lncRNA"]),
             lnc_features        = sum(bt == "lncRNA", na.rm = TRUE),
             unknown_genes       = uniqueN(gtab$gene_id[is.na(sym)]),
             unknown_features    = sum(is.na(sym)),
             cons_unknown_genes  = uniqueN(intersect(gtab$gene_id, cons_unknown)),
             cons_unknown_features = sum(gtab$gene_id %in% cons_unknown))
}
feat_summary <- rbind(feat_counts(dmp_genes, "DMPs"),
                      if (nrow(dmr_genes)) feat_counts(dmr_genes, "DMRs"))
fwrite(feat_summary, file.path(DAT, "biotype_annotation_feature_counts.tsv"), sep = "\t")
cu_dmr <- dmr_genes[gene_id %in% cons_unknown, .(n_dmr = .N), by = gene_id]
cu_dmp <- dmp_genes[gene_id %in% cons_unknown, .(n_dmp = .N), by = gene_id]
cu <- merge(cu_dmp, cu_dmr, by = "gene_id", all = TRUE)
cu[is.na(n_dmp), n_dmp := 0L][is.na(n_dmr), n_dmr := 0L]; setorder(cu, -n_dmr, -n_dmp)
fwrite(cu, file.path(DAT, "conserved_unknown_dmp_dmr.tsv"), sep = "\t")
cat(sprintf("  conserved-unknown (OMARK, %d loci): %d with >=1 DMP, %d with >=1 DMR; top: %s (%d DMRs, %d DMPs)\n",
            length(cons_unknown), nrow(cu_dmp), nrow(cu_dmr),
            cu$gene_id[1], cu$n_dmr[1], cu$n_dmp[1]))
print(feat_summary)

# Step 8 - fig5a DMP volcano: hex-density background of tested CpGs plus DMPs
suppressPackageStartupMessages({ library(ggrepel); library(hexbin) })
dml <- as.data.table(dml_test)
# Per-chromosome DMLtests are BH-adjusted separately; a genome-wide FDR is recomputed for
# the plot. x = delta beta (amputated - control) = -diff; y = -log10(FDR), capped at 50
dml[, fdr := p.adjust(pval, "BH")]
dml[, x_db := -diff]
dml[, neglog := pmin(-log10(pmax(fdr, 1e-300)), 50)]
hb <- hexbin(dml$x_db, dml$neglog, xbins = 250, IDs = FALSE)
hb_dt <- data.table(x = hcell2xy(hb)$x, y = hcell2xy(hb)$y, count = hb@count)

# Overlay the called DMPs on the same genome-wide FDR scale
dmps[, x_db := -diff]
dmps[dml, on = c("chr", "pos"), fdr_gw := i.fdr]
dmps[, neglog := pmin(-log10(pmax(fdr_gw, 1e-300)), 50)]
dmps[, gene_name := ifelse(!is.na(gene_id) & !is.na(sym_all[gene_id]),
                           unname(sym_all[gene_id]), "")]
# One label per gene: the 15 top-ranked distinct symbols
top_lab <- dmps[gene_name != ""][order(-abs(x_db) * neglog)][!duplicated(gene_name)][seq_len(min(15, .N))]
pa <- ggplot() +
  geom_hex(data = hb_dt, aes(x, y, fill = count), stat = "identity") +
  scale_fill_viridis_c(name = "CpGs / bin", trans = "log10", labels = comma) +
  geom_vline(xintercept = c(-0.10, 0.10), linetype = 2, colour = "grey40") +
  # No FDR reference line: DSS calls DMPs on the posterior probability, not on BH FDR
  geom_point(data = dmps, aes(x_db, neglog, colour = direction), alpha = 0.6, size = 0.5) +
  scale_colour_manual(values = COL_DIR[c("Hyper", "Hypo")], name = "DMP direction") +
  geom_text_repel(data = top_lab, aes(x_db, neglog, label = gene_name),
                  size = 2.6, max.overlaps = 20, min.segment.length = 0,
                  box.padding = 0.4, colour = "black", seed = 20260426) +
  labs(x = expression(Delta*beta~"(Amputated - Control)"),
       y = expression(-log[10]*"(BH-adjusted "*italic(p)*")"),
       title = "DMP volcano",
       subtitle = sprintf("%s DMPs of %s CpGs", comma(nrow(dmps)), comma(nrow(dml)))) +
  theme_pub()
save_fig(pa, "fig5a_dmp_volcano", 4.0, 2.67)
rm(hb, hb_dt); gc(verbose = FALSE)   # dml is reused in Step 11

# Step 9 - fig5b DMR volcano: delta beta vs |areaStat| for the called DMRs
# callDMR gives no per-region p-value; |areaStat| is the DSS region statistic used on y
if (nrow(dmrs)) {
  dmrs[, x_db := -diff.Methy]
  dmrs[, y_area := abs(areaStat)]
  dmrs[, gene_name := ifelse(!is.na(gene_id) & !is.na(sym_all[gene_id]),
                             unname(sym_all[gene_id]), "")]
  top_lab_dmr <- dmrs[gene_name != ""][order(-abs(x_db) * y_area)][!duplicated(gene_name)][seq_len(min(15, .N))]
  pb <- ggplot(dmrs, aes(x_db, y_area, colour = direction)) +
    geom_vline(xintercept = c(-0.10, 0.10), linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.6, size = 0.9) +
    scale_colour_manual(values = COL_DIR[c("Hyper", "Hypo")], name = "DMR direction") +
    geom_text_repel(data = top_lab_dmr, aes(x_db, y_area, label = gene_name),
                    size = 2.6, max.overlaps = 20, min.segment.length = 0,
                    box.padding = 0.4, colour = "black", seed = 20260426) +
    labs(x = expression(Delta*beta~"(Amputated - Control)"),
         y = "DMR area statistic (|areaStat|)",
         title = "DMR volcano",
         subtitle = sprintf("%s DMRs (%s Hyper / %s Hypo); DSS callDMR, |Delta-beta| >= 0.10",
                            comma(nrow(dmrs)),
                            comma(sum(dmrs$direction == "Hyper")),
                            comma(sum(dmrs$direction == "Hypo")))) +
    theme_pub()
  save_fig(pb, "fig5b_dmr_volcano", 4.0, 2.67)
}

# Step 10 - fig5c and fig5d region-distribution pies (DMP, DMR)
mk_pie <- function(dt, title) {
  tab <- dt[!is.na(region), .N, by = region][, frac := N/sum(N)][]
  tab[, region := factor(region, levels = names(COL_REGION))]
  ggplot(tab, aes("", frac, fill = region)) +
    geom_col(width = 1, colour = "white") + coord_polar(theta = "y") +
    geom_text(aes(label = sprintf("%s\n%.0f%%", region, 100*frac)),
              position = position_stack(vjust = 0.5), size = 2.6) +
    scale_fill_manual(values = COL_REGION, guide = "none") +
    labs(title = title, x = NULL, y = NULL) + theme_pub() +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank())
}
save_fig(mk_pie(dmps, "DMP region distribution"), "fig5c_dmp_region_pie", 4.2, 4.0)
if (nrow(dmrs)) save_fig(mk_pie(dmrs, "DMR region distribution"), "fig5d_dmr_region_pie", 2.9, 2.76)

# Step 11 - Region enrichment vs genomic background (Fisher); fig5e and fig5e2
set.seed(20260426)   # reseed so the sample is independent of earlier RNG use
# DMP background: 200k random tested CpGs annotated by position; four separate regions
bgsamp <- dml[sample(.N, min(.N, 200000))]
bgsamp[, region := annotate_region(chr, pos)]
region_enrichment_fig <- function(feat, bg, kind, name) {
  enr <- rbindlist(lapply(names(COL_REGION), function(rg) {
    a <- sum(feat$region == rg, na.rm = TRUE); b <- nrow(feat) - a
    c_ <- sum(bg$region == rg, na.rm = TRUE); d <- nrow(bg) - c_
    # STAT TEST: two-sided Fisher's exact test on the 2x2 (feature vs genomic CpG background)
    ft <- fisher.test(matrix(c(a, b, c_, d), 2))
    data.table(region = rg, OR = unname(ft$estimate), lo = ft$conf.int[1],
               hi = ft$conf.int[2], p = ft$p.value)
  }))
  enr[, fdr := p.adjust(p, "BH")]; enr[, region := factor(region, levels = names(COL_REGION))]
  fwrite(enr, file.path(DAT, sprintf("%s_region_enrichment.tsv", kind)), sep = "\t")
  # Bar per region: odds ratio with 95% CI and a significance label
  enr[, sig := ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", "")))]
  pe <- ggplot(enr, aes(region, OR, fill = region)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2fx%s", OR, sig), y = hi), vjust = -0.4, size = 2.8) +
    scale_fill_manual(values = COL_REGION, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Odds ratio vs matched background",
         title = sprintf("%s region enrichment", toupper(kind))) + theme_pub()
  save_fig(pe, name, 4.3, 2.87)
}
region_enrichment_fig(dmps, bgsamp, "dmp", "fig5e_dmp_region_enrichment")
if (nrow(dmrs)) {
  set.seed(20260426)
  # DMR background: random intervals with the DMR width distribution, anchored at tested CpGs
  n_bg   <- 200000L
  dmr_w  <- dmrs$end - dmrs$start + 1L
  anc    <- dml[sample(.N, n_bg, replace = TRUE)]
  bg_dmr <- data.table(chr = anc$chr, start = anc$pos,
                       end = anc$pos + sample(dmr_w, n_bg, replace = TRUE) - 1L)
  bg_dmr[, region := annotate_region_range(chr, start, end)]
  region_enrichment_fig(dmrs, bg_dmr, "dmr", "fig5e2_dmr_region_enrichment")
}

# Step 12 - Per-gene DMP and DMR burden tables; fig5f top DMP-burden genes
# Bar fill = fraction of features hypermethylated in amputated; labels = symbol else LOC id
top_burden_fig <- function(dt, kind, name, xlab, draw = TRUE) {
  b <- dt[!is.na(gene_id), .(n = .N, hyper = sum(direction == "Hyper")),
          by = gene_id][order(-n)]
  b[, symbol := disp_name(gene_id)]
  # The burden table is written even when no figure is drawn
  fwrite(b, file.path(DAT, sprintf("%s_gene_burden.tsv", kind)), sep = "\t")
  if (!draw) return(invisible(b))
  top <- head(b, 20)
  top[, lab := factor(make.unique(symbol), levels = rev(make.unique(symbol)))]
  p <- ggplot(top, aes(n, lab, fill = hyper / n)) + geom_col() +
    scale_fill_gradient2(low = COL_DIR["Hypo"], mid = "grey85", high = COL_DIR["Hyper"],
                         midpoint = 0.5, name = "Frac Hyper", limits = c(0, 1)) +
    labs(x = xlab, y = NULL, title = sprintf("Top %s-burden genes", toupper(kind))) +
    theme_pub() + theme(axis.text.y = element_text(size = 7))
  save_fig(p, name, 4.3, 4.3)
}
top_burden_fig(dmp_genes, "dmp", "fig5f_top_dmp_burden_genes", "DMPs per gene")
# DMR burden: table only (draw = FALSE)
if (nrow(dmr_genes)) top_burden_fig(dmr_genes, "dmr", "fig5g_top_dmr_burden_genes",
                                    "DMRs per gene", draw = FALSE)

# Step 13 - GO/KEGG over-representation (STRING v12 terms); fig5h and fig5i dotplots
cat("[7b] GO/KEGG enrichment (STRING v12)\n")
suppressPackageStartupMessages(library(clusterProfiler))
# STRING protein id -> gene id; universe per category = genes carrying that category's terms
sterms <- fread(STRING_ENR, col.names = c("string_id", "category", "term", "description"))
sterms[, gene_id := sub("^[^.]+\\.", "", string_id)]
GO_CATS <- c(BP = "Biological Process (Gene Ontology)",
             MF = "Molecular Function (Gene Ontology)",
             CC = "Cellular Component (Gene Ontology)",
             KEGG = "KEGG (Kyoto Encyclopedia of Genes and Genomes)")

run_enr <- function(sig, cat_label) {   # one category ORA
  sub <- sterms[category == cat_label]
  # Universe restricted to the analysis gene set (chr1-31 + mito), not the whole STRING table
  sub <- sub[gene_id %in% gid_of_gene]
  uni <- unique(sub$gene_id)
  hit <- intersect(sig, uni)
  if (length(hit) < 5) return(data.table())
  # STAT TEST: hypergeometric over-representation (clusterProfiler::enricher), BH within category
  res <- enricher(gene = hit, universe = uni,
                  TERM2GENE = sub[, .(term, gene_id)],
                  TERM2NAME = unique(sub[, .(term, description)]),
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = 5, maxGSSize = 500)
  if (is.null(res)) return(data.table())
  as.data.table(as.data.frame(res))
}

# Enrichment across the four categories; writes <tag>_go_kegg_enrichment.tsv
enrich_all <- function(sig, kind) {
  tag <- gsub("[ /]", "_", tolower(kind))
  all_enr <- rbindlist(lapply(names(GO_CATS), function(o) {
    d <- run_enr(sig, GO_CATS[[o]]); if (nrow(d)) d[, ontology := o]; d
  }), fill = TRUE)
  all_enr[, n_genes := length(sig)]
  fwrite(all_enr, file.path(DAT, sprintf("%s_go_kegg_enrichment.tsv", tag)), sep = "\t")
  all_enr
}

# Dotplot: x = fold enrichment, colour = BH FDR, size = gene count; all ontologies together
go_dotplot_main <- function(all_enr, kind, name, n_top = 18, supp = FALSE) {
  go <- all_enr[ontology %in% c("BP", "MF", "CC", "KEGG")]
  if (!nrow(go)) { cat(sprintf("  %s: no terms for dotplot\n", kind)); return(invisible()) }
  # Reserve up to a quarter of the panel for significant KEGG terms, then fill by FDR
  n_kegg <- min(sum(go$ontology == "KEGG" & go$p.adjust < 0.05), floor(n_top / 4))
  keg <- if (n_kegg > 0) go[ontology == "KEGG" & p.adjust < 0.05][order(p.adjust)][seq_len(n_kegg)]
         else go[0]
  rest <- go[!(ID %in% keg$ID)][order(p.adjust)][seq_len(min(n_top - nrow(keg), .N))]
  top  <- rbind(keg, rest)[order(p.adjust)]
  cat(sprintf("  %s dotplot: %d terms plotted (%d KEGG)\n",
              kind, nrow(top), sum(top$ontology == "KEGG")))
  top[, lab := factor(make.unique(sprintf("%s (%s)", Description, ontology)),
                      levels = rev(make.unique(sprintf("%s (%s)", Description, ontology))))]
  p <- ggplot(top, aes(FoldEnrichment, lab, colour = p.adjust, size = Count)) +
    geom_point() +
    scale_colour_gradient(low = "#C0392B", high = "#2471A3", name = "BH FDR", limits = c(0, 1), oob = scales::squish,
                          guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(2, 8), name = "Genes enriched") +
    labs(x = "Fold enrichment (observed / expected)", y = NULL,
         title = sprintf("%s gene GO and KEGG enrichment", kind),
         subtitle = sprintf("%d of %d terms at FDR < 0.05", sum(go$p.adjust < 0.05), nrow(go))) +
    theme_pub() + theme(axis.text.y = element_text(size = 9))
  (if (supp) save_supp else save_fig)(p, name, 7.0, 4.4)
}

# Main panel: all DMP genes (fig5h); supplementary: hyper/hypo splits and the DMR panels
dmp_enr <- enrich_all(unique(dmp_genes$gene_id), "DMP")
go_dotplot_main(dmp_enr, "DMP", "fig5h_dmp_go_dotplot")
dmp_hyper_enr <- enrich_all(unique(dmp_genes[direction == "Hyper", gene_id]), "DMP hyper")
go_dotplot_main(dmp_hyper_enr, "DMP hyper", "fig5h_dmp_hyper_go", supp = TRUE)
dmp_hypo_enr <- enrich_all(unique(dmp_genes[direction == "Hypo",  gene_id]), "DMP hypo")
go_dotplot_main(dmp_hypo_enr,  "DMP hypo",  "fig5h_dmp_hypo_go",  supp = TRUE)
if (nrow(dmr_genes)) {
  dmr_enr <- enrich_all(unique(dmr_genes$gene_id), "DMR")
  go_dotplot_main(dmr_enr, "DMR", "fig5i_dmr_go_dotplot", supp = TRUE)
  go_dotplot_main(enrich_all(unique(dmr_genes[direction == "Hyper", gene_id]), "DMR hyper"), "DMR hyper", "fig5i_dmr_hyper_go", supp = TRUE)
  go_dotplot_main(enrich_all(unique(dmr_genes[direction == "Hypo",  gene_id]), "DMR hypo"),  "DMR hypo",  "fig5i_dmr_hypo_go",  supp = TRUE)
}

# Step 14 - Target-size adjusted GO/KEGG (goseq Wallenius; bias = CpGs per gene)
cat("[7c] GOmeth target-size adjusted enrichment (goseq Wallenius; bias = analysed CpGs per gene target)\n")
# Longer genes carry more analysed CpGs and collect more DMPs; goseq weights each gene by
# the number of analysed CpGs over its target (and, as a sensitivity check, by gene length)
# Share of genes with a DMP per gene-length quintile
lq <- data.table(gene_id = gid_of_gene, len = width(gene_gr))[, lenq := cut(len, quantile(len, seq(0, 1, 0.2)), include.lowest = TRUE, labels = FALSE)]
lq[, has_dmp := gene_id %in% dmp_genes$gene_id]
fwrite(lq[, .(n_genes = .N, median_len = as.numeric(median(len)), pct_with_dmp = 100 * mean(has_dmp)), by = lenq][order(lenq)],
       file.path(DAT, "dmp_share_by_length_quintile.tsv"), sep = "\t")
suppressPackageStartupMessages(library(goseq))
tgt_gr <- GRanges(seqnames(gene_gr),
                  IRanges(pmax(1L, pmin(start(gene_gr), start(prom_gr))), pmax(end(gene_gr), end(prom_gr))))
bias_ncpg <- setNames(countOverlaps(tgt_gr, granges(bs)), gid_of_gene)   # analysed CpGs per target
bias_len  <- setNames(width(gene_gr), gid_of_gene)   # gene length, sensitivity covariate
run_goseq <- function(sig, cat_label, bias) {
  sub <- sterms[category == cat_label][gene_id %in% gid_of_gene]
  uni <- unique(sub$gene_id); hit <- intersect(sig, uni)
  if (length(hit) < 5) return(data.table())
  tsize <- sub[, .N, by = term]; sub <- sub[term %in% tsize[N >= 5 & N <= 500, term]]
  vec <- setNames(as.integer(uni %in% hit), uni)
  # STAT TEST: GOseq Wallenius non-central hypergeometric with a probability weighting
  # function fitted on the per-gene bias covariate; BH within the category.
  pwf <- nullp(vec, bias.data = pmax(bias[uni], 1), plot.fit = FALSE)
  res <- goseq(pwf, gene2cat = as.data.frame(sub[, .(gene_id, term)]),
               method = "Wallenius", use_genes_without_cat = TRUE)   # universe = every annotated gene, as in enricher
  dt <- as.data.table(res)[, .(term = category, p = over_represented_pvalue,
                               n_hit = numDEInCat, n_term = numInCat)]
  dt <- dt[n_hit >= 1]   # BH over terms with at least one hit gene, as in enricher
  dt[, fdr := p.adjust(p, "BH")]
  merge(dt, unique(sub[, .(term, description)]), by = "term")
}
# Per gene set: goseq per category with both bias covariates, merged with the enricher FDR
goseq_set <- function(sig, base, tag) {
  out <- rbindlist(lapply(names(GO_CATS), function(o) {
    a <- run_goseq(sig, GO_CATS[[o]], bias_ncpg); if (!nrow(a)) return(data.table())
    b <- run_goseq(sig, GO_CATS[[o]], bias_len)
    d <- merge(a, b[, .(term, goseq_fdr_length = fdr)], by = "term", all.x = TRUE)
    d[, ontology := o]; d
  }), fill = TRUE)
  if (!nrow(out)) return(out)
  setnames(out, c("p", "fdr"), c("goseq_p_ncpg", "goseq_fdr_ncpg"))
  if (nrow(base)) out <- merge(out, base[, .(term = ID, ontology, enricher_fdr = p.adjust)],
                               by = c("term", "ontology"), all.x = TRUE) else out[, enricher_fdr := NA_real_]
  setorder(out, goseq_fdr_ncpg)
  fwrite(out, file.path(DAT, sprintf("goseq_%s.tsv", tag)), sep = "\t")
  out
}
goseq_sets <- list(
  dmp = goseq_set(unique(dmp_genes$gene_id), dmp_enr, "dmp"),
  dmp_hyper = goseq_set(unique(dmp_genes[direction == "Hyper", gene_id]), dmp_hyper_enr, "dmp_hyper"),
  dmp_hypo  = goseq_set(unique(dmp_genes[direction == "Hypo",  gene_id]), dmp_hypo_enr, "dmp_hypo"))
if (nrow(dmr_genes)) goseq_sets$dmr <- goseq_set(unique(dmr_genes$gene_id), dmr_enr, "dmr")
# Summary: significant terms before and after the adjustment, per gene set
goseq_summary <- rbindlist(lapply(names(goseq_sets), function(nm) {
  d <- goseq_sets[[nm]]
  data.table(set = nm, n_terms_tested = nrow(d),
             enricher_sig = sum(d$enricher_fdr < 0.05, na.rm = TRUE),
             goseq_sig_ncpg = sum(d$goseq_fdr_ncpg < 0.05, na.rm = TRUE),
             goseq_sig_length = sum(d$goseq_fdr_length < 0.05, na.rm = TRUE),
             sig_in_both = sum(d$enricher_fdr < 0.05 & d$goseq_fdr_ncpg < 0.05, na.rm = TRUE),
             lost_after_adjust = sum(d$enricher_fdr < 0.05 & d$goseq_fdr_ncpg >= 0.05, na.rm = TRUE),
             best_goseq_fdr_ncpg = if (nrow(d)) min(d$goseq_fdr_ncpg, na.rm = TRUE) else NA_real_)
}))
fwrite(goseq_summary, file.path(DAT, "goseq_summary.tsv"), sep = "\t"); print(goseq_summary)

# Step 15 - fig5j three-way overlap of DMP, DMR and DE genes (supplementary)
cat("[9] DMP / DMR / DE Venn\n")
# DE set = FDR < 0.05 and |log2FC| >= 1 from 01_genome_toolkit/data/gene_de_tail.tsv
de_path <- file.path(PIPE, "01_genome_toolkit/data/gene_de_tail.tsv")
if (nrow(dmr_genes) > 0 && file.exists(de_path)) {
  suppressPackageStartupMessages(library(VennDiagram))
  de <- fread(de_path)
  de_set  <- unique(de[!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1, gene_id])
  dmp_set <- unique(dmp_genes$gene_id); dmr_set <- unique(dmr_genes$gene_id)
  sets <- list(DMP = dmp_set, DMR = dmr_set, DE = de_set)
  all3 <- Reduce(intersect, sets)
  fwrite(data.table(gene_id = all3, symbol = disp_name(all3)),
         file.path(DAT, "dmp_dmr_de_intersection.tsv"), sep = "\t")
  vp <- venn.diagram(sets, filename = NULL, disable.logging = TRUE,
    fill = c("#2C7FB8", "#6A51A3", "#B2182B"), alpha = 0.45, lwd = 1,
    cex = 1.2, fontface = "bold", cat.cex = 1.2, cat.fontface = "bold",
    margin = 0.10, main = "DMP / DMR / DE gene overlap", main.cex = 1.3)
  dv <- function() { grid::grid.newpage(); grid::grid.draw(vp) }
  pdf(file.path(FIGS, "fig5j_dmp_dmr_de_venn.pdf"), width = 6, height = 5.5); dv(); dev.off()
  png(file.path(FIGS, "fig5j_dmp_dmr_de_venn.png"), width = 6, height = 5.5, units = "in", res = 150); dv(); dev.off()
  svglite::svglite(file.path(FIGS, "fig5j_dmp_dmr_de_venn.svg"), width = 6, height = 5.5); dv(); dev.off()
  cat(sprintf("  saved fig5j | DMP %d, DMR %d, DE %d, all-three %d\n",
              length(dmp_set), length(dmr_set), length(de_set), length(all3)))
}

# Step 16 - fig5k DMP/DMR genes by biotype (supplementary)
cat("[10] gene biotype breakdown (protein_coding / lncRNA; pseudogenes excluded)\n")
# Pseudogenes are excluded from the panel; protein_coding and lncRNA are shown
kn <- data.table::melt(gene_summary[, .(set, protein_coding, lncRNA)], id.vars = "set",
           variable.name = "biotype", value.name = "n")
kn[, biotype := factor(as.character(biotype), levels = BIOTYPES)]
kn[, frac := n / sum(n), by = set]
kn[, set := factor(set, levels = c("DMP genes", "DMR genes"))]
# Dot at each biotype's percentage; classes below 1% keep one decimal
kn[, pct := 100 * frac]
bt_ord <- kn[, .(tot = sum(n)), by = biotype][order(tot), as.character(biotype)]
kn[, biotype := factor(as.character(biotype), levels = bt_ord)]
kn[, lab := sprintf("%s  (%s)", fifelse(pct < 1, sprintf("%.1f%%", pct), sprintf("%.0f%%", pct)), comma(n))]
pk <- ggplot(kn, aes(pct, biotype, colour = biotype)) +
  geom_segment(aes(xend = 0, yend = biotype), colour = "grey85", linewidth = 0.5) +
  geom_point(size = 4) +
  geom_text(aes(label = lab), hjust = -0.12, size = 2.7, colour = "grey20") +
  facet_wrap(~ set) +
  scale_colour_manual(values = COL_BT, guide = "none") +
  scale_y_discrete(labels = LAB_BT) +
  scale_x_continuous(limits = c(0, 175), expand = expansion(mult = c(0, 0))) +
  labs(title = "DMP / DMR genes by gene biotype", x = "Percentage of genes", y = NULL) +
  theme_pub() + theme(strip.background = element_blank(), panel.grid.major.y = element_blank())
save_supp(pk, "fig5k_gene_coding_pie", 6.6, 2.2)

# Step 17 - Biotype selectivity: raw vs gene-length-adjusted odds ratio (table only)
cat("[10b] biotype selectivity, raw vs gene-length-adjusted\n")
# Protein-coding genes are longer than lncRNAs; the CMH test stratifies by length quintile
bt_ann <- data.table(gene_id = gid_of_gene, biotype = biotype_all[gid_of_gene],
                     len = width(gene_gr))[biotype %in% c("protein_coding", "lncRNA")]
bt_rows <- rbindlist(lapply(c("DMP", "DMR"), function(k) {
  hit <- if (k == "DMP") unique(dmp_genes$gene_id) else unique(dmr_genes$gene_id)
  d <- copy(bt_ann)[, has := gene_id %in% hit]
  d[, lenq := cut(len, quantile(len, seq(0, 1, 0.2)), include.lowest = TRUE, labels = FALSE)]
  # STAT TEST: two-sided Fisher (raw) and Cochran-Mantel-Haenszel common OR stratified by gene-length quintile
  raw <- fisher.test(table(d$biotype, d$has))
  cmh <- mantelhaen.test(table(d$biotype, d$has, d$lenq))
  data.table(set = k, n_protein_coding = d[biotype == "protein_coding" & has, .N],
             n_lncRNA = d[biotype == "lncRNA" & has, .N],
             pct_coding = 100 * d[(has), mean(biotype == "protein_coding")],
             raw_OR = raw$estimate, raw_p = raw$p.value,
             lenadj_OR = cmh$estimate, lenadj_lo = cmh$conf.int[1],
             lenadj_hi = cmh$conf.int[2], lenadj_p = cmh$p.value)
}))
fwrite(bt_rows, file.path(DAT, "biotype_selectivity_length_adjusted.tsv"), sep = "\t")
print(bt_rows)

# Step 18 - fig5l known vs unknown genes among DMP/DMR genes (supplementary)
cat("[11] gene known vs unknown pies\n")
ku <- data.table::melt(gene_summary[, .(set, known, unknown)], id.vars = "set",
           variable.name = "status", value.name = "n")
ku[, status := factor(fifelse(status == "known", "Known", "Unknown"),
                      levels = c("Known", "Unknown"))]
ku[, frac := n / sum(n), by = set]
ku[, set := factor(set, levels = c("DMP genes", "DMR genes"))]
pl <- ggplot(ku, aes("", frac, fill = status)) +
  geom_col(width = 1, colour = "white") + coord_polar(theta = "y") +
  geom_text(aes(label = sprintf("%.0f%%\n(%s)", 100 * frac, format(n, big.mark = ","))),
            position = position_stack(vjust = 0.5), size = 2.8) +
  facet_wrap(~ set) +
  scale_fill_manual(values = c("Known" = "#1B9E9E", "Unknown" = "#BDBDBD"), name = NULL) +
  labs(title = "DMP / DMR genes: known vs unknown", x = NULL, y = NULL) +
  theme_pub() + theme(axis.line = element_blank(), axis.text = element_blank(),
                      axis.ticks = element_blank(), panel.grid = element_blank(),
                      strip.background = element_blank())
save_supp(pl, "fig5l_gene_known_unknown_pie", 6.0, 3.6)

# Step 19 - DMRs in gene-body TEs: permutation enrichment test
cat("[12] DMRs in gene-body TEs (permutation test)\n")
# Null: the same DMRs (widths kept) placed at random inside the gene-body universe.
# GRanges set operations are strand-aware, so genes are unstranded before intersecting with TEs
te <- fread(TE)[chrom %in% keep_chr]
tcf <- as.character(te$class_family)
te[, class := ifelse(grepl("^LINE", tcf), "LINE", ifelse(grepl("^SINE", tcf), "SINE",
              ifelse(grepl("^LTR", tcf), "LTR", ifelse(grepl("^DNA", tcf), "DNA",
              ifelse(grepl("^RC", tcf), "RC", "Unknown")))))]
te <- te[class %in% c("LINE", "SINE", "LTR", "DNA", "RC")]   # five scored TE classes; Unknown copies dropped
te_gr <- GRanges(te$chrom, IRanges(te$start, te$end), class = te$class, te_name = te$te_name)

gene_us <- gene_gr; strand(gene_us) <- "*"   # unstrand before reduce
gene_uni <- reduce(gene_us)   # gene-body universe
te_in_gene <- intersect(reduce(te_gr), gene_uni)   # TE bases inside gene bodies
dmr_gr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
in_gene <- overlapsAny(dmr_gr, gene_uni)

# Random placement: map a uniform draw onto the concatenated gene-body universe
gb <- as.data.table(gene_uni)[, .(seqnames = as.character(seqnames), start, w = width)]
gb_tot <- sum(as.numeric(gb$w)); gb_cum <- cumsum(as.numeric(gb$w))
place <- function(widths) {
  p <- runif(length(widths), 0, gb_tot)
  i <- findInterval(p, c(0, gb_cum[-length(gb_cum)]))
  stopifnot(!anyNA(i), all(i >= 1L), all(i <= nrow(gb)))   # stop on an invalid interval index
  GRanges(gb$seqnames[i], IRanges(gb$start[i] + as.integer(p - c(0, gb_cum)[i]), width = widths))
}
set.seed(20260426)   # reseed so the null is independent of earlier RNG use
NPERM <- 1000
# STAT TEST: permutation test. Null = NPERM random placements of same-width regions;
# p = (1 + #{null overlaps >= observed}) / (NPERM + 1), a one-sided enrichment p-value.
perm_test <- function(dmr_sub, target) {
  if (!length(dmr_sub)) return(NULL)
  obs <- sum(overlapsAny(dmr_sub, target)); w <- width(dmr_sub)
  null <- replicate(NPERM, sum(overlapsAny(place(w), target)))
  list(obs = obs, null = null, n = length(dmr_sub), exp = mean(null),
       fold = obs / max(mean(null), 1e-9), p = (1 + sum(null >= obs)) / (NPERM + 1))
}
ov_all <- perm_test(dmr_gr[in_gene], te_in_gene)
cat(sprintf("  %d/%d gene-body DMRs hit a TE (%.1f%%); TE = %.1f%% of gene-body bp\n",
            ov_all$obs, ov_all$n, 100*ov_all$obs/ov_all$n,
            100*sum(width(te_in_gene))/sum(width(gene_uni))))
cat(sprintf("  expected %.1f | observed %d | fold %.2fx | p = %.4f\n",
            ov_all$exp, ov_all$obs, ov_all$fold, ov_all$p))

# Per TE class: is the signal driven by one class?
cls <- c("DNA", "LINE", "LTR", "SINE", "RC")
enr <- rbindlist(lapply(c("All TEs", cls), function(k) {
  tgt <- if (k == "All TEs") te_in_gene else intersect(reduce(te_gr[te_gr$class == k]), gene_uni)
  r <- if (k == "All TEs") ov_all else perm_test(dmr_gr[in_gene], tgt)   # reuse the all-TE result so table and figure share one null
  data.table(class = k, n_dmr = r$n, observed = r$obs, expected = round(r$exp, 1),
             fold = round(r$fold, 2), p = r$p,
             pct_gene_bp = round(100*sum(width(tgt))/sum(width(gene_uni)), 2))
}))
enr[, fdr := p.adjust(p, "BH")]
fwrite(enr, file.path(DAT, "dmr_te_gene_body_enrichment.tsv"), sep = "\t")
print(enr)

# Every DMR x gene-body TE pair with its host gene
tg <- te_gr[overlapsAny(te_gr, gene_uni)]
h  <- findOverlaps(tg, dmr_gr)
pairs <- data.table(
  te_name = tg$te_name[queryHits(h)], te_class = tg$class[queryHits(h)],
  te_chr = as.character(seqnames(tg))[queryHits(h)],
  te_start = start(tg)[queryHits(h)], te_end = end(tg)[queryHits(h)],
  dmr_chr = dmrs$chr[subjectHits(h)], dmr_start = dmrs$start[subjectHits(h)],
  dmr_end = dmrs$end[subjectHits(h)], direction = dmrs$direction[subjectHits(h)],
  diff.Methy = dmrs$diff.Methy[subjectHits(h)], region = dmrs$region[subjectHits(h)],
  host_gene = dmrs$gene_id[subjectHits(h)])
pairs[, host_symbol := disp_name(host_gene)]
setorder(pairs, -"diff.Methy")   # descending diff.Methy: strongest Hypo pairs first
fwrite(pairs, file.path(DAT, "dmr_te_gene_body_pairs.tsv"), sep = "\t")
cat(sprintf("  %d TE-DMR pairs across %d TE copies -> dmr_te_gene_body_pairs.tsv\n",
            nrow(pairs), uniqueN(pairs[, .(te_chr, te_start)])))

# Step 20 - Controls for the TE enrichment (direction, function); fig5m
# Control 1: is the Hyper/Hypo split different for on-TE vs off-TE gene-body DMRs?
is_in_gene <- in_gene
is_on_te   <- overlapsAny(dmr_gr, te_in_gene)
dg <- copy(dmrs)[, `:=`(on_te = is_on_te, in_gene = is_in_gene)][in_gene == TRUE]
dir_tab <- rbindlist(lapply(c("All TEs", cls), function(k) {
  tgt <- if (k == "All TEs") te_in_gene else intersect(reduce(te_gr[te_gr$class == k]), gene_uni)
  d <- dg[overlapsAny(GRanges(dg$chr, IRanges(dg$start, dg$end)), tgt)]
  data.table(set = k, DMRs = nrow(d), Hyper = sum(d$direction == "Hyper"),
             Hypo = sum(d$direction == "Hypo"),
             pct_hyper = round(100*mean(d$direction == "Hyper"), 1),
             median_abs_delta = round(median(abs(d$diff.Methy)), 3))
}))
dir_tab <- rbind(dir_tab, data.table(set = "Off-TE (gene body)", DMRs = sum(!dg$on_te),
  Hyper = sum(!dg$on_te & dg$direction == "Hyper"), Hypo = sum(!dg$on_te & dg$direction == "Hypo"),
  pct_hyper = round(100*mean(dg[on_te == FALSE, direction == "Hyper"]), 1),
  median_abs_delta = round(median(abs(dg[on_te == FALSE, diff.Methy])), 3)))
fwrite(dir_tab, file.path(DAT, "dmr_te_direction.tsv"), sep = "\t")
print(dir_tab)
# STAT TEST: two-sided Fisher's exact test — is the Hyper/Hypo split associated with TE overlap?
dir_p <- fisher.test(table(dg$on_te, dg$direction))$p.value
cat(sprintf("  on-TE vs off-TE Hyper/Hypo: Fisher p = %.3g  (%s)\n", dir_p,
            if (dir_p < 0.05) "DIRECTIONAL BIAS DETECTED" else "no directional bias"))

# Control 2: GO of TE-DMR host genes against the DMR-gene universe (not all genes)
te_host <- unique(na.omit(dmr_genes[row %in% dg[on_te == TRUE, row], gene_id]))
dmr_universe <- unique(na.omit(dmr_genes$gene_id))
cat(sprintf("[12b] GO: %d TE-DMR host genes vs %d DMR genes (universe = DMR genes)\n",
            length(te_host), length(dmr_universe)))
run_enr_uni <- function(sig, cat_label, universe) {   # ORA with a custom universe
  sub <- sterms[category == cat_label]
  uni <- intersect(universe, unique(sub$gene_id))
  hit <- intersect(sig, uni)
  if (length(hit) < 5) return(data.table())
  # STAT TEST: hypergeometric over-representation (clusterProfiler::enricher) on a custom universe
  res <- suppressWarnings(enricher(gene = hit, universe = uni,
                  TERM2GENE = sub[, .(term, gene_id)],
                  TERM2NAME = unique(sub[, .(term, description)]),
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = 5, maxGSSize = 500))
  if (is.null(res)) return(data.table())
  as.data.table(as.data.frame(res))
}
te_go <- rbindlist(lapply(names(GO_CATS), function(o) {
  d <- run_enr_uni(te_host, GO_CATS[[o]], dmr_universe); if (nrow(d)) d[, ontology := o]; d
}), fill = TRUE)
fwrite(te_go, file.path(DAT, "dmr_te_go_vs_dmr_universe.tsv"), sep = "\t")
n_sig <- if (nrow(te_go)) sum(te_go$p.adjust < 0.05, na.rm = TRUE) else 0L
cat(sprintf("  %d terms tested, %d significant at FDR<0.05 -> TE-DMR genes are %s\n",
            nrow(te_go), n_sig,
            if (n_sig == 0L) "NOT functionally distinct from other DMR genes; the TE effect is POSITIONAL."
            else "FUNCTIONALLY DISTINCT from other DMR genes at FDR<0.05 -- re-read before quoting."))

# fig5m (a): permutation null vs observed overlap
nulldt <- data.table(x = ov_all$null)
p_null <- ggplot(nulldt, aes(x)) +
  geom_histogram(bins = 40, fill = "grey80", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = ov_all$obs, colour = COL_DIR["Hyper"], linewidth = 0.6) +
  annotate("text", x = ov_all$obs, y = Inf, label = sprintf(" observed = %d", ov_all$obs),
           hjust = 1.05, vjust = 1.8, size = 2.6, colour = COL_DIR["Hyper"]) +
  labs(x = "Gene-body DMRs overlapping a TE (1,000 random placements)", y = "Permutations",
       title = sprintf("%s (%.2fx, p = %.3f)",
                       if (ov_all$p < 0.05 & ov_all$fold > 1) "DMRs favour gene-body TEs"
                       else "Gene-body DMRs vs TE placement", ov_all$fold, ov_all$p)) +
  theme_pub()
# fig5m (b): fold enrichment per TE class
enr[, lab := factor(class, levels = rev(c("All TEs", cls)))]
p_cls <- ggplot(enr, aes(fold, lab, fill = fdr < 0.05)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = COL_DIR["Hyper"], "FALSE" = "grey75"),
                    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "n.s."), name = NULL) +
  labs(x = "Fold enrichment vs random placement", y = NULL,
       title = "By TE class") + theme_pub()
save_supp(p_null / p_cls + plot_layout(heights = c(1.3, 1)),
          "fig5m_dmr_te_gene_body", 6.0, 6.0)

# Step 21 - Label-swap null: the same DSS test on the two mislabelled 2-vs-2 splits
cat("[2b] label-swap null (same DSS settings on the two mislabelled 2-vs-2 splits)\n")
# Swaps pair one control with one amputated animal on each side; same CpG universe, same
# per-chromosome smoothed DMLtest and thresholds; cached per chromosome like the real contrast
run_dss <- function(g1, g2, tag) {
  rds <- file.path(OBJ, sprintf("dmltest_%s_chrmt.rds", tag))
  if (file.exists(rds) && file.mtime(rds) < BS_MTIME) stop(sprintf("%s is older than the bsseq object: delete it and rerun", basename(rds)), call. = FALSE)
  if (file.exists(rds)) return(readRDS(rds))
  chrs <- intersect(keep_chr, as.character(unique(seqnames(bs))))
  dl <- lapply(chrs, function(ch) {
    cf <- file.path(OBJ, sprintf("dmltest_%s_chrmt_%s.rds", tag, ch))
    if (file.exists(cf)) return(readRDS(cf))
    cat(sprintf("  DMLtest %s %s\n", tag, ch))
    r <- DMLtest(bs[as.character(seqnames(bs)) == ch, ], group1 = g1, group2 = g2, smoothing = TRUE)
    saveRDS(r, cf); r
  })
  d <- do.call(rbind, dl); saveRDS(d, rds); d
}
swaps <- list(real  = list(g1 = c("C1", "C2"), g2 = c("A1", "A2")),
              swapA = list(g1 = c("C1", "A1"), g2 = c("C2", "A2")),
              swapB = list(g1 = c("C1", "A2"), g2 = c("C2", "A1")))
swap_keys <- list(); swap_rows <- list()
for (nm in names(swaps)) {
  d <- if (nm == "real") dml_test else run_dss(swaps[[nm]]$g1, swaps[[nm]]$g2, nm)
  # STAT TEST: identical to [2] -- DSS callDML (p < 0.05, |delta| >= 0.10) and callDMR
  # (p < 0.05, delta 0.10, minlen 50, minCG 3) on the smoothed per-chromosome DMLtest.
  p <- as.data.table(callDML(d, p.threshold = 0.05, delta = 0.10))
  r <- callDMR(d, p.threshold = 0.05, delta = 0.10, minlen = 50, minCG = 3)
  r <- if (is.null(r)) data.table() else as.data.table(r)
  swap_keys[[nm]] <- p[, paste(chr, pos)]
  swap_rows[[nm]] <- data.table(
    contrast = nm, group1 = paste(swaps[[nm]]$g1, collapse = "+"), group2 = paste(swaps[[nm]]$g2, collapse = "+"),
    n_dmp = nrow(p), n_dmp_group2_higher = sum(p$diff < 0), n_dmp_group2_lower = sum(p$diff > 0),
    n_dmr = nrow(r), n_dmr_group2_higher = if (nrow(r)) sum(r$diff.Methy < 0) else 0L,
    n_dmr_group2_lower = if (nrow(r)) sum(r$diff.Methy > 0) else 0L)
  if (nm != "real") rm(d); invisible(gc())
}
swap_tab <- rbindlist(swap_rows)
real_key <- swap_keys$real
swap_tab[, real_dmps_recovered := vapply(contrast, function(nm) sum(real_key %in% swap_keys[[nm]]), numeric(1))]
swap_tab[, frac_real_dmps_recovered := real_dmps_recovered / length(real_key)]
swap_tab[, frac_real_dmps_in_either_swap := mean(real_key %in% c(swap_keys$swapA, swap_keys$swapB))]
fwrite(swap_tab, file.path(DAT, "label_swap_counts.tsv"), sep = "\t")
cat("  contrast  DMPs (g2 higher / lower)   DMRs (g2 higher / lower)   real DMPs recovered\n")
for (i in seq_len(nrow(swap_tab))) with(swap_tab[i], cat(sprintf(
  "  %-7s %8s (%s / %s)   %6s (%s / %s)   %s (%.1f%%)\n", contrast, format(n_dmp, big.mark = ","),
  format(n_dmp_group2_higher, big.mark = ","), format(n_dmp_group2_lower, big.mark = ","),
  format(n_dmr, big.mark = ","), n_dmr_group2_higher, n_dmr_group2_lower,
  format(real_dmps_recovered, big.mark = ","), 100 * frac_real_dmps_recovered)))
cat(sprintf("  real DMPs called by either swap: %.1f%%\n", 100 * swap_tab$frac_real_dmps_in_either_swap[1]))
# fig5n (supplementary): DMP and DMR counts per contrast, split by direction
sw_long <- data.table::melt(swap_tab[, .(contrast, `DMPs, group 2 higher` = n_dmp_group2_higher, `DMPs, group 2 lower` = n_dmp_group2_lower,
                                         `DMRs, group 2 higher` = n_dmr_group2_higher, `DMRs, group 2 lower` = n_dmr_group2_lower)],
                            id.vars = "contrast", variable.name = "what", value.name = "n")
sw_long[, `:=`(feature = ifelse(grepl("^DMP", what), "DMPs", "DMRs"),
               direction = ifelse(grepl("higher", what), "Group 2 higher", "Group 2 lower"))]
sw_long[, contrast := factor(contrast, levels = c("real", "swapA", "swapB"),
                             labels = c("C1+C2 vs A1+A2\n(true labels)", "C1+A1 vs C2+A2\n(swap A)", "C1+A2 vs C2+A1\n(swap B)"))]
p_sw <- ggplot(sw_long, aes(contrast, n, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(data = function(d) d[n > 0], aes(label = scales::comma(n)), position = position_stack(vjust = 0.5), size = 2.3, colour = "white") +
  facet_wrap(~ feature, scales = "free_y") +
  scale_fill_manual(values = c(`Group 2 higher` = COL_DIR[["Hyper"]], `Group 2 lower` = COL_DIR[["Hypo"]]), name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Count", title = "Label-swap null: the same DSS test on mislabelled 2-vs-2 splits") +
  theme_pub() + theme(strip.background = element_blank(), axis.text.x = element_text(size = 7),
                      legend.position = "bottom", plot.title = element_text(size = rel(1)))
save_supp(p_sw, "fig5n_label_swap_null", 7.2, 3.4)

# Step 22 - Argonaute loci: DMP and DMR counts per locus (table only)
ago <- fread(file.path(PIPE, "01_genome_toolkit/data/argonaute_census.tsv"))
ago[, n_dmp := vapply(gene_id, function(g) dmp_genes[gene_id == g, .N], integer(1))]
ago[, n_dmr := vapply(gene_id, function(g) if (nrow(dmr_genes)) dmr_genes[gene_id == g, .N] else 0L, integer(1))]
fwrite(ago[, .(gene_id, clade, gff_symbol, kept, baseMean, log2FoldChange, padj, n_dmp, n_dmr)],
       file.path(DAT, "argonaute_methylation.tsv"), sep = "\t")
cat(sprintf("  Argonaute loci: %s\n",
            paste(ago[kept == TRUE, sprintf("%s (%s) DMPs %d DMRs %d", gff_symbol, clade, n_dmp, n_dmr)], collapse = "; ")))

# Step 23 - Record package versions (sessionInfo)
writeLines(capture.output(sessionInfo()), file.path(BATCH, "sessionInfo_05_differential.txt"))
cat("[05_differential] done\n")
