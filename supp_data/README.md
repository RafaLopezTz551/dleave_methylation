# Supplementary data (Data S1–S5)

Assembled 2026-10-08, updated 2026-10-09 (Data S2 and S3 sheets added) from the methylation pipeline tables on Fenix (meth_paper/main/methylation_pipeline); files copied unchanged.

Data_S1_human_promoter_class_GO_KEGG.tsv
      GO (BP, MF, CC) and KEGG overrepresentation for human promoter CpG classes (HCP, ICP, LCP), complete ranked term list.
      Source: 03_promoters/data/human_weber_class_go_enrichment.tsv

Data_S2_WGCNA_sweep_and_length_adjusted_enrichment.xlsx
      sheet variance_quantile_sweep: variance-filter quantile, genes kept, best scale-free R2 and soft-threshold power (NA = R2 >= 0.85 not reached at power <= 20).
        Source: 07_wgcna/data/variance_filter_sweep.tsv
      sheet module_DMP_length_adjusted: per module Fisher (OR, lo, hi, p, fdr) and gene-length-stratified Cochran-Mantel-Haenszel (OR_adj ... fdr_adj; sensitivity analysis).
        Source: 07_wgcna/data/module_dmp_enrichment_length_adjusted.tsv
      sheet module_DMP_CMH_strata: Cochran-Mantel-Haenszel module tests stratified by gene length, by gene length x baseline gene-body methylation (5 x 5), and by movable CpGs per gene.
        Source: 07_wgcna/data/module_dmp_enrichment_length_meth_adjusted.tsv
      sheet module_rate_per_callable_CpG: negative binomial rate ratio of DMP events (DMPs within 500 bp merged) per callable CpG (control beta 0.10-0.80), module vs rest, true contrast and swap A (swapA_ columns); paper Fisher OR alongside.
        Source: analysis/dmp_rate_callable/data/module_dmp_rate_ratio_callable.tsv; analysis/swapA_calls/data/swapA_module_summary_14.tsv (exploratory analysis folder, not yet a pipeline module)

Data_S3_DMP_DMR_GO_KEGG_by_direction.xlsx
      GO and KEGG overrepresentation of DMP- and DMR-associated genes: all, hypermethylated, hypomethylated (6 sheets).
      Source: 05_differential/data/{dmp,dmp_hyper,dmp_hypo,dmr,dmr_hyper,dmr_hypo}_go_kegg_enrichment.tsv
      sheet dmp_goseq_bias_corrected: goseq (Wallenius) re-test of DMP genes with tested CpGs per gene (goseq_fdr_ncpg) and gene length (goseq_fdr_length) as bias covariates.
      Source: 05_differential/data/goseq_dmp.tsv

Data_S4_JASPAR_motif_to_Dlaeve_assignment.xlsx
      sheet jaspar_ortholog_bridge: JASPAR motif to D. laeve locus mapping (01_genome_toolkit/data/jaspar_ortholog_bridge.tsv)
      sheet tf_dlaeve_guide: per motif evidence type, DBD identity and locus (01_genome_toolkit/data/tf_dlaeve_guide.tsv)
      sheet dbd_identity: DNA binding domain identity per bait/locus pair (01_genome_toolkit/data/dbd_identity.tsv)
      sheet orthology_funnel_gates: alternative gates and library sizes (01_genome_toolkit/data/tf_orthology_funnel.tsv)
      sheet tree_coortholog_clades: tree placements of co-orthologs (01_genome_toolkit/data/tree_coortholog_clades.tsv)
      sheet cisbp_family_thresholds: per family Cis-BP transfer threshold (01_genome_toolkit/objects/cisbp_family_thresholds.tsv)

Data_S5_promoter_class_methylation.tsv
      Per gene promoter CpG class, CpG O/E, covered CpGs, mean beta, methylation state and DESeq2 baseMean (23,216 protein-coding promoters).
      Source: doc/overleaf_package/supp_data/Data_S5_promoter_class_methylation.tsv (same file as in Overleaf)
