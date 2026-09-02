# CLAUDE.md — rules for anyone (human or assistant) editing this repository

This is the analysis code of *DNA methylation landscape during tail regeneration in the land
slug Deroceras laeve*. It is a linear pipeline of eight R scripts, one per folder, each a
self-contained module that reads raw data or the outputs of LOWER numbered modules only and
writes only into its own folder (`data/`, `objects/`, `figures/`).

1. **One script per module.** `NN_name/code/NN_name.R` produces every table, object and figure of
   that module; `NN_name.slurm` only loads modules, sets `R_LIBS` and calls that one script.
   External tools (HOMER, HMMER, DIAMOND, MAFFT, trimAl, IQ-TREE, BLAST+) are invoked from inside the
   R script with `system2()`, so the exact command line is in the script.
2. **Downstream only.** Module N may read modules 1..N-1, never a higher one.
3. **Chromosome universe.** Every data type is restricted to `chr1..chr31` plus the mitochondrial
   scaffold `HiC_scaffold_1563` (`keep_chr` at the top of each script).
4. **Fixed seed** `20260426`, set first in every script. No hidden regular expressions.
5. **Two DE definitions exist**: the paper's set is `padj < 0.05 & |log2FC| >= 1` (143 genes); the
   discovery set `padj < 0.05` (1,776). Never mix them.
6. **Gene sets come from the genome annotation** (protein coding + lncRNA), never from an expression
   filtered set, unless the script says so in a comment.
7. **Gene length is the default confound.** Any gene set enrichment is reported with its
   length-stratified (Cochran–Mantel–Haenszel) version, and that is the one interpreted.
8. **Honest statistics**: effect size and FDR together; per animal points over bars; two animals per
   condition is stated wherever a differential result is quoted.
9. **Paths.** Each script has one block of ALL-CAPS path constants at its top; change those lines to
   your storage layout and nothing else.
10. **Never edit a script while its job runs**, and never patch an output by hand: fix the script and
    rerun the module.
