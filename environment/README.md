# environment — what the scripts were run with

* **Cluster**: Fénix (LAVIS/UNAM), Slurm, one job per module (`NN_name/code/NN_name.slurm`);
  every launcher loads its modules explicitly and sets `R_LIBS` to the project library.
* **R**: 4.4.1 with Bioconductor 3.20, called by absolute path (`/opt/apps/r/4.4.1-studio/bin/Rscript`).
  Every module ends by writing `sessionInfo_<module>.txt` next to its outputs; the copies of
  the final run are in this folder (`sessionInfo_01_genome_toolkit.txt` … `sessionInfo_08_motifs.txt`)
  and are the authoritative list of package versions (`renv` was not used).
* **Environment modules** loaded by the launchers: `modules_loaded_by_launchers.txt`
  (diamond 2.1.0, mafft 7.299, trimal 1.2, iqtree3 3.0.1, samtools 1.22.1).
* **Project-local tools** (not modules): HMMER 3.4 (`tools/hmmer`), HOMER 5.1 (`tools/homer`),
  Pfam-A pressed with HMMER 3.4 (`tools/pfam`), the JASPAR2024 CORE SQLite (`tools/jaspar`),
  cached KEGG tables (`tools/kegg`, fetched 2026-07-28), STRING v12 per species term and alias
  files for the cross species promoter comparison (`tools/string`, fetched 2026-09-01 by
  `tools/string/fetch_species.sh`), and a conda environment with fibertools 0.13.0 for the 6mA
  re call. The `tools/` folder is not in this repository because of its size; each script names
  the exact file it expects, and the fetch scripts in `tools/string` and `01_genome_toolkit/dataset`
  show how they were staged.
* **No network on compute nodes**: every external resource (UniProt sequences, Cis-BP thresholds,
  NCBI assemblies and annotations, STRING files) is staged once from a login node into `dataset/` or
  `tools/`, and the scripts read it by path and stop with a clear message if it is missing.
* **Python** is not used by the analysis. The only Python in the project is the one-off that built the
  Cis-BP per family thresholds table (`01_genome_toolkit/dataset/cisbp/01_build_thresholds.py`, with
  its raw downloads and a provenance note); the analysis reads the resulting TSV by path.
