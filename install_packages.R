## ============================================================
## Multi-omics R package installer
## Uses Tsinghua / USTC mirrors to avoid network timeouts in mainland China.
## ============================================================

# --- 1. Configure mirrors ---
options(
  repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
  BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor",
  timeout = 600
)

# --- 2. Install BiocManager ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

# --- 3. Package lists ---

# Proteomics / proteogenomics core
proteomics_pkgs <- c(
  "SummarizedExperiment",
  "GenomicRanges",
  "DESeq2",
  "vsn",
  "preprocessCore",
  "impute",
  "pcaMethods",
  "sva"
)

# Pathway enrichment
pathway_pkgs <- c(
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "ReactomePA",
  "pathview",
  "GO.db",
  "KEGGREST",
  "msigdbr"
)

# Heatmaps and visualisation
viz_pkgs <- c(
  "ComplexHeatmap",
  "circlize",
  "pheatmap",
  "ggpubr",
  "ggrepel",
  "survminer",
  "RColorBrewer",
  "viridis",
  "patchwork",
  "cowplot",
  "VennDiagram",
  "UpSetR",
  "corrplot"
)

# Metabolomics
metabo_pkgs <- c(
  "xcms",
  "MSnbase",
  "CAMERA",
  "MetaboCoreUtils",
  "Spectra"
)

# Multi-omics integration
integration_pkgs <- c(
  "MOFA2",
  "mixOmics",
  "WGCNA",
  "STRINGdb"
)

# Survival and statistics
stat_pkgs <- c(
  "survival",
  "survminer",
  "maxstat",
  "timeROC",
  "pROC",
  "glmnet",
  "randomForest",
  "ConsensusClusterPlus"
)

# --- 4. Combine and install ---
all_pkgs <- unique(c(
  proteomics_pkgs,
  pathway_pkgs,
  viz_pkgs,
  metabo_pkgs,
  integration_pkgs,
  stat_pkgs
))

cat("\n========================================\n")
cat("Installing", length(all_pkgs), "R packages\n")
cat("========================================\n\n")

BiocManager::install(all_pkgs, ask = FALSE, update = FALSE, Ncpus = 4)

# --- 5. Verify ---
cat("\n\n========== Verification ==========\n")
installed <- rownames(installed.packages())
success <- all_pkgs[all_pkgs %in% installed]
failed  <- all_pkgs[!all_pkgs %in% installed]

cat("Installed:", length(success), "/", length(all_pkgs), "\n")
if (length(failed) > 0) {
  cat("Failed packages:\n")
  for (p in failed) cat("  -", p, "\n")

  cat("\n--- Retrying failed packages ---\n")
  for (p in failed) {
    tryCatch({
      BiocManager::install(p, ask = FALSE, update = FALSE)
    }, error = function(e) {
      cat("  Final failure:", p, "-", conditionMessage(e), "\n")
    })
  }

  installed2 <- rownames(installed.packages())
  still_failed <- all_pkgs[!all_pkgs %in% installed2]
  if (length(still_failed) > 0) {
    cat("\nPackages still missing:\n")
    for (p in still_failed) cat("  -", p, "\n")
  } else {
    cat("\nAll packages installed successfully.\n")
  }
} else {
  cat("All packages installed successfully.\n")
}
