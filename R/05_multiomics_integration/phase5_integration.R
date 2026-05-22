#!/usr/bin/env Rscript
# ============================================================================
# Phase 5: Multi-omics Integration (MOFA2 + DIABLO)
# HAE Multi-omics Integration Study
# ============================================================================
# Strategy:
#   1. MOFA2: Unsupervised latent factor discovery across 3 omics
#   2. DIABLO (mixOmics): Supervised multi-omics integration (N vs Adjacent)
#   3. Cross-omics sample clustering (consensus)
#   4. Factor-clinical association analysis
# ============================================================================

# Configure reticulate to use the multiomics conda Python BEFORE loading MOFA2
library(reticulate)
use_python("/Users/rishat/miniforge3/envs/multiomics/bin/python", required = TRUE)

suppressPackageStartupMessages({
  library(MOFA2)
  library(mixOmics)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
OUT_DIR  <- file.path(PROJECT, "analysis/results/phase5_integration")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 5: Multi-omics Integration\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD & PREPARE MATCHED DATA
# ============================================================================
cat(">>> 0. Loading and matching multi-omics data\n")

# Load sample correspondence (12 patients with all 3 omics)
sample_corr <- read.csv(file.path(PROC_DIR, "sample_correspondence.csv"), check.names = FALSE)
cat(sprintf("  Complete 3-omics patients: %d\n", nrow(sample_corr)))

# Load differential results for ID->symbol mapping
degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# ---- Transcriptomics ----
tc_vst_raw <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                       check.names = FALSE, row.names = 1)
id2sym_tc <- setNames(degs$gene_name, degs$gene_id)
id2sym_tc <- id2sym_tc[!is.na(id2sym_tc) & id2sym_tc != ""]
tc_sym <- id2sym_tc[rownames(tc_vst_raw)]
keep_tc <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst_raw[keep_tc, ])
rownames(tc_mat) <- tc_sym[keep_tc]

# ---- Proteomics ----
pr_raw <- read.csv(file.path(PROC_DIR, "proteomics_log2_norm.csv"),
                   check.names = FALSE, row.names = 1)
id2sym_pr <- setNames(deps$gene_name, deps$Protein)
id2sym_pr <- id2sym_pr[!is.na(id2sym_pr) & id2sym_pr != "" & id2sym_pr != "_"]
pr_sym <- id2sym_pr[rownames(pr_raw)]
keep_pr <- !is.na(pr_sym) & pr_sym != "" & !duplicated(pr_sym)
pr_mat <- as.matrix(pr_raw[keep_pr, ])
rownames(pr_mat) <- pr_sym[keep_pr]

# ---- Metabolomics ----
met_raw <- read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"),
                    check.names = FALSE, row.names = 1)
met_mat <- as.matrix(met_raw)

cat(sprintf("  Transcriptomics: %d genes x %d samples\n", nrow(tc_mat), ncol(tc_mat)))
cat(sprintf("  Proteomics: %d proteins x %d samples\n", nrow(pr_mat), ncol(pr_mat)))
cat(sprintf("  Metabolomics: %d metabolites x %d samples\n", nrow(met_mat), ncol(met_mat)))

# ---- Harmonize sample names across omics ----
# All matrices use: Normal1..14, Adjacent1..14
shared_samples <- Reduce(intersect, list(colnames(tc_mat), colnames(pr_mat), colnames(met_mat)))
cat(sprintf("  Shared samples across all 3 omics: %d\n", length(shared_samples)))

tc_matched <- tc_mat[, shared_samples]
pr_matched <- pr_mat[, shared_samples]
met_matched <- met_mat[, shared_samples]

# Group info
group_vec <- ifelse(grepl("^Normal", shared_samples), "Normal", "Adjacent")
names(group_vec) <- shared_samples
cat(sprintf("  Groups: Normal=%d, Adjacent=%d\n",
            sum(group_vec == "Normal"), sum(group_vec == "Adjacent")))


# ============================================================================
# 1. FEATURE SELECTION FOR INTEGRATION
# ============================================================================
cat("\n>>> 1. Feature selection for integration\n")

# Select top variable features to reduce dimensionality
select_top_var <- function(mat, n = 2000) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  top <- head(order(vars, decreasing = TRUE), min(n, nrow(mat)))
  mat[top, ]
}

tc_top <- select_top_var(tc_matched, 2000)
pr_top <- select_top_var(pr_matched, 2000)
met_top <- select_top_var(met_matched, 500)

cat(sprintf("  Selected features: TC=%d, PR=%d, MET=%d\n",
            nrow(tc_top), nrow(pr_top), nrow(met_top)))


# ============================================================================
# 2. MOFA2: Multi-Omics Factor Analysis
# ============================================================================
cat("\n>>> 2. MOFA2 Analysis\n")

# Create MOFA object
mofa_data <- list(
  transcriptomics = tc_top,
  proteomics = pr_top,
  metabolomics = met_top
)

cat("  Creating MOFA object...\n")
mofa_obj <- create_mofa(mofa_data)

# Configure MOFA
data_opts <- get_default_data_options(mofa_obj)
model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 5

train_opts <- get_default_training_options(mofa_obj)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$maxiter <- 1000
train_opts$verbose <- FALSE

mofa_obj <- prepare_mofa(mofa_obj,
                          data_options = data_opts,
                          model_options = model_opts,
                          training_options = train_opts)

cat("  Training MOFA model...\n")
mofa_model <- run_mofa(mofa_obj, use_basilisk = FALSE)

cat(sprintf("  MOFA model: %d factors\n", get_dimensions(mofa_model)$K))

# Save model
saveRDS(mofa_model, file.path(OUT_DIR, "mofa2_model.rds"))

# ---- Variance explained ----
r2 <- get_variance_explained(mofa_model)
cat("  Variance explained per factor (%):\n")
for (v in names(r2$r2_per_factor)) {
  vals <- round(r2$r2_per_factor[[v]][1, ], 1)
  cat(sprintf("    %s: %s\n", v, paste(vals, collapse = ", ")))
}

# Plot variance explained
pdf(file.path(FIG_DIR, "MOFA2_variance_explained.pdf"), width = 8, height = 5)
plot_variance_explained(mofa_model, plot_total = TRUE)[[2]]
dev.off()
cat("  Variance explained plot saved.\n")

# ---- Factor values & group association ----
factors <- get_factors(mofa_model)[[1]]

# Add group info
factor_df <- as.data.frame(factors)
factor_df$sample <- rownames(factor_df)
factor_df$group <- group_vec[factor_df$sample]

# Test factor-group association (paired Wilcoxon for matched samples)
cat("\n  Factor-group associations:\n")
factor_group <- data.frame()

# Build paired sample vectors
extract_id_mofa <- function(x) as.integer(gsub("Normal|Adjacent", "", x))
adj_fac <- rownames(factors)[grepl("^Adjacent", rownames(factors))]
nor_fac <- rownames(factors)[grepl("^Normal", rownames(factors))]
adj_fac_ids <- extract_id_mofa(adj_fac)
nor_fac_ids <- extract_id_mofa(nor_fac)
paired_fac_ids <- intersect(adj_fac_ids, nor_fac_ids)
adj_fac_paired <- paste0("Adjacent", paired_fac_ids)
nor_fac_paired <- paste0("Normal", paired_fac_ids)

for (f in colnames(factors)) {
  wt <- wilcox.test(factors[adj_fac_paired, f],
                     factors[nor_fac_paired, f], paired = TRUE)
  fc <- mean(factors[group_vec == "Adjacent", f], na.rm = TRUE) -
        mean(factors[group_vec == "Normal", f], na.rm = TRUE)
  cat(sprintf("    %s: diff=%.3f, P=%.4f\n", f, fc, wt$p.value))
  factor_group <- rbind(factor_group, data.frame(
    factor = f, mean_diff = fc, pvalue = wt$p.value
  ))
}
write.csv(factor_group, file.path(OUT_DIR, "MOFA2_factor_group_assoc.csv"), row.names = FALSE)

# Factor scatter plot (Factor 1 vs Factor 2)
if (ncol(factors) >= 2) {
  p_mofa <- ggplot(factor_df, aes(Factor1, Factor2, color = group)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = sample), size = 2.5) +
    scale_color_manual(values = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")) +
    labs(title = "MOFA2 Factor Space", x = "Factor 1", y = "Factor 2") +
    theme_bw(base_size = 12)
  ggsave(file.path(FIG_DIR, "MOFA2_factor1vs2.pdf"), p_mofa, width = 8, height = 6)
  cat("  Factor scatter plot saved.\n")
}

# ---- Top weights per factor ----
for (fi in 1:min(3, ncol(factors))) {
  fname <- paste0("Factor", fi)
  for (view in c("transcriptomics", "proteomics", "metabolomics")) {
    tryCatch({
      w <- get_weights(mofa_model, views = view, factors = fi)[[1]]
      w_sorted <- sort(abs(w[,1]), decreasing = TRUE)
      top_features <- head(names(w_sorted), 20)
      cat(sprintf("    %s %s top-5: %s\n", fname, view,
                  paste(head(top_features, 5), collapse = ", ")))
    }, error = function(e) NULL)
  }
}


# ============================================================================
# 3. DIABLO: Supervised Multi-omics Integration
# ============================================================================
cat("\n>>> 3. DIABLO Analysis\n")

# Prepare data for DIABLO (samples in rows, features in columns)
X <- list(
  transcriptomics = t(tc_top),
  proteomics = t(pr_top),
  metabolomics = t(met_top)
)
Y <- factor(group_vec[shared_samples], levels = c("Normal", "Adjacent"))
names(Y) <- shared_samples

# Paired design: extract patient IDs for multilevel decomposition
extract_id <- function(x) as.integer(gsub("Normal|Adjacent", "", x))
patient_ids <- factor(extract_id(shared_samples))
names(patient_ids) <- shared_samples

cat(sprintf("  DIABLO input: %d samples, Y: Normal=%d, Adjacent=%d\n",
            length(Y), sum(Y == "Normal"), sum(Y == "Adjacent")))
cat(sprintf("  Paired patients: %d\n", nlevels(patient_ids)))

# Design matrix (moderate correlation between views)
design <- matrix(0.1, nrow = 3, ncol = 3)
diag(design) <- 0
rownames(design) <- colnames(design) <- names(X)

# ---- DIABLO keepX parameter tuning via tune.block.splsda() ----
# NOTE: multilevel parameter handles within-subject centering for paired design
cat("  Tuning DIABLO keepX parameters via repeated cross-validation...\n")
ncomp <- 2

# Grid of candidate keepX values per omics layer
test_keepX <- list(
  transcriptomics = c(5, 10, 15, 20, 30),
  proteomics      = c(5, 10, 15, 20, 30),
  metabolomics    = c(5, 10, 15, 20)
)

# With n=12 paired patients (24 samples), use leave-one-out CV for
# maximum effective test set utilisation and repeat 10 times
# for stability. folds = length(Y) triggers LOO.
set.seed(42)
tune_diablo <- tryCatch({
  tune.block.splsda(
    X, Y,
    ncomp       = ncomp,
    test.keepX  = test_keepX,
    design      = design,
    multilevel  = patient_ids,
    validation  = "Mfold",
    folds       = length(Y),          # LOO-CV
    nrepeat     = 10,
    dist        = "centroids.dist",
    progressBar = TRUE,
    seed        = 42
  )
}, error = function(e) {
  cat(sprintf("  [WARN] tune.block.splsda failed: %s\n", e$message))
  cat("  Falling back to literature-informed keepX defaults.\n")
  NULL
})

# Apply tuned or default keepX
if (!is.null(tune_diablo)) {
  keepX <- tune_diablo$choice.keepX
  cat("  Tuned keepX per component:\n")
  for (vn in names(keepX)) {
    cat(sprintf("    %s: %s\n", vn, paste(keepX[[vn]], collapse = ", ")))
  }
  # Save tuning error rate summary
  tune_summary <- data.frame(
    comp1_error = tune_diablo$error.rate[1, 1],
    comp2_error = ifelse(ncol(tune_diablo$error.rate) >= 2,
                         tune_diablo$error.rate[1, 2], NA)
  )
  write.csv(tune_summary, file.path(OUT_DIR, "DIABLO_tune_error_rate.csv"),
            row.names = FALSE)
  # Save the full tuning object for reproducibility
  saveRDS(tune_diablo, file.path(OUT_DIR, "DIABLO_tune_result.rds"))
  cat("  Tuning results saved.\n")
} else {
  # Fallback: literature-informed defaults for rare disease multi-omics
  # Rationale: 15-20 features per block per component is a conservative
  # choice for n=12 paired observations (Rohart et al., PLoS Comput Biol 2017)
  keepX <- list(
    transcriptomics = c(20, 20),
    proteomics      = c(20, 20),
    metabolomics    = c(15, 15)
  )
  cat("  Using default keepX (based on Rohart et al. 2017 guidelines):\n")
  for (vn in names(keepX)) {
    cat(sprintf("    %s: %s\n", vn, paste(keepX[[vn]], collapse = ", ")))
  }
}

# ---- Fit final DIABLO model ----
cat("  Running DIABLO (block.splsda with multilevel for paired design)...\n")

diablo_model <- block.splsda(X, Y, ncomp = ncomp, keepX = keepX,
                              design = design, multilevel = patient_ids)

cat("  DIABLO model fitted.\n")
saveRDS(diablo_model, file.path(OUT_DIR, "diablo_model.rds"))

# ---- DIABLO sample plot ----
pdf(file.path(FIG_DIR, "DIABLO_sample_plot.pdf"), width = 8, height = 6)
plotIndiv(diablo_model, comp = c(1, 2), ind.names = TRUE,
          group = Y, legend = TRUE, title = "DIABLO Sample Plot",
          ellipse = TRUE)
dev.off()
cat("  DIABLO sample plot saved.\n")

# ---- Selected features per component ----
extract_diablo_features <- function(model, comp) {
  result <- data.frame()
  for (view in names(model$X)) {
    tryCatch({
      sel <- selectVar(model, comp = comp, block = view)
      if (is.list(sel) && !is.null(sel[[view]])) {
        vals <- sel[[view]]$value
      } else if (is.data.frame(sel$value)) {
        vals <- sel$value
      } else {
        vals <- sel[[1]]$value
      }
      if (!is.null(vals) && nrow(vals) > 0) {
        vals$feature <- rownames(vals)
        vals$view <- view
        vals$component <- comp
        result <- rbind(result, vals)
      }
    }, error = function(e) {
      cat(sprintf("    [skip] %s comp%d: %s\n", view, comp, e$message))
    })
  }
  return(result)
}

diablo_features <- rbind(
  extract_diablo_features(diablo_model, 1),
  extract_diablo_features(diablo_model, 2)
)

write.csv(diablo_features, file.path(OUT_DIR, "DIABLO_selected_features.csv"), row.names = FALSE)
cat(sprintf("  DIABLO selected features: %d total across %d components\n",
            nrow(diablo_features), ncomp))

# Print top features per view
for (view in names(X)) {
  feat <- diablo_features[diablo_features$view == view & diablo_features$component == 1, ]
  if (nrow(feat) > 0) {
    val_col <- grep("value", colnames(feat), value = TRUE)[1]
    if (!is.na(val_col)) feat <- feat[order(-abs(feat[[val_col]])), ]
    cat(sprintf("    %s comp1 top-5: %s\n", view,
                paste(head(feat$feature, 5), collapse = ", ")))
  }
}

# ---- Circos plot ----
tryCatch({
  pdf(file.path(FIG_DIR, "DIABLO_circos.pdf"), width = 8, height = 8)
  circosPlot(diablo_model, cutoff = 0.5, comp = 1,
             color.blocks = c("#3C5488", "#E64B35", "#00A087"),
             color.cor = c("#B2182B", "#2166AC"))
  dev.off()
  cat("  DIABLO circos plot saved.\n")
}, error = function(e) cat("  Circos plot skipped:", e$message, "\n"))

# ---- Loading plot ----
tryCatch({
  pdf(file.path(FIG_DIR, "DIABLO_loadings_comp1.pdf"), width = 10, height = 8)
  plotLoadings(diablo_model, comp = 1, contrib = "max", method = "median")
  dev.off()
  cat("  DIABLO loading plot saved.\n")
}, error = function(e) cat("  Loading plot skipped:", e$message, "\n"))


# ============================================================================
# 4. CROSS-OMICS CONSENSUS CLUSTERING
# ============================================================================
cat("\n>>> 4. Cross-omics Sample Clustering\n")

# Use MOFA factors for consensus clustering
if (ncol(factors) >= 2) {
  # Hierarchical clustering on MOFA factors
  d <- dist(factors[, 1:min(5, ncol(factors))])
  hc <- hclust(d, method = "ward.D2")
  
  # Cut into 2 clusters
  clusters <- cutree(hc, k = 2)
  cluster_df <- data.frame(
    sample = names(clusters),
    cluster = paste0("C", clusters),
    group = group_vec[names(clusters)]
  )
  
  write.csv(cluster_df, file.path(OUT_DIR, "MOFA2_sample_clusters.csv"), row.names = FALSE)
  
  # Cross-tabulation
  ct <- table(cluster_df$cluster, cluster_df$group)
  cat("  MOFA factor clustering vs Group:\n")
  print(ct)
  
  # Test
  ft <- fisher.test(ct)
  cat(sprintf("  Fisher's exact P = %.4f\n", ft$p.value))
}


# ============================================================================
# 5. MULTI-OMICS FEATURE OVERLAP
# ============================================================================
cat("\n>>> 5. Multi-omics Feature Overlap\n")

# Overlap between DIABLO-selected genes and MOFA top-weight genes
mofa_top_tc <- tryCatch({
  w <- get_weights(mofa_model, views = "transcriptomics", factors = 1)[[1]]
  head(names(sort(abs(w[,1]), decreasing = TRUE)), 50)
}, error = function(e) character(0))

mofa_top_pr <- tryCatch({
  w <- get_weights(mofa_model, views = "proteomics", factors = 1)[[1]]
  head(names(sort(abs(w[,1]), decreasing = TRUE)), 50)
}, error = function(e) character(0))

diablo_tc <- diablo_features$feature[diablo_features$view == "transcriptomics"]
diablo_pr <- diablo_features$feature[diablo_features$view == "proteomics"]

# Cross-method overlap
if (length(mofa_top_tc) > 0) {
  overlap_tc <- intersect(mofa_top_tc, diablo_tc)
  cat(sprintf("  MOFA-DIABLO TC overlap: %d genes (MOFA top50: %d, DIABLO: %d)\n",
              length(overlap_tc), length(mofa_top_tc), length(diablo_tc)))
  if (length(overlap_tc) > 0) cat(sprintf("    Shared: %s\n", paste(overlap_tc, collapse = ", ")))
}

if (length(mofa_top_pr) > 0) {
  overlap_pr <- intersect(mofa_top_pr, diablo_pr)
  cat(sprintf("  MOFA-DIABLO PR overlap: %d proteins\n", length(overlap_pr)))
  if (length(overlap_pr) > 0) cat(sprintf("    Shared: %s\n", paste(overlap_pr, collapse = ", ")))
}

# Cross-omics: genes appearing in both TC and PR selections
tc_genes_all <- union(mofa_top_tc, diablo_tc)
pr_genes_all <- union(mofa_top_pr, diablo_pr)
cross_omics_genes <- intersect(tc_genes_all, pr_genes_all)
cat(sprintf("  Cross-omics key genes (TC & PR): %d\n", length(cross_omics_genes)))
if (length(cross_omics_genes) > 0) {
  cat(sprintf("    %s\n", paste(cross_omics_genes, collapse = ", ")))
  write.csv(data.frame(gene = cross_omics_genes), 
            file.path(OUT_DIR, "cross_omics_key_genes.csv"), row.names = FALSE)
}


# ============================================================================
# 6. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 5 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Matched samples: %d (%d Normal, %d Adjacent)\n",
            length(shared_samples), sum(group_vec == "Normal"), sum(group_vec == "Adjacent")))
cat(sprintf("\nMOFA2: %d factors extracted\n", ncol(factors)))
cat(sprintf("  Sig factor-group associations (P<0.05): %d\n",
            sum(factor_group$pvalue < 0.05)))
cat(sprintf("\nDIABLO: %d components, %d selected features\n",
            ncomp, nrow(diablo_features)))
cat(sprintf("  Cross-omics key genes: %d\n", length(cross_omics_genes)))
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 5 COMPLETE.\n")
cat("Phase 5 COMPLETE.\n")
