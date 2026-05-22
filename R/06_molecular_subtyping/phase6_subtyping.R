#!/usr/bin/env Rscript
# ============================================================================
# Phase 6: Molecular Subtyping via Consensus Clustering
# HAE Multi-omics Integration Study
# ============================================================================
# Strategy:
#   1. Use Adjacent (lesion) tissue samples (n=12) for subtyping
#   2. Input: MOFA2 factors + concatenated multi-omics top-variable features
#   3. ConsensusClusterPlus for K=2-6 with multiple algorithms
#   4. Optimal K selection (CDF, delta-area, PAC, silhouette)
#   5. Clinical phenotype association with identified subtypes
# ============================================================================

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(cluster)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(RColorBrewer)
  library(dplyr)
  library(readxl)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
P5_DIR   <- file.path(PROJECT, "analysis/results/phase5_integration")
OUT_DIR  <- file.path(PROJECT, "analysis/results/phase6_subtyping")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 6: Molecular Subtyping\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD & PREPARE DATA (Adjacent samples only)
# ============================================================================
cat(">>> 0. Loading multi-omics data (Adjacent samples)\n")

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

# ---- Match samples across omics ----
shared_samples <- Reduce(intersect, list(colnames(tc_mat), colnames(pr_mat), colnames(met_mat)))

# Select Adjacent samples ONLY for subtyping
adj_samples <- grep("^Adjacent", shared_samples, value = TRUE)
nor_samples <- grep("^Normal", shared_samples, value = TRUE)
cat(sprintf("  Total shared samples: %d (Normal=%d, Adjacent=%d)\n",
            length(shared_samples), length(nor_samples), length(adj_samples)))

# Extract patient IDs from Adjacent sample names
patient_ids <- gsub("^Adjacent", "", adj_samples)
cat(sprintf("  Patients for subtyping: %s\n", paste(patient_ids, collapse = ", ")))

tc_adj <- tc_mat[, adj_samples]
pr_adj <- pr_mat[, adj_samples]
met_adj <- met_mat[, adj_samples]

cat(sprintf("  Transcriptomics: %d genes x %d samples\n", nrow(tc_adj), ncol(tc_adj)))
cat(sprintf("  Proteomics: %d proteins x %d samples\n", nrow(pr_adj), ncol(pr_adj)))
cat(sprintf("  Metabolomics: %d metabolites x %d samples\n", nrow(met_adj), ncol(met_adj)))


# ============================================================================
# 1. FEATURE SELECTION FOR CLUSTERING
# ============================================================================
cat("\n>>> 1. Feature selection for clustering\n")

# Select top variable features per omics
select_top_var <- function(mat, n) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  top <- head(order(vars, decreasing = TRUE), min(n, nrow(mat)))
  mat[top, ]
}

# Use top variable features - scaled
tc_sel <- select_top_var(tc_adj, 1000)
pr_sel <- select_top_var(pr_adj, 1000)
met_sel <- select_top_var(met_adj, 300)

cat(sprintf("  Selected: TC=%d, PR=%d, MET=%d\n",
            nrow(tc_sel), nrow(pr_sel), nrow(met_sel)))

# Scale each omics independently
tc_scaled <- t(scale(t(tc_sel)))
pr_scaled <- t(scale(t(pr_sel)))
met_scaled <- t(scale(t(met_sel)))

# Handle NAs from zero-variance features
tc_scaled <- tc_scaled[complete.cases(tc_scaled), ]
pr_scaled <- pr_scaled[complete.cases(pr_scaled), ]
met_scaled <- met_scaled[complete.cases(met_scaled), ]

cat(sprintf("  After scaling/filtering: TC=%d, PR=%d, MET=%d\n",
            nrow(tc_scaled), nrow(pr_scaled), nrow(met_scaled)))


# ============================================================================
# 2. MULTI-OMICS CONSENSUS CLUSTERING
# ============================================================================
cat("\n>>> 2. Consensus Clustering\n")

# ---- Strategy A: Concatenated multi-omics matrix ----
cat("  Strategy A: Concatenated multi-omics features\n")

# Prefix feature names to distinguish omics
rownames(tc_scaled) <- paste0("TC_", rownames(tc_scaled))
rownames(pr_scaled) <- paste0("PR_", rownames(pr_scaled))
rownames(met_scaled) <- paste0("MET_", rownames(met_scaled))

concat_mat <- rbind(tc_scaled, pr_scaled, met_scaled)
cat(sprintf("  Concatenated matrix: %d features x %d samples\n", nrow(concat_mat), ncol(concat_mat)))

# Run ConsensusClusterPlus
cc_dir_A <- file.path(FIG_DIR, "consensus_concat")
dir.create(cc_dir_A, recursive = TRUE, showWarnings = FALSE)

cc_A <- ConsensusClusterPlus(
  d = concat_mat,
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "pearson",
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2",
  seed = 42,
  title = cc_dir_A,
  plot = "pdf"
)

# Calculate cluster metrics
cat("  Computing cluster quality metrics (Strategy A)...\n")
icl_A <- calcICL(cc_A, title = cc_dir_A, plot = "pdf")

# ---- Strategy B: MOFA2 factor-based clustering ----
cat("\n  Strategy B: MOFA2 factor-based clustering\n")

# Load MOFA2 model and extract factors for Adjacent samples
mofa_model <- readRDS(file.path(P5_DIR, "mofa2_model.rds"))
library(MOFA2)
all_factors <- get_factors(mofa_model)[[1]]
adj_factor_names <- intersect(adj_samples, rownames(all_factors))

if (length(adj_factor_names) >= 6) {
  mofa_factors <- all_factors[adj_factor_names, , drop = FALSE]
  cat(sprintf("  MOFA factors for Adjacent: %d samples x %d factors\n",
              nrow(mofa_factors), ncol(mofa_factors)))
  
  # Scale MOFA factors
  mofa_scaled <- t(scale(mofa_factors))
  
  cc_dir_B <- file.path(FIG_DIR, "consensus_mofa")
  dir.create(cc_dir_B, recursive = TRUE, showWarnings = FALSE)
  
  cc_B <- ConsensusClusterPlus(
    d = mofa_scaled,
    maxK = 5,
    reps = 1000,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "euclidean",
    innerLinkage = "ward.D2",
    finalLinkage = "ward.D2",
    seed = 42,
    title = cc_dir_B,
    plot = "pdf"
  )
  
  icl_B <- calcICL(cc_B, title = cc_dir_B, plot = "pdf")
} else {
  cat("  [SKIP] Not enough Adjacent samples in MOFA model\n")
  cc_B <- NULL
}

# ---- Strategy C: DIABLO component scores ----
cat("\n  Strategy C: DIABLO component-based clustering\n")

diablo_model <- readRDS(file.path(P5_DIR, "diablo_model.rds"))
library(mixOmics)

# Extract DIABLO variates for Adjacent samples
diablo_variates <- list()
for (view in names(diablo_model$variates)) {
  if (view == "Y") next
  v <- diablo_model$variates[[view]]
  # Match Adjacent samples
  adj_idx <- grep("^Adjacent", rownames(v))
  if (length(adj_idx) > 0) {
    diablo_variates[[view]] <- v[adj_idx, , drop = FALSE]
  }
}

if (length(diablo_variates) > 0) {
  # Concatenate DIABLO variates
  diablo_concat <- do.call(cbind, diablo_variates)
  cat(sprintf("  DIABLO variates for Adjacent: %d samples x %d components\n",
              nrow(diablo_concat), ncol(diablo_concat)))
  
  diablo_scaled <- t(scale(diablo_concat))
  
  cc_dir_C <- file.path(FIG_DIR, "consensus_diablo")
  dir.create(cc_dir_C, recursive = TRUE, showWarnings = FALSE)
  
  cc_C <- ConsensusClusterPlus(
    d = diablo_scaled,
    maxK = 5,
    reps = 1000,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "euclidean",
    innerLinkage = "ward.D2",
    finalLinkage = "ward.D2",
    seed = 42,
    title = cc_dir_C,
    plot = "pdf"
  )
  
  icl_C <- calcICL(cc_C, title = cc_dir_C, plot = "pdf")
} else {
  cat("  [SKIP] Could not extract DIABLO variates\n")
  cc_C <- NULL
}


# ============================================================================
# 3. OPTIMAL K SELECTION
# ============================================================================
cat("\n>>> 3. Optimal K selection\n")

# Function to compute PAC (Proportion of Ambiguous Clustering)
compute_PAC <- function(cc_result, maxK) {
  pac_vals <- numeric(maxK - 1)
  for (k in 2:maxK) {
    M <- cc_result[[k]]$consensusMatrix
    m_lower <- M[lower.tri(M)]
    pac_vals[k - 1] <- mean(m_lower > 0.1 & m_lower < 0.9)
  }
  names(pac_vals) <- paste0("K=", 2:maxK)
  return(pac_vals)
}

# Function to compute silhouette widths
compute_silhouette <- function(cc_result, maxK) {
  sil_vals <- numeric(maxK - 1)
  for (k in 2:maxK) {
    M <- cc_result[[k]]$consensusMatrix
    d <- as.dist(1 - M)
    cl <- cc_result[[k]]$consensusClass
    if (length(unique(cl)) >= 2) {
      sil <- silhouette(cl, d)
      sil_vals[k - 1] <- mean(sil[, "sil_width"])
    } else {
      sil_vals[k - 1] <- NA
    }
  }
  names(sil_vals) <- paste0("K=", 2:maxK)
  return(sil_vals)
}

# Evaluate Strategy A (concatenated)
pac_A <- compute_PAC(cc_A, 6)
sil_A <- compute_silhouette(cc_A, 6)

cat("  Strategy A (Concatenated multi-omics):\n")
cat(sprintf("    PAC: %s\n", paste(sprintf("%s=%.3f", names(pac_A), pac_A), collapse=", ")))
cat(sprintf("    Silhouette: %s\n", paste(sprintf("%s=%.3f", names(sil_A), sil_A), collapse=", ")))

# Optimal K for A: lowest PAC or highest silhouette
best_k_pac_A <- which.min(pac_A) + 1
best_k_sil_A <- which.max(sil_A) + 1
cat(sprintf("    Best K by PAC: %d (PAC=%.3f)\n", best_k_pac_A, pac_A[best_k_pac_A - 1]))
cat(sprintf("    Best K by silhouette: %d (sil=%.3f)\n", best_k_sil_A, sil_A[best_k_sil_A - 1]))

# Evaluate Strategy B (MOFA)
if (!is.null(cc_B)) {
  pac_B <- compute_PAC(cc_B, 5)
  sil_B <- compute_silhouette(cc_B, 5)
  cat("\n  Strategy B (MOFA2 factors):\n")
  cat(sprintf("    PAC: %s\n", paste(sprintf("%s=%.3f", names(pac_B), pac_B), collapse=", ")))
  cat(sprintf("    Silhouette: %s\n", paste(sprintf("%s=%.3f", names(sil_B), sil_B), collapse=", ")))
  best_k_pac_B <- which.min(pac_B) + 1
  best_k_sil_B <- which.max(sil_B) + 1
  cat(sprintf("    Best K by PAC: %d\n", best_k_pac_B))
  cat(sprintf("    Best K by silhouette: %d\n", best_k_sil_B))
}

# Evaluate Strategy C (DIABLO)
if (!is.null(cc_C)) {
  pac_C <- compute_PAC(cc_C, 5)
  sil_C <- compute_silhouette(cc_C, 5)
  cat("\n  Strategy C (DIABLO variates):\n")
  cat(sprintf("    PAC: %s\n", paste(sprintf("%s=%.3f", names(pac_C), pac_C), collapse=", ")))
  cat(sprintf("    Silhouette: %s\n", paste(sprintf("%s=%.3f", names(sil_C), sil_C), collapse=", ")))
}

# ---- Additional stability metrics: Cophenetic & Dispersion ----
cat("\n  Computing cophenetic correlation and dispersion...\n")

compute_cophenetic <- function(cc_result, maxK) {
  coph_vals <- numeric(maxK - 1)
  for (k in 2:maxK) {
    M <- cc_result[[k]]$consensusMatrix
    d_orig <- as.dist(1 - M)
    hc_cons <- hclust(d_orig, method = "average")
    d_coph <- cophenetic(hc_cons)
    coph_vals[k - 1] <- cor(d_orig, d_coph)
  }
  names(coph_vals) <- paste0("K=", 2:maxK)
  return(coph_vals)
}

compute_dispersion <- function(cc_result, maxK) {
  disp_vals <- numeric(maxK - 1)
  for (k in 2:maxK) {
    M <- cc_result[[k]]$consensusMatrix
    m_lower <- M[lower.tri(M)]
    disp_vals[k - 1] <- sum(4 * (m_lower - 0.5)^2) / length(m_lower)
  }
  names(disp_vals) <- paste0("K=", 2:maxK)
  return(disp_vals)
}

coph_A <- compute_cophenetic(cc_A, 6)
disp_A <- compute_dispersion(cc_A, 6)
cat(sprintf("  Strategy A cophenetic: %s\n",
            paste(sprintf("%s=%.3f", names(coph_A), coph_A), collapse = ", ")))
cat(sprintf("  Strategy A dispersion: %s\n",
            paste(sprintf("%s=%.3f", names(disp_A), disp_A), collapse = ", ")))

coph_B <- disp_B <- NULL
if (!is.null(cc_B)) {
  coph_B <- compute_cophenetic(cc_B, 5)
  disp_B <- compute_dispersion(cc_B, 5)
}
coph_C <- disp_C <- NULL
if (!is.null(cc_C)) {
  coph_C <- compute_cophenetic(cc_C, 5)
  disp_C <- compute_dispersion(cc_C, 5)
}

# ---- Per-sample silhouette widths for selected K values ----
cat("  Saving per-sample silhouette widths...\n")
for (k in 2:3) {
  M <- cc_A[[k]]$consensusMatrix
  d <- as.dist(1 - M)
  cl <- cc_A[[k]]$consensusClass
  if (length(unique(cl)) >= 2) {
    sil <- silhouette(cl, d)
    sil_df <- data.frame(
      sample = names(cl),
      patient_id = gsub("^Adjacent", "", names(cl)),
      cluster = sil[, "cluster"],
      neighbor = sil[, "neighbor"],
      sil_width = round(sil[, "sil_width"], 4)
    )
    sil_df <- sil_df[order(sil_df$cluster, -sil_df$sil_width), ]
    write.csv(sil_df, file.path(OUT_DIR, sprintf("silhouette_per_sample_K%d.csv", k)),
              row.names = FALSE)
    cat(sprintf("    K=%d: mean=%.3f, min=%.3f, negative=%d/%d\n",
                k, mean(sil_df$sil_width), min(sil_df$sil_width),
                sum(sil_df$sil_width < 0), nrow(sil_df)))
  }
}

# ---- Jaccard bootstrap stability (subsampling 80% x 500 repeats) ----
cat("  Computing Jaccard bootstrap cluster stability (Strategy A, K=2)...\n")
set.seed(42)
n_boot <- 500
n_samp <- ncol(concat_mat)
ref_cl_k2 <- cc_A[[2]]$consensusClass

jaccard_boot <- function(ref_cl, data_mat, k, n_boot, subsample_frac = 0.8) {
  n <- ncol(data_mat)
  n_sub <- floor(n * subsample_frac)
  ref_labels <- unique(ref_cl)
  jaccard_per_cluster <- matrix(NA, nrow = n_boot, ncol = length(ref_labels))
  colnames(jaccard_per_cluster) <- paste0("CS", sort(ref_labels))
  
  for (b in 1:n_boot) {
    idx <- sort(sample(1:n, n_sub))
    sub_mat <- data_mat[, idx]
    d_sub <- as.dist(1 - cor(sub_mat, method = "pearson"))
    hc_sub <- hclust(d_sub, method = "ward.D2")
    cl_sub <- cutree(hc_sub, k = k)
    
    # Align cluster labels to reference via max overlap
    cl_aligned <- cl_sub
    used_ref <- c()
    for (ci in sort(unique(cl_sub))) {
      overlap <- sapply(sort(ref_labels), function(ri) {
        shared_samps <- intersect(names(ref_cl), names(cl_sub))
        sum(cl_sub[shared_samps] == ci & ref_cl[shared_samps] == ri)
      })
      best_ref <- sort(ref_labels)[which.max(overlap)]
      cl_aligned[cl_sub == ci] <- best_ref
    }
    
    # Compute Jaccard per cluster
    shared_samps <- intersect(names(ref_cl), names(cl_aligned))
    for (ri in sort(ref_labels)) {
      set_ref <- shared_samps[ref_cl[shared_samps] == ri]
      set_boot <- shared_samps[cl_aligned[shared_samps] == ri]
      inter <- length(intersect(set_ref, set_boot))
      union_n <- length(union(set_ref, set_boot))
      col_idx <- which(sort(ref_labels) == ri)
      jaccard_per_cluster[b, col_idx] <- ifelse(union_n > 0, inter / union_n, NA)
    }
  }
  return(jaccard_per_cluster)
}

jaccard_result <- jaccard_boot(ref_cl_k2, concat_mat, k = 2, n_boot = n_boot)
jaccard_means <- colMeans(jaccard_result, na.rm = TRUE)
jaccard_sds   <- apply(jaccard_result, 2, sd, na.rm = TRUE)

cat(sprintf("  Jaccard bootstrap stability (K=2, %d iterations):\n", n_boot))
for (i in 1:ncol(jaccard_result)) {
  cat(sprintf("    %s: %.3f +/- %.3f\n", colnames(jaccard_result)[i],
              jaccard_means[i], jaccard_sds[i]))
}
cat(sprintf("    Overall mean Jaccard: %.3f\n", mean(jaccard_means)))

# ---- Cross-strategy ARI (comprehensive) ----
cat("\n  Cross-strategy Adjusted Rand Index:\n")
ari_results <- data.frame()

compute_ari_safe <- function(cl1, cl2) {
  shared <- intersect(names(cl1), names(cl2))
  if (length(shared) < 4) return(NA)
  if (requireNamespace("mclust", quietly = TRUE)) {
    return(mclust::adjustedRandIndex(cl1[shared], cl2[shared]))
  }
  return(NA)
}

for (k in 2:3) {
  cl_A <- cc_A[[k]]$consensusClass
  if (!is.null(cc_B) && k <= 5) {
    cl_B <- cc_B[[k]]$consensusClass
    ari_ab <- compute_ari_safe(cl_A, cl_B)
    cat(sprintf("    K=%d Concat vs MOFA: ARI=%.3f\n", k, ari_ab))
    ari_results <- rbind(ari_results, data.frame(K = k, comparison = "Concat_vs_MOFA", ARI = ari_ab))
  }
  if (!is.null(cc_C) && k <= 5) {
    cl_C <- cc_C[[k]]$consensusClass
    ari_ac <- compute_ari_safe(cl_A, cl_C)
    cat(sprintf("    K=%d Concat vs DIABLO: ARI=%.3f\n", k, ari_ac))
    ari_results <- rbind(ari_results, data.frame(K = k, comparison = "Concat_vs_DIABLO", ARI = ari_ac))
  }
  if (!is.null(cc_B) && !is.null(cc_C) && k <= 5) {
    cl_B <- cc_B[[k]]$consensusClass
    cl_C <- cc_C[[k]]$consensusClass
    ari_bc <- compute_ari_safe(cl_B, cl_C)
    cat(sprintf("    K=%d MOFA vs DIABLO: ARI=%.3f\n", k, ari_bc))
    ari_results <- rbind(ari_results, data.frame(K = k, comparison = "MOFA_vs_DIABLO", ARI = ari_bc))
  }
}
if (nrow(ari_results) > 0) {
  write.csv(ari_results, file.path(OUT_DIR, "cross_strategy_ARI.csv"), row.names = FALSE)
}

# ---- Comprehensive K selection metrics table ----
k_metrics <- data.frame(
  K = 2:6,
  PAC_concat      = pac_A,
  Silhouette_concat = sil_A,
  Cophenetic_concat = coph_A,
  Dispersion_concat = disp_A
)
if (!is.null(cc_B)) {
  k_metrics$PAC_mofa        <- c(pac_B, NA)
  k_metrics$Silhouette_mofa <- c(sil_B, NA)
  k_metrics$Cophenetic_mofa <- c(coph_B, NA)
  k_metrics$Dispersion_mofa <- c(disp_B, NA)
}
if (!is.null(cc_C)) {
  k_metrics$PAC_diablo        <- c(pac_C, NA)
  k_metrics$Silhouette_diablo <- c(sil_C, NA)
  k_metrics$Cophenetic_diablo <- c(coph_C, NA)
  k_metrics$Dispersion_diablo <- c(disp_C, NA)
}

# Add Jaccard stability for K=2 (Strategy A)
k_metrics$Jaccard_mean_concat <- NA
k_metrics$Jaccard_mean_concat[1] <- mean(jaccard_means)

write.csv(k_metrics, file.path(OUT_DIR, "K_selection_metrics_comprehensive.csv"), row.names = FALSE)
cat("\n  Comprehensive metrics table saved.\n")

# ---- K selection summary plots (4-panel) ----
pdf(file.path(FIG_DIR, "K_selection_metrics.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

# PAC plot
plot(2:6, pac_A, type = "b", pch = 16, col = "#3C5488", lwd = 2,
     xlab = "K", ylab = "PAC",
     main = "Proportion of Ambiguous Clustering",
     ylim = c(0, max(pac_A, na.rm = TRUE) * 1.2))
if (!is.null(cc_B)) {
  lines(2:5, pac_B, type = "b", pch = 17, col = "#E64B35", lwd = 2)
}
if (!is.null(cc_C)) {
  lines(2:5, pac_C, type = "b", pch = 15, col = "#00A087", lwd = 2)
}
legend("topright", legend = c("Multi-omics", "MOFA2", "DIABLO")[1:(1 + !is.null(cc_B) + !is.null(cc_C))],
       col = c("#3C5488", "#E64B35", "#00A087")[1:(1 + !is.null(cc_B) + !is.null(cc_C))],
       pch = c(16, 17, 15)[1:(1 + !is.null(cc_B) + !is.null(cc_C))], lwd = 2, cex = 0.8)

# Silhouette plot
plot(2:6, sil_A, type = "b", pch = 16, col = "#3C5488", lwd = 2,
     xlab = "K", ylab = "Mean Silhouette Width",
     main = "Silhouette Analysis",
     ylim = c(min(sil_A, na.rm = TRUE) * 0.8, 1))
if (!is.null(cc_B)) lines(2:5, sil_B, type = "b", pch = 17, col = "#E64B35", lwd = 2)
if (!is.null(cc_C)) lines(2:5, sil_C, type = "b", pch = 15, col = "#00A087", lwd = 2)

# Cophenetic plot
plot(2:6, coph_A, type = "b", pch = 16, col = "#3C5488", lwd = 2,
     xlab = "K", ylab = "Cophenetic Correlation",
     main = "Cophenetic Correlation",
     ylim = c(min(coph_A, na.rm = TRUE) * 0.8, 1))
if (!is.null(cc_B)) lines(2:5, coph_B, type = "b", pch = 17, col = "#E64B35", lwd = 2)
if (!is.null(cc_C)) lines(2:5, coph_C, type = "b", pch = 15, col = "#00A087", lwd = 2)

# Dispersion plot
plot(2:6, disp_A, type = "b", pch = 16, col = "#3C5488", lwd = 2,
     xlab = "K", ylab = "Dispersion Coefficient",
     main = "Dispersion Coefficient",
     ylim = c(0, 1))
if (!is.null(cc_B)) lines(2:5, disp_B, type = "b", pch = 17, col = "#E64B35", lwd = 2)
if (!is.null(cc_C)) lines(2:5, disp_C, type = "b", pch = 15, col = "#00A087", lwd = 2)

dev.off()
cat("  K selection plots (4-panel) saved.\n")


# ============================================================================
# 4. ASSIGN FINAL SUBTYPES
# ============================================================================
cat("\n>>> 4. Final subtype assignment\n")

# Use the multi-omics concatenated strategy (most information)
# Select K based on PAC: prefer smallest PAC with biological interpretability
# With n=12, K=2 or K=3 are most reliable

# Use K=2 as primary (most stable with small sample)
# Also evaluate K=3 for potential finer resolution
final_K <- best_k_pac_A
if (final_K > 3) final_K <- 2  # Constrain to K<=3 for n=12
cat(sprintf("  Selected final K=%d (from Strategy A - concatenated)\n", final_K))

# Assign subtypes
subtype_assignment <- cc_A[[final_K]]$consensusClass
subtype_labels <- paste0("CS", subtype_assignment)
names(subtype_labels) <- names(subtype_assignment)

# Create subtype dataframe
subtype_df <- data.frame(
  sample = names(subtype_labels),
  patient_id = gsub("^Adjacent", "", names(subtype_labels)),
  subtype = subtype_labels,
  stringsAsFactors = FALSE
)

cat("  Subtype distribution:\n")
print(table(subtype_df$subtype))

# Also compute for K=2 and K=3 for comparison
for (k in 2:min(3, 6)) {
  cl <- cc_A[[k]]$consensusClass
  cat(sprintf("  K=%d distribution: %s\n", k,
              paste(sprintf("CS%d: %d", 1:k, table(factor(cl, levels = 1:k))), collapse = ", ")))
}

# Save assignments
write.csv(subtype_df, file.path(OUT_DIR, "subtype_assignments.csv"), row.names = FALSE)

# Also save K=2 and K=3 assignments
for (k in 2:3) {
  cl <- cc_A[[k]]$consensusClass
  df_k <- data.frame(
    sample = names(cl),
    patient_id = gsub("^Adjacent", "", names(cl)),
    subtype = paste0("CS", cl)
  )
  write.csv(df_k, file.path(OUT_DIR, sprintf("subtype_K%d.csv", k)), row.names = FALSE)
}


# ============================================================================
# 5. CONSENSUS HEATMAP WITH ANNOTATIONS
# ============================================================================
cat("\n>>> 5. Consensus heatmap visualization\n")

# Consensus matrix heatmap for selected K
for (k in 2:3) {
  cons_mat <- cc_A[[k]]$consensusMatrix
  cl <- cc_A[[k]]$consensusClass
  
  # Create display names: Adjacent1 -> P1, etc.
  display_names <- gsub("Adjacent", "P", names(cl))
  rownames(cons_mat) <- colnames(cons_mat) <- display_names
  
  subtype_colors <- c("CS1" = "#3C5488", "CS2" = "#E64B35", "CS3" = "#00A087")[1:k]
  
  # Build annotation using original sample names (keys of cl) mapped to display names
  anno_df <- data.frame(
    Subtype = paste0("CS", cl),
    row.names = display_names
  )
  
  anno_colors <- list(Subtype = subtype_colors)
  
  pdf(file.path(FIG_DIR, sprintf("consensus_heatmap_K%d.pdf", k)), width = 8, height = 7)
  pheatmap(cons_mat,
           color = colorRampPalette(c("white", "#3C5488"))(100),
           annotation_col = anno_df,
           annotation_row = anno_df,
           annotation_colors = anno_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           clustering_method = "ward.D2",
           main = sprintf("Consensus Matrix (K=%d, 1000 iterations)", k),
           fontsize = 10)
  dev.off()
  cat(sprintf("  Consensus heatmap K=%d saved.\n", k))
}


# ============================================================================
# 6. SILHOUETTE PLOT FOR SELECTED K
# ============================================================================
cat("\n>>> 6. Silhouette analysis\n")

for (k in 2:3) {
  cons_mat <- cc_A[[k]]$consensusMatrix
  d <- as.dist(1 - cons_mat)
  cl <- cc_A[[k]]$consensusClass
  
  if (length(unique(cl)) >= 2) {
    sil <- silhouette(cl, d)
    
    pdf(file.path(FIG_DIR, sprintf("silhouette_K%d.pdf", k)), width = 8, height = 6)
    par(mar = c(5, 6, 4, 2))
    plot(sil, col = c("#3C5488", "#E64B35", "#00A087")[1:k],
         main = sprintf("Silhouette Plot (K=%d, mean=%.3f)", k, mean(sil[, "sil_width"])),
         border = NA)
    dev.off()
    cat(sprintf("  Silhouette K=%d: mean width=%.3f\n", k, mean(sil[, "sil_width"])))
  }
}


# ============================================================================
# 7. CLINICAL DATA PARSING & ASSOCIATION
# ============================================================================
cat("\n>>> 7. Clinical association analysis\n")

# Parse clinical data
clin_raw <- read_excel(file.path(PROJECT, "14 例肝泡型棘球蚴病患者临床信息.xlsx"), col_names = FALSE)

# Row 3 has actual column names, rows 4-17 have data
col_names <- as.character(clin_raw[3, ])
clin_data <- clin_raw[4:17, ]
colnames(clin_data) <- col_names

# Create clean clinical dataframe
clinical <- data.frame(
  patient_id = as.character(clin_data[["编号"]]),
  age = as.numeric(clin_data[["年龄"]]),
  sex = as.character(clin_data[["性别"]]),
  ethnicity = as.character(clin_data[["民族"]]),
  endemic_area = as.character(clin_data[["居住地是否为疫区"]]),
  BMI = as.numeric(clin_data[["BMI"]]),
  disease_duration = as.numeric(clin_data[["确诊时间（年）"]]),
  surgery_history = as.character(clin_data[["是否接受过手术"]]),
  ABZ_treatment = as.character(clin_data[["阿苯达唑用药史"]]),
  lesion_size = as.numeric(clin_data[["大小（cm2）"]]),
  lesion_count = as.numeric(clin_data[["数目"]]),
  lesion_location = as.character(clin_data[["位置"]]),
  stringsAsFactors = FALSE
)

# PNM staging (columns 22-24 based on row 3)
pnm_cols <- which(col_names %in% c("病灶", "临近器官", "转移病灶"))
if (length(pnm_cols) >= 3) {
  # PNM columns: P (col 22), N (col 23), M (col 24)
  clinical$PNM_P <- as.character(clin_data[[pnm_cols[1]]])
  clinical$PNM_N <- as.character(clin_data[[pnm_cols[2]]])
  clinical$PNM_M <- as.character(clin_data[[pnm_cols[3]]])
}

# PIVM staging (columns 25-28)
pivm_cols <- which(col_names %in% c("病灶", "侵犯胆道", "血管侵犯", "转移病灶"))
if (length(pivm_cols) >= 4) {
  clinical$PIVM_stage <- as.character(clin_data[[pivm_cols[2] - 1]])  # PIVM stage
  clinical$bile_duct_invasion <- as.character(clin_data[[pivm_cols[2]]])
  clinical$vascular_invasion <- as.character(clin_data[[pivm_cols[3]]])
  clinical$PIVM_M <- as.character(clin_data[[pivm_cols[4]]])
}

# Lab values
clinical$WBC <- as.numeric(clin_data[["白细胞"]])
clinical$neutrophil_pct <- as.numeric(clin_data[["中性粒细胞（%）"]])
clinical$lymphocyte_pct <- as.numeric(clin_data[["淋巴细胞（%）"]])
clinical$ALT <- as.numeric(clin_data[["ALT(U/L)"]])
clinical$AST <- as.numeric(clin_data[["AST(U/L)"]])
clinical$ALP <- as.numeric(clin_data[["ALP(U)"]])
clinical$GGT <- as.numeric(clin_data[["GGT(U/L)"]])
clinical$total_bilirubin <- as.numeric(clin_data[["总胆红素（umol/L）"]])
clinical$albumin <- as.numeric(clin_data[["白蛋白(g/L)"]])
clinical$glucose <- as.numeric(clin_data[["葡萄糖(mmol/L)"]])
clinical$lactate <- as.numeric(clin_data[["乳酸(mmol/L)"]])
clinical$cholesterol <- as.numeric(clin_data[["胆固醇（mmol/L）"]])

# CRP and IL-6 (may contain text like "<1.5")
crp_raw <- as.character(clin_data[["CRP（mg/L)"]])
clinical$CRP <- suppressWarnings(as.numeric(gsub("[<>]", "", crp_raw)))
il6_raw <- as.character(clin_data[["IL-6（pg/ml）"]])
clinical$IL6 <- suppressWarnings(as.numeric(gsub("[<>]", "", il6_raw)))

cat(sprintf("  Parsed clinical data: %d patients x %d variables\n", nrow(clinical), ncol(clinical)))

# Save parsed clinical
write.csv(clinical, file.path(OUT_DIR, "clinical_parsed.csv"), row.names = FALSE)

# ---- Merge subtypes with clinical ----
# Use final K assignment
subtype_k2 <- read.csv(file.path(OUT_DIR, "subtype_K2.csv"), stringsAsFactors = FALSE)
subtype_k3 <- read.csv(file.path(OUT_DIR, "subtype_K3.csv"), stringsAsFactors = FALSE)

clin_sub <- merge(clinical, subtype_k2[, c("patient_id", "subtype")],
                  by = "patient_id", all.x = TRUE)
colnames(clin_sub)[colnames(clin_sub) == "subtype"] <- "subtype_K2"

clin_sub <- merge(clin_sub, subtype_k3[, c("patient_id", "subtype")],
                  by = "patient_id", all.x = TRUE)
colnames(clin_sub)[colnames(clin_sub) == "subtype"] <- "subtype_K3"

# Only keep patients with subtype assignment
clin_typed <- clin_sub[!is.na(clin_sub$subtype_K2), ]
cat(sprintf("  Patients with subtype + clinical: %d\n", nrow(clin_typed)))
write.csv(clin_typed, file.path(OUT_DIR, "clinical_with_subtypes.csv"), row.names = FALSE)


# ============================================================================
# 8. CLINICAL-SUBTYPE ASSOCIATION TESTS
# ============================================================================
cat("\n>>> 8. Clinical-subtype association tests\n")

# For each K, test association with clinical variables
run_association_tests <- function(clin_df, subtype_col) {
  results <- data.frame()
  st <- clin_df[[subtype_col]]
  
  if (length(unique(st)) < 2) {
    cat("    Only one subtype - skipping tests\n")
    return(results)
  }
  
  # Continuous variables: Wilcoxon rank-sum / Kruskal-Wallis
  cont_vars <- c("age", "BMI", "disease_duration", "lesion_size", "lesion_count",
                  "WBC", "neutrophil_pct", "lymphocyte_pct", "ALT", "AST", "ALP",
                  "GGT", "total_bilirubin", "albumin", "CRP", "IL6",
                  "glucose", "lactate", "cholesterol")
  
  for (v in cont_vars) {
    if (!v %in% colnames(clin_df)) next
    vals <- as.numeric(clin_df[[v]])
    if (sum(!is.na(vals)) < 4) next
    
    tryCatch({
      if (length(unique(st)) == 2) {
        test <- wilcox.test(vals ~ st)
        test_name <- "Wilcoxon"
      } else {
        test <- kruskal.test(vals ~ factor(st))
        test_name <- "Kruskal-Wallis"
      }
      
      # Compute per-subtype means
      means <- tapply(vals, st, mean, na.rm = TRUE)
      
      results <- rbind(results, data.frame(
        variable = v,
        type = "continuous",
        test = test_name,
        pvalue = test$p.value,
        subtype_means = paste(sprintf("%s=%.2f", names(means), means), collapse = "; "),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  
  # Categorical variables: Fisher's exact test
  cat_vars <- c("sex", "ethnicity", "endemic_area", "surgery_history",
                "ABZ_treatment", "lesion_location",
                "PNM_P", "PNM_N", "PNM_M",
                "bile_duct_invasion", "vascular_invasion", "PIVM_stage", "PIVM_M")
  
  for (v in cat_vars) {
    if (!v %in% colnames(clin_df)) next
    vals <- clin_df[[v]]
    if (sum(!is.na(vals)) < 4 || length(unique(vals[!is.na(vals)])) < 2) next
    
    tryCatch({
      tab <- table(st, vals)
      ft <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
      
      results <- rbind(results, data.frame(
        variable = v,
        type = "categorical",
        test = "Fisher",
        pvalue = ft$p.value,
        subtype_means = paste(capture.output(print(tab)), collapse = " | "),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  
  results <- results[order(results$pvalue), ]
  return(results)
}

# Test K=2
cat("  Testing K=2 associations:\n")
assoc_k2 <- run_association_tests(clin_typed, "subtype_K2")
if (nrow(assoc_k2) > 0) {
  sig_k2 <- assoc_k2[assoc_k2$pvalue < 0.05, ]
  cat(sprintf("    Total tests: %d, Significant (P<0.05): %d\n", nrow(assoc_k2), nrow(sig_k2)))
  if (nrow(sig_k2) > 0) {
    for (i in 1:nrow(sig_k2)) {
      cat(sprintf("    * %s: P=%.4f (%s)\n", sig_k2$variable[i], sig_k2$pvalue[i], sig_k2$subtype_means[i]))
    }
  }
  # Show top 10
  cat("    Top 10 by P-value:\n")
  top10 <- head(assoc_k2, 10)
  for (i in 1:nrow(top10)) {
    cat(sprintf("      %s: P=%.4f [%s]\n", top10$variable[i], top10$pvalue[i], top10$type[i]))
  }
  write.csv(assoc_k2, file.path(OUT_DIR, "clinical_association_K2.csv"), row.names = FALSE)
}

# Test K=3
cat("\n  Testing K=3 associations:\n")
assoc_k3 <- run_association_tests(clin_typed, "subtype_K3")
if (nrow(assoc_k3) > 0) {
  sig_k3 <- assoc_k3[assoc_k3$pvalue < 0.05, ]
  cat(sprintf("    Total tests: %d, Significant (P<0.05): %d\n", nrow(assoc_k3), nrow(sig_k3)))
  if (nrow(sig_k3) > 0) {
    for (i in 1:nrow(sig_k3)) {
      cat(sprintf("    * %s: P=%.4f (%s)\n", sig_k3$variable[i], sig_k3$pvalue[i], sig_k3$subtype_means[i]))
    }
  }
  write.csv(assoc_k3, file.path(OUT_DIR, "clinical_association_K3.csv"), row.names = FALSE)
}


# ============================================================================
# 9. SUBTYPE MOLECULAR CHARACTERIZATION OVERVIEW
# ============================================================================
cat("\n>>> 9. Subtype molecular overview\n")

# For each K, show per-subtype differentially expressed features
# Use the selected K
selected_K <- final_K
cl <- cc_A[[selected_K]]$consensusClass
subtype_vec <- paste0("CS", cl)
names(subtype_vec) <- names(cl)

# ---- Per-subtype mean expression heatmap ----
# Use the top-variable concatenated matrix
# Compute per-subtype means
subtypes_unique <- sort(unique(subtype_vec))
mean_mat <- sapply(subtypes_unique, function(s) {
  samps <- names(subtype_vec)[subtype_vec == s]
  rowMeans(concat_mat[, samps, drop = FALSE], na.rm = TRUE)
})

# Select top features by variance across subtypes
feat_var <- apply(mean_mat, 1, var)
top_feats <- head(names(sort(feat_var, decreasing = TRUE)), 50)

# Heatmap of subtype profiles
pdf(file.path(FIG_DIR, sprintf("subtype_profile_K%d.pdf", selected_K)), width = 10, height = 12)
mat_plot <- concat_mat[top_feats, ]
# Order samples by subtype
samp_order <- names(sort(subtype_vec))
mat_plot <- mat_plot[, samp_order]

anno_col <- data.frame(
  Subtype = subtype_vec[samp_order]
)
rownames(anno_col) <- samp_order

# Feature omics annotation
feat_omics <- ifelse(grepl("^TC_", top_feats), "Transcriptomics",
               ifelse(grepl("^PR_", top_feats), "Proteomics", "Metabolomics"))
anno_row <- data.frame(Omics = feat_omics)
rownames(anno_row) <- top_feats

subtype_colors <- c("CS1" = "#3C5488", "CS2" = "#E64B35", "CS3" = "#00A087")[1:selected_K]
omics_colors <- c("Transcriptomics" = "#4DBBD5", "Proteomics" = "#E64B35", "Metabolomics" = "#00A087")

pheatmap(mat_plot,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(Subtype = subtype_colors, Omics = omics_colors),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         main = sprintf("Top 50 Discriminating Features (K=%d)", selected_K),
         breaks = seq(-3, 3, length.out = 101))
dev.off()
cat("  Subtype profile heatmap saved.\n")


# ============================================================================
# 10. SUBTYPE-SPECIFIC DIFFERENTIAL FEATURES
# ============================================================================
cat("\n>>> 10. Subtype-specific differential features\n")

# For K=2: one-vs-other comparison using Wilcoxon on the original data
for (k in 2:min(3, 6)) {
  cl <- cc_A[[k]]$consensusClass
  st <- paste0("CS", cl)
  names(st) <- names(cl)
  
  cat(sprintf("\n  K=%d subtype-specific features:\n", k))
  
  all_diff <- data.frame()
  
  for (s in sort(unique(st))) {
    is_s <- st == s
    other <- !is_s
    
    if (sum(is_s) < 2 || sum(other) < 2) next
    
    # Test on each omics separately
    for (omics_name in c("TC", "PR", "MET")) {
      if (omics_name == "TC") mat_test <- tc_adj
      else if (omics_name == "PR") mat_test <- pr_adj
      else mat_test <- met_adj
      
      # Ensure sample alignment
      samps_s <- names(st)[is_s]
      samps_o <- names(st)[other]
      samps_s <- intersect(samps_s, colnames(mat_test))
      samps_o <- intersect(samps_o, colnames(mat_test))
      
      if (length(samps_s) < 2 || length(samps_o) < 2) next
      
      pvals <- apply(mat_test, 1, function(x) {
        tryCatch(wilcox.test(x[samps_s], x[samps_o])$p.value, error = function(e) NA)
      })
      
      fc <- rowMeans(mat_test[, samps_s, drop = FALSE], na.rm = TRUE) -
            rowMeans(mat_test[, samps_o, drop = FALSE], na.rm = TRUE)
      
      diff_df <- data.frame(
        feature = names(pvals),
        omics = omics_name,
        subtype = s,
        log2FC = fc,
        pvalue = pvals,
        padj = p.adjust(pvals, method = "BH"),
        stringsAsFactors = FALSE
      )
      diff_df <- diff_df[!is.na(diff_df$pvalue), ]
      all_diff <- rbind(all_diff, diff_df)
    }
    
    sig_feats <- all_diff[all_diff$subtype == s & all_diff$pvalue < 0.05, ]
    cat(sprintf("    %s: %d features P<0.05 (TC=%d, PR=%d, MET=%d)\n",
                s, nrow(sig_feats),
                sum(sig_feats$omics == "TC"),
                sum(sig_feats$omics == "PR"),
                sum(sig_feats$omics == "MET")))
  }
  
  if (nrow(all_diff) > 0) {
    write.csv(all_diff, file.path(OUT_DIR, sprintf("subtype_diff_features_K%d.csv", k)),
              row.names = FALSE)
  }
}


# ============================================================================
# 11. CROSS-STRATEGY CONCORDANCE
# ============================================================================
cat("\n>>> 11. Cross-strategy concordance\n")

# Compare subtype assignments across strategies
concordance <- data.frame(sample = adj_samples, stringsAsFactors = FALSE)
concordance$concat_K2 <- paste0("CS", cc_A[[2]]$consensusClass[adj_samples])

if (!is.null(cc_B)) {
  mofa_cl <- cc_B[[2]]$consensusClass
  concordance$mofa_K2 <- paste0("CS", mofa_cl[intersect(adj_samples, names(mofa_cl))])
}

if (!is.null(cc_C)) {
  diablo_cl <- cc_C[[2]]$consensusClass
  concordance$diablo_K2 <- paste0("CS", diablo_cl[intersect(adj_samples, names(diablo_cl))])
}

write.csv(concordance, file.path(OUT_DIR, "subtype_cross_strategy.csv"), row.names = FALSE)

# ARI (Adjusted Rand Index) between strategies
if (!is.null(cc_B) && requireNamespace("mclust", quietly = TRUE)) {
  tryCatch({
    library(mclust)
    shared <- intersect(names(cc_A[[2]]$consensusClass), names(cc_B[[2]]$consensusClass))
    if (length(shared) > 4) {
      ari_AB <- adjustedRandIndex(cc_A[[2]]$consensusClass[shared], cc_B[[2]]$consensusClass[shared])
      cat(sprintf("  ARI (Concat vs MOFA): %.3f\n", ari_AB))
    }
  }, error = function(e) cat("  ARI computation skipped\n"))
}


# ============================================================================
# 12. COMPREHENSIVE SUBTYPE VISUALIZATION
# ============================================================================
cat("\n>>> 12. Comprehensive subtype visualization\n")

# ---- PCA colored by subtype ----
pca_mat <- t(concat_mat)
pca_res <- prcomp(pca_mat, scale. = FALSE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample <- rownames(pca_df)
pca_df$patient <- gsub("^Adjacent", "P", pca_df$sample)

# Add subtypes for K=2 and K=3
cl2 <- cc_A[[2]]$consensusClass
pca_df$subtype_K2 <- paste0("CS", cl2[pca_df$sample])
cl3 <- cc_A[[3]]$consensusClass
pca_df$subtype_K3 <- paste0("CS", cl3[pca_df$sample])

var_explained <- summary(pca_res)$importance[2, 1:2] * 100

p1 <- ggplot(pca_df, aes(PC1, PC2, color = subtype_K2)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = patient), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("CS1" = "#3C5488", "CS2" = "#E64B35")) +
  labs(title = "Multi-omics PCA (K=2)",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       color = "Subtype") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

p2 <- ggplot(pca_df, aes(PC1, PC2, color = subtype_K3)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = patient), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("CS1" = "#3C5488", "CS2" = "#E64B35", "CS3" = "#00A087")) +
  labs(title = "Multi-omics PCA (K=3)",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       color = "Subtype") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

pdf(file.path(FIG_DIR, "subtype_PCA.pdf"), width = 14, height = 6)
ggarrange(p1, p2, ncol = 2, common.legend = FALSE)
dev.off()
cat("  Subtype PCA plot saved.\n")

# ---- Clinical association barplot (for K=2) ----
# Visualize key clinical variables by subtype
if (nrow(clin_typed) > 0 && "subtype_K2" %in% colnames(clin_typed)) {
  
  # Continuous variables boxplots
  cont_to_plot <- c("age", "BMI", "lesion_size", "ALT", "AST", "ALP", "GGT",
                     "WBC", "neutrophil_pct", "lymphocyte_pct", "CRP", "albumin")
  cont_to_plot <- intersect(cont_to_plot, colnames(clin_typed))
  
  plot_list <- list()
  for (v in cont_to_plot) {
    df_plot <- clin_typed[, c(v, "subtype_K2")]
    df_plot[[v]] <- as.numeric(df_plot[[v]])
    df_plot <- df_plot[!is.na(df_plot[[v]]), ]
    
    if (nrow(df_plot) < 4) next
    
    p <- ggplot(df_plot, aes(x = subtype_K2, y = .data[[v]], fill = subtype_K2)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      scale_fill_manual(values = c("CS1" = "#3C5488", "CS2" = "#E64B35", "CS3" = "#00A087")) +
      labs(title = v, x = "", y = v) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none")
    
    # Add p-value
    tryCatch({
      pval <- wilcox.test(as.numeric(df_plot[[v]]) ~ df_plot$subtype_K2)$p.value
      p <- p + annotate("text", x = 1.5, y = max(df_plot[[v]], na.rm = TRUE),
                        label = sprintf("P=%.3f", pval), size = 3)
    }, error = function(e) NULL)
    
    plot_list[[v]] <- p
  }
  
  if (length(plot_list) >= 4) {
    pdf(file.path(FIG_DIR, "clinical_boxplots_K2.pdf"), width = 16, height = 12)
    print(ggarrange(plotlist = plot_list, ncol = 4, nrow = 3))
    dev.off()
    cat("  Clinical boxplots saved.\n")
  }
  
  # Categorical variables: stacked barplot
  cat_to_plot <- c("sex", "endemic_area", "vascular_invasion", "bile_duct_invasion")
  cat_to_plot <- intersect(cat_to_plot, colnames(clin_typed))
  
  if (length(cat_to_plot) > 0) {
    pdf(file.path(FIG_DIR, "clinical_categorical_K2.pdf"), width = 12, height = 8)
    par(mfrow = c(2, 2))
    for (v in cat_to_plot) {
      vals <- clin_typed[[v]]
      st <- clin_typed$subtype_K2
      if (sum(!is.na(vals)) < 4) next
      tab <- table(st, vals)
      barplot(t(tab), beside = TRUE, col = brewer.pal(max(3, ncol(tab)), "Set2"),
              main = v, xlab = "Subtype", ylab = "Count",
              legend.text = colnames(tab), args.legend = list(cex = 0.8))
    }
    dev.off()
    cat("  Clinical categorical plots saved.\n")
  }
}


# ============================================================================
# 13. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 6 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Samples for subtyping: %d Adjacent tissue samples\n", length(adj_samples)))
cat(sprintf("Features used: TC=%d + PR=%d + MET=%d = %d total\n",
            nrow(tc_scaled), nrow(pr_scaled), nrow(met_scaled), nrow(concat_mat)))
cat(sprintf("\nConsensus clustering strategies:\n"))
cat(sprintf("  A. Concatenated multi-omics (K=2-6)\n"))
cat(sprintf("  B. MOFA2 factors (K=2-5)\n"))
cat(sprintf("  C. DIABLO variates (K=2-5)\n"))
cat(sprintf("\nOptimal K selection (Strategy A):\n"))
cat(sprintf("  Best by PAC: K=%d (PAC=%.3f)\n", best_k_pac_A, pac_A[best_k_pac_A - 1]))
cat(sprintf("  Best by silhouette: K=%d (sil=%.3f)\n", best_k_sil_A, sil_A[best_k_sil_A - 1]))
cat(sprintf("  Cophenetic (K=2): %.3f\n", coph_A[1]))
cat(sprintf("  Dispersion (K=2): %.3f\n", disp_A[1]))
cat(sprintf("  Jaccard bootstrap stability (K=2): %.3f\n", mean(jaccard_means)))

# Final subtype distribution
for (k in 2:3) {
  cl <- cc_A[[k]]$consensusClass
  cat(sprintf("\nK=%d distribution: %s\n", k,
              paste(sprintf("CS%d=%d", 1:k, table(factor(cl, levels = 1:k))), collapse = ", ")))
}

cat(sprintf("\nClinical associations (K=2): %d tests, %d significant (P<0.05)\n",
            nrow(assoc_k2), sum(assoc_k2$pvalue < 0.05)))
if (nrow(assoc_k2) > 0 && sum(assoc_k2$pvalue < 0.05) > 0) {
  sig <- assoc_k2[assoc_k2$pvalue < 0.05, ]
  for (i in 1:nrow(sig)) {
    cat(sprintf("  * %s: P=%.4f\n", sig$variable[i], sig$pvalue[i]))
  }
}

cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 6 COMPLETE.\n")
