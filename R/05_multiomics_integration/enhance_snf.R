#!/usr/bin/env Rscript
# ============================================================================
# Enhancement: Similarity Network Fusion (SNF) Independent Integration
# ============================================================================
# PURPOSE: Validate multi-omics integration findings using a method
#   completely independent of MOFA2 and DIABLO. SNF constructs patient
#   similarity networks from each omics layer, fuses them, and clusters
#   the fused network. Agreement with MOFA2/DIABLO/ConsensusClusterPlus
#   results demonstrates method-independence of conclusions.
#
# Why this matters for reviewers:
#   "Results depend on MOFA2/DIABLO choice" → refuted by SNF concordance
# ============================================================================

suppressPackageStartupMessages({
  library(SNFtool)
  library(cluster)
  library(mclust)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
P5_DIR   <- file.path(PROJECT, "analysis/results/phase5_integration")
P6_DIR   <- file.path(PROJECT, "analysis/results/phase6_subtyping")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhance_snf")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("SNF Independent Multi-omics Integration\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD & PREPARE DATA
# ============================================================================
cat(">>> 0. Loading data\n")

degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# TC
tc_vst <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                    check.names = FALSE, row.names = 1)
id2sym <- setNames(degs$gene_name, degs$gene_id)
id2sym <- id2sym[!is.na(id2sym) & id2sym != ""]
tc_sym <- id2sym[rownames(tc_vst)]
keep <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst[keep, ]); rownames(tc_mat) <- tc_sym[keep]

# PR
pr_raw <- read.csv(file.path(PROC_DIR, "proteomics_log2_norm.csv"),
                    check.names = FALSE, row.names = 1)
id2sym_pr <- setNames(deps$gene_name, deps$Protein)
id2sym_pr <- id2sym_pr[!is.na(id2sym_pr) & id2sym_pr != "" & id2sym_pr != "_"]
pr_sym <- id2sym_pr[rownames(pr_raw)]
keep_pr <- !is.na(pr_sym) & pr_sym != "" & !duplicated(pr_sym)
pr_mat <- as.matrix(pr_raw[keep_pr, ]); rownames(pr_mat) <- pr_sym[keep_pr]

# MET
met_mat <- as.matrix(read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"),
                               check.names = FALSE, row.names = 1))

# Match samples
shared <- Reduce(intersect, list(colnames(tc_mat), colnames(pr_mat), colnames(met_mat)))
group_vec <- ifelse(grepl("^Normal", shared), "Normal", "Adjacent")
names(group_vec) <- shared

cat(sprintf("  Shared samples: %d (Normal=%d, Adjacent=%d)\n",
            length(shared), sum(group_vec == "Normal"), sum(group_vec == "Adjacent")))

# Select top variable features per omics
select_top <- function(mat, samps, n) {
  m <- mat[, samps]
  v <- apply(m, 1, var, na.rm = TRUE)
  m[head(order(v, decreasing = TRUE), min(n, nrow(m))), ]
}

tc_top  <- select_top(tc_mat, shared, 1000)
pr_top  <- select_top(pr_mat, shared, 1000)
met_top <- select_top(met_mat, shared, 300)


# ============================================================================
# 1. SNF: ALL SAMPLES (Normal + Adjacent)
# ============================================================================
cat("\n>>> 1. SNF on all samples (N vs A classification)\n")

# SNF requires: samples in rows, features in columns
# Standardize each omics independently
tc_snf <- t(scale(t(tc_top))); tc_snf <- tc_snf[complete.cases(tc_snf), ]
pr_snf <- t(scale(t(pr_top))); pr_snf <- pr_snf[complete.cases(pr_snf), ]
met_snf <- t(scale(t(met_top))); met_snf <- met_snf[complete.cases(met_snf), ]

# Transpose: samples in rows
Data1 <- t(tc_snf)
Data2 <- t(pr_snf)
Data3 <- t(met_snf)

# SNF parameters
K_nn <- 10    # number of nearest neighbors
alpha_snf <- 0.5  # hyperparameter for local model
T_iter <- 20  # number of iterations

# Construct distance matrices
Dist1 <- (dist2(Data1, Data1))^(1/2)
Dist2 <- (dist2(Data2, Data2))^(1/2)
Dist3 <- (dist2(Data3, Data3))^(1/2)

# Construct affinity matrices
W1 <- affinityMatrix(Dist1, K = K_nn, sigma = alpha_snf)
W2 <- affinityMatrix(Dist2, K = K_nn, sigma = alpha_snf)
W3 <- affinityMatrix(Dist3, K = K_nn, sigma = alpha_snf)

# Fuse networks
W_fused <- SNF(list(W1, W2, W3), K = K_nn, t = T_iter)
rownames(W_fused) <- colnames(W_fused) <- shared

cat("  SNF fusion complete.\n")

# Spectral clustering on fused network
snf_labels_k2 <- spectralClustering(W_fused, K = 2)
names(snf_labels_k2) <- shared

# Compare with true group labels
group_numeric <- ifelse(group_vec == "Normal", 1, 2)
ari_snf_group <- adjustedRandIndex(snf_labels_k2, group_numeric)
cat(sprintf("  SNF K=2 vs True Group (Normal/Adjacent): ARI = %.3f\n", ari_snf_group))

# Classification accuracy
tab <- table(SNF = snf_labels_k2, Group = group_vec)
cat("  SNF cluster vs Group:\n")
print(tab)

# Determine which SNF cluster corresponds to which group
acc1 <- sum(diag(tab)) / sum(tab)
acc2 <- sum(diag(tab[2:1, ])) / sum(tab)
snf_accuracy <- max(acc1, acc2)
cat(sprintf("  SNF classification accuracy: %.1f%%\n", snf_accuracy * 100))


# ============================================================================
# 2. SNF: ADJACENT SAMPLES ONLY (Subtype discovery)
# ============================================================================
cat("\n>>> 2. SNF subtyping (Adjacent samples only)\n")

adj_samples <- shared[grepl("^Adjacent", shared)]

# Subset affinity matrices for Adjacent samples
adj_idx <- which(shared %in% adj_samples)

Data1_adj <- Data1[adj_idx, ]
Data2_adj <- Data2[adj_idx, ]
Data3_adj <- Data3[adj_idx, ]

Dist1_adj <- (dist2(Data1_adj, Data1_adj))^(1/2)
Dist2_adj <- (dist2(Data2_adj, Data2_adj))^(1/2)
Dist3_adj <- (dist2(Data3_adj, Data3_adj))^(1/2)

K_adj <- min(K_nn, length(adj_samples) - 1)
W1_adj <- affinityMatrix(Dist1_adj, K = K_adj, sigma = alpha_snf)
W2_adj <- affinityMatrix(Dist2_adj, K = K_adj, sigma = alpha_snf)
W3_adj <- affinityMatrix(Dist3_adj, K = K_adj, sigma = alpha_snf)

W_fused_adj <- SNF(list(W1_adj, W2_adj, W3_adj), K = K_adj, t = T_iter)
rownames(W_fused_adj) <- colnames(W_fused_adj) <- adj_samples

# Cluster for K=2,3
for (k in 2:3) {
  snf_sub <- spectralClustering(W_fused_adj, K = k)
  names(snf_sub) <- adj_samples
  cat(sprintf("  SNF K=%d: %s\n", k,
              paste(sprintf("C%d=%d", 1:k, table(factor(snf_sub, 1:k))), collapse = ", ")))
  
  # Save assignments
  snf_df <- data.frame(
    sample = adj_samples,
    patient_id = gsub("^Adjacent", "", adj_samples),
    snf_cluster = paste0("SNF_C", snf_sub)
  )
  write.csv(snf_df, file.path(OUT_DIR, sprintf("snf_subtype_K%d.csv", k)), row.names = FALSE)
}

snf_k2_adj <- spectralClustering(W_fused_adj, K = 2)
names(snf_k2_adj) <- adj_samples


# ============================================================================
# 3. CROSS-METHOD CONCORDANCE
# ============================================================================
cat("\n>>> 3. Cross-method concordance analysis\n")

# Load ConsensusClusterPlus results (Phase 6)
cc_k2_file <- file.path(P6_DIR, "subtype_K2.csv")
if (file.exists(cc_k2_file)) {
  cc_k2 <- read.csv(cc_k2_file, stringsAsFactors = FALSE)
  cc_labels <- setNames(cc_k2$subtype, cc_k2$sample)
  
  # Match samples
  matched_samps <- intersect(names(snf_k2_adj), names(cc_labels))
  if (length(matched_samps) >= 4) {
    ari_snf_cc <- adjustedRandIndex(snf_k2_adj[matched_samps],
                                     as.integer(factor(cc_labels[matched_samps])))
    cat(sprintf("  SNF vs ConsensusClusterPlus (K=2): ARI = %.3f\n", ari_snf_cc))
  }
}

# Load MOFA2 clusters (Phase 5)
mofa_cl_file <- file.path(P5_DIR, "MOFA2_sample_clusters.csv")
if (file.exists(mofa_cl_file)) {
  mofa_cl <- read.csv(mofa_cl_file, stringsAsFactors = FALSE)
  mofa_adj <- mofa_cl[mofa_cl$sample %in% adj_samples, ]
  if (nrow(mofa_adj) >= 4) {
    mofa_labels <- setNames(mofa_adj$cluster, mofa_adj$sample)
    matched2 <- intersect(names(snf_k2_adj), names(mofa_labels))
    if (length(matched2) >= 4) {
      ari_snf_mofa <- adjustedRandIndex(snf_k2_adj[matched2],
                                         as.integer(factor(mofa_labels[matched2])))
      cat(sprintf("  SNF vs MOFA2 clusters: ARI = %.3f\n", ari_snf_mofa))
    }
  }
}

# Concordance summary
concordance_df <- data.frame(
  comparison = character(),
  ARI = numeric(),
  stringsAsFactors = FALSE
)

concordance_df <- rbind(concordance_df, data.frame(
  comparison = "SNF_vs_TrueGroup (all samples)", ARI = ari_snf_group))

if (exists("ari_snf_cc")) {
  concordance_df <- rbind(concordance_df, data.frame(
    comparison = "SNF_vs_ConsensusClusterPlus (adj K=2)", ARI = ari_snf_cc))
}
if (exists("ari_snf_mofa")) {
  concordance_df <- rbind(concordance_df, data.frame(
    comparison = "SNF_vs_MOFA2 (adj)", ARI = ari_snf_mofa))
}

write.csv(concordance_df, file.path(OUT_DIR, "cross_method_concordance.csv"), row.names = FALSE)
cat("\n  Concordance summary:\n")
print(concordance_df)


# ============================================================================
# 4. VISUALIZATION
# ============================================================================
cat("\n>>> 4. Generating visualizations\n")

# ---- Fused network heatmap (all samples) ----
display_names <- gsub("Normal", "N", gsub("Adjacent", "A", shared))
anno_df <- data.frame(
  Group = group_vec,
  SNF_Cluster = paste0("C", snf_labels_k2),
  row.names = display_names
)
W_display <- W_fused
rownames(W_display) <- colnames(W_display) <- display_names

pdf(file.path(FIG_DIR, "SNF_fused_network_all.pdf"), width = 10, height = 8)
pheatmap(W_display,
         color = colorRampPalette(c("white", "#3C5488", "#00008B"))(100),
         annotation_col = anno_df, annotation_row = anno_df,
         annotation_colors = list(
           Group = c(Normal = "#4DBBD5", Adjacent = "#E64B35"),
           SNF_Cluster = c(C1 = "#3C5488", C2 = "#E64B35")
         ),
         main = "SNF Fused Network (All Samples)",
         fontsize = 8)
dev.off()
cat("  Fused network heatmap (all) saved.\n")

# ---- Fused network heatmap (adjacent only) ----
adj_display <- gsub("Adjacent", "P", adj_samples)
W_adj_display <- W_fused_adj
rownames(W_adj_display) <- colnames(W_adj_display) <- adj_display

anno_adj <- data.frame(
  SNF_Subtype = paste0("C", snf_k2_adj),
  row.names = adj_display
)
if (exists("cc_labels")) {
  anno_adj$ConsensusCluster <- cc_labels[adj_samples]
}

pdf(file.path(FIG_DIR, "SNF_fused_network_adjacent.pdf"), width = 8, height = 7)
pheatmap(W_adj_display,
         color = colorRampPalette(c("white", "#3C5488", "#00008B"))(100),
         annotation_col = anno_adj, annotation_row = anno_adj,
         annotation_colors = list(
           SNF_Subtype = c(C1 = "#3C5488", C2 = "#E64B35"),
           ConsensusCluster = c(CS1 = "#3C5488", CS2 = "#E64B35")
         ),
         main = "SNF Fused Network (Adjacent Samples, K=2)",
         fontsize = 10)
dev.off()
cat("  Fused network heatmap (adjacent) saved.\n")

# ---- Per-omics similarity comparison ----
pdf(file.path(FIG_DIR, "SNF_per_omics_vs_fused.pdf"), width = 16, height = 4)
par(mfrow = c(1, 4), mar = c(2, 2, 3, 1))

for (nm in list(list("Transcriptomics", W1), list("Proteomics", W2),
                 list("Metabolomics", W3), list("Fused", W_fused))) {
  m <- nm[[2]]; rownames(m) <- colnames(m) <- display_names
  image(m[nrow(m):1, ], col = colorRampPalette(c("white", "#3C5488"))(50),
        axes = FALSE, main = nm[[1]])
}
dev.off()
cat("  Per-omics comparison saved.\n")


# ============================================================================
# 5. SNF EIGEN-GAP FOR OPTIMAL K
# ============================================================================
cat("\n>>> 5. Eigen-gap analysis for optimal K\n")

estimate_k <- estimateNumberOfClustersGivenGraph(W_fused_adj, NUMC = 2:5)
cat("  SNF eigen-gap estimates:\n")
for (nm in names(estimate_k)) {
  cat(sprintf("    %s: K = %s\n", nm, paste(estimate_k[[nm]], collapse = ", ")))
}

# Save eigen-gap
eigen_df <- data.frame(
  method = names(estimate_k),
  optimal_K = sapply(estimate_k, function(x) x[1])
)
write.csv(eigen_df, file.path(OUT_DIR, "snf_optimal_K_estimates.csv"), row.names = FALSE)


# ============================================================================
# 6. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("SNF Integration SUMMARY\n")
cat("========================================\n")
cat(sprintf("All samples: SNF K=2 vs True Group ARI=%.3f, Accuracy=%.1f%%\n",
            ari_snf_group, snf_accuracy * 100))
cat(sprintf("Adjacent subtyping: K=2 → %s\n",
            paste(sprintf("C%d=%d", 1:2, table(factor(snf_k2_adj, 1:2))), collapse = ", ")))
if (exists("ari_snf_cc")) {
  cat(sprintf("Cross-method: SNF vs ConsensusClusterPlus ARI=%.3f\n", ari_snf_cc))
}
cat("\nInterpretation:\n")
cat("  ARI > 0.6 : Strong concordance (method-independent conclusion)\n")
cat("  ARI 0.3-0.6 : Moderate concordance\n")
cat("  ARI < 0.3 : Weak concordance (results may be method-dependent)\n")
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("SNF Integration COMPLETE.\n")
