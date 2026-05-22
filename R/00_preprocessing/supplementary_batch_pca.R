#!/usr/bin/env Rscript
# ============================================================================
# Supplementary: PCA Batch Effect Assessment
# HAE Multi-omics Integration Study
# ============================================================================
# Purpose: Generate PCA plots for each omics layer and the combined matrix
#   to evaluate potential batch effects and ensure that biological group
#   (Normal vs Adjacent) is the dominant source of variance, not technical
#   confounders (platform, run order, patient identity).
#
# Output: supplementary PCA figures for batch effect transparency
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(RColorBrewer)
  library(dplyr)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
OUT_DIR  <- file.path(PROJECT, "analysis/results/supplementary_batch_pca")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Supplementary: PCA Batch Effect Assessment\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA
# ============================================================================
cat(">>> 0. Loading multi-omics data\n")

# Load differential results for ID -> symbol mapping
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

cat(sprintf("  TC: %d genes x %d samples\n", nrow(tc_mat), ncol(tc_mat)))
cat(sprintf("  PR: %d proteins x %d samples\n", nrow(pr_mat), ncol(pr_mat)))
cat(sprintf("  MET: %d metabolites x %d samples\n", nrow(met_mat), ncol(met_mat)))


# ============================================================================
# 1. BUILD SAMPLE METADATA
# ============================================================================
cat("\n>>> 1. Building sample metadata\n")

build_meta <- function(sample_names, omics_label) {
  group <- ifelse(grepl("^Normal", sample_names), "Normal", "Adjacent")
  patient_id <- gsub("Normal|Adjacent", "", sample_names)
  data.frame(
    sample = sample_names,
    group = group,
    patient_id = patient_id,
    omics = omics_label,
    stringsAsFactors = FALSE
  )
}

meta_tc  <- build_meta(colnames(tc_mat), "Transcriptomics")
meta_pr  <- build_meta(colnames(pr_mat), "Proteomics")
meta_met <- build_meta(colnames(met_mat), "Metabolomics")


# ============================================================================
# 2. PER-OMICS PCA (colored by Group, shaped by Patient)
# ============================================================================
cat("\n>>> 2. Per-omics PCA\n")

run_pca_plot <- function(mat, meta, title_prefix) {
  # Remove zero-variance features
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat_filt <- mat[vars > 0 & !is.na(vars), ]
  
  # PCA
  pca <- prcomp(t(mat_filt), scale. = TRUE, center = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    sample = rownames(pca$x)
  )
  pca_df <- merge(pca_df, meta, by = "sample")
  
  # --- Plot A: Colored by Group ---
  p_group <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_text_repel(aes(label = sample), size = 2.2, max.overlaps = 20,
                    show.legend = FALSE) +
    scale_color_manual(values = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")) +
    labs(
      title = sprintf("%s - Colored by Group", title_prefix),
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
  
  # --- Plot B: Colored by Patient ID (to detect patient-driven batch effects) ---
  n_patients <- length(unique(pca_df$patient_id))
  if (n_patients <= 14) {
    pal <- colorRampPalette(brewer.pal(min(n_patients, 12), "Set3"))(n_patients)
  } else {
    pal <- rainbow(n_patients)
  }
  
  p_patient <- ggplot(pca_df, aes(PC1, PC2, color = patient_id)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_line(aes(group = patient_id), linetype = "dashed", alpha = 0.4, linewidth = 0.4) +
    geom_text_repel(aes(label = sample), size = 2.2, max.overlaps = 20,
                    show.legend = FALSE) +
    scale_color_manual(values = setNames(pal, sort(unique(pca_df$patient_id)))) +
    labs(
      title = sprintf("%s - Paired by Patient", title_prefix),
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Patient"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
  
  # --- Compute PERMANOVA-like test (adonis2 or manual) ---
  # Simple metric: proportion of variance on PC1-2 separating groups
  group_r2 <- tryCatch({
    fit <- aov(PC1 ~ group, data = pca_df)
    r2_pc1 <- summary(fit)[[1]]["group", "Sum Sq"] /
              sum(summary(fit)[[1]][, "Sum Sq"])
    fit2 <- aov(PC2 ~ group, data = pca_df)
    r2_pc2 <- summary(fit2)[[1]]["group", "Sum Sq"] /
              sum(summary(fit2)[[1]][, "Sum Sq"])
    c(PC1_R2 = round(r2_pc1, 3), PC2_R2 = round(r2_pc2, 3))
  }, error = function(e) c(PC1_R2 = NA, PC2_R2 = NA))
  
  patient_r2 <- tryCatch({
    fit <- aov(PC1 ~ patient_id, data = pca_df)
    r2_pc1 <- summary(fit)[[1]]["patient_id", "Sum Sq"] /
              sum(summary(fit)[[1]][, "Sum Sq"])
    fit2 <- aov(PC2 ~ patient_id, data = pca_df)
    r2_pc2 <- summary(fit2)[[1]]["patient_id", "Sum Sq"] /
              sum(summary(fit2)[[1]][, "Sum Sq"])
    c(PC1_R2 = round(r2_pc1, 3), PC2_R2 = round(r2_pc2, 3))
  }, error = function(e) c(PC1_R2 = NA, PC2_R2 = NA))
  
  cat(sprintf("  %s: Group R2 on PC1=%.3f, PC2=%.3f | Patient R2 on PC1=%.3f, PC2=%.3f\n",
              title_prefix, group_r2[1], group_r2[2], patient_r2[1], patient_r2[2]))
  
  return(list(p_group = p_group, p_patient = p_patient,
              var_exp = var_exp, group_r2 = group_r2, patient_r2 = patient_r2))
}

pca_tc  <- run_pca_plot(tc_mat, meta_tc, "Transcriptomics")
pca_pr  <- run_pca_plot(pr_mat, meta_pr, "Proteomics")
pca_met <- run_pca_plot(met_mat, meta_met, "Metabolomics")


# ============================================================================
# 3. COMBINED MULTI-OMICS PCA (shared samples only)
# ============================================================================
cat("\n>>> 3. Combined multi-omics PCA\n")

shared_samples <- Reduce(intersect, list(colnames(tc_mat), colnames(pr_mat), colnames(met_mat)))
cat(sprintf("  Shared samples: %d\n", length(shared_samples)))

# Select top variable features and concatenate
select_top <- function(mat, n) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  mat[head(order(v, decreasing = TRUE), min(n, nrow(mat))), ]
}

tc_top  <- select_top(tc_mat[, shared_samples], 1000)
pr_top  <- select_top(pr_mat[, shared_samples], 1000)
met_top <- select_top(met_mat[, shared_samples], 300)

# Z-score scale each layer independently before concatenation
tc_z  <- t(scale(t(tc_top)));  tc_z <- tc_z[complete.cases(tc_z), ]
pr_z  <- t(scale(t(pr_top)));  pr_z <- pr_z[complete.cases(pr_z), ]
met_z <- t(scale(t(met_top))); met_z <- met_z[complete.cases(met_z), ]

rownames(tc_z)  <- paste0("TC_", rownames(tc_z))
rownames(pr_z)  <- paste0("PR_", rownames(pr_z))
rownames(met_z) <- paste0("MET_", rownames(met_z))

combined <- rbind(tc_z, pr_z, met_z)
cat(sprintf("  Combined matrix: %d features x %d samples\n", nrow(combined), ncol(combined)))

meta_combined <- build_meta(shared_samples, "Combined")
pca_combined <- run_pca_plot(combined, meta_combined, "Combined Multi-omics")


# ============================================================================
# 4. FEATURE-ORIGIN PCA (color by omics layer)
# ============================================================================
cat("\n>>> 4. Feature-origin PCA (loading space)\n")

# PCA on combined matrix - show loadings colored by omics origin
pca_load <- prcomp(t(combined), scale. = FALSE, center = TRUE)
load_df <- data.frame(
  PC1 = pca_load$rotation[, 1],
  PC2 = pca_load$rotation[, 2],
  feature = rownames(pca_load$rotation)
)
load_df$omics <- ifelse(grepl("^TC_", load_df$feature), "Transcriptomics",
                  ifelse(grepl("^PR_", load_df$feature), "Proteomics", "Metabolomics"))

p_load <- ggplot(load_df, aes(PC1, PC2, color = omics)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("Transcriptomics" = "#4DBBD5",
                                "Proteomics" = "#E64B35",
                                "Metabolomics" = "#00A087")) +
  labs(title = "PCA Loading Space - Feature Origin",
       x = "PC1 loading", y = "PC2 loading", color = "Omics Layer") +
  theme_bw(base_size = 11) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  theme(legend.position = "right")


# ============================================================================
# 5. SAVE FIGURES
# ============================================================================
cat("\n>>> 5. Saving figures\n")

# --- Figure S_batch_A: Per-omics PCA colored by group (3 panels) ---
pdf(file.path(FIG_DIR, "FigS_batch_PCA_by_group.pdf"), width = 18, height = 5.5)
print(ggarrange(pca_tc$p_group, pca_pr$p_group, pca_met$p_group,
                ncol = 3, labels = c("A", "B", "C")))
dev.off()
cat("  FigS_batch_PCA_by_group.pdf saved.\n")

# --- Figure S_batch_B: Per-omics PCA colored by patient (3 panels) ---
pdf(file.path(FIG_DIR, "FigS_batch_PCA_by_patient.pdf"), width = 18, height = 5.5)
print(ggarrange(pca_tc$p_patient, pca_pr$p_patient, pca_met$p_patient,
                ncol = 3, labels = c("A", "B", "C")))
dev.off()
cat("  FigS_batch_PCA_by_patient.pdf saved.\n")

# --- Figure S_batch_C: Combined PCA + Loading space (2 panels) ---
pdf(file.path(FIG_DIR, "FigS_batch_combined_PCA.pdf"), width = 14, height = 5.5)
print(ggarrange(pca_combined$p_group, p_load, ncol = 2, labels = c("A", "B")))
dev.off()
cat("  FigS_batch_combined_PCA.pdf saved.\n")

# --- Figure S_batch_D: All-in-one 8-panel overview ---
pdf(file.path(FIG_DIR, "FigS_batch_overview.pdf"), width = 16, height = 14)
print(ggarrange(
  pca_tc$p_group, pca_tc$p_patient,
  pca_pr$p_group, pca_pr$p_patient,
  pca_met$p_group, pca_met$p_patient,
  pca_combined$p_group, p_load,
  ncol = 2, nrow = 4,
  labels = LETTERS[1:8]
))
dev.off()
cat("  FigS_batch_overview.pdf saved.\n")


# ============================================================================
# 6. VARIANCE DECOMPOSITION SUMMARY TABLE
# ============================================================================
cat("\n>>> 6. Variance decomposition summary\n")

var_summary <- data.frame(
  Omics = c("Transcriptomics", "Proteomics", "Metabolomics", "Combined"),
  PC1_var_pct = c(pca_tc$var_exp[1], pca_pr$var_exp[1], pca_met$var_exp[1], pca_combined$var_exp[1]),
  PC2_var_pct = c(pca_tc$var_exp[2], pca_pr$var_exp[2], pca_met$var_exp[2], pca_combined$var_exp[2]),
  Group_R2_PC1 = c(pca_tc$group_r2[1], pca_pr$group_r2[1], pca_met$group_r2[1], pca_combined$group_r2[1]),
  Group_R2_PC2 = c(pca_tc$group_r2[2], pca_pr$group_r2[2], pca_met$group_r2[2], pca_combined$group_r2[2]),
  Patient_R2_PC1 = c(pca_tc$patient_r2[1], pca_pr$patient_r2[1], pca_met$patient_r2[1], pca_combined$patient_r2[1]),
  Patient_R2_PC2 = c(pca_tc$patient_r2[2], pca_pr$patient_r2[2], pca_met$patient_r2[2], pca_combined$patient_r2[2])
)

write.csv(var_summary, file.path(OUT_DIR, "variance_decomposition_summary.csv"), row.names = FALSE)

cat("\n  Variance decomposition:\n")
for (i in 1:nrow(var_summary)) {
  cat(sprintf("    %s: PC1=%.1f%%, PC2=%.1f%% | Group R2: PC1=%.3f, PC2=%.3f | Patient R2: PC1=%.3f, PC2=%.3f\n",
              var_summary$Omics[i],
              var_summary$PC1_var_pct[i], var_summary$PC2_var_pct[i],
              var_summary$Group_R2_PC1[i], var_summary$Group_R2_PC2[i],
              var_summary$Patient_R2_PC1[i], var_summary$Patient_R2_PC2[i]))
}


# ============================================================================
# 7. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Supplementary Batch Effect PCA SUMMARY\n")
cat("========================================\n")
cat("Interpretation guide:\n")
cat("  - If Group R2 >> Patient R2 on PC1: biological signal dominates\n")
cat("  - If Patient R2 >> Group R2 on PC1: patient-specific batch effects present\n")
cat("  - Loading-space PCA: if omics layers form separate clusters,\n")
cat("    platform effects may dominate the combined analysis\n")
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Supplementary Batch PCA COMPLETE.\n")
