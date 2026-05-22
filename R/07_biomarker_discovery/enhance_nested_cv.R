#!/usr/bin/env Rscript
# ============================================================================
# Enhancement: Nested Cross-Validation for Biomarker Discovery
# ============================================================================
# PURPOSE: Eliminate data leakage by placing feature selection INSIDE the
#   cross-validation loop. The original phase8 selects top-variance features
#   on the FULL dataset, then performs CV -- creating information leakage.
#   This script implements proper nested CV where:
#     Outer loop: evaluates model performance (honest AUC)
#     Inner loop: selects features + tunes model (per fold)
#
# OUTPUT: Nested CV AUC vs. Original AUC comparison table
#   If gap is small → original results are valid
#   If gap is large → overfitting was present
# ============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhance_nested_cv")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Nested Cross-Validation Biomarker Analysis\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA (same as phase8)
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

# Group labels
build_group <- function(samps) {
  g <- ifelse(grepl("^Normal", samps), 0, 1)
  names(g) <- samps
  g
}

cat(sprintf("  TC: %d x %d, PR: %d x %d, MET: %d x %d\n",
            nrow(tc_mat), ncol(tc_mat), nrow(pr_mat), ncol(pr_mat),
            nrow(met_mat), ncol(met_mat)))


# ============================================================================
# 1. NESTED CV FUNCTION
# ============================================================================
cat("\n>>> 1. Defining nested CV procedure\n")

nested_cv <- function(mat, group_vec, omics_name,
                      n_outer_repeats = 50,
                      n_outer_folds = 5,
                      n_top_features = 500,
                      alpha = 0.5,
                      seed = 42) {
  
  set.seed(seed)
  n <- length(group_vec)
  Y <- factor(group_vec, levels = c(0, 1))
  
  # --- Original approach (feature selection on FULL data → CV) ---
  vars_full <- apply(mat, 1, var, na.rm = TRUE)
  top_idx_full <- head(order(vars_full, decreasing = TRUE),
                       min(n_top_features, nrow(mat)))
  X_full <- t(mat[top_idx_full, ])
  
  original_cv_aucs <- rep(NA_real_, n_outer_repeats)
  for (r in 1:n_outer_repeats) {
    folds <- sample(rep(1:n_outer_folds, length.out = n))
    pred_cv <- numeric(n)
    for (k in 1:n_outer_folds) {
      test_idx <- which(folds == k)
      train_idx <- which(folds != k)
      if (length(unique(Y[train_idx])) < 2) next
      fit <- tryCatch(
        cv.glmnet(X_full[train_idx, , drop = FALSE], Y[train_idx],
                  family = "binomial", alpha = alpha, nfolds = min(5, length(train_idx))),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        pred_cv[test_idx] <- as.numeric(predict(fit, X_full[test_idx, , drop = FALSE],
                                                 s = "lambda.min", type = "response"))
      }
    }
    if (all(pred_cv != 0) && length(unique(pred_cv)) > 1) {
      original_cv_aucs[r] <- tryCatch(
        as.numeric(auc(roc(Y, pred_cv, quiet = TRUE))),
        error = function(e) NA
      )
    }
  }
  
  # --- Nested approach (feature selection INSIDE each fold) ---
  nested_cv_aucs <- rep(NA_real_, n_outer_repeats)
  for (r in 1:n_outer_repeats) {
    folds <- sample(rep(1:n_outer_folds, length.out = n))
    pred_nested <- numeric(n)
    for (k in 1:n_outer_folds) {
      test_idx <- which(folds == k)
      train_idx <- which(folds != k)
      if (length(unique(Y[train_idx])) < 2) next
      
      # Feature selection on TRAINING data only
      train_mat <- mat[, train_idx, drop = FALSE]
      vars_train <- apply(train_mat, 1, var, na.rm = TRUE)
      top_train <- head(order(vars_train, decreasing = TRUE),
                        min(n_top_features, nrow(train_mat)))
      
      X_train <- t(mat[top_train, train_idx, drop = FALSE])
      X_test  <- t(mat[top_train, test_idx, drop = FALSE])
      
      fit <- tryCatch(
        cv.glmnet(X_train, Y[train_idx],
                  family = "binomial", alpha = alpha,
                  nfolds = min(5, length(train_idx))),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        pred_nested[test_idx] <- as.numeric(predict(fit, X_test,
                                                     s = "lambda.min", type = "response"))
      }
    }
    if (all(pred_nested != 0) && length(unique(pred_nested)) > 1) {
      nested_cv_aucs[r] <- tryCatch(
        as.numeric(auc(roc(Y, pred_nested, quiet = TRUE))),
        error = function(e) NA
      )
    }
  }
  
  original_mean <- mean(original_cv_aucs, na.rm = TRUE)
  nested_mean   <- mean(nested_cv_aucs, na.rm = TRUE)
  gap <- original_mean - nested_mean
  
  cat(sprintf("  %s: Original CV AUC=%.3f, Nested CV AUC=%.3f, Gap=%.3f\n",
              omics_name, original_mean, nested_mean, gap))
  
  return(data.frame(
    omics = omics_name,
    original_cv_auc_mean = round(original_mean, 4),
    original_cv_auc_sd   = round(sd(original_cv_aucs, na.rm = TRUE), 4),
    nested_cv_auc_mean   = round(nested_mean, 4),
    nested_cv_auc_sd     = round(sd(nested_cv_aucs, na.rm = TRUE), 4),
    auc_gap              = round(gap, 4),
    leakage_pct          = round(gap / original_mean * 100, 1),
    n_valid_original     = sum(!is.na(original_cv_aucs)),
    n_valid_nested       = sum(!is.na(nested_cv_aucs)),
    original_aucs        = I(list(original_cv_aucs)),
    nested_aucs          = I(list(nested_cv_aucs))
  ))
}


# ============================================================================
# 2. RUN NESTED CV FOR EACH OMICS
# ============================================================================
cat("\n>>> 2. Running nested CV per omics layer\n")

results <- list()

# TC (n=24: 12 Normal + 12 Adjacent)
tc_samps <- colnames(tc_mat)
tc_group <- build_group(tc_samps)
cat("  Running Transcriptomics (n=%d)...\n", length(tc_samps))
results$TC <- nested_cv(tc_mat, tc_group, "Transcriptomics",
                         n_top_features = 500, n_outer_repeats = 50)

# PR (n=28: 14 Normal + 14 Adjacent)
pr_samps <- colnames(pr_mat)
pr_group <- build_group(pr_samps)
cat("  Running Proteomics (n=%d)...\n", length(pr_samps))
results$PR <- nested_cv(pr_mat, pr_group, "Proteomics",
                         n_top_features = 500, n_outer_repeats = 50)

# MET (n=28: 14 Normal + 14 Adjacent)
met_samps <- colnames(met_mat)
met_group <- build_group(met_samps)
cat("  Running Metabolomics (n=%d)...\n", length(met_samps))
results$MET <- nested_cv(met_mat, met_group, "Metabolomics",
                          n_top_features = 300, n_outer_repeats = 50)

# Multi-omics combined (shared samples)
shared <- Reduce(intersect, list(colnames(tc_mat), colnames(pr_mat), colnames(met_mat)))
tc_top_shared <- tc_mat[head(order(apply(tc_mat[, shared], 1, var), decreasing = TRUE), 300), shared]
pr_top_shared <- pr_mat[head(order(apply(pr_mat[, shared], 1, var), decreasing = TRUE), 300), shared]
met_top_shared <- met_mat[head(order(apply(met_mat[, shared], 1, var), decreasing = TRUE), 200), shared]

# Scale and concatenate
sc <- function(m) { s <- t(scale(t(m))); s[complete.cases(s), ] }
combined <- rbind(sc(tc_top_shared), sc(pr_top_shared), sc(met_top_shared))
combined_group <- build_group(shared)
cat("  Running Combined multi-omics (n=%d)...\n", length(shared))
results$Combined <- nested_cv(combined, combined_group, "Combined",
                               n_top_features = 500, n_outer_repeats = 50)


# ============================================================================
# 3. RESULTS TABLE
# ============================================================================
cat("\n>>> 3. Compiling results\n")

summary_df <- do.call(rbind, lapply(results, function(r) {
  r[, !names(r) %in% c("original_aucs", "nested_aucs")]
}))
rownames(summary_df) <- NULL

write.csv(summary_df, file.path(OUT_DIR, "nested_cv_comparison.csv"), row.names = FALSE)
cat("\n  Summary:\n")
print(summary_df[, 1:7])

cat("\n  Interpretation:\n")
cat("    Gap < 0.05 : minimal leakage, original results reliable\n")
cat("    Gap 0.05-0.10 : moderate leakage, nested AUC is the honest estimate\n")
cat("    Gap > 0.10 : substantial leakage, must report nested AUC\n")


# ============================================================================
# 4. VISUALIZATION
# ============================================================================
cat("\n>>> 4. Generating comparison plots\n")

# Collect all AUC distributions
all_aucs <- data.frame()
for (nm in names(results)) {
  r <- results[[nm]]
  orig <- r$original_aucs[[1]]
  nest <- r$nested_aucs[[1]]
  all_aucs <- rbind(all_aucs,
    data.frame(omics = nm, method = "Original CV", auc = orig[!is.na(orig)]),
    data.frame(omics = nm, method = "Nested CV", auc = nest[!is.na(nest)])
  )
}
all_aucs$method <- factor(all_aucs$method, levels = c("Original CV", "Nested CV"))
all_aucs$omics  <- factor(all_aucs$omics, levels = c("Transcriptomics", "Proteomics",
                                                       "Metabolomics", "Combined"))

p_box <- ggplot(all_aucs, aes(x = omics, y = auc, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Original CV" = "#E64B35", "Nested CV" = "#3C5488")) +
  labs(title = "Nested CV vs Original CV: Data Leakage Assessment",
       subtitle = "Gap between methods indicates degree of information leakage",
       x = "", y = "Cross-Validated AUC", fill = "Method") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.5, y = 0.52, label = "Random", hjust = 0, size = 3, color = "grey50") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 15, hjust = 1))

# Add gap annotations
for (i in 1:nrow(summary_df)) {
  p_box <- p_box +
    annotate("text",
             x = i, y = max(all_aucs$auc, na.rm = TRUE) + 0.02,
             label = sprintf("Gap=%.3f", summary_df$auc_gap[i]),
             size = 3, fontface = "bold")
}

ggsave(file.path(FIG_DIR, "nested_cv_comparison.pdf"), p_box, width = 10, height = 6)
cat("  Comparison plot saved.\n")

# ---- Bar plot of leakage percentage ----
p_leak <- ggplot(summary_df, aes(x = omics, y = leakage_pct, fill = omics)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", leakage_pct)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Transcriptomics" = "#4DBBD5", "Proteomics" = "#E64B35",
                                "Metabolomics" = "#00A087", "Combined" = "#3C5488")) +
  labs(title = "Information Leakage Assessment",
       x = "", y = "AUC Inflation (%)", fill = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  ylim(0, max(summary_df$leakage_pct, na.rm = TRUE) * 1.3)

ggsave(file.path(FIG_DIR, "leakage_assessment.pdf"), p_leak, width = 7, height = 5)
cat("  Leakage plot saved.\n")


# ============================================================================
# 5. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Nested CV Analysis SUMMARY\n")
cat("========================================\n")
for (i in 1:nrow(summary_df)) {
  status <- ifelse(summary_df$auc_gap[i] < 0.05, "MINIMAL LEAKAGE",
             ifelse(summary_df$auc_gap[i] < 0.10, "MODERATE LEAKAGE", "SUBSTANTIAL LEAKAGE"))
  cat(sprintf("  %s: Original=%.3f, Nested=%.3f, Gap=%.3f [%s]\n",
              summary_df$omics[i], summary_df$original_cv_auc_mean[i],
              summary_df$nested_cv_auc_mean[i], summary_df$auc_gap[i], status))
}
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Nested CV Analysis COMPLETE.\n")
