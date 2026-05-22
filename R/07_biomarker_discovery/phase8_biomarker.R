#!/usr/bin/env Rscript
# ============================================================================
# Phase 8-9: Biomarker Discovery & Figure Preparation
# HAE Multi-omics Integration Study
# ============================================================================
# Strategy:
#   1. LASSO/Elastic-net feature selection for subtype classification
#   2. Multi-omics biomarker panel (protein + metabolite)
#   3. ROC curve evaluation
#   4. Cross-omics concordance analysis (mRNA-protein correlation)
#   5. Publication-quality summary figures
# ============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
P6_DIR   <- file.path(PROJECT, "analysis/results/phase6_subtyping")
P7_DIR   <- file.path(PROJECT, "analysis/results/phase7_characterization")
OUT_DIR  <- file.path(PROJECT, "analysis/results/phase8_biomarker")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 8-9: Biomarker Discovery\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA
# ============================================================================
cat(">>> 0. Loading data\n")

# Subtype assignments
subtypes <- read.csv(file.path(P6_DIR, "subtype_K2.csv"), stringsAsFactors = FALSE)
subtype_vec <- setNames(subtypes$subtype, subtypes$sample)

# Differential results for ID mapping
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
met_mat <- as.matrix(read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"),
                               check.names = FALSE, row.names = 1))

# Match samples
adj_samples <- names(subtype_vec)
nor_samples <- gsub("^Adjacent", "Normal", adj_samples)

# All samples
all_samples <- intersect(c(adj_samples, nor_samples), colnames(tc_mat))
all_samples <- intersect(all_samples, colnames(pr_mat))
all_samples <- intersect(all_samples, colnames(met_mat))

# Create group vector (Adjacent vs Normal for disease biomarkers)
group_all <- ifelse(grepl("^Normal", all_samples), "Normal", "Adjacent")
names(group_all) <- all_samples

cat(sprintf("  Total samples: %d (Normal=%d, Adjacent=%d)\n",
            length(all_samples), sum(group_all == "Normal"), sum(group_all == "Adjacent")))
cat(sprintf("  Subtypes: CS1=%d, CS2=%d\n",
            sum(subtype_vec == "CS1"), sum(subtype_vec == "CS2")))


# ============================================================================
# 1. DISEASE BIOMARKERS (Adjacent vs Normal)
# ============================================================================
cat("\n>>> 1. Disease biomarker discovery (Adjacent vs Normal)\n")

# ---- LASSO on each omics ----
run_lasso_biomarker <- function(mat, samples, labels, omics_name, alpha = 0.5) {
  cat(sprintf("  Running elastic-net on %s...\n", omics_name))
  
  # Pre-filter: select top variable features
  vars <- apply(mat[, samples], 1, var, na.rm = TRUE)
  top_n <- min(500, sum(vars > 0))
  top_feats <- names(head(sort(vars, decreasing = TRUE), top_n))
  
  X <- t(mat[top_feats, samples])
  Y <- factor(labels[samples], levels = c("Normal", "Adjacent"))
  
  # Remove NA columns
  na_cols <- colSums(is.na(X)) > 0
  X <- X[, !na_cols]
  
  # LOOCV for lambda selection
  set.seed(42)
  cv_fit <- tryCatch({
    cv.glmnet(X, Y, family = "binomial", alpha = alpha, nfolds = min(10, nrow(X)),
              type.measure = "auc")
  }, error = function(e) {
    cat(sprintf("    CV failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(cv_fit)) return(list(features = character(0), auc = NA, auc_resub = NA,
                                    model = NULL, roc = NULL, X = NULL, Y = NULL, cv_aucs = NULL))
  
  # Extract selected features at lambda.min
  coefs <- coef(cv_fit, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  selected <- setdiff(selected, "(Intercept)")
  
  # Resubstitution prediction (optimistic - predicting on training data)
  pred_resub <- predict(cv_fit, newx = X, s = "lambda.min", type = "response")
  roc_resub <- roc(Y, pred_resub[, 1], quiet = TRUE)
  
  # Also report resubstitution AUC clearly labeled
  cat(sprintf("    Selected features: %d, Resubstitution AUC: %.3f (optimistic estimate)\n",
              length(selected), auc(roc_resub)))
  
  # Use repeated 5-fold CV for honest AUC estimate
  set.seed(42)
  n_repeats <- 20
  cv_aucs <- rep(NA_real_, n_repeats)
  for (rep_i in 1:n_repeats) {
    folds <- sample(rep(1:5, length.out = nrow(X)))
    pred_cv <- numeric(nrow(X))
    for (fi in 1:5) {
      train_idx <- which(folds != fi)
      test_idx <- which(folds == fi)
      if (length(unique(Y[train_idx])) < 2) next
      fit_cv <- tryCatch(
        glmnet(X[train_idx, , drop = FALSE], Y[train_idx], family = "binomial",
               alpha = alpha, lambda = cv_fit$lambda.min),
        error = function(e) NULL)
      if (!is.null(fit_cv)) {
        pred_cv[test_idx] <- predict(fit_cv, X[test_idx, , drop = FALSE], type = "response")[, 1]
      }
    }
    if (all(pred_cv != 0)) {
      cv_aucs[rep_i] <- tryCatch(as.numeric(auc(roc(Y, pred_cv, quiet = TRUE))), error = function(e) NA)
    }
  }
  honest_auc <- mean(cv_aucs, na.rm = TRUE)
  cat(sprintf("    Repeated 5-fold CV AUC (x%d): %.3f (SD=%.3f)\n",
              n_repeats, honest_auc, sd(cv_aucs, na.rm = TRUE)))
  
  return(list(features = selected, auc = honest_auc, auc_resub = auc(roc_resub),
              model = cv_fit, roc = roc_resub, X = X, Y = Y,
              cv_aucs = cv_aucs))
}

lasso_tc <- run_lasso_biomarker(tc_mat, all_samples, group_all, "Transcriptomics")
lasso_pr <- run_lasso_biomarker(pr_mat, all_samples, group_all, "Proteomics")
lasso_met <- run_lasso_biomarker(met_mat, all_samples, group_all, "Metabolomics")

# Combine biomarker features
all_biomarkers <- data.frame()
for (info in list(list(lasso_tc, "TC"), list(lasso_pr, "PR"), list(lasso_met, "MET"))) {
  result <- info[[1]]
  omics <- info[[2]]
  if (length(result$features) > 0) {
    all_biomarkers <- rbind(all_biomarkers, data.frame(
      feature = result$features,
      omics = omics,
      AUC = as.numeric(result$auc),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(all_biomarkers, file.path(OUT_DIR, "disease_biomarkers.csv"), row.names = FALSE)
cat(sprintf("\n  Total disease biomarkers: %d (TC=%d, PR=%d, MET=%d)\n",
            nrow(all_biomarkers),
            sum(all_biomarkers$omics == "TC"),
            sum(all_biomarkers$omics == "PR"),
            sum(all_biomarkers$omics == "MET")))

# ---- ROC curves ----
pdf(file.path(FIG_DIR, "ROC_disease_biomarkers.pdf"), width = 8, height = 7)
plot(1, type = "n", xlim = c(1, 0), ylim = c(0, 1),
     xlab = "Specificity", ylab = "Sensitivity",
     main = "Disease Biomarker ROC Curves (Adjacent vs Normal, Resubstitution)")
abline(a = 0, b = 1, lty = 2, col = "grey50")

colors <- c("#3C5488", "#E64B35", "#00A087")
labels_roc <- c()
idx <- 1

for (info in list(list(lasso_tc, "Transcriptomics"), list(lasso_pr, "Proteomics"), list(lasso_met, "Metabolomics"))) {
  result <- info[[1]]
  name <- info[[2]]
  if (!is.null(result$roc)) {
    lines(result$roc, col = colors[idx], lwd = 2)
    labels_roc <- c(labels_roc, sprintf("%s (AUC=%.3f, n=%d)",
                                         name, auc(result$roc), length(result$features)))
  }
  idx <- idx + 1
}

legend("bottomright", legend = labels_roc, col = colors[1:length(labels_roc)],
       lwd = 2, cex = 0.9)
dev.off()
cat("  ROC curve saved.\n")


# ============================================================================
# 2. SUBTYPE BIOMARKERS (CS1 vs CS2)
# ============================================================================
cat("\n>>> 2. Subtype biomarker discovery (CS1 vs CS2)\n")

# Use only Adjacent samples
subtype_labels <- subtype_vec[adj_samples]

run_lasso_subtype <- function(mat, samples, labels, omics_name, alpha = 0.5) {
  cat(sprintf("  Running elastic-net on %s (subtypes)...\n", omics_name))
  
  samps <- intersect(samples, colnames(mat))
  vars <- apply(mat[, samps, drop = FALSE], 1, var, na.rm = TRUE)
  top_n <- min(300, sum(vars > 0))
  top_feats <- names(head(sort(vars, decreasing = TRUE), top_n))
  
  X <- t(mat[top_feats, samps])
  Y <- factor(labels[samps], levels = c("CS1", "CS2"))
  
  na_cols <- colSums(is.na(X)) > 0
  X <- X[, !na_cols]
  
  set.seed(42)
  cv_fit <- tryCatch({
    # Use repeated 5-fold CV instead of LOOCV (more stable for n=12)
    cv.glmnet(X, Y, family = "binomial", alpha = alpha, nfolds = 5,
              type.measure = "auc")
  }, error = function(e) {
    cat(sprintf("    CV failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(cv_fit)) return(list(features = character(0), auc = NA, auc_resub = NA,
                                    model = NULL, roc = NULL))
  
  coefs <- coef(cv_fit, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  selected <- setdiff(selected, "(Intercept)")
  
  # Resubstitution AUC (for reference only)
  pred_resub <- predict(cv_fit, X, s = "lambda.min", type = "response")
  roc_resub <- roc(Y, pred_resub[, 1], quiet = TRUE)
  
  # Honest AUC via repeated 5-fold CV
  n_repeats <- 20
  cv_aucs <- rep(NA_real_, n_repeats)
  for (rep_i in 1:n_repeats) {
    folds <- sample(rep(1:5, length.out = nrow(X)))
    pred_cv <- numeric(nrow(X))
    for (fi in 1:5) {
      train_idx <- which(folds != fi)
      test_idx <- which(folds == fi)
      if (length(unique(Y[train_idx])) < 2) next
      fit_cv <- tryCatch(
        glmnet(X[train_idx, , drop = FALSE], Y[train_idx], family = "binomial",
               alpha = alpha, lambda = cv_fit$lambda.min),
        error = function(e) NULL)
      if (!is.null(fit_cv)) {
        pred_cv[test_idx] <- predict(fit_cv, X[test_idx, , drop = FALSE], type = "response")[, 1]
      }
    }
    if (all(pred_cv != 0)) {
      cv_aucs[rep_i] <- tryCatch(as.numeric(auc(roc(Y, pred_cv, quiet = TRUE))), error = function(e) NA)
    }
  }
  honest_auc <- mean(cv_aucs, na.rm = TRUE)
  
  cat(sprintf("    Selected features: %d, Resub AUC: %.3f, CV AUC(x%d): %.3f (SD=%.3f)\n",
              length(selected), auc(roc_resub), n_repeats, honest_auc, sd(cv_aucs, na.rm = TRUE)))
  
  return(list(features = selected, auc = honest_auc, auc_resub = auc(roc_resub),
              model = cv_fit, roc = roc_resub))
}

sub_tc <- run_lasso_subtype(tc_mat, adj_samples, subtype_labels, "Transcriptomics")
sub_pr <- run_lasso_subtype(pr_mat, adj_samples, subtype_labels, "Proteomics")
sub_met <- run_lasso_subtype(met_mat, adj_samples, subtype_labels, "Metabolomics")

subtype_biomarkers <- data.frame()
for (info in list(list(sub_tc, "TC"), list(sub_pr, "PR"), list(sub_met, "MET"))) {
  result <- info[[1]]
  omics <- info[[2]]
  if (length(result$features) > 0) {
    subtype_biomarkers <- rbind(subtype_biomarkers, data.frame(
      feature = result$features,
      omics = omics,
      AUC = as.numeric(result$auc),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(subtype_biomarkers, file.path(OUT_DIR, "subtype_biomarkers.csv"), row.names = FALSE)
cat(sprintf("  Subtype biomarkers: %d (TC=%d, PR=%d, MET=%d)\n",
            nrow(subtype_biomarkers),
            sum(subtype_biomarkers$omics == "TC"),
            sum(subtype_biomarkers$omics == "PR"),
            sum(subtype_biomarkers$omics == "MET")))

# ---- ROC for subtype classification ----
pdf(file.path(FIG_DIR, "ROC_subtype_biomarkers.pdf"), width = 8, height = 7)
plot(1, type = "n", xlim = c(1, 0), ylim = c(0, 1),
     xlab = "Specificity", ylab = "Sensitivity",
     main = "Subtype Classification ROC (CS1 vs CS2, Resubstitution)")
abline(a = 0, b = 1, lty = 2, col = "grey50")

colors <- c("#3C5488", "#E64B35", "#00A087")
labels_sub <- c()
idx <- 1
for (info in list(list(sub_tc, "Transcriptomics"), list(sub_pr, "Proteomics"), list(sub_met, "Metabolomics"))) {
  result <- info[[1]]
  name <- info[[2]]
  if (!is.null(result$roc)) {
    lines(result$roc, col = colors[idx], lwd = 2)
    labels_sub <- c(labels_sub, sprintf("%s (AUC=%.3f, n=%d)",
                                         name, auc(result$roc), length(result$features)))
  }
  idx <- idx + 1
}
legend("bottomright", legend = labels_sub, col = colors[1:length(labels_sub)],
       lwd = 2, cex = 0.9)
dev.off()
cat("  Subtype ROC curve saved.\n")


# ============================================================================
# 3. mRNA-PROTEIN CORRELATION ANALYSIS
# ============================================================================
cat("\n>>> 3. mRNA-Protein correlation analysis\n")

# Find common genes between TC and PR
common_genes <- intersect(rownames(tc_mat), rownames(pr_mat))
cat(sprintf("  Common genes (TC & PR): %d\n", length(common_genes)))

# Compute gene-wise correlation across all shared samples
shared_for_corr <- intersect(colnames(tc_mat), colnames(pr_mat))
cat(sprintf("  Shared samples for correlation: %d\n", length(shared_for_corr)))

gene_corr <- data.frame()
for (g in common_genes) {
  tc_vals <- tc_mat[g, shared_for_corr]
  pr_vals <- pr_mat[g, shared_for_corr]
  
  if (sd(tc_vals) == 0 || sd(pr_vals) == 0) next
  
  ct <- cor.test(tc_vals, pr_vals, method = "spearman")
  gene_corr <- rbind(gene_corr, data.frame(
    gene = g,
    rho = ct$estimate,
    pvalue = ct$p.value,
    stringsAsFactors = FALSE
  ))
}

gene_corr$padj <- p.adjust(gene_corr$pvalue, method = "BH")
gene_corr <- gene_corr[order(-abs(gene_corr$rho)), ]

cat(sprintf("  Gene-wise correlations computed: %d\n", nrow(gene_corr)))
cat(sprintf("  Significant (padj<0.05): %d (%.1f%%)\n",
            sum(gene_corr$padj < 0.05),
            100 * sum(gene_corr$padj < 0.05) / nrow(gene_corr)))
cat(sprintf("  Positive: %d, Negative: %d\n",
            sum(gene_corr$rho > 0), sum(gene_corr$rho < 0)))
cat(sprintf("  Median rho: %.3f\n", median(gene_corr$rho)))

write.csv(gene_corr, file.path(OUT_DIR, "mRNA_protein_correlation.csv"), row.names = FALSE)

# ---- Correlation density plot ----
p_corr <- ggplot(gene_corr, aes(x = rho)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#3C5488", alpha = 0.7) +
  geom_density(color = "#E64B35", linewidth = 1) +
  geom_vline(xintercept = median(gene_corr$rho), linetype = "dashed", color = "black") +
  annotate("text", x = median(gene_corr$rho) + 0.05, y = Inf,
           label = sprintf("Median=%.3f", median(gene_corr$rho)),
           vjust = 2, size = 4) +
  labs(title = "mRNA-Protein Correlation Distribution",
       subtitle = sprintf("n=%d genes, %d significant (padj<0.05)",
                           nrow(gene_corr), sum(gene_corr$padj < 0.05)),
       x = "Spearman rho", y = "Density") +
  theme_bw(base_size = 12)

ggsave(file.path(FIG_DIR, "mRNA_protein_correlation_dist.pdf"), p_corr, width = 8, height = 5)
cat("  Correlation distribution plot saved.\n")

# ---- Top correlated genes scatter plots ----
top_pos <- head(gene_corr[gene_corr$rho > 0, ], 6)
top_neg <- head(gene_corr[gene_corr$rho < 0, ], 3)
top_examples <- rbind(top_pos, top_neg)

if (nrow(top_examples) >= 4) {
  plot_list <- list()
  for (i in 1:min(9, nrow(top_examples))) {
    g <- top_examples$gene[i]
    df_scatter <- data.frame(
      mRNA = tc_mat[g, shared_for_corr],
      Protein = pr_mat[g, shared_for_corr],
      Group = group_all[shared_for_corr]
    )
    
    p_s <- ggplot(df_scatter, aes(mRNA, Protein, color = Group)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.5) +
      scale_color_manual(values = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")) +
      labs(title = sprintf("%s (rho=%.2f)", g, top_examples$rho[i]),
           x = "mRNA (VST)", y = "Protein (log2)") +
      theme_bw(base_size = 9) +
      theme(legend.position = "none")
    plot_list[[i]] <- p_s
  }
  
  pdf(file.path(FIG_DIR, "mRNA_protein_scatter_top.pdf"), width = 12, height = 10)
  print(ggarrange(plotlist = plot_list, ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom"))
  dev.off()
  cat("  mRNA-protein scatter plots saved.\n")
}


# ============================================================================
# 4. INTEGRATED BIOMARKER PANEL
# ============================================================================
cat("\n>>> 4. Integrated multi-omics biomarker panel\n")

# Combine disease biomarkers from all omics
# Select the most important features from each omics
combine_features <- function(lasso_tc, lasso_pr, lasso_met, mat_tc, mat_pr, mat_met, samples, labels) {
  # Get all features
  all_feats <- c()
  feat_mats <- list()
  
  if (length(lasso_tc$features) > 0) {
    tc_feats <- intersect(lasso_tc$features, rownames(mat_tc))
    if (length(tc_feats) > 0) {
      feat_mats[["TC"]] <- t(mat_tc[tc_feats, samples, drop = FALSE])
      colnames(feat_mats[["TC"]]) <- paste0("TC_", colnames(feat_mats[["TC"]]))
    }
  }
  if (length(lasso_pr$features) > 0) {
    pr_feats <- intersect(lasso_pr$features, rownames(mat_pr))
    if (length(pr_feats) > 0) {
      feat_mats[["PR"]] <- t(mat_pr[pr_feats, samples, drop = FALSE])
      colnames(feat_mats[["PR"]]) <- paste0("PR_", colnames(feat_mats[["PR"]]))
    }
  }
  if (length(lasso_met$features) > 0) {
    met_feats <- intersect(lasso_met$features, rownames(mat_met))
    if (length(met_feats) > 0) {
      feat_mats[["MET"]] <- t(mat_met[met_feats, samples, drop = FALSE])
      colnames(feat_mats[["MET"]]) <- paste0("MET_", colnames(feat_mats[["MET"]]))
    }
  }
  
  if (length(feat_mats) == 0) return(NULL)
  
  X_combined <- do.call(cbind, feat_mats)
  Y <- factor(labels[samples])
  
  # Remove NA
  na_cols <- colSums(is.na(X_combined)) > 0
  X_combined <- X_combined[, !na_cols, drop = FALSE]
  
  if (ncol(X_combined) < 2) return(NULL)
  
  # Fit combined model
  set.seed(42)
  cv_combined <- tryCatch({
    cv.glmnet(X_combined, Y, family = "binomial", alpha = 0.5,
              nfolds = min(10, nrow(X_combined)), type.measure = "auc")
  }, error = function(e) NULL)
  
  if (is.null(cv_combined)) return(NULL)
  
  # Resubstitution AUC (for reference only - optimistic)
  pred_resub <- predict(cv_combined, X_combined, s = "lambda.min", type = "response")
  roc_resub <- roc(Y, pred_resub[, 1], quiet = TRUE)
  
  # Honest repeated 5-fold CV AUC
  n_repeats <- 20
  cv_aucs <- rep(NA_real_, n_repeats)
  for (rep_i in 1:n_repeats) {
    folds <- sample(rep(1:5, length.out = nrow(X_combined)))
    pred_cv <- numeric(nrow(X_combined))
    for (fi in 1:5) {
      train_idx <- which(folds != fi)
      test_idx <- which(folds == fi)
      if (length(unique(Y[train_idx])) < 2) next
      fit_cv <- tryCatch(
        glmnet(X_combined[train_idx, , drop = FALSE], Y[train_idx], family = "binomial",
               alpha = 0.5, lambda = cv_combined$lambda.min),
        error = function(e) NULL)
      if (!is.null(fit_cv)) {
        pred_cv[test_idx] <- predict(fit_cv, X_combined[test_idx, , drop = FALSE], type = "response")[, 1]
      }
    }
    if (all(pred_cv != 0)) {
      cv_aucs[rep_i] <- tryCatch(as.numeric(auc(roc(Y, pred_cv, quiet = TRUE))), error = function(e) NA)
    }
  }
  honest_auc <- mean(cv_aucs, na.rm = TRUE)
  
  coefs <- coef(cv_combined, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  selected <- setdiff(selected, "(Intercept)")
  
  return(list(auc = honest_auc, auc_resub = auc(roc_resub), roc = roc_resub,
              features = selected, model = cv_combined))
}

# Combined disease panel
combined_disease <- combine_features(lasso_tc, lasso_pr, lasso_met,
                                      tc_mat, pr_mat, met_mat,
                                      all_samples, group_all)

if (!is.null(combined_disease)) {
  cat(sprintf("  Combined disease panel: AUC=%.3f, %d features\n",
              combined_disease$auc, length(combined_disease$features)))
  
  # Add to ROC plot
  pdf(file.path(FIG_DIR, "ROC_combined_disease.pdf"), width = 8, height = 7)
  plot(1, type = "n", xlim = c(1, 0), ylim = c(0, 1),
       xlab = "Specificity", ylab = "Sensitivity",
       main = "Disease Biomarker Panels (Adjacent vs Normal, Resubstitution)")
  abline(a = 0, b = 1, lty = 2, col = "grey50")
  
  colors <- c("#3C5488", "#E64B35", "#00A087", "#F39B7F")
  legend_labels <- c()
  idx <- 1
  
  for (info in list(list(lasso_tc, "TC"), list(lasso_pr, "PR"), list(lasso_met, "MET"))) {
    r <- info[[1]]
    n <- info[[2]]
    if (!is.null(r$roc)) {
      lines(r$roc, col = colors[idx], lwd = 2, lty = 2)
      legend_labels <- c(legend_labels, sprintf("%s (AUC=%.3f)", n, auc(r$roc)))
    }
    idx <- idx + 1
  }
  
  lines(combined_disease$roc, col = colors[4], lwd = 3)
  legend_labels <- c(legend_labels, sprintf("Combined (AUC=%.3f)", combined_disease$auc))
  
  legend("bottomright", legend = legend_labels,
         col = colors[1:length(legend_labels)],
         lwd = c(rep(2, length(legend_labels) - 1), 3),
         lty = c(rep(2, length(legend_labels) - 1), 1),
         cex = 0.9)
  dev.off()
  cat("  Combined ROC saved.\n")
}

# Combined subtype panel
combined_subtype <- combine_features(sub_tc, sub_pr, sub_met,
                                      tc_mat, pr_mat, met_mat,
                                      adj_samples, subtype_labels)

if (!is.null(combined_subtype)) {
  cat(sprintf("  Combined subtype panel: AUC=%.3f, %d features\n",
              combined_subtype$auc, length(combined_subtype$features)))
}


# ============================================================================
# 5. BIOMARKER PANEL HEATMAP
# ============================================================================
cat("\n>>> 5. Biomarker panel visualization\n")

# Disease biomarker panel heatmap
if (nrow(all_biomarkers) > 0) {
  biomarker_genes_tc <- all_biomarkers$feature[all_biomarkers$omics == "TC"]
  biomarker_genes_pr <- all_biomarkers$feature[all_biomarkers$omics == "PR"]
  biomarker_mets <- all_biomarkers$feature[all_biomarkers$omics == "MET"]
  
  # Build heatmap matrix
  samp_order <- c(
    sort(intersect(nor_samples, all_samples)),
    sort(intersect(adj_samples, all_samples))
  )
  
  hm_parts <- list()
  
  # TC layer
  tc_bio <- intersect(biomarker_genes_tc, rownames(tc_mat))
  if (length(tc_bio) >= 2) {
    m <- t(scale(t(tc_mat[tc_bio, samp_order])))
    m[m > 3] <- 3; m[m < -3] <- -3
    rownames(m) <- paste0("[TC] ", rownames(m))
    hm_parts[["TC"]] <- m
  }
  
  # PR layer
  pr_bio <- intersect(biomarker_genes_pr, rownames(pr_mat))
  if (length(pr_bio) >= 2) {
    m <- t(scale(t(pr_mat[pr_bio, samp_order])))
    m[m > 3] <- 3; m[m < -3] <- -3
    rownames(m) <- paste0("[PR] ", rownames(m))
    hm_parts[["PR"]] <- m
  }
  
  # MET layer
  met_bio <- intersect(biomarker_mets, rownames(met_mat))
  if (length(met_bio) >= 2) {
    m <- t(scale(t(met_mat[met_bio, samp_order])))
    m[m > 3] <- 3; m[m < -3] <- -3
    rownames(m) <- paste0("[MET] ", rownames(m))
    hm_parts[["MET"]] <- m
  }
  
  if (length(hm_parts) > 0) {
    bio_mat <- do.call(rbind, hm_parts)
    
    anno_col <- data.frame(
      Group = ifelse(grepl("^Normal", samp_order), "Normal", "Adjacent")
    )
    rownames(anno_col) <- samp_order
    
    # Add subtype info for Adjacent samples
    anno_col$Subtype <- NA
    adj_in <- intersect(names(subtype_vec), samp_order)
    anno_col[adj_in, "Subtype"] <- subtype_vec[adj_in]
    
    pdf(file.path(FIG_DIR, "biomarker_panel_heatmap.pdf"), width = 12, height = max(6, nrow(bio_mat) * 0.3 + 3))
    pheatmap(bio_mat,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             annotation_col = anno_col,
             annotation_colors = list(
               Group = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35"),
               Subtype = c("CS1" = "#3C5488", "CS2" = "#F39B7F")
             ),
             cluster_cols = FALSE,
             cluster_rows = TRUE,
             show_colnames = TRUE,
             fontsize_row = 7,
             fontsize_col = 7,
             main = "Multi-omics Biomarker Panel",
             breaks = seq(-3, 3, length.out = 101))
    dev.off()
    cat("  Biomarker panel heatmap saved.\n")
  }
}


# ============================================================================
# 6. CLINICAL VARIABLE PREDICTION
# ============================================================================
cat("\n>>> 6. Clinical variable prediction\n")

# Load clinical data
clinical <- read.csv(file.path(P6_DIR, "clinical_with_subtypes.csv"), stringsAsFactors = FALSE)

# Test association: biomarker expression vs key clinical variables
# Focus on invasion status (bile duct & vascular)
cat("  Testing biomarker-clinical associations...\n")

clin_assoc <- data.frame()

# Get Adjacent samples with clinical data
adj_with_clin <- intersect(adj_samples, paste0("Adjacent", clinical$patient_id))

for (feat in head(all_biomarkers$feature, 30)) {
  omics <- all_biomarkers$omics[all_biomarkers$feature == feat][1]
  
  if (omics == "TC" && feat %in% rownames(tc_mat)) {
    vals <- tc_mat[feat, adj_with_clin]
  } else if (omics == "PR" && feat %in% rownames(pr_mat)) {
    vals <- pr_mat[feat, adj_with_clin]
  } else if (omics == "MET" && feat %in% rownames(met_mat)) {
    vals <- met_mat[feat, adj_with_clin]
  } else {
    next
  }
  
  # Map to patient IDs
  pids <- gsub("^Adjacent", "", names(vals))
  
  # Test vs bile duct invasion
  if ("bile_duct_invasion" %in% colnames(clinical)) {
    bd <- clinical$bile_duct_invasion[match(pids, clinical$patient_id)]
    if (length(unique(bd[!is.na(bd)])) >= 2) {
      tryCatch({
        wt <- wilcox.test(vals ~ bd)
        clin_assoc <- rbind(clin_assoc, data.frame(
          feature = feat, omics = omics, clinical_var = "bile_duct_invasion",
          pvalue = wt$p.value, stringsAsFactors = FALSE
        ))
      }, error = function(e) NULL)
    }
  }
  
  # Test vs vascular invasion
  if ("vascular_invasion" %in% colnames(clinical)) {
    vi <- clinical$vascular_invasion[match(pids, clinical$patient_id)]
    # Note: vascular_invasion is complex (liver segment codes), skip if too many categories
  }
  
  # Test vs PNM staging
  for (pnm in c("PNM_P", "PNM_N", "PNM_M")) {
    if (pnm %in% colnames(clinical)) {
      stage <- clinical[[pnm]][match(pids, clinical$patient_id)]
      if (length(unique(stage[!is.na(stage)])) >= 2) {
        tryCatch({
          ct <- cor.test(vals, as.numeric(stage), method = "spearman")
          clin_assoc <- rbind(clin_assoc, data.frame(
            feature = feat, omics = omics, clinical_var = pnm,
            pvalue = ct$p.value, stringsAsFactors = FALSE
          ))
        }, error = function(e) NULL)
      }
    }
  }
}

if (nrow(clin_assoc) > 0) {
  clin_assoc$padj <- p.adjust(clin_assoc$pvalue, method = "BH")
  clin_assoc <- clin_assoc[order(clin_assoc$pvalue), ]
  write.csv(clin_assoc, file.path(OUT_DIR, "biomarker_clinical_association.csv"), row.names = FALSE)
  
  sig_clin <- clin_assoc[clin_assoc$pvalue < 0.05, ]
  cat(sprintf("  Biomarker-clinical associations: %d tests, %d sig (P<0.05)\n",
              nrow(clin_assoc), nrow(sig_clin)))
  
  if (nrow(sig_clin) > 0) {
    cat("  Top associations:\n")
    for (i in 1:min(10, nrow(sig_clin))) {
      cat(sprintf("    %s (%s) ~ %s: P=%.4f\n",
                  sig_clin$feature[i], sig_clin$omics[i],
                  sig_clin$clinical_var[i], sig_clin$pvalue[i]))
    }
  }
}


# ============================================================================
# 7. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 8-9 SUMMARY\n")
cat("========================================\n")

cat("1. Disease biomarkers (Adjacent vs Normal) [Repeated 5-fold CV AUC]:\n")
cat(sprintf("   TC: %d features, CV-AUC=%.3f (resub=%.3f)\n",
            length(lasso_tc$features), lasso_tc$auc, lasso_tc$auc_resub))
cat(sprintf("   PR: %d features, CV-AUC=%.3f (resub=%.3f)\n",
            length(lasso_pr$features), lasso_pr$auc, lasso_pr$auc_resub))
cat(sprintf("   MET: %d features, CV-AUC=%.3f (resub=%.3f)\n",
            length(lasso_met$features), lasso_met$auc, lasso_met$auc_resub))
if (!is.null(combined_disease)) {
  cat(sprintf("   Combined: %d features, CV-AUC=%.3f (resub=%.3f)\n",
              length(combined_disease$features), combined_disease$auc, combined_disease$auc_resub))
}

cat("\n2. Subtype biomarkers (CS1 vs CS2) [Repeated 5-fold CV AUC]:\n")
cat(sprintf("   TC: %d features, CV-AUC=%.3f (resub=%.3f)\n",
            length(sub_tc$features), sub_tc$auc, sub_tc$auc_resub))
cat(sprintf("   PR: %d features, CV-AUC=%.3f (resub=%.3f)\n",
            length(sub_pr$features), sub_pr$auc, sub_pr$auc_resub))
cat(sprintf("   MET: %d features, CV-AUC=%.3f (resub=%.3f)\n",
            length(sub_met$features), sub_met$auc, sub_met$auc_resub))
if (!is.null(combined_subtype)) {
  cat(sprintf("   Combined: %d features, CV-AUC=%.3f (resub=%.3f)\n",
              length(combined_subtype$features), combined_subtype$auc, combined_subtype$auc_resub))
}

cat(sprintf("\n3. mRNA-Protein correlation:\n"))
cat(sprintf("   %d gene-wise correlations, median rho=%.3f\n",
            nrow(gene_corr), median(gene_corr$rho)))
cat(sprintf("   %d significant (padj<0.05)\n", sum(gene_corr$padj < 0.05)))

cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 8-9 COMPLETE.\n")
