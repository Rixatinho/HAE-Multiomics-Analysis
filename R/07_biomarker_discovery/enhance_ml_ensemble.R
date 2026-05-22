# ==============================================================================
# HAE Multi-omics Machine Learning Ensemble Analysis
# ==============================================================================
# Purpose: Implement 3 advanced ML methods to complement existing RF/LASSO analysis
#   1. XGBoost subtype classification (compare with RF)
#   2. SVM-RFE feature selection (compare with LASSO)
#   3. Stacking ensemble learning (RF + SVM + LASSO + meta-learner)
# 
# Author: HAE Multi-omics Project
# Date: 2026-03-23
# ==============================================================================

# Set seed for reproducibility
set.seed(42)

# ============================================================================
# 1. ENVIRONMENT SETUP
# ============================================================================
cat("=============================================================================\n")
cat(" HAE Multi-omics Machine Learning Ensemble Analysis\n")
cat("=============================================================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Load/install required packages
required_pkgs <- c("xgboost", "caret", "randomForest", "glmnet", "e1071", 
                   "pROC", "ggplot2", "dplyr", "tidyr", "VennDiagram",
                   "RColorBrewer", "gridExtra", "scales", "pheatmap")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("[OK] All packages loaded successfully.\n\n")

# ============================================================================
# 2. DATA LOADING AND PREPARATION
# ============================================================================
cat("[Step 1] Loading and preparing data...\n")

# Set paths
base_dir <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
data_dir <- file.path(base_dir, "analysis/data/processed")
results_dir <- file.path(base_dir, "analysis/results/enhancement_ml_ensemble")

# Create output directory
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
cat("  - Output directory:", results_dir, "\n")

# Load proteomics data
prot_data <- read.csv(file.path(data_dir, "proteomics_log2_norm.csv"), row.names = 1)
cat("  - Proteomics:", nrow(prot_data), "features x", ncol(prot_data), "samples\n")

# Load metabolomics data
metab_data <- read.csv(file.path(data_dir, "metabolomics_log2_merged.csv"), row.names = 1)
cat("  - Metabolomics:", nrow(metab_data), "features x", ncol(metab_data), "samples\n")

# Load subtype labels
subtype_file <- file.path(base_dir, "analysis/results/phase6_subtyping/subtype_K2.csv")
subtype_df <- read.csv(subtype_file, stringsAsFactors = FALSE)
cat("  - Subtype labels loaded:", nrow(subtype_df), "samples\n")

# Create subtype mapping
subtype_map <- setNames(subtype_df$subtype, subtype_df$sample)

# Get Adjacent samples with subtype labels
adj_samples <- grep("^Adjacent", colnames(prot_data), value = TRUE)
adj_samples <- adj_samples[adj_samples %in% names(subtype_map)]
cat("  - Adjacent samples with subtypes:", length(adj_samples), "\n")

# Extract labels for Adjacent samples
y_subtype <- factor(subtype_map[adj_samples])
cat("  - Subtype distribution:", table(y_subtype), "\n")

# ============================================================================
# 2.1 Feature Selection and Matrix Preparation
# ============================================================================
cat("\n  - Preparing feature matrix...\n")

# Variance-based feature filtering for proteomics (top 500)
prot_adj <- prot_data[, adj_samples]
prot_vars <- apply(prot_adj, 1, var, na.rm = TRUE)
prot_top500 <- names(sort(prot_vars, decreasing = TRUE))[1:min(500, length(prot_vars))]
prot_filtered <- prot_adj[prot_top500, ]

# Variance-based feature filtering for metabolomics (top 200)
metab_adj <- metab_data[, adj_samples]
metab_vars <- apply(metab_adj, 1, var, na.rm = TRUE)
metab_top200 <- names(sort(metab_vars, decreasing = TRUE))[1:min(200, length(metab_vars))]
metab_filtered <- metab_adj[metab_top200, ]

# Add suffix to distinguish omics
rownames(prot_filtered) <- paste0(rownames(prot_filtered), "_PROT")
rownames(metab_filtered) <- paste0(rownames(metab_filtered), "_METAB")

# Combine features
combined_mat <- rbind(prot_filtered, metab_filtered)
combined_mat <- combined_mat[complete.cases(combined_mat), ]

# Transpose: samples x features
X_all <- t(as.matrix(combined_mat))
X_all <- scale(X_all)  # Z-score normalize

cat("    - Combined feature matrix:", nrow(X_all), "samples x", ncol(X_all), "features\n")

# Verify sample alignment
stopifnot(all(rownames(X_all) == adj_samples))

# ============================================================================
# 3. XGBOOST SUBTYPE CLASSIFICATION WITH LOOCV
# ============================================================================
cat("\n[Step 2] XGBoost Subtype Classification (LOOCV)...\n")

n_samples <- nrow(X_all)
y_numeric <- as.numeric(y_subtype) - 1  # Convert to 0/1 for xgboost

# XGBoost parameters (conservative for small sample)
xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  max_depth = 3,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  lambda = 1,
  alpha = 0
)
nrounds <- 100

# LOOCV for XGBoost
xgb_preds <- rep(NA, n_samples)
xgb_probs <- rep(NA, n_samples)
xgb_importance_list <- list()

cat("  - Running LOOCV (n=", n_samples, ")...\n")
pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (i in 1:n_samples) {
  # Train/test split
  X_train <- X_all[-i, , drop = FALSE]
  y_train <- y_numeric[-i]
  X_test <- X_all[i, , drop = FALSE]
  
  # Create DMatrix
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgb.DMatrix(data = X_test)
  
  # Train model
  suppressWarnings({
    xgb_model <- xgb.train(
      params = xgb_params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )
  })
  
  # Predict
  prob <- predict(xgb_model, dtest)
  xgb_probs[i] <- prob
  xgb_preds[i] <- ifelse(prob > 0.5, 1, 0)
  
  # Feature importance
  imp <- xgb.importance(model = xgb_model)
  xgb_importance_list[[i]] <- imp
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Calculate metrics
xgb_accuracy <- mean(xgb_preds == y_numeric)
xgb_roc <- roc(y_numeric, xgb_probs, quiet = TRUE)
xgb_auc <- auc(xgb_roc)

# Confusion matrix
xgb_cm <- confusionMatrix(factor(xgb_preds, levels = c(0, 1)), 
                          factor(y_numeric, levels = c(0, 1)))
xgb_f1 <- xgb_cm$byClass["F1"]
if (is.na(xgb_f1)) xgb_f1 <- 0

cat("\n  - XGBoost Results:\n")
cat("    Accuracy:", sprintf("%.1f%%", xgb_accuracy * 100), "\n")
cat("    AUC:", sprintf("%.3f", xgb_auc), "\n")
cat("    F1 Score:", sprintf("%.3f", xgb_f1), "\n")

# Aggregate feature importance
xgb_imp_all <- do.call(rbind, xgb_importance_list)
xgb_imp_mean <- xgb_imp_all %>%
  group_by(Feature) %>%
  summarise(Gain = mean(Gain, na.rm = TRUE)) %>%
  arrange(desc(Gain)) %>%
  head(20)

cat("  - Top 5 XGBoost features:\n")
print(head(xgb_imp_mean, 5))

# ============================================================================
# 4. RANDOM FOREST WITH LOOCV (for comparison)
# ============================================================================
cat("\n[Step 3] Random Forest Subtype Classification (LOOCV)...\n")

rf_preds <- rep(NA, n_samples)
rf_probs <- rep(NA, n_samples)
rf_importance_list <- list()

pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (i in 1:n_samples) {
  X_train <- X_all[-i, , drop = FALSE]
  y_train <- y_subtype[-i]
  X_test <- X_all[i, , drop = FALSE]
  
  rf_model <- randomForest(X_train, y_train, ntree = 500, importance = TRUE)
  
  rf_preds[i] <- as.character(predict(rf_model, X_test))
  rf_probs[i] <- predict(rf_model, X_test, type = "prob")[, "CS2"]
  
  # Feature importance
  imp <- importance(rf_model, type = 1)
  rf_importance_list[[i]] <- data.frame(
    Feature = rownames(imp),
    MeanDecreaseAccuracy = imp[, 1]
  )
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Calculate metrics
rf_accuracy <- mean(rf_preds == as.character(y_subtype))
rf_roc <- roc(as.numeric(y_subtype) - 1, rf_probs, quiet = TRUE)
rf_auc <- auc(rf_roc)

rf_cm <- confusionMatrix(factor(rf_preds, levels = levels(y_subtype)), y_subtype)
rf_f1 <- rf_cm$byClass["F1"]
if (is.na(rf_f1)) rf_f1 <- 0

cat("\n  - Random Forest Results:\n")
cat("    Accuracy:", sprintf("%.1f%%", rf_accuracy * 100), "\n")
cat("    AUC:", sprintf("%.3f", rf_auc), "\n")
cat("    F1 Score:", sprintf("%.3f", rf_f1), "\n")

# Aggregate RF importance
rf_imp_all <- do.call(rbind, rf_importance_list)
rf_imp_mean <- rf_imp_all %>%
  group_by(Feature) %>%
  summarise(Importance = mean(MeanDecreaseAccuracy, na.rm = TRUE)) %>%
  arrange(desc(Importance)) %>%
  head(20)

# ============================================================================
# 5. SVM-RFE FEATURE SELECTION
# ============================================================================
cat("\n[Step 4] SVM-RFE Feature Selection...\n")

# Use a subset of features for RFE (computational constraint)
# Select top 100 by variance
top_var_idx <- order(apply(X_all, 2, var), decreasing = TRUE)[1:100]
X_rfe <- X_all[, top_var_idx]

# Define RFE control with LOOCV
rfe_ctrl <- rfeControl(
  functions = caretFuncs,
  method = "LOOCV",
  verbose = FALSE
)

# Define subset sizes to evaluate
subset_sizes <- c(5, 10, 15, 20, 30, 50, 75, 100)

cat("  - Running SVM-RFE with LOOCV...\n")
cat("  - Feature subset sizes:", paste(subset_sizes, collapse = ", "), "\n")

# Run RFE with SVM
suppressWarnings({
  rfe_result <- rfe(
    x = X_rfe,
    y = y_subtype,
    sizes = subset_sizes,
    rfeControl = rfe_ctrl,
    method = "svmRadial"
  )
})

cat("  - SVM-RFE Results:\n")
cat("    Optimal features:", rfe_result$optsize, "\n")
cat("    Best accuracy:", sprintf("%.1f%%", max(rfe_result$results$Accuracy) * 100), "\n")

# Get optimal features
svm_rfe_features <- predictors(rfe_result)
cat("    Selected features:", length(svm_rfe_features), "\n")

# RFE accuracy curve data
rfe_curve_data <- rfe_result$results[, c("Variables", "Accuracy")]
colnames(rfe_curve_data) <- c("NumFeatures", "Accuracy")

# ============================================================================
# 6. SVM CLASSIFICATION WITH LOOCV (using optimal RFE features)
# ============================================================================
cat("\n[Step 5] SVM Classification with optimal features (LOOCV)...\n")

# Use RFE selected features
X_svm <- X_all[, svm_rfe_features[1:min(20, length(svm_rfe_features))]]

svm_preds <- rep(NA, n_samples)
svm_probs <- rep(NA, n_samples)

pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (i in 1:n_samples) {
  X_train <- X_svm[-i, , drop = FALSE]
  y_train <- y_subtype[-i]
  X_test <- X_svm[i, , drop = FALSE]
  
  svm_model <- svm(X_train, y_train, kernel = "radial", probability = TRUE)
  
  pred <- predict(svm_model, X_test, probability = TRUE)
  svm_preds[i] <- as.character(pred)
  svm_probs[i] <- attr(pred, "probabilities")[, "CS2"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Calculate SVM metrics
svm_accuracy <- mean(svm_preds == as.character(y_subtype))
svm_roc <- roc(as.numeric(y_subtype) - 1, svm_probs, quiet = TRUE)
svm_auc <- auc(svm_roc)

svm_cm <- confusionMatrix(factor(svm_preds, levels = levels(y_subtype)), y_subtype)
svm_f1 <- svm_cm$byClass["F1"]
if (is.na(svm_f1)) svm_f1 <- 0

cat("\n  - SVM Results:\n")
cat("    Accuracy:", sprintf("%.1f%%", svm_accuracy * 100), "\n")
cat("    AUC:", sprintf("%.3f", svm_auc), "\n")
cat("    F1 Score:", sprintf("%.3f", svm_f1), "\n")

# ============================================================================
# 7. LASSO FEATURE SELECTION (for comparison)
# ============================================================================
cat("\n[Step 6] LASSO Feature Selection with LOOCV...\n")

lasso_preds <- rep(NA, n_samples)
lasso_probs <- rep(NA, n_samples)
lasso_coef_list <- list()

pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (i in 1:n_samples) {
  X_train <- X_all[-i, , drop = FALSE]
  y_train <- as.numeric(y_subtype[-i]) - 1
  X_test <- X_all[i, , drop = FALSE]
  
  # Fit LASSO with internal CV for lambda
  cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
  
  # Predict
  prob <- predict(cv_lasso, newx = X_test, s = "lambda.min", type = "response")[1]
  lasso_probs[i] <- prob
  lasso_preds[i] <- ifelse(prob > 0.5, 1, 0)
  
  # Get coefficients
  coef_vec <- as.vector(coef(cv_lasso, s = "lambda.min"))[-1]
  names(coef_vec) <- colnames(X_train)
  lasso_coef_list[[i]] <- coef_vec
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Calculate LASSO metrics
lasso_accuracy <- mean(lasso_preds == y_numeric)
lasso_roc <- roc(y_numeric, lasso_probs, quiet = TRUE)
lasso_auc <- auc(lasso_roc)

lasso_cm <- confusionMatrix(factor(lasso_preds, levels = c(0, 1)), 
                            factor(y_numeric, levels = c(0, 1)))
lasso_f1 <- lasso_cm$byClass["F1"]
if (is.na(lasso_f1)) lasso_f1 <- 0

cat("\n  - LASSO Results:\n")
cat("    Accuracy:", sprintf("%.1f%%", lasso_accuracy * 100), "\n")
cat("    AUC:", sprintf("%.3f", lasso_auc), "\n")
cat("    F1 Score:", sprintf("%.3f", lasso_f1), "\n")

# Aggregate LASSO importance
lasso_coef_mat <- do.call(rbind, lasso_coef_list)
lasso_imp <- data.frame(
  Feature = colnames(lasso_coef_mat),
  AbsCoef = apply(abs(lasso_coef_mat), 2, mean)
) %>%
  arrange(desc(AbsCoef)) %>%
  filter(AbsCoef > 0) %>%
  head(20)

lasso_features <- lasso_imp$Feature

# ============================================================================
# 8. STACKING ENSEMBLE LEARNING
# ============================================================================
cat("\n[Step 7] Stacking Ensemble Learning...\n")

# Generate out-of-fold predictions from all base learners
stack_meta_X <- matrix(NA, nrow = n_samples, ncol = 3)
colnames(stack_meta_X) <- c("RF", "SVM", "LASSO")

# We already have LOOCV predictions from each model
stack_meta_X[, "RF"] <- rf_probs
stack_meta_X[, "SVM"] <- svm_probs
stack_meta_X[, "LASSO"] <- lasso_probs

# Meta-learner: Logistic Regression with LOOCV
stack_preds <- rep(NA, n_samples)
stack_probs <- rep(NA, n_samples)

pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (i in 1:n_samples) {
  X_train <- stack_meta_X[-i, , drop = FALSE]
  y_train <- y_numeric[-i]
  X_test <- stack_meta_X[i, , drop = FALSE]
  
  # Fit logistic regression as meta-learner
  meta_data <- data.frame(y = y_train, X_train)
  meta_model <- glm(y ~ ., data = meta_data, family = binomial)
  
  # Predict
  test_data <- data.frame(X_test)
  prob <- predict(meta_model, newdata = test_data, type = "response")
  stack_probs[i] <- prob
  stack_preds[i] <- ifelse(prob > 0.5, 1, 0)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Calculate Stacking metrics
stack_accuracy <- mean(stack_preds == y_numeric)
stack_roc <- roc(y_numeric, stack_probs, quiet = TRUE)
stack_auc <- auc(stack_roc)

stack_cm <- confusionMatrix(factor(stack_preds, levels = c(0, 1)), 
                            factor(y_numeric, levels = c(0, 1)))
stack_f1 <- stack_cm$byClass["F1"]
if (is.na(stack_f1)) stack_f1 <- 0

cat("\n  - Stacking Ensemble Results:\n")
cat("    Accuracy:", sprintf("%.1f%%", stack_accuracy * 100), "\n")
cat("    AUC:", sprintf("%.3f", stack_auc), "\n")
cat("    F1 Score:", sprintf("%.3f", stack_f1), "\n")

# ============================================================================
# 9. RESULTS SUMMARY AND COMPARISON
# ============================================================================
cat("\n[Step 8] Generating results summary...\n")

# Create comparison table
comparison_df <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM", "LASSO", "Stacking Ensemble"),
  Accuracy = c(rf_accuracy, xgb_accuracy, svm_accuracy, lasso_accuracy, stack_accuracy),
  AUC = c(as.numeric(rf_auc), as.numeric(xgb_auc), as.numeric(svm_auc), 
          as.numeric(lasso_auc), as.numeric(stack_auc)),
  F1 = c(rf_f1, xgb_f1, svm_f1, lasso_f1, stack_f1)
)

cat("\n  Model Performance Comparison:\n")
print(comparison_df)

# Save comparison table
write.csv(comparison_df, file.path(results_dir, "model_comparison.csv"), row.names = FALSE)

# ============================================================================
# 10. FEATURE OVERLAP ANALYSIS
# ============================================================================
cat("\n[Step 9] Feature overlap analysis...\n")

# Get top 20 features from each method
rf_top20 <- rf_imp_mean$Feature[1:20]
xgb_top20 <- xgb_imp_mean$Feature[1:20]
svm_top20 <- svm_rfe_features[1:min(20, length(svm_rfe_features))]
lasso_top20 <- lasso_features[1:min(20, length(lasso_features))]

# Calculate Jaccard indices
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_matrix <- matrix(NA, 4, 4)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- c("RF", "XGBoost", "SVM-RFE", "LASSO")
feature_lists <- list(RF = rf_top20, XGBoost = xgb_top20, `SVM-RFE` = svm_top20, LASSO = lasso_top20)

for (i in 1:4) {
  for (j in 1:4) {
    jaccard_matrix[i, j] <- jaccard(feature_lists[[i]], feature_lists[[j]])
  }
}

cat("\n  Feature Overlap (Jaccard Index):\n")
print(round(jaccard_matrix, 3))

# Save Jaccard matrix
write.csv(jaccard_matrix, file.path(results_dir, "feature_jaccard_overlap.csv"))

# ============================================================================
# 11. VISUALIZATIONS
# ============================================================================
cat("\n[Step 10] Generating visualizations...\n")

# Theme for publication
pub_theme <- theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# --- 11.1 Multi-classifier Performance Bar Plot ---
cat("  - Creating performance comparison bar plot...\n")

perf_long <- comparison_df %>%
  pivot_longer(cols = c(Accuracy, AUC, F1), names_to = "Metric", values_to = "Value") %>%
  mutate(Model = factor(Model, levels = c("Random Forest", "XGBoost", "SVM", "LASSO", "Stacking Ensemble")))

p_perf <- ggplot(perf_long, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Value)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(title = "HAE Subtype Classification: Model Performance Comparison",
       x = "", y = "Score", fill = "Metric") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(results_dir, "fig_performance_comparison.pdf"), p_perf, 
       width = 10, height = 6)

# --- 11.2 ROC Curves Overlay ---
cat("  - Creating ROC curves overlay...\n")

roc_data <- rbind(
  data.frame(Model = "Random Forest", 
             Sensitivity = rf_roc$sensitivities, 
             Specificity = 1 - rf_roc$specificities),
  data.frame(Model = "XGBoost", 
             Sensitivity = xgb_roc$sensitivities, 
             Specificity = 1 - xgb_roc$specificities),
  data.frame(Model = "SVM", 
             Sensitivity = svm_roc$sensitivities, 
             Specificity = 1 - svm_roc$specificities),
  data.frame(Model = "LASSO", 
             Sensitivity = lasso_roc$sensitivities, 
             Specificity = 1 - lasso_roc$specificities),
  data.frame(Model = "Stacking", 
             Sensitivity = stack_roc$sensitivities, 
             Specificity = 1 - stack_roc$specificities)
)

# Create AUC labels
auc_labels <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM", "LASSO", "Stacking"),
  AUC = c(rf_auc, xgb_auc, svm_auc, lasso_auc, stack_auc)
) %>%
  mutate(Label = sprintf("%s (AUC=%.3f)", Model, AUC))

roc_data <- roc_data %>%
  left_join(auc_labels[, c("Model", "Label")], by = "Model") %>%
  mutate(Model = Label)

p_roc <- ggplot(roc_data, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "ROC Curves: HAE Subtype Classification (LOOCV)",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)",
       color = "Model") +
  pub_theme +
  theme(legend.position = c(0.7, 0.25),
        legend.background = element_rect(fill = "white", color = "gray80"))

ggsave(file.path(results_dir, "fig_roc_curves.pdf"), p_roc, 
       width = 8, height = 7)

# --- 11.3 RFE Accuracy Curve ---
cat("  - Creating RFE accuracy curve...\n")

p_rfe <- ggplot(rfe_curve_data, aes(x = NumFeatures, y = Accuracy)) +
  geom_line(color = "#2E86AB", linewidth = 1.2) +
  geom_point(color = "#2E86AB", size = 3) +
  geom_vline(xintercept = rfe_result$optsize, linetype = "dashed", color = "red") +
  annotate("text", x = rfe_result$optsize + 5, y = max(rfe_curve_data$Accuracy) - 0.05,
           label = sprintf("Optimal: %d features", rfe_result$optsize), color = "red") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "SVM-RFE Feature Selection: Accuracy vs Number of Features",
       x = "Number of Features", y = "LOOCV Accuracy") +
  pub_theme

ggsave(file.path(results_dir, "fig_rfe_curve.pdf"), p_rfe, 
       width = 8, height = 6)

# --- 11.4 Feature Importance Heatmap ---
cat("  - Creating feature importance heatmap...\n")

# Get all unique top features
all_top_features <- unique(c(rf_top20, xgb_top20, svm_top20, lasso_top20))

# Create importance matrix
imp_matrix <- matrix(0, nrow = length(all_top_features), ncol = 4)
rownames(imp_matrix) <- all_top_features
colnames(imp_matrix) <- c("RF", "XGBoost", "SVM-RFE", "LASSO")

# Fill RF importance (rank-based)
rf_ranks <- setNames(1:20, rf_top20)
for (f in names(rf_ranks)) {
  if (f %in% rownames(imp_matrix)) {
    imp_matrix[f, "RF"] <- (21 - rf_ranks[f]) / 20
  }
}

# Fill XGBoost importance
xgb_ranks <- setNames(1:20, xgb_top20)
for (f in names(xgb_ranks)) {
  if (f %in% rownames(imp_matrix)) {
    imp_matrix[f, "XGBoost"] <- (21 - xgb_ranks[f]) / 20
  }
}

# Fill SVM-RFE importance
svm_ranks <- setNames(1:length(svm_top20), svm_top20)
for (f in names(svm_ranks)) {
  if (f %in% rownames(imp_matrix)) {
    imp_matrix[f, "SVM-RFE"] <- (21 - svm_ranks[f]) / 20
  }
}

# Fill LASSO importance
lasso_ranks <- setNames(1:length(lasso_top20), lasso_top20)
for (f in names(lasso_ranks)) {
  if (f %in% rownames(imp_matrix)) {
    imp_matrix[f, "LASSO"] <- (21 - lasso_ranks[f]) / 20
  }
}

# Select features that appear in at least 2 methods
feature_counts <- rowSums(imp_matrix > 0)
top_features <- names(sort(feature_counts, decreasing = TRUE))[1:min(30, sum(feature_counts >= 2))]

if (length(top_features) < 10) {
  top_features <- names(sort(feature_counts, decreasing = TRUE))[1:min(30, length(feature_counts))]
}

imp_matrix_subset <- imp_matrix[top_features, , drop = FALSE]

# Create heatmap
pdf(file.path(results_dir, "fig_feature_importance_heatmap.pdf"), width = 10, height = 12)
pheatmap(imp_matrix_subset,
         main = "Feature Importance Comparison Across Methods",
         color = colorRampPalette(c("white", "#FEE0D2", "#FC9272", "#DE2D26"))(50),
         cluster_cols = FALSE,
         fontsize_row = 8,
         fontsize_col = 12,
         border_color = "gray90")
dev.off()

# --- 11.5 Venn Diagram ---
cat("  - Creating Venn diagram...\n")

# Prepare feature lists (clean names for display)
venn_lists <- list(
  RF = rf_top20,
  XGBoost = xgb_top20,
  `SVM-RFE` = svm_top20,
  LASSO = lasso_top20
)

pdf(file.path(results_dir, "fig_feature_venn.pdf"), width = 10, height = 10)
venn.plot <- venn.diagram(
  x = venn_lists,
  filename = NULL,
  category.names = names(venn_lists),
  fill = brewer.pal(4, "Set2"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold",
  main = "Top 20 Feature Overlap Across Methods",
  main.cex = 1.5
)
grid.draw(venn.plot)
dev.off()

# ============================================================================
# 12. SAVE DETAILED RESULTS
# ============================================================================
cat("\n[Step 11] Saving detailed results...\n")

# Save XGBoost importance
write.csv(xgb_imp_mean, file.path(results_dir, "xgboost_top20_features.csv"), row.names = FALSE)

# Save RF importance
write.csv(rf_imp_mean, file.path(results_dir, "rf_top20_features.csv"), row.names = FALSE)

# Save SVM-RFE features
svm_rfe_df <- data.frame(
  Rank = 1:length(svm_rfe_features),
  Feature = svm_rfe_features
)
write.csv(svm_rfe_df, file.path(results_dir, "svm_rfe_selected_features.csv"), row.names = FALSE)

# Save LASSO features
write.csv(lasso_imp, file.path(results_dir, "lasso_top_features.csv"), row.names = FALSE)

# Save RFE curve data
write.csv(rfe_curve_data, file.path(results_dir, "rfe_accuracy_curve.csv"), row.names = FALSE)

# Save all predictions
predictions_df <- data.frame(
  Sample = adj_samples,
  Actual = as.character(y_subtype),
  RF_Pred = rf_preds,
  RF_Prob = rf_probs,
  XGBoost_Pred = ifelse(xgb_preds == 1, "CS2", "CS1"),
  XGBoost_Prob = xgb_probs,
  SVM_Pred = svm_preds,
  SVM_Prob = svm_probs,
  LASSO_Pred = ifelse(lasso_preds == 1, "CS2", "CS1"),
  LASSO_Prob = lasso_probs,
  Stacking_Pred = ifelse(stack_preds == 1, "CS2", "CS1"),
  Stacking_Prob = stack_probs
)
write.csv(predictions_df, file.path(results_dir, "all_predictions.csv"), row.names = FALSE)

# Save confusion matrices
cm_summary <- data.frame(
  Model = c("RF", "XGBoost", "SVM", "LASSO", "Stacking"),
  Sensitivity = c(rf_cm$byClass["Sensitivity"], 
                  xgb_cm$byClass["Sensitivity"],
                  svm_cm$byClass["Sensitivity"],
                  lasso_cm$byClass["Sensitivity"],
                  stack_cm$byClass["Sensitivity"]),
  Specificity = c(rf_cm$byClass["Specificity"],
                  xgb_cm$byClass["Specificity"],
                  svm_cm$byClass["Specificity"],
                  lasso_cm$byClass["Specificity"],
                  stack_cm$byClass["Specificity"]),
  PPV = c(rf_cm$byClass["Pos Pred Value"],
          xgb_cm$byClass["Pos Pred Value"],
          svm_cm$byClass["Pos Pred Value"],
          lasso_cm$byClass["Pos Pred Value"],
          stack_cm$byClass["Pos Pred Value"]),
  NPV = c(rf_cm$byClass["Neg Pred Value"],
          xgb_cm$byClass["Neg Pred Value"],
          svm_cm$byClass["Neg Pred Value"],
          lasso_cm$byClass["Neg Pred Value"],
          stack_cm$byClass["Neg Pred Value"])
)
write.csv(cm_summary, file.path(results_dir, "confusion_matrix_metrics.csv"), row.names = FALSE)

# ============================================================================
# 13. FINAL SUMMARY
# ============================================================================
cat("\n=============================================================================\n")
cat(" ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")

cat("Results saved to:", results_dir, "\n\n")

cat("OUTPUT FILES:\n")
cat("  Tables:\n")
cat("    - model_comparison.csv             : Overall performance metrics\n")
cat("    - feature_jaccard_overlap.csv      : Feature overlap analysis\n")
cat("    - xgboost_top20_features.csv       : XGBoost feature importance\n")
cat("    - rf_top20_features.csv            : Random Forest feature importance\n")
cat("    - svm_rfe_selected_features.csv    : SVM-RFE selected features\n")
cat("    - lasso_top_features.csv           : LASSO coefficients\n")
cat("    - rfe_accuracy_curve.csv           : RFE accuracy vs features\n")
cat("    - all_predictions.csv              : Sample-level predictions\n")
cat("    - confusion_matrix_metrics.csv     : Detailed classifier metrics\n")
cat("\n  Figures:\n")
cat("    - fig_performance_comparison.pdf   : Multi-model performance bar plot\n")
cat("    - fig_roc_curves.pdf               : ROC curves overlay\n")
cat("    - fig_rfe_curve.pdf                : SVM-RFE accuracy curve\n")
cat("    - fig_feature_importance_heatmap.pdf : Cross-method importance heatmap\n")
cat("    - fig_feature_venn.pdf             : Top 20 feature Venn diagram\n")

cat("\n\nKEY FINDINGS:\n")
cat(sprintf("  Best single model: %s (Accuracy=%.1f%%, AUC=%.3f)\n",
            comparison_df$Model[which.max(comparison_df$AUC)],
            comparison_df$Accuracy[which.max(comparison_df$AUC)] * 100,
            max(comparison_df$AUC)))
cat(sprintf("  Stacking ensemble: Accuracy=%.1f%%, AUC=%.3f\n",
            stack_accuracy * 100, stack_auc))
cat(sprintf("  Optimal SVM-RFE features: %d (Best accuracy=%.1f%%)\n",
            rfe_result$optsize, max(rfe_curve_data$Accuracy) * 100))

cat("\n[DONE] ML Ensemble Analysis completed successfully.\n")
