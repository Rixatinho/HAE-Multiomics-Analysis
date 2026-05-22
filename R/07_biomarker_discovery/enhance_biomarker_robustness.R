################################################################################
# enhance_biomarker_robustness.R
# 
# Purpose: Biomarker robustness analysis with LOO-CV, parsimonious models,
#          and external validation for HAE multi-omics study
#
# Author: Multi-omics Analysis Pipeline
# Date: March 2026
################################################################################

# Set seed for reproducibility
set.seed(42)

# ==============================================================================
# SECTION 1: Library Loading and Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# Define project paths
project_root <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
data_dir <- file.path(project_root, "analysis/data/processed")
results_dir <- file.path(project_root, "analysis/results/phase8_biomarker")
output_dir <- file.path(project_root, "analysis/results/biomarker_robustness")
external_dir <- file.path(project_root, "analysis/external_data")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# NC theme style (using Helvetica as fallback for Arial)
theme_nc <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 9)
    )
}

message("=== Starting Biomarker Robustness Analysis ===\n")

# ==============================================================================
# SECTION 2: Load Data
# ==============================================================================

message("Loading data files...")

# Load transcriptomics data
tc_data <- read.csv(file.path(data_dir, "transcriptomics_logcpm_paired.csv"), 
                    row.names = 1, check.names = FALSE)
message(sprintf("Transcriptomics: %d features x %d samples", nrow(tc_data), ncol(tc_data)))

# Load proteomics data
pr_data <- read.csv(file.path(data_dir, "proteomics_log2_norm.csv"),
                    row.names = 1, check.names = FALSE)
message(sprintf("Proteomics: %d features x %d samples", nrow(pr_data), ncol(pr_data)))

# Load metabolomics data
met_data <- read.csv(file.path(data_dir, "metabolomics_log2_merged.csv"),
                     row.names = 1, check.names = FALSE)
message(sprintf("Metabolomics: %d features x %d samples", nrow(met_data), ncol(met_data)))

# Load biomarker information
biomarkers <- read.csv(file.path(results_dir, "disease_biomarkers.csv"))
panel <- read.csv(file.path(results_dir, "multiomics_diagnostic_panel.csv"))

message(sprintf("Biomarkers loaded: %d top features, %d panel features", 
                nrow(biomarkers), nrow(panel)))

# ==============================================================================
# SECTION 3: Identify Sample Structure
# ==============================================================================

message("\nAnalyzing sample structure...")

# Sample naming: Normal = Control, Adjacent = Disease (Tumor-adjacent)
extract_patient_id <- function(sample_name) {
  # Extract numeric ID from sample names like "Normal1", "Adjacent14"
  as.numeric(gsub("Normal|Adjacent", "", sample_name))
}

# Get patient IDs for each omics
tc_samples <- colnames(tc_data)
tc_normal <- tc_samples[grepl("^Normal", tc_samples)]
tc_adjacent <- tc_samples[grepl("^Adjacent", tc_samples)]
tc_normal_ids <- sapply(tc_normal, extract_patient_id)
tc_adjacent_ids <- sapply(tc_adjacent, extract_patient_id)
tc_paired_ids <- intersect(tc_normal_ids, tc_adjacent_ids)

pr_samples <- colnames(pr_data)
pr_normal <- pr_samples[grepl("^Normal", pr_samples)]
pr_adjacent <- pr_samples[grepl("^Adjacent", pr_samples)]
pr_normal_ids <- sapply(pr_normal, extract_patient_id)
pr_adjacent_ids <- sapply(pr_adjacent, extract_patient_id)
pr_paired_ids <- intersect(pr_normal_ids, pr_adjacent_ids)

met_samples <- colnames(met_data)
met_normal <- met_samples[grepl("^Normal", met_samples)]
met_adjacent <- met_samples[grepl("^Adjacent", met_samples)]
met_normal_ids <- sapply(met_normal, extract_patient_id)
met_adjacent_ids <- sapply(met_adjacent, extract_patient_id)
met_paired_ids <- intersect(met_normal_ids, met_adjacent_ids)

# Find patients with all three omics
complete_patients <- Reduce(intersect, list(tc_paired_ids, pr_paired_ids, met_paired_ids))

message(sprintf("Transcriptomics paired samples: %d patients", length(tc_paired_ids)))
message(sprintf("Proteomics paired samples: %d patients", length(pr_paired_ids)))
message(sprintf("Metabolomics paired samples: %d patients", length(met_paired_ids)))
message(sprintf("Complete 3-omics patients: %d patients", length(complete_patients)))

# ==============================================================================
# SECTION 4: Map Biomarker Names to Data Row Names
# ==============================================================================

message("\nMapping biomarker names to data row names...")

# Get top 10 biomarkers for each omics
tc_biomarkers <- biomarkers %>% filter(omics == "TC") %>% head(10)
pr_biomarkers <- biomarkers %>% filter(omics == "PR") %>% head(10)
met_biomarkers <- biomarkers %>% filter(omics == "MET") %>% head(10)

# Function to find matching row in data
find_feature_row <- function(feature_name, data_df, omics_type) {
  rownames_data <- rownames(data_df)
  
  if (omics_type == "MET") {
    # Direct match for metabolomics
    if (feature_name %in% rownames_data) {
      return(feature_name)
    }
  } else if (omics_type == "TC" || omics_type == "PR") {
    # For transcriptomics/proteomics, search by gene symbol in ENSG/ENSP IDs
    # Pattern: ENSG00000xxx_n or ENSP00000xxx.n (gene symbol may be in feature_name)
    pattern <- paste0("_", feature_name, "$|\\.", feature_name, "$")
    matches <- grep(pattern, rownames_data, value = TRUE, ignore.case = TRUE)
    
    if (length(matches) > 0) {
      return(matches[1])
    }
    
    # Also try gene symbol directly
    pattern2 <- paste0("^", feature_name, "$|_", feature_name, "_")
    matches2 <- grep(pattern2, rownames_data, value = TRUE, ignore.case = TRUE)
    if (length(matches2) > 0) {
      return(matches2[1])
    }
  }
  return(NA)
}

# Need to get actual row mappings from panel data (which has feature_id)
# Extract gene symbol to feature_id mapping from panel
panel_mapping <- panel %>%
  select(feature_id, feature_name) %>%
  mutate(omics_type = case_when(
    grepl("^ENSG", feature_id) ~ "TC",
    grepl("^ENSP", feature_id) ~ "PR",
    grepl("^Com_", feature_id) ~ "MET"
  ))

# Build comprehensive mapping by searching the data files
build_feature_mapping <- function(gene_symbol, omics_type, tc_data, pr_data, met_data) {
  if (omics_type == "MET") {
    # Metabolites: direct match with Com_xxx_pos/neg format
    if (gene_symbol %in% rownames(met_data)) {
      return(gene_symbol)
    }
    return(NA)
  }
  
  if (omics_type == "TC") {
    # Search transcriptomics data for gene symbol
    # Format: ENSG00000177807_10 where 10 might relate to gene symbol
    rownames_tc <- rownames(tc_data)
    
    # Check panel for this gene
    panel_row <- panel %>% filter(feature_name == gene_symbol & grepl("^ENSG", feature_id))
    if (nrow(panel_row) > 0) {
      fid <- panel_row$feature_id[1]
      if (fid %in% rownames_tc) return(fid)
    }
    return(NA)
  }
  
  if (omics_type == "PR") {
    # Search proteomics data
    rownames_pr <- rownames(pr_data)
    
    # Check panel for this gene
    panel_row <- panel %>% filter(feature_name == gene_symbol & grepl("^ENSP", feature_id))
    if (nrow(panel_row) > 0) {
      fid <- panel_row$feature_id[1]
      if (fid %in% rownames_pr) return(fid)
    }
    return(NA)
  }
  return(NA)
}

# For biomarkers in disease_biomarkers.csv, we need to find their IDs
# Since we don't have a direct mapping, we need to search

# Load the phase1_differential results which should have gene symbol mappings
search_gene_in_data <- function(gene_symbol, data_df) {
  # Gene symbols in our data may be embedded or we need external annotation
  # For now, return NA and handle later
  return(NA)
}

message("Note: Gene symbol to feature ID mapping requires annotation files.")
message("Using available panel feature IDs for analysis.\n")

# ==============================================================================
# SECTION 5: Prepare Analysis Data Matrices
# ==============================================================================

message("Preparing analysis data matrices...")

# Use the 16-feature LASSO panel for analysis
# Extract feature IDs from panel
panel_tc_features <- panel %>% filter(grepl("^ENSG", feature_id)) %>% pull(feature_id)
panel_pr_features <- panel %>% filter(grepl("^ENSP", feature_id)) %>% pull(feature_id)
panel_met_features <- panel %>% filter(grepl("^Com_", feature_id)) %>% pull(feature_id)

message(sprintf("Panel features - TC: %d, PR: %d, MET: %d",
                length(panel_tc_features), length(panel_pr_features), length(panel_met_features)))

# Check which features are available in data
available_tc <- panel_tc_features[panel_tc_features %in% rownames(tc_data)]
available_pr <- panel_pr_features[panel_pr_features %in% rownames(pr_data)]
available_met <- panel_met_features[panel_met_features %in% rownames(met_data)]

message(sprintf("Available features - TC: %d/%d, PR: %d/%d, MET: %d/%d",
                length(available_tc), length(panel_tc_features),
                length(available_pr), length(panel_pr_features),
                length(available_met), length(panel_met_features)))

# Create combined data matrix for complete patients only
# Use proteomics/metabolomics patients (14) since transcriptomics has only 12

# For LOO-CV analysis, use each omics separately first
# Then combine for multi-omics analysis

# Prepare single-omics data for LOO-CV
prepare_omics_data <- function(data_df, paired_ids, feature_list = NULL) {
  samples_normal <- paste0("Normal", paired_ids)
  samples_adjacent <- paste0("Adjacent", paired_ids)
  
  # Filter to existing columns
  samples_normal <- samples_normal[samples_normal %in% colnames(data_df)]
  samples_adjacent <- samples_adjacent[samples_adjacent %in% colnames(data_df)]
  
  # Get common patient IDs
  normal_ids <- as.numeric(gsub("Normal", "", samples_normal))
  adjacent_ids <- as.numeric(gsub("Adjacent", "", samples_adjacent))
  common_ids <- intersect(normal_ids, adjacent_ids)
  
  samples_normal <- paste0("Normal", common_ids)
  samples_adjacent <- paste0("Adjacent", common_ids)
  
  if (!is.null(feature_list)) {
    available_features <- feature_list[feature_list %in% rownames(data_df)]
    if (length(available_features) == 0) {
      return(NULL)
    }
    data_df <- data_df[available_features, , drop = FALSE]
  }
  
  # Combine Normal (0) and Adjacent (1) samples
  X <- cbind(data_df[, samples_normal, drop = FALSE], 
             data_df[, samples_adjacent, drop = FALSE])
  y <- c(rep(0, length(samples_normal)), rep(1, length(samples_adjacent)))
  patient_ids <- c(common_ids, common_ids)
  sample_type <- c(rep("Normal", length(samples_normal)), 
                   rep("Adjacent", length(samples_adjacent)))
  
  list(X = X, y = y, patient_ids = patient_ids, sample_type = sample_type,
       paired_ids = common_ids)
}

# Prepare data for each omics
tc_analysis <- prepare_omics_data(tc_data, 1:14)
pr_analysis <- prepare_omics_data(pr_data, 1:14)
met_analysis <- prepare_omics_data(met_data, 1:14)

message(sprintf("\nAnalysis datasets prepared:"))
message(sprintf("  Transcriptomics: %d samples from %d paired patients", 
                length(tc_analysis$y), length(tc_analysis$paired_ids)))
message(sprintf("  Proteomics: %d samples from %d paired patients",
                length(pr_analysis$y), length(pr_analysis$paired_ids)))
message(sprintf("  Metabolomics: %d samples from %d paired patients",
                length(met_analysis$y), length(met_analysis$paired_ids)))

# ==============================================================================
# SECTION 6: Single-Feature LOO-CV Analysis
# ==============================================================================

message("\n=== Single-Feature LOO-CV Analysis ===\n")

# Function to perform LOO-CV for a single feature using logistic regression
loo_cv_single_feature <- function(feature_values, y, patient_ids, feature_name) {
  n_patients <- length(unique(patient_ids))
  unique_patients <- unique(patient_ids)
  
  predictions <- numeric(length(y))
  
  for (i in seq_along(unique_patients)) {
    # Leave out one patient (both Normal and Adjacent samples)
    test_patient <- unique_patients[i]
    test_idx <- which(patient_ids == test_patient)
    train_idx <- which(patient_ids != test_patient)
    
    # Training data
    X_train <- feature_values[train_idx]
    y_train <- y[train_idx]
    
    # Fit logistic regression
    df_train <- data.frame(x = X_train, y = y_train)
    model <- tryCatch(
      glm(y ~ x, data = df_train, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    
    if (is.null(model)) {
      predictions[test_idx] <- 0.5
    } else {
      # Predict on test
      X_test <- feature_values[test_idx]
      df_test <- data.frame(x = X_test)
      predictions[test_idx] <- predict(model, newdata = df_test, type = "response")
    }
  }
  
  # Calculate ROC and AUC
  roc_obj <- tryCatch(
    roc(y, predictions, quiet = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(roc_obj)) {
    return(list(auc = NA, auc_ci_low = NA, auc_ci_high = NA, predictions = predictions))
  }
  
  auc_val <- as.numeric(auc(roc_obj))
  
  # Bootstrap CI
  ci_obj <- tryCatch(
    ci.auc(roc_obj, conf.level = 0.95, method = "bootstrap", boot.n = 2000, quiet = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(ci_obj)) {
    ci_low <- NA
    ci_high <- NA
  } else {
    ci_low <- ci_obj[1]
    ci_high <- ci_obj[3]
  }
  
  list(auc = auc_val, auc_ci_low = ci_low, auc_ci_high = ci_high, 
       predictions = predictions, roc = roc_obj)
}

# Analyze top features from each omics
# Using panel features since we have their IDs

# Get top features by coefficient magnitude from panel
panel_sorted <- panel %>%
  mutate(abs_coef = abs(coefficient)) %>%
  arrange(desc(abs_coef))

# Prepare features for single-feature analysis
# Top metabolomics: Com_288_pos (1-Methylxanthine)
# We'll use all panel features

single_feature_results <- list()

message("Analyzing Transcriptomics panel features...")
for (fid in available_tc) {
  fname <- panel %>% filter(feature_id == fid) %>% pull(feature_name)
  if (length(fname) == 0) fname <- fid
  
  # Get feature values
  samples <- colnames(tc_analysis$X)
  fvals <- as.numeric(tc_analysis$X[fid, ])
  
  result <- loo_cv_single_feature(fvals, tc_analysis$y, tc_analysis$patient_ids, fname)
  
  # Also calculate resubstitution AUC
  df_full <- data.frame(x = fvals, y = tc_analysis$y)
  model_full <- tryCatch(
    glm(y ~ x, data = df_full, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  if (!is.null(model_full)) {
    pred_full <- predict(model_full, type = "response")
    roc_full <- tryCatch(roc(tc_analysis$y, pred_full, quiet = TRUE), error = function(e) NULL)
    resub_auc <- if (!is.null(roc_full)) as.numeric(auc(roc_full)) else NA
  } else {
    resub_auc <- NA
  }
  
  single_feature_results[[fid]] <- data.frame(
    feature_id = fid,
    feature_name = fname,
    omics = "Transcriptomics",
    resub_auc = resub_auc,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  %s (%s): LOO-CV AUC = %.3f [%.3f-%.3f]", 
                  fid, fname, result$auc, result$auc_ci_low, result$auc_ci_high))
}

message("\nAnalyzing Proteomics panel features...")
for (fid in available_pr) {
  fname <- panel %>% filter(feature_id == fid) %>% pull(feature_name)
  if (length(fname) == 0) fname <- fid
  
  fvals <- as.numeric(pr_analysis$X[fid, ])
  result <- loo_cv_single_feature(fvals, pr_analysis$y, pr_analysis$patient_ids, fname)
  
  df_full <- data.frame(x = fvals, y = pr_analysis$y)
  model_full <- tryCatch(
    glm(y ~ x, data = df_full, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  resub_auc <- NA
  if (!is.null(model_full)) {
    pred_full <- predict(model_full, type = "response")
    roc_full <- tryCatch(roc(pr_analysis$y, pred_full, quiet = TRUE), error = function(e) NULL)
    resub_auc <- if (!is.null(roc_full)) as.numeric(auc(roc_full)) else NA
  }
  
  single_feature_results[[fid]] <- data.frame(
    feature_id = fid,
    feature_name = fname,
    omics = "Proteomics",
    resub_auc = resub_auc,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  %s (%s): LOO-CV AUC = %.3f [%.3f-%.3f]",
                  fid, fname, result$auc, result$auc_ci_low, result$auc_ci_high))
}

message("\nAnalyzing Metabolomics panel features...")
for (fid in available_met) {
  fname <- panel %>% filter(feature_id == fid) %>% pull(feature_name)
  if (length(fname) == 0) fname <- fid
  
  fvals <- as.numeric(met_analysis$X[fid, ])
  result <- loo_cv_single_feature(fvals, met_analysis$y, met_analysis$patient_ids, fname)
  
  df_full <- data.frame(x = fvals, y = met_analysis$y)
  model_full <- tryCatch(
    glm(y ~ x, data = df_full, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  resub_auc <- NA
  if (!is.null(model_full)) {
    pred_full <- predict(model_full, type = "response")
    roc_full <- tryCatch(roc(met_analysis$y, pred_full, quiet = TRUE), error = function(e) NULL)
    resub_auc <- if (!is.null(roc_full)) as.numeric(auc(roc_full)) else NA
  }
  
  single_feature_results[[fid]] <- data.frame(
    feature_id = fid,
    feature_name = fname,
    omics = "Metabolomics",
    resub_auc = resub_auc,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  %s (%s): LOO-CV AUC = %.3f [%.3f-%.3f]",
                  fid, fname, result$auc, result$auc_ci_low, result$auc_ci_high))
}

# Combine single feature results
single_feature_df <- bind_rows(single_feature_results)

# Save single feature results
write.csv(single_feature_df, file.path(output_dir, "single_feature_loo_auc.csv"), row.names = FALSE)
message("\nSingle feature results saved.")

# ==============================================================================
# SECTION 7: Multi-Feature LOO-CV Analysis
# ==============================================================================

message("\n=== Multi-Feature LOO-CV Analysis ===\n")

# Function to perform LOO-CV for multiple features
loo_cv_multi_feature <- function(X_matrix, y, patient_ids, method = "logistic") {
  n_patients <- length(unique(patient_ids))
  unique_patients <- unique(patient_ids)
  
  predictions <- numeric(length(y))
  
  for (i in seq_along(unique_patients)) {
    test_patient <- unique_patients[i]
    test_idx <- which(patient_ids == test_patient)
    train_idx <- which(patient_ids != test_patient)
    
    X_train <- as.matrix(t(X_matrix[, train_idx]))
    y_train <- y[train_idx]
    X_test <- as.matrix(t(X_matrix[, test_idx]))
    
    if (method == "logistic") {
      df_train <- data.frame(y = y_train, X_train)
      model <- tryCatch(
        glm(y ~ ., data = df_train, family = binomial(link = "logit")),
        error = function(e) NULL
      )
      if (is.null(model)) {
        predictions[test_idx] <- 0.5
      } else {
        df_test <- data.frame(X_test)
        colnames(df_test) <- colnames(df_train)[-1]
        predictions[test_idx] <- predict(model, newdata = df_test, type = "response")
      }
    } else if (method == "lasso") {
      # LASSO with internal CV for lambda selection
      cv_fit <- tryCatch(
        cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = min(5, nrow(X_train)-1)),
        error = function(e) NULL
      )
      if (is.null(cv_fit)) {
        predictions[test_idx] <- 0.5
      } else {
        predictions[test_idx] <- as.numeric(predict(cv_fit, newx = X_test, s = "lambda.min", type = "response"))
      }
    }
  }
  
  # Calculate ROC and AUC
  roc_obj <- tryCatch(
    roc(y, predictions, quiet = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(roc_obj)) {
    return(list(auc = NA, auc_ci_low = NA, auc_ci_high = NA, predictions = predictions))
  }
  
  auc_val <- as.numeric(auc(roc_obj))
  ci_obj <- tryCatch(
    ci.auc(roc_obj, conf.level = 0.95, method = "bootstrap", boot.n = 2000, quiet = TRUE),
    error = function(e) NULL
  )
  
  ci_low <- if (!is.null(ci_obj)) ci_obj[1] else NA
  ci_high <- if (!is.null(ci_obj)) ci_obj[3] else NA
  
  list(auc = auc_val, auc_ci_low = ci_low, auc_ci_high = ci_high, 
       predictions = predictions, roc = roc_obj)
}

# Calculate resubstitution AUC for multi-feature model
resub_auc_multi <- function(X_matrix, y, method = "logistic") {
  X_mat <- as.matrix(t(X_matrix))
  
  if (method == "logistic") {
    df <- data.frame(y = y, X_mat)
    model <- tryCatch(
      glm(y ~ ., data = df, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    if (is.null(model)) return(NA)
    pred <- predict(model, type = "response")
  } else if (method == "lasso") {
    cv_fit <- tryCatch(
      cv.glmnet(X_mat, y, family = "binomial", alpha = 1, nfolds = 5),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) return(NA)
    pred <- as.numeric(predict(cv_fit, newx = X_mat, s = "lambda.min", type = "response"))
  }
  
  roc_obj <- tryCatch(roc(y, pred, quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(NA)
  as.numeric(auc(roc_obj))
}

# Prepare combined data matrix for complete patients (those with all omics)
# Use 12 patients that have transcriptomics (which is the limiting factor)
complete_ids <- tc_analysis$paired_ids

# Build combined multi-omics matrix
build_combined_matrix <- function(tc_data, pr_data, met_data, patient_ids, 
                                  tc_features, pr_features, met_features) {
  samples_normal <- paste0("Normal", patient_ids)
  samples_adjacent <- paste0("Adjacent", patient_ids)
  
  # Filter features
  tc_available <- tc_features[tc_features %in% rownames(tc_data)]
  pr_available <- pr_features[pr_features %in% rownames(pr_data)]
  met_available <- met_features[met_features %in% rownames(met_data)]
  
  # Extract data for each omics
  tc_subset <- if(length(tc_available) > 0) tc_data[tc_available, c(samples_normal, samples_adjacent), drop=FALSE] else NULL
  pr_subset <- if(length(pr_available) > 0) pr_data[pr_available, c(samples_normal, samples_adjacent), drop=FALSE] else NULL
  met_subset <- if(length(met_available) > 0) met_data[met_available, c(samples_normal, samples_adjacent), drop=FALSE] else NULL
  
  # Combine matrices
  combined_list <- list()
  if (!is.null(tc_subset)) combined_list$tc <- tc_subset
  if (!is.null(pr_subset)) combined_list$pr <- pr_subset
  if (!is.null(met_subset)) combined_list$met <- met_subset
  
  if (length(combined_list) == 0) return(NULL)
  
  X <- do.call(rbind, combined_list)
  y <- c(rep(0, length(samples_normal)), rep(1, length(samples_adjacent)))
  patient_vec <- c(patient_ids, patient_ids)
  
  list(X = X, y = y, patient_ids = patient_vec)
}

# Build analysis matrices for different model configurations
model_results <- list()

# Model 1: Best single metabolite (Com_288_pos)
message("Model 1: Best single metabolite (1-Methylxanthine)...")
best_met_id <- "Com_288_pos"
if (best_met_id %in% rownames(met_data)) {
  met_single_data <- build_combined_matrix(tc_data, pr_data, met_data, met_analysis$paired_ids,
                                            character(0), character(0), best_met_id)
  result <- loo_cv_multi_feature(met_single_data$X, met_single_data$y, met_single_data$patient_ids, "logistic")
  resub <- resub_auc_multi(met_single_data$X, met_single_data$y, "logistic")
  model_results[["1_Met_single"]] <- data.frame(
    model = "1-feature (MET: 1-Methylxanthine)",
    n_features = 1,
    resub_auc = resub,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  Resub AUC: %.3f, LOO-CV AUC: %.3f [%.3f-%.3f]", 
                  resub, result$auc, result$auc_ci_low, result$auc_ci_high))
}

# Model 2: Best single proteomics (from panel with highest abs coef)
message("Model 2: Best single proteomics feature...")
best_pr_id <- panel %>% filter(grepl("^ENSP", feature_id)) %>% 
  arrange(desc(abs(coefficient))) %>% slice(1) %>% pull(feature_id)
if (length(best_pr_id) > 0 && best_pr_id %in% rownames(pr_data)) {
  pr_single_data <- build_combined_matrix(tc_data, pr_data, met_data, pr_analysis$paired_ids,
                                           character(0), best_pr_id, character(0))
  result <- loo_cv_multi_feature(pr_single_data$X, pr_single_data$y, pr_single_data$patient_ids, "logistic")
  resub <- resub_auc_multi(pr_single_data$X, pr_single_data$y, "logistic")
  fname <- panel %>% filter(feature_id == best_pr_id) %>% pull(feature_name)
  model_results[["2_PR_single"]] <- data.frame(
    model = paste0("1-feature (PR: ", fname, ")"),
    n_features = 1,
    resub_auc = resub,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  %s: Resub AUC: %.3f, LOO-CV AUC: %.3f [%.3f-%.3f]", 
                  fname, resub, result$auc, result$auc_ci_low, result$auc_ci_high))
}

# Model 3: Best single transcriptomics
message("Model 3: Best single transcriptomics feature...")
best_tc_id <- panel %>% filter(grepl("^ENSG", feature_id)) %>%
  arrange(desc(abs(coefficient))) %>% slice(1) %>% pull(feature_id)
if (length(best_tc_id) > 0 && best_tc_id %in% rownames(tc_data)) {
  tc_single_data <- build_combined_matrix(tc_data, pr_data, met_data, tc_analysis$paired_ids,
                                           best_tc_id, character(0), character(0))
  result <- loo_cv_multi_feature(tc_single_data$X, tc_single_data$y, tc_single_data$patient_ids, "logistic")
  resub <- resub_auc_multi(tc_single_data$X, tc_single_data$y, "logistic")
  fname <- panel %>% filter(feature_id == best_tc_id) %>% pull(feature_name)
  model_results[["3_TC_single"]] <- data.frame(
    model = paste0("1-feature (TC: ", fname, ")"),
    n_features = 1,
    resub_auc = resub,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  %s: Resub AUC: %.3f, LOO-CV AUC: %.3f [%.3f-%.3f]",
                  fname, resub, result$auc, result$auc_ci_low, result$auc_ci_high))
}

# Model 4: 2-feature model (best Met + best PR)
message("Model 4: 2-feature cross-omics (Met + PR)...")
combo2_data <- build_combined_matrix(tc_data, pr_data, met_data, pr_analysis$paired_ids,
                                      character(0), best_pr_id, best_met_id)
if (!is.null(combo2_data)) {
  result <- loo_cv_multi_feature(combo2_data$X, combo2_data$y, combo2_data$patient_ids, "logistic")
  resub <- resub_auc_multi(combo2_data$X, combo2_data$y, "logistic")
  model_results[["4_Met_PR_2feat"]] <- data.frame(
    model = "2-feature (MET + PR)",
    n_features = 2,
    resub_auc = resub,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  Resub AUC: %.3f, LOO-CV AUC: %.3f [%.3f-%.3f]",
                  resub, result$auc, result$auc_ci_low, result$auc_ci_high))
}

# Model 5: 3-feature model (1 from each omics)
message("Model 5: 3-feature cross-omics (1 per omics)...")
combo3_data <- build_combined_matrix(tc_data, pr_data, met_data, tc_analysis$paired_ids,
                                      best_tc_id, best_pr_id, best_met_id)
if (!is.null(combo3_data)) {
  result <- loo_cv_multi_feature(combo3_data$X, combo3_data$y, combo3_data$patient_ids, "logistic")
  resub <- resub_auc_multi(combo3_data$X, combo3_data$y, "logistic")
  model_results[["5_3omics_3feat"]] <- data.frame(
    model = "3-feature (1 per omics)",
    n_features = 3,
    resub_auc = resub,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  Resub AUC: %.3f, LOO-CV AUC: %.3f [%.3f-%.3f]",
                  resub, result$auc, result$auc_ci_low, result$auc_ci_high))
}

# Model 6: Full 16-feature LASSO panel
message("Model 6: Full 16-feature LASSO panel...")
full_panel_data <- build_combined_matrix(tc_data, pr_data, met_data, tc_analysis$paired_ids,
                                          available_tc, available_pr, available_met)
if (!is.null(full_panel_data) && nrow(full_panel_data$X) >= 3) {
  result <- loo_cv_multi_feature(full_panel_data$X, full_panel_data$y, 
                                  full_panel_data$patient_ids, "lasso")
  resub <- resub_auc_multi(full_panel_data$X, full_panel_data$y, "lasso")
  model_results[["6_Full_LASSO"]] <- data.frame(
    model = "16-feature LASSO",
    n_features = nrow(full_panel_data$X),
    resub_auc = resub,
    loo_cv_auc = result$auc,
    ci_low = result$auc_ci_low,
    ci_high = result$auc_ci_high
  )
  message(sprintf("  %d features: Resub AUC: %.3f, LOO-CV AUC: %.3f [%.3f-%.3f]",
                  nrow(full_panel_data$X), resub, result$auc, result$auc_ci_low, result$auc_ci_high))
}

# Combine all model results
model_comparison_df <- bind_rows(model_results)
write.csv(model_comparison_df, file.path(output_dir, "parsimonious_model_comparison.csv"), row.names = FALSE)
message("\nModel comparison results saved.")

# ==============================================================================
# SECTION 8: GSE124362 External Validation
# ==============================================================================

message("\n=== GSE124362 External Validation ===\n")

# Read GSE124362 series matrix
gse_file <- file.path(external_dir, "GSE124362_series_matrix.txt.gz")

if (file.exists(gse_file)) {
  message("Reading GSE124362 series matrix...")
  
  # Read the file
  gse_lines <- readLines(gzfile(gse_file))
  
  # Find table start
  table_start <- grep("^!series_matrix_table_begin", gse_lines)
  table_end <- grep("^!series_matrix_table_end", gse_lines)
  
  if (length(table_start) > 0 && length(table_end) > 0) {
    # Extract data section
    data_lines <- gse_lines[(table_start + 1):(table_end - 1)]
    
    # Parse header and data
    header <- strsplit(data_lines[1], "\t")[[1]]
    header <- gsub('"', '', header)
    
    # Parse expression data
    expr_data <- lapply(data_lines[-1], function(line) {
      vals <- strsplit(line, "\t")[[1]]
      vals <- gsub('"', '', vals)
      return(vals)
    })
    
    expr_matrix <- do.call(rbind, expr_data)
    rownames(expr_matrix) <- expr_matrix[, 1]
    expr_matrix <- expr_matrix[, -1, drop = FALSE]
    colnames(expr_matrix) <- header[-1]
    
    # Convert to numeric
    expr_numeric <- apply(expr_matrix, 2, as.numeric)
    rownames(expr_numeric) <- rownames(expr_matrix)
    
    message(sprintf("GSE124362 expression matrix: %d probes x %d samples", 
                    nrow(expr_numeric), ncol(expr_numeric)))
    
    # Parse sample information
    # From metadata: rep1-1 = periparasitic (disease), rep1-2 = distal (control)
    # Samples ending in -1 are disease, -2 are control
    sample_info <- data.frame(
      sample_id = colnames(expr_numeric),
      group = ifelse(grepl("-1$", sapply(strsplit(colnames(expr_numeric), "_"), function(x) x[1])) |
                       grepl("Periparasitic", colnames(expr_numeric)),
                     "Disease", "Control")
    )
    
    # More robust: check source name from metadata
    source_line <- gse_lines[grep("^!Sample_source_name", gse_lines)[1]]
    if (!is.null(source_line)) {
      source_parts <- strsplit(source_line, "\t")[[1]]
      source_parts <- gsub('"', '', source_parts)[-1]
      sample_info$source <- source_parts
      sample_info$group <- ifelse(grepl("Periparasitic", sample_info$source), "Disease", "Control")
    }
    
    message("Sample groups:")
    print(table(sample_info$group))
    
    # Note: GSE124362 uses Arraystar lncRNA chip (GPL16956)
    # This chip may not have direct gene symbol mappings for our mRNA markers
    # We'll attempt to find matches based on probe annotations
    
    message("\nNote: GSE124362 uses Arraystar lncRNA microarray (GPL16956).")
    message("This platform primarily measures lncRNAs, not mRNAs.")
    message("Direct validation of our mRNA biomarkers may not be possible.\n")
    
    # Save validation attempt results
    gse_validation <- data.frame(
      dataset = "GSE124362",
      platform = "GPL16956 (Arraystar lncRNA)",
      n_samples = ncol(expr_numeric),
      n_disease = sum(sample_info$group == "Disease"),
      n_control = sum(sample_info$group == "Control"),
      validation_status = "Platform mismatch - lncRNA chip cannot validate mRNA markers",
      notes = "GSE124362 measures lncRNAs, our biomarkers are mRNAs/proteins/metabolites"
    )
    
  } else {
    message("Could not parse GSE124362 series matrix.")
    gse_validation <- data.frame(
      dataset = "GSE124362",
      validation_status = "Parse error",
      notes = "Could not locate expression matrix in series_matrix file"
    )
  }
} else {
  message("GSE124362 file not found.")
  gse_validation <- data.frame(
    dataset = "GSE124362",
    validation_status = "File not found",
    notes = paste("Expected file:", gse_file)
  )
}

write.csv(gse_validation, file.path(output_dir, "gse124362_independent_validation.csv"), row.names = FALSE)
message("GSE124362 validation results saved.")

# ==============================================================================
# SECTION 9: Compile LOO-CV Results
# ==============================================================================

message("\n=== Compiling LOO-CV Results ===\n")

# Combine all LOO-CV results into single table
single_feature_for_merge <- single_feature_df %>% 
  mutate(model_type = "Single Feature", 
         model = paste0(omics, ": ", feature_name),
         n_features = 1) %>%
  select(model, model_type, n_features, resub_auc, loo_cv_auc, ci_low, ci_high)

multi_feature_for_merge <- model_comparison_df %>%
  mutate(model_type = "Multi-Feature") %>%
  select(model, model_type, n_features, resub_auc, loo_cv_auc, ci_low, ci_high)

loo_results <- bind_rows(single_feature_for_merge, multi_feature_for_merge)

write.csv(loo_results, file.path(output_dir, "loo_cv_results.csv"), row.names = FALSE)
message("Complete LOO-CV results saved.")

# ==============================================================================
# SECTION 10: Visualization
# ==============================================================================

message("\n=== Generating Visualizations ===\n")

# Figure 1: LOO-CV ROC curves
message("Creating LOO-CV ROC curves...")

# Recalculate ROC objects for plotting
roc_objects <- list()

# Best single metabolite
if (best_met_id %in% rownames(met_data)) {
  met_single_data <- build_combined_matrix(tc_data, pr_data, met_data, met_analysis$paired_ids,
                                            character(0), character(0), best_met_id)
  result <- loo_cv_multi_feature(met_single_data$X, met_single_data$y, met_single_data$patient_ids, "logistic")
  if (!is.null(result$roc)) {
    roc_objects[["1-Methylxanthine"]] <- result$roc
  }
}

# Best proteomics
if (length(best_pr_id) > 0 && best_pr_id %in% rownames(pr_data)) {
  pr_single_data <- build_combined_matrix(tc_data, pr_data, met_data, pr_analysis$paired_ids,
                                           character(0), best_pr_id, character(0))
  result <- loo_cv_multi_feature(pr_single_data$X, pr_single_data$y, pr_single_data$patient_ids, "logistic")
  fname <- panel %>% filter(feature_id == best_pr_id) %>% pull(feature_name)
  if (!is.null(result$roc)) {
    roc_objects[[fname]] <- result$roc
  }
}

# Best transcriptomics
if (length(best_tc_id) > 0 && best_tc_id %in% rownames(tc_data)) {
  tc_single_data <- build_combined_matrix(tc_data, pr_data, met_data, tc_analysis$paired_ids,
                                           best_tc_id, character(0), character(0))
  result <- loo_cv_multi_feature(tc_single_data$X, tc_single_data$y, tc_single_data$patient_ids, "logistic")
  fname <- panel %>% filter(feature_id == best_tc_id) %>% pull(feature_name)
  if (!is.null(result$roc)) {
    roc_objects[[fname]] <- result$roc
  }
}

# 3-feature model
if (!is.null(combo3_data)) {
  result <- loo_cv_multi_feature(combo3_data$X, combo3_data$y, combo3_data$patient_ids, "logistic")
  if (!is.null(result$roc)) {
    roc_objects[["3-feature"]] <- result$roc
  }
}

# Full LASSO model
if (!is.null(full_panel_data) && nrow(full_panel_data$X) >= 3) {
  result <- loo_cv_multi_feature(full_panel_data$X, full_panel_data$y, 
                                  full_panel_data$patient_ids, "lasso")
  if (!is.null(result$roc)) {
    roc_objects[["16-feature LASSO"]] <- result$roc
  }
}

# Create ROC plot
if (length(roc_objects) > 0) {
  pdf(file.path(output_dir, "Fig_biomarker_loo_roc.pdf"), width = 6, height = 6)
  
  # Color palette
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
  
  # Plot first ROC
  plot(roc_objects[[1]], col = colors[1], lwd = 2,
       main = "LOO-CV ROC Curves", 
       legacy.axes = TRUE,
       xlab = "False Positive Rate",
       ylab = "True Positive Rate")
  
  # Add remaining ROCs
  for (i in 2:length(roc_objects)) {
    lines(roc_objects[[i]], col = colors[i], lwd = 2)
  }
  
  # Add legend
  legend_text <- sapply(names(roc_objects), function(name) {
    auc_val <- as.numeric(auc(roc_objects[[name]]))
    paste0(name, " (AUC = ", sprintf("%.3f", auc_val), ")")
  })
  
  legend("bottomright", legend = legend_text, col = colors[1:length(roc_objects)],
         lwd = 2, cex = 0.8, bty = "n")
  
  dev.off()
  message("  Fig_biomarker_loo_roc.pdf saved.")
}

# Figure 2: Model comparison bar plot
message("Creating model comparison plot...")

if (nrow(model_comparison_df) > 0) {
  # Prepare data for plotting
  plot_data <- model_comparison_df %>%
    select(model, resub_auc, loo_cv_auc) %>%
    pivot_longer(cols = c(resub_auc, loo_cv_auc), 
                 names_to = "metric", 
                 values_to = "AUC") %>%
    mutate(metric = factor(metric, 
                           levels = c("resub_auc", "loo_cv_auc"),
                           labels = c("Resubstitution", "LOO-CV")))
  
  p_comparison <- ggplot(plot_data, aes(x = model, y = AUC, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Resubstitution" = "#E41A1C", "LOO-CV" = "#377EB8")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = "", y = "AUC", fill = "Metric",
         title = "Model Performance: Resubstitution vs LOO-CV") +
    theme_nc() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "top")
  
  ggsave(file.path(output_dir, "Fig_biomarker_model_comparison.pdf"), 
         p_comparison, width = 8, height = 6)
  message("  Fig_biomarker_model_comparison.pdf saved.")
}

# Figure 3: GSE124362 validation - REMOVED (platform mismatch, text-only figure not suitable)
# The GSE124362 dataset uses Arraystar lncRNA microarray (GPL16956) which is incompatible
# with our mRNA/protein/metabolite markers. Instead of a text-only figure, this panel
# has been removed from the supplementary figure layout.
message("Skipping GSE124362 validation figure (platform mismatch - lncRNA chip)...")
message("  Note: Panel c removed from SuppFig_21 layout due to platform incompatibility.")

# ==============================================================================
# SECTION 11: Summary Report
# ==============================================================================

message("\n=== Analysis Summary ===\n")

# Print summary statistics
message("Key Findings:")
message("--------------")

if (nrow(model_comparison_df) > 0) {
  # Best parsimonious model by LOO-CV
  best_model <- model_comparison_df %>% 
    arrange(desc(loo_cv_auc)) %>% 
    slice(1)
  
  message(sprintf("\nBest performing model (by LOO-CV AUC):"))
  message(sprintf("  Model: %s", best_model$model))
  message(sprintf("  LOO-CV AUC: %.3f [%.3f-%.3f]", 
                  best_model$loo_cv_auc, best_model$ci_low, best_model$ci_high))
  message(sprintf("  Resubstitution AUC: %.3f", best_model$resub_auc))
  message(sprintf("  Optimism (Resub - LOO): %.3f", 
                  best_model$resub_auc - best_model$loo_cv_auc))
  
  # Overfitting assessment
  message("\nOverfitting Assessment:")
  for (i in 1:nrow(model_comparison_df)) {
    row <- model_comparison_df[i, ]
    optimism <- row$resub_auc - row$loo_cv_auc
    message(sprintf("  %s: Optimism = %.3f %s", 
                    row$model, optimism,
                    ifelse(optimism > 0.3, "(HIGH)", 
                           ifelse(optimism > 0.15, "(MODERATE)", "(LOW)"))))
  }
}

# Print output files
message("\n\nOutput Files:")
message("-------------")
output_files <- list.files(output_dir, full.names = FALSE)
for (f in output_files) {
  message(sprintf("  %s", f))
}

message("\n=== Biomarker Robustness Analysis Complete ===\n")
