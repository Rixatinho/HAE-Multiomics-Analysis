#!/usr/bin/env Rscript
# ============================================================================
# Phase 0: Data Preprocessing & Quality Control
# HAE Multi-omics Integration Study
# ============================================================================
# Input:  Raw omics data from Novogene sequencing
# Output: Normalized matrices + QC report + PCA plots
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(vsn)
  library(impute)
  library(sva)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(pheatmap)
})

# ---- Paths ----
DATA_ROOT <- "/Volumes/阿热-实验数据备份/工作/1.肝包虫/1.泡型包虫/1.2024-11-泡球蚴"
PROJECT   <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
OUT_DIR   <- file.path(PROJECT, "analysis/results/phase0_qc")
DATA_DIR  <- file.path(PROJECT, "analysis/data/processed")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Color scheme (Nature style) ----
col_group <- c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")

cat("========================================\n")
cat("Phase 0: Data Preprocessing & QC\n")
cat("========================================\n\n")

# ============================================================================
# 1. TRANSCRIPTOMICS PREPROCESSING
# ============================================================================
cat(">>> 1. Transcriptomics Preprocessing\n")

# Read raw counts
tc_file <- file.path(DATA_ROOT, "1.转录组学测序/03.Result_X101SC24101695-Z02-F001_homo_sapiens/Result_X101SC24101695-Z02-F001_homo_sapiens/4.Quant/1.Count/gene_count.xls")
tc_raw <- read.delim(tc_file, check.names = FALSE)
rownames(tc_raw) <- tc_raw$gene_id

# Extract sample columns only (Normal*/Adjacent*)
meta_cols <- c("gene_id", "gene_name", "gene_chr", "gene_start", "gene_end",
               "gene_strand", "gene_length", "gene_biotype", "gene_description", "Family")
sample_cols <- grep("^(Normal|Adjacent)", colnames(tc_raw), value = TRUE)
tc_counts <- tc_raw[, sample_cols]

# Identify Normal and Adjacent samples
normal_samples_tc  <- sort(grep("^Normal", colnames(tc_counts), value = TRUE))
adjacent_samples_tc <- sort(grep("^Adjacent", colnames(tc_counts), value = TRUE))
cat(sprintf("  Normal samples: %d (%s)\n", length(normal_samples_tc),
            paste(normal_samples_tc, collapse = ", ")))
cat(sprintf("  Adjacent samples: %d (%s)\n", length(adjacent_samples_tc),
            paste(adjacent_samples_tc, collapse = ", ")))

# Note missing samples
all_expected <- paste0("Adjacent", 1:14)
missing_adj <- setdiff(all_expected, adjacent_samples_tc)
if (length(missing_adj) > 0) {
  cat(sprintf("  WARNING: Missing Adjacent samples in transcriptomics: %s\n",
              paste(missing_adj, collapse = ", ")))
}

# Build sample metadata for DESeq2
extract_id <- function(x) as.integer(gsub("Normal|Adjacent", "", x))
sample_info_tc <- data.frame(
  sample = colnames(tc_counts),
  group  = ifelse(grepl("^Normal", colnames(tc_counts)), "Normal", "Adjacent"),
  patient_id = extract_id(colnames(tc_counts)),
  stringsAsFactors = FALSE
)
rownames(sample_info_tc) <- sample_info_tc$sample
sample_info_tc$group <- factor(sample_info_tc$group, levels = c("Normal", "Adjacent"))
sample_info_tc$patient_id <- factor(sample_info_tc$patient_id)

# Filter: keep protein-coding genes with reasonable expression
if ("gene_biotype" %in% colnames(tc_raw)) {
  protein_coding <- rownames(tc_raw)[tc_raw$gene_biotype == "protein_coding"]
  cat(sprintf("  Protein-coding genes: %d / %d total\n",
              length(protein_coding), nrow(tc_counts)))
  tc_counts_pc <- tc_counts[intersect(protein_coding, rownames(tc_counts)), ]
} else {
  tc_counts_pc <- tc_counts
}

# Filter low-expression genes: at least 5 counts in >= 50% of samples per group
keep_genes <- apply(tc_counts_pc, 1, function(x) {
  n_normal   <- sum(x[normal_samples_tc] >= 5)
  n_adjacent <- sum(x[adjacent_samples_tc] >= 5)
  (n_normal >= length(normal_samples_tc) * 0.5) |
    (n_adjacent >= length(adjacent_samples_tc) * 0.5)
})
tc_counts_filt <- tc_counts_pc[keep_genes, ]
cat(sprintf("  After filtering: %d genes retained\n", nrow(tc_counts_filt)))

# Ensure integer counts
tc_counts_filt <- round(tc_counts_filt)

# DESeq2 normalization (paired design)
# Only use paired samples for DESeq2
paired_ids_tc <- intersect(
  extract_id(normal_samples_tc),
  extract_id(adjacent_samples_tc)
)
paired_samples <- c(paste0("Normal", paired_ids_tc), paste0("Adjacent", paired_ids_tc))
tc_paired <- tc_counts_filt[, paired_samples]
si_paired <- sample_info_tc[paired_samples, ]

dds <- DESeqDataSetFromMatrix(
  countData = tc_paired,
  colData   = si_paired,
  design    = ~ patient_id + group
)
dds <- DESeq(dds)

# VST transformation for downstream analysis
vst_data <- vst(dds, blind = FALSE)
tc_vst <- assay(vst_data)
cat(sprintf("  VST-normalized matrix: %d genes x %d samples\n",
            nrow(tc_vst), ncol(tc_vst)))

# Also save full normalized counts (not just paired)
dds_full <- DESeqDataSetFromMatrix(
  countData = tc_counts_filt,
  colData   = sample_info_tc[colnames(tc_counts_filt), ],
  design    = ~ group
)
dds_full <- estimateSizeFactors(dds_full)
tc_norm_counts <- counts(dds_full, normalized = TRUE)

# Save
write.csv(tc_vst, file.path(DATA_DIR, "transcriptomics_vst_paired.csv"))
write.csv(tc_norm_counts, file.path(DATA_DIR, "transcriptomics_norm_counts.csv"))
saveRDS(dds, file.path(DATA_DIR, "transcriptomics_dds_paired.rds"))

# PCA plot
pca_tc <- prcomp(t(tc_vst), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_tc$x[, 1], PC2 = pca_tc$x[, 2],
  group = si_paired$group, patient = si_paired$patient_id,
  sample = rownames(si_paired)
)
var_explained <- round(100 * summary(pca_tc)$importance[2, 1:2], 1)

p_pca_tc <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.5, max.overlaps = 20) +
  scale_color_manual(values = col_group) +
  labs(
    title = "Transcriptomics PCA (VST, paired samples)",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "PCA_transcriptomics.pdf"), p_pca_tc,
       width = 7, height = 6)
cat("  PCA plot saved.\n\n")

# ============================================================================
# 2. PROTEOMICS PREPROCESSING
# ============================================================================
cat(">>> 2. Proteomics Preprocessing\n")

pr_file <- file.path(DATA_ROOT, "2.蛋白质组学/Result-X101SC24101695_Z01_J001_B1_43/3.ProtExprQuantification/all_sample.xls")
pr_raw <- read.delim(pr_file, check.names = FALSE)

# Extract annotation and expression columns
pr_anno_cols <- c("Protein", "Description", "Gene")
pr_expr <- pr_raw[, setdiff(colnames(pr_raw), pr_anno_cols)]
pr_anno <- pr_raw[, pr_anno_cols]
rownames(pr_expr) <- pr_raw$Protein

normal_samples_pr  <- sort(grep("^Normal", colnames(pr_expr), value = TRUE))
adjacent_samples_pr <- sort(grep("^Adjacent", colnames(pr_expr), value = TRUE))
cat(sprintf("  Normal samples: %d, Adjacent samples: %d\n",
            length(normal_samples_pr), length(adjacent_samples_pr)))
cat(sprintf("  Raw proteins: %d\n", nrow(pr_expr)))

# Convert to numeric, handle any non-numeric
pr_expr <- apply(pr_expr, 2, as.numeric)
rownames(pr_expr) <- pr_raw$Protein

# Missing value analysis
missing_frac <- rowMeans(is.na(pr_expr) | pr_expr == 0)
cat(sprintf("  Proteins with >50%% missing: %d\n", sum(missing_frac > 0.5)))

# Filter: keep proteins with <50% missing across all samples
pr_filt <- pr_expr[missing_frac <= 0.5, ]
cat(sprintf("  After missing filter: %d proteins\n", nrow(pr_filt)))

# Log2 transformation (replace 0 with NA first)
pr_filt[pr_filt == 0] <- NA
pr_log2 <- log2(pr_filt)

# Imputation: KNN for proteins with few missing values
# For proteins with many missing, use minimum value
na_frac_row <- rowMeans(is.na(pr_log2))
# Proteins with low missing rate (<30%): KNN imputation
low_missing <- na_frac_row < 0.3 & na_frac_row > 0
high_missing <- na_frac_row >= 0.3

cat(sprintf("  KNN imputation for %d proteins (<%d%% missing)\n",
            sum(low_missing), 30))
cat(sprintf("  MinProb imputation for %d proteins (>=%d%% missing)\n",
            sum(high_missing), 30))

if (sum(low_missing) > 0) {
  knn_result <- impute.knn(as.matrix(pr_log2), k = 3)
  pr_imputed <- knn_result$data
} else {
  pr_imputed <- as.matrix(pr_log2)
}

# For remaining NAs (high missing rate), impute with column-wise minimum - 1.8 SD
for (j in 1:ncol(pr_imputed)) {
  col_vals <- pr_imputed[!is.na(pr_imputed[, j]), j]
  min_val <- mean(col_vals) - 1.8 * sd(col_vals)
  pr_imputed[is.na(pr_imputed[, j]), j] <- min_val
}

# Median centering normalization
col_medians <- apply(pr_imputed, 2, median, na.rm = TRUE)
global_median <- median(col_medians)
pr_norm <- sweep(pr_imputed, 2, col_medians - global_median)

cat(sprintf("  Normalized matrix: %d proteins x %d samples\n",
            nrow(pr_norm), ncol(pr_norm)))

# Build sample info for proteomics
sample_info_pr <- data.frame(
  sample = colnames(pr_norm),
  group  = ifelse(grepl("^Normal", colnames(pr_norm)), "Normal", "Adjacent"),
  patient_id = extract_id(colnames(pr_norm)),
  stringsAsFactors = FALSE
)
rownames(sample_info_pr) <- sample_info_pr$sample

# Map gene names
gene_map <- setNames(pr_raw$Gene, pr_raw$Protein)

# Save
write.csv(pr_norm, file.path(DATA_DIR, "proteomics_log2_norm.csv"))
write.csv(pr_raw[rownames(pr_norm), pr_anno_cols],
          file.path(DATA_DIR, "proteomics_annotation.csv"))

# PCA
pca_pr <- prcomp(t(pr_norm), scale. = TRUE)
pca_pr_df <- data.frame(
  PC1 = pca_pr$x[, 1], PC2 = pca_pr$x[, 2],
  group = sample_info_pr$group, patient = sample_info_pr$patient_id,
  sample = sample_info_pr$sample
)
var_pr <- round(100 * summary(pca_pr)$importance[2, 1:2], 1)

p_pca_pr <- ggplot(pca_pr_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.5, max.overlaps = 20) +
  scale_color_manual(values = col_group) +
  labs(
    title = "Proteomics PCA (log2, median-centered)",
    x = paste0("PC1 (", var_pr[1], "%)"),
    y = paste0("PC2 (", var_pr[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "PCA_proteomics.pdf"), p_pca_pr,
       width = 7, height = 6)
cat("  PCA plot saved.\n\n")

# ============================================================================
# 3. METABOLOMICS PREPROCESSING
# ============================================================================
cat(">>> 3. Metabolomics Preprocessing\n")

met_root <- file.path(DATA_ROOT, "3.非靶标代谢组学/Result-X101SC24101695-Z03-J001-B1-42/1.MetQuant-QC")

# Read positive and negative ion mode data
met_pos <- as.data.frame(readxl::read_excel(file.path(met_root, "meta_intensity_pos.xlsx")))
met_neg <- as.data.frame(readxl::read_excel(file.path(met_root, "meta_intensity_neg.xlsx")))

cat(sprintf("  Positive ion: %d metabolites\n", nrow(met_pos)))
cat(sprintf("  Negative ion: %d metabolites\n", nrow(met_neg)))

# Annotation columns
met_anno_cols <- c("Compound_ID", "Name", "ChineseName", "IonMode", "Formula",
                   "MolecularWeight", "m/z", "MassError", "Adduct", "RT (min)",
                   "Score", "Level", "Column", "ClassI", "ClassI (Chinese)",
                   "ClassII", "ClassII (Chinese)", "ClassIII", "ClassIII (Chinese)",
                   "CAS", "HMDB_ID", "KEGG_ID", "KEGG_MapID", "Lipidmaps_ID",
                   "PubChemID", "SMILES", "InChIKey")

# Extract sample columns
pos_sample_cols <- grep("^pos_(Adjacent|Normal)", colnames(met_pos), value = TRUE)
neg_sample_cols <- grep("^neg_(Adjacent|Normal)", colnames(met_neg), value = TRUE)
pos_qc_cols <- grep("^pos_QC", colnames(met_pos), value = TRUE)
neg_qc_cols <- grep("^neg_QC", colnames(met_neg), value = TRUE)

cat(sprintf("  POS sample cols: %d, QC cols: %d\n", length(pos_sample_cols), length(pos_qc_cols)))
cat(sprintf("  NEG sample cols: %d, QC cols: %d\n", length(neg_sample_cols), length(neg_qc_cols)))

# QC check: coefficient of variation in QC samples
if (length(pos_qc_cols) > 0) {
  qc_pos <- met_pos[, pos_qc_cols]
  qc_pos <- apply(qc_pos, 2, as.numeric)
  cv_pos <- apply(qc_pos, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)
  cat(sprintf("  POS QC CV: median=%.1f%%, <30%%: %d/%d (%.1f%%)\n",
              median(cv_pos, na.rm = TRUE),
              sum(cv_pos < 30, na.rm = TRUE), length(cv_pos),
              100 * sum(cv_pos < 30, na.rm = TRUE) / length(cv_pos)))
}

if (length(neg_qc_cols) > 0) {
  qc_neg <- met_neg[, neg_qc_cols]
  qc_neg <- apply(qc_neg, 2, as.numeric)
  cv_neg <- apply(qc_neg, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)
  cat(sprintf("  NEG QC CV: median=%.1f%%, <30%%: %d/%d (%.1f%%)\n",
              median(cv_neg, na.rm = TRUE),
              sum(cv_neg < 30, na.rm = TRUE), length(cv_neg),
              100 * sum(cv_neg < 30, na.rm = TRUE) / length(cv_neg)))
}

# Process positive ion mode
pos_expr <- apply(met_pos[, pos_sample_cols], 2, as.numeric)
rownames(pos_expr) <- met_pos$Compound_ID
# Standardize column names: pos_Adjacent1 -> Adjacent1
colnames(pos_expr) <- gsub("^pos_", "", colnames(pos_expr))

# Process negative ion mode
neg_expr <- apply(met_neg[, neg_sample_cols], 2, as.numeric)
rownames(neg_expr) <- met_neg$Compound_ID
colnames(neg_expr) <- gsub("^neg_", "", colnames(neg_expr))

# Merge positive and negative (ensure same column order)
common_samples <- intersect(colnames(pos_expr), colnames(neg_expr))
cat(sprintf("  Common samples across pos/neg: %d\n", length(common_samples)))
met_merged <- rbind(pos_expr[, common_samples], neg_expr[, common_samples])

# Merge annotations
pos_anno <- met_pos[, intersect(met_anno_cols, colnames(met_pos))]
neg_anno <- met_neg[, intersect(met_anno_cols, colnames(met_neg))]
met_anno <- rbind(pos_anno, neg_anno)
rownames(met_anno) <- c(met_pos$Compound_ID, met_neg$Compound_ID)

cat(sprintf("  Merged metabolites: %d\n", nrow(met_merged)))

# Filter: remove metabolites with >50% missing
met_missing <- rowMeans(is.na(met_merged) | met_merged == 0)
met_filt <- met_merged[met_missing <= 0.5, ]
cat(sprintf("  After missing filter: %d metabolites\n", nrow(met_filt)))

# Replace 0 with NA, then log2 transform
met_filt[met_filt == 0] <- NA

# Impute metabolomics: unified KNN(k=3) + MinProb(mean-1.8SD) strategy
# Consistent with proteomics imputation (QC-verified)
met_filt_pre <- met_filt  # still in raw scale with NAs
met_log2_pre <- log2(met_filt_pre)
na_frac_met <- rowMeans(is.na(met_log2_pre))
low_miss_met <- na_frac_met < 0.3 & na_frac_met > 0
high_miss_met <- na_frac_met >= 0.3

cat(sprintf("  Metabolomics KNN imputation: %d metabolites (<30%% missing)\n", sum(low_miss_met)))
cat(sprintf("  Metabolomics MinProb imputation: %d metabolites (>=30%% missing)\n", sum(high_miss_met)))

if (any(!is.na(met_log2_pre) & is.finite(met_log2_pre))) {
  knn_result_met <- impute.knn(as.matrix(met_log2_pre), k = 3)
  met_log2 <- knn_result_met$data
} else {
  met_log2 <- as.matrix(met_log2_pre)
}

# For remaining NAs (high missing rate), impute with MinProb: column-wise mean - 1.8*SD
for (j in 1:ncol(met_log2)) {
  na_idx_j <- is.na(met_log2[, j])
  if (any(na_idx_j)) {
    col_vals <- met_log2[!na_idx_j, j]
    min_val <- mean(col_vals) - 1.8 * sd(col_vals)
    met_log2[na_idx_j, j] <- min_val
  }
}

# Pareto scaling (column-wise centering, row-wise sqrt(sd) scaling)
met_centered <- sweep(met_log2, 1, rowMeans(met_log2))
row_sds <- apply(met_log2, 1, sd)
met_pareto <- sweep(met_centered, 1, sqrt(row_sds), "/")

# Also keep log2 version (more commonly used for DE)
# For differential analysis, we use log2; for clustering, we use Pareto
cat(sprintf("  Final normalized matrix (log2): %d x %d\n",
            nrow(met_log2), ncol(met_log2)))

# Sample info
normal_samples_met  <- sort(grep("^Normal", colnames(met_log2), value = TRUE))
adjacent_samples_met <- sort(grep("^Adjacent", colnames(met_log2), value = TRUE))

sample_info_met <- data.frame(
  sample = colnames(met_log2),
  group  = ifelse(grepl("^Normal", colnames(met_log2)), "Normal", "Adjacent"),
  patient_id = extract_id(colnames(met_log2)),
  stringsAsFactors = FALSE
)
rownames(sample_info_met) <- sample_info_met$sample

# Save
write.csv(met_log2, file.path(DATA_DIR, "metabolomics_log2_merged.csv"))
write.csv(met_pareto, file.path(DATA_DIR, "metabolomics_pareto_merged.csv"))
write.csv(met_anno[rownames(met_log2), ], file.path(DATA_DIR, "metabolomics_annotation.csv"))

# PCA
pca_met <- prcomp(t(met_log2), scale. = TRUE)
pca_met_df <- data.frame(
  PC1 = pca_met$x[, 1], PC2 = pca_met$x[, 2],
  group = sample_info_met$group, patient = sample_info_met$patient_id,
  sample = sample_info_met$sample
)
var_met <- round(100 * summary(pca_met)$importance[2, 1:2], 1)

p_pca_met <- ggplot(pca_met_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.5, max.overlaps = 20) +
  scale_color_manual(values = col_group) +
  labs(
    title = "Metabolomics PCA (log2, merged pos+neg)",
    x = paste0("PC1 (", var_met[1], "%)"),
    y = paste0("PC2 (", var_met[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "PCA_metabolomics.pdf"), p_pca_met,
       width = 7, height = 6)
cat("  PCA plot saved.\n\n")

# ============================================================================
# 4. SAMPLE CORRESPONDENCE TABLE
# ============================================================================
cat(">>> 4. Sample Correspondence Summary\n")

# Identify samples available in each omics
all_patients <- 1:14
sample_table <- data.frame(
  Patient = all_patients,
  Transcriptomics_Normal = paste0("Normal", all_patients) %in% normal_samples_tc,
  Transcriptomics_Adjacent = paste0("Adjacent", all_patients) %in% adjacent_samples_tc,
  Proteomics_Normal = paste0("Normal", all_patients) %in% normal_samples_pr,
  Proteomics_Adjacent = paste0("Adjacent", all_patients) %in% adjacent_samples_pr,
  Metabolomics_Normal = paste0("Normal", all_patients) %in% normal_samples_met,
  Metabolomics_Adjacent = paste0("Adjacent", all_patients) %in% adjacent_samples_met
)
sample_table$Complete_3omics <- sample_table$Transcriptomics_Adjacent &
  sample_table$Proteomics_Adjacent &
  sample_table$Metabolomics_Adjacent

cat("Sample availability matrix:\n")
print(sample_table)
cat(sprintf("\nPatients with complete 3-omics pairs: %d / 14\n",
            sum(sample_table$Complete_3omics)))

write.csv(sample_table, file.path(DATA_DIR, "sample_correspondence.csv"),
          row.names = FALSE)

# ============================================================================
# 5. COMBINED QC VISUALIZATION
# ============================================================================
cat("\n>>> 5. Combined QC Visualization\n")

# Sample correlation heatmap (per omics)
for (omics_name in c("Transcriptomics", "Proteomics", "Metabolomics")) {
  if (omics_name == "Transcriptomics") {
    mat <- tc_vst
  } else if (omics_name == "Proteomics") {
    mat <- pr_norm
  } else {
    mat <- as.matrix(met_log2)
  }

  # Sample-sample correlation
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")

  # Annotation
  ann_df <- data.frame(
    Group = ifelse(grepl("^Normal", colnames(mat)), "Normal", "Adjacent"),
    row.names = colnames(mat)
  )
  ann_colors <- list(Group = col_group)

  pdf(file.path(FIG_DIR, paste0("correlation_heatmap_", omics_name, ".pdf")),
      width = 8, height = 7)
  pheatmap(
    cor_mat,
    annotation_col = ann_df,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = paste0(omics_name, " - Sample Correlation (Spearman)"),
    fontsize = 8,
    display_numbers = FALSE
  )
  dev.off()
  cat(sprintf("  %s correlation heatmap saved.\n", omics_name))
}

# ============================================================================
# 6. SUMMARY REPORT
# ============================================================================
cat("\n========================================\n")
cat("Phase 0 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Transcriptomics:\n"))
cat(sprintf("  Input:  %d genes x %d samples\n", nrow(tc_raw), length(sample_cols)))
cat(sprintf("  Output: %d genes x %d paired samples (VST)\n", nrow(tc_vst), ncol(tc_vst)))
cat(sprintf("  Paired patients: %s\n", paste(paired_ids_tc, collapse = ", ")))
cat(sprintf("\nProteomics:\n"))
cat(sprintf("  Input:  %d proteins x %d samples\n", nrow(pr_raw), ncol(pr_expr)))
cat(sprintf("  Output: %d proteins x %d samples (log2, median-centered)\n",
            nrow(pr_norm), ncol(pr_norm)))
cat(sprintf("\nMetabolomics:\n"))
cat(sprintf("  Input:  %d (pos) + %d (neg) = %d metabolites\n",
            nrow(met_pos), nrow(met_neg), nrow(met_pos) + nrow(met_neg)))
cat(sprintf("  Output: %d metabolites x %d samples (log2)\n",
            nrow(met_log2), ncol(met_log2)))
cat(sprintf("\nMulti-omics integration: %d patients with complete 3-omics pairs\n",
            sum(sample_table$Complete_3omics)))
cat(sprintf("Output directory: %s\n", DATA_DIR))
cat("========================================\n")
cat("Phase 0 COMPLETE.\n")
