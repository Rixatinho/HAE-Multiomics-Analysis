#!/usr/bin/env Rscript
# ============================================================================
# Enhancement: Pathway-Level External Validation
# ============================================================================
# PURPOSE: Validate multi-omics findings at the pathway level rather than
#   individual gene level. Pathway signatures are more robust across
#   platforms and cohorts. This addresses the reviewer concern that
#   individual gene validation is circular.
#
# Strategy:
#   1. Derive pathway activity signatures from our HAE data (ssGSEA/GSVA)
#   2. Apply same pathway scoring to external dataset (GSE124362)
#   3. Correlate pathway activity patterns between cohorts
#   4. Test whether HAE-derived pathway signatures discriminate in external data
#   5. Compare Hallmark, KEGG, and Reactome pathway-level patterns
# ============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(msigdbr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
EXT_DIR  <- file.path(PROJECT, "analysis/results/enhancement9_external_validation")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhance_pathway_validation")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Pathway-Level External Validation\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA
# ============================================================================
cat(">>> 0. Loading data\n")

# Our transcriptomic data
degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
tc_vst <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                    check.names = FALSE, row.names = 1)
id2sym <- setNames(degs$gene_name, degs$gene_id)
id2sym <- id2sym[!is.na(id2sym) & id2sym != ""]
tc_sym <- id2sym[rownames(tc_vst)]
keep <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst[keep, ]); rownames(tc_mat) <- tc_sym[keep]

tc_group <- ifelse(grepl("^Normal", colnames(tc_mat)), "Normal", "Adjacent")
names(tc_group) <- colnames(tc_mat)

cat(sprintf("  Our TC data: %d genes x %d samples\n", nrow(tc_mat), ncol(tc_mat)))

# External data (if available)
ext_file <- file.path(EXT_DIR, "GSE124362_expr_matrix.csv")
ext_available <- file.exists(ext_file)
if (ext_available) {
  ext_mat <- as.matrix(read.csv(ext_file, check.names = FALSE, row.names = 1))
  ext_group_file <- file.path(EXT_DIR, "GSE124362_sample_groups.csv")
  if (file.exists(ext_group_file)) {
    ext_groups <- read.csv(ext_group_file, stringsAsFactors = FALSE)
    ext_group <- setNames(ext_groups$group, ext_groups$sample)
  } else {
    ext_group <- ifelse(grepl("normal|control|healthy", colnames(ext_mat), ignore.case = TRUE),
                        "Normal", "Disease")
    names(ext_group) <- colnames(ext_mat)
  }
  cat(sprintf("  External data: %d genes x %d samples\n", nrow(ext_mat), ncol(ext_mat)))
} else {
  cat("  [INFO] External dataset not found. Running internal pathway validation only.\n")
}


# ============================================================================
# 1. LOAD GENE SETS
# ============================================================================
cat("\n>>> 1. Loading gene set collections\n")

# Hallmark gene sets
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat(sprintf("  Hallmark: %d gene sets\n", length(hallmark_list)))

# KEGG pathways
kegg_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_list <- split(kegg_df$gene_symbol, kegg_df$gs_name)
cat(sprintf("  KEGG: %d gene sets\n", length(kegg_list)))

# Reactome
reactome_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_list <- split(reactome_df$gene_symbol, reactome_df$gs_name)
# Filter to reasonable size
reactome_list <- reactome_list[sapply(reactome_list, length) >= 15 &
                                sapply(reactome_list, length) <= 500]
cat(sprintf("  Reactome (filtered): %d gene sets\n", length(reactome_list)))


# ============================================================================
# 2. COMPUTE PATHWAY ACTIVITY SCORES (ssGSEA)
# ============================================================================
cat("\n>>> 2. Computing pathway activity scores\n")

run_ssgsea <- function(expr_mat, gene_sets, label) {
  cat(sprintf("  Running ssGSEA for %s...\n", label))
  param <- ssgseaParam(expr_mat, gene_sets, normalize = TRUE)
  scores <- gsva(param, verbose = FALSE)
  cat(sprintf("    %s: %d pathways x %d samples\n", label, nrow(scores), ncol(scores)))
  return(scores)
}

# Our data
hallmark_scores <- run_ssgsea(tc_mat, hallmark_list, "Hallmark (our data)")
kegg_scores     <- run_ssgsea(tc_mat, kegg_list, "KEGG (our data)")

# External data
if (ext_available) {
  hallmark_ext <- run_ssgsea(ext_mat, hallmark_list, "Hallmark (external)")
  kegg_ext     <- run_ssgsea(ext_mat, kegg_list, "KEGG (external)")
}


# ============================================================================
# 3. DIFFERENTIAL PATHWAY ANALYSIS (OUR DATA)
# ============================================================================
cat("\n>>> 3. Differential pathway analysis (our data)\n")

diff_pathways <- function(scores, group_vec, name) {
  adj_idx <- which(group_vec == "Adjacent")
  nor_idx <- which(group_vec == "Normal")
  
  # Paired patient IDs
  adj_names <- names(group_vec)[adj_idx]
  nor_names <- names(group_vec)[nor_idx]
  adj_ids <- as.integer(gsub("Adjacent|Normal", "", adj_names))
  nor_ids <- as.integer(gsub("Adjacent|Normal", "", nor_names))
  paired_ids <- intersect(adj_ids, nor_ids)
  adj_paired <- paste0("Adjacent", paired_ids)
  nor_paired <- paste0("Normal", paired_ids)
  
  results <- data.frame()
  for (pw in rownames(scores)) {
    adj_vals <- scores[pw, adj_paired]
    nor_vals <- scores[pw, nor_paired]
    wt <- tryCatch(wilcox.test(adj_vals, nor_vals, paired = TRUE), error = function(e) NULL)
    if (!is.null(wt)) {
      results <- rbind(results, data.frame(
        pathway = pw,
        mean_diff = mean(adj_vals - nor_vals),
        pvalue = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(results$pvalue), ]
  cat(sprintf("  %s: %d pathways tested, %d significant (padj<0.05)\n",
              name, nrow(results), sum(results$padj < 0.05)))
  return(results)
}

hallmark_diff <- diff_pathways(hallmark_scores, tc_group, "Hallmark")
kegg_diff     <- diff_pathways(kegg_scores, tc_group, "KEGG")

write.csv(hallmark_diff, file.path(OUT_DIR, "hallmark_diff_pathways.csv"), row.names = FALSE)
write.csv(kegg_diff, file.path(OUT_DIR, "kegg_diff_pathways.csv"), row.names = FALSE)


# ============================================================================
# 4. PATHWAY-LEVEL CROSS-COHORT CORRELATION
# ============================================================================
if (ext_available) {
  cat("\n>>> 4. Cross-cohort pathway correlation\n")
  
  # Compute mean pathway difference (Disease - Normal) in external data
  diff_ext <- function(scores, group_vec) {
    dis_idx <- which(group_vec != "Normal")
    nor_idx <- which(group_vec == "Normal")
    if (length(dis_idx) < 2 || length(nor_idx) < 2) return(NULL)
    rowMeans(scores[, dis_idx, drop = FALSE]) - rowMeans(scores[, nor_idx, drop = FALSE])
  }
  
  ext_hallmark_diff <- diff_ext(hallmark_ext, ext_group)
  our_hallmark_diff <- setNames(hallmark_diff$mean_diff, hallmark_diff$pathway)
  
  if (!is.null(ext_hallmark_diff)) {
    common_pw <- intersect(names(our_hallmark_diff), names(ext_hallmark_diff))
    if (length(common_pw) >= 10) {
      cor_test <- cor.test(our_hallmark_diff[common_pw], ext_hallmark_diff[common_pw],
                            method = "spearman")
      cat(sprintf("  Hallmark pathway correlation: rho=%.3f, P=%.4f (%d pathways)\n",
                  cor_test$estimate, cor_test$p.value, length(common_pw)))
      
      # Scatter plot
      cor_df <- data.frame(
        pathway = common_pw,
        our_diff = our_hallmark_diff[common_pw],
        ext_diff = ext_hallmark_diff[common_pw]
      )
      # Flag significant in our data
      sig_pw <- hallmark_diff$pathway[hallmark_diff$padj < 0.05]
      cor_df$significant <- cor_df$pathway %in% sig_pw
      cor_df$label <- ifelse(cor_df$significant, gsub("HALLMARK_", "", cor_df$pathway), "")
      
      p_cor <- ggplot(cor_df, aes(our_diff, ext_diff)) +
        geom_point(aes(color = significant), size = 2.5, alpha = 0.8) +
        geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.8) +
        geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
        scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#E64B35")) +
        labs(title = sprintf("Pathway-Level Cross-Cohort Concordance (rho=%.3f, P=%.2e)",
                              cor_test$estimate, cor_test$p.value),
             x = "Our cohort: mean pathway difference (Adjacent - Normal)",
             y = "External cohort: mean pathway difference (Disease - Normal)",
             color = "Significant\nin our data") +
        theme_bw(base_size = 11) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3)
      
      ggsave(file.path(FIG_DIR, "pathway_cross_cohort_correlation.pdf"), p_cor, width = 10, height = 8)
      cat("  Cross-cohort correlation plot saved.\n")
      
      # Save correlation result
      write.csv(cor_df, file.path(OUT_DIR, "pathway_cross_cohort_comparison.csv"), row.names = FALSE)
    }
  }
  
  # KEGG correlation
  ext_kegg_diff <- diff_ext(kegg_ext, ext_group)
  our_kegg_diff <- setNames(kegg_diff$mean_diff, kegg_diff$pathway)
  
  if (!is.null(ext_kegg_diff)) {
    common_kegg <- intersect(names(our_kegg_diff), names(ext_kegg_diff))
    if (length(common_kegg) >= 10) {
      cor_kegg <- cor.test(our_kegg_diff[common_kegg], ext_kegg_diff[common_kegg], method = "spearman")
      cat(sprintf("  KEGG pathway correlation: rho=%.3f, P=%.4f (%d pathways)\n",
                  cor_kegg$estimate, cor_kegg$p.value, length(common_kegg)))
    }
  }
} else {
  cat("\n>>> 4. [SKIP] External data not available for cross-cohort correlation\n")
}


# ============================================================================
# 5. INTERNAL SPLIT-HALF PATHWAY VALIDATION
# ============================================================================
cat("\n>>> 5. Internal split-half pathway validation\n")

# Split patients into discovery (60%) and validation (40%)
set.seed(42)
adj_idx <- which(tc_group == "Adjacent")
nor_idx <- which(tc_group == "Normal")
adj_ids <- as.integer(gsub("Adjacent", "", names(tc_group)[adj_idx]))
n_disc <- ceiling(length(adj_ids) * 0.6)

n_splits <- 100
concordance_rates <- numeric(n_splits)

for (s in 1:n_splits) {
  disc_ids <- sort(sample(adj_ids, n_disc))
  val_ids  <- setdiff(adj_ids, disc_ids)
  
  disc_samps <- c(paste0("Normal", disc_ids), paste0("Adjacent", disc_ids))
  val_samps  <- c(paste0("Normal", val_ids), paste0("Adjacent", val_ids))
  
  disc_samps <- intersect(disc_samps, colnames(hallmark_scores))
  val_samps  <- intersect(val_samps, colnames(hallmark_scores))
  
  if (length(disc_samps) < 4 || length(val_samps) < 4) next
  
  # Discovery: find significant pathways
  disc_group <- tc_group[disc_samps]
  disc_adj <- disc_samps[disc_group == "Adjacent"]
  disc_nor <- disc_samps[disc_group == "Normal"]
  
  disc_pvals <- apply(hallmark_scores[, disc_samps], 1, function(x) {
    tryCatch(wilcox.test(x[disc_adj], x[disc_nor])$p.value, error = function(e) NA)
  })
  disc_sig <- names(which(disc_pvals < 0.05))
  
  # Validation: check direction consistency
  if (length(disc_sig) >= 3) {
    val_group <- tc_group[val_samps]
    val_adj <- val_samps[val_group == "Adjacent"]
    val_nor <- val_samps[val_group == "Normal"]
    
    disc_dir <- sign(rowMeans(hallmark_scores[disc_sig, disc_adj, drop = FALSE]) -
                      rowMeans(hallmark_scores[disc_sig, disc_nor, drop = FALSE]))
    val_dir  <- sign(rowMeans(hallmark_scores[disc_sig, val_adj, drop = FALSE]) -
                      rowMeans(hallmark_scores[disc_sig, val_nor, drop = FALSE]))
    
    concordance_rates[s] <- mean(disc_dir == val_dir, na.rm = TRUE)
  }
}

concordance_rates <- concordance_rates[concordance_rates > 0]
mean_concordance <- mean(concordance_rates, na.rm = TRUE)
cat(sprintf("  Split-half concordance (100 splits): mean=%.3f, sd=%.3f\n",
            mean_concordance, sd(concordance_rates, na.rm = TRUE)))
cat(sprintf("  Interpretation: >0.7 = good reproducibility, >0.8 = excellent\n"))

write.csv(data.frame(concordance = concordance_rates),
          file.path(OUT_DIR, "split_half_concordance.csv"), row.names = FALSE)


# ============================================================================
# 6. VISUALIZATION
# ============================================================================
cat("\n>>> 6. Generating pathway validation figures\n")

# ---- Top differential pathways heatmap ----
top_hallmark <- head(hallmark_diff, 20)
top_pw <- top_hallmark$pathway

mat_plot <- hallmark_scores[top_pw, ]
# Order: Normal first, then Adjacent
samp_order <- c(names(tc_group)[tc_group == "Normal"], names(tc_group)[tc_group == "Adjacent"])
samp_order <- intersect(samp_order, colnames(mat_plot))
mat_plot <- mat_plot[, samp_order]

anno_col <- data.frame(Group = tc_group[samp_order], row.names = samp_order)
anno_colors <- list(Group = c(Normal = "#4DBBD5", Adjacent = "#E64B35"))

# Clean pathway names
clean_names <- gsub("HALLMARK_", "", rownames(mat_plot))
clean_names <- gsub("_", " ", clean_names)
rownames(mat_plot) <- clean_names

pdf(file.path(FIG_DIR, "top_hallmark_pathways_heatmap.pdf"), width = 12, height = 8)
pheatmap(mat_plot,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         cluster_cols = FALSE,
         fontsize_row = 8,
         main = "Top 20 Differential Hallmark Pathways (ssGSEA)")
dev.off()
cat("  Pathway heatmap saved.\n")

# ---- Split-half concordance histogram ----
p_conc <- ggplot(data.frame(x = concordance_rates), aes(x)) +
  geom_histogram(bins = 25, fill = "#3C5488", alpha = 0.8) +
  geom_vline(xintercept = mean_concordance, color = "#E64B35", linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = mean_concordance, y = Inf, vjust = 2,
           label = sprintf("Mean=%.3f", mean_concordance), color = "#E64B35", fontface = "bold") +
  labs(title = "Split-Half Pathway Direction Concordance",
       subtitle = "100 random 60/40 splits of patients",
       x = "Concordance rate (fraction of pathways with same direction)",
       y = "Frequency") +
  theme_bw(base_size = 12)

ggsave(file.path(FIG_DIR, "split_half_concordance.pdf"), p_conc, width = 8, height = 5)
cat("  Concordance histogram saved.\n")


# ============================================================================
# 7. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Pathway Validation SUMMARY\n")
cat("========================================\n")
cat(sprintf("Hallmark: %d/%d pathways significant (padj<0.05)\n",
            sum(hallmark_diff$padj < 0.05), nrow(hallmark_diff)))
cat(sprintf("KEGG: %d/%d pathways significant (padj<0.05)\n",
            sum(kegg_diff$padj < 0.05), nrow(kegg_diff)))
if (ext_available && exists("cor_test")) {
  cat(sprintf("Cross-cohort Hallmark correlation: rho=%.3f (P=%.2e)\n",
              cor_test$estimate, cor_test$p.value))
}
cat(sprintf("Split-half concordance: %.3f +/- %.3f\n",
            mean_concordance, sd(concordance_rates, na.rm = TRUE)))
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Pathway Validation COMPLETE.\n")
