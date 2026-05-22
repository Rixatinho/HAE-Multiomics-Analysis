#!/usr/bin/env Rscript
# ============================================================================
# Phase 7: Multi-dimensional Subtype Characterization
# HAE Multi-omics Integration Study
# ============================================================================
# Strategy:
#   1. Subtype-specific pathway profiling (GSVA scoring per subtype)
#   2. Immune microenvironment differences between subtypes
#   3. Metabolic landscape per subtype
#   4. Key molecular markers for each subtype
#   5. Comprehensive oncoplot-style summary heatmap
# ============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(fgsea)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
P2_DIR   <- file.path(PROJECT, "analysis/results/phase2_enrichment")
P4_DIR   <- file.path(PROJECT, "analysis/results/phase4_immune")
P6_DIR   <- file.path(PROJECT, "analysis/results/phase6_subtyping")
GMT_DIR  <- file.path(PROJECT, "analysis/data/gmt")
OUT_DIR  <- file.path(PROJECT, "analysis/results/phase7_characterization")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 7: Subtype Characterization\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA & SUBTYPE ASSIGNMENTS
# ============================================================================
cat(">>> 0. Loading data and subtype assignments\n")

# Subtype assignments (K=2)
subtypes <- read.csv(file.path(P6_DIR, "subtype_K2.csv"), stringsAsFactors = FALSE)
subtype_vec <- setNames(subtypes$subtype, subtypes$sample)
patient_subtypes <- setNames(subtypes$subtype, subtypes$patient_id)
cat(sprintf("  Subtypes: CS1=%d, CS2=%d\n",
            sum(subtype_vec == "CS1"), sum(subtype_vec == "CS2")))

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
met_mat <- as.matrix(read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"),
                               check.names = FALSE, row.names = 1))

# Match to Adjacent samples
adj_samples <- names(subtype_vec)
tc_adj <- tc_mat[, intersect(adj_samples, colnames(tc_mat))]
pr_adj <- pr_mat[, intersect(adj_samples, colnames(pr_mat))]
met_adj <- met_mat[, intersect(adj_samples, colnames(met_mat))]

cat(sprintf("  TC: %d genes x %d samples\n", nrow(tc_adj), ncol(tc_adj)))
cat(sprintf("  PR: %d proteins x %d samples\n", nrow(pr_adj), ncol(pr_adj)))
cat(sprintf("  MET: %d metabolites x %d samples\n", nrow(met_adj), ncol(met_adj)))


# ============================================================================
# 1. PATHWAY GSVA SCORING PER SUBTYPE
# ============================================================================
cat("\n>>> 1. Pathway GSVA scoring per subtype\n")

# Load Hallmark gene sets
hallmark_gmt <- file.path(GMT_DIR, "h.all.v2024.1.Hs.symbols.gmt")
if (file.exists(hallmark_gmt)) {
  hallmark_sets <- gmtPathways(hallmark_gmt)
  cat(sprintf("  Hallmark gene sets loaded: %d\n", length(hallmark_sets)))
} else {
  cat("  [WARN] Hallmark GMT not found, downloading...\n")
  BASE_URL <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs"
  download.file(paste0(BASE_URL, "/h.all.v2024.1.Hs.symbols.gmt"),
                hallmark_gmt, quiet = TRUE)
  hallmark_sets <- gmtPathways(hallmark_gmt)
}

# Run ssGSEA on transcriptomics (Adjacent samples only)
cat("  Running ssGSEA on transcriptomics...\n")
tc_param <- ssgseaParam(tc_adj, hallmark_sets, normalize = TRUE)
tc_gsva <- gsva(tc_param, verbose = FALSE)

# Run ssGSEA on proteomics
cat("  Running ssGSEA on proteomics...\n")
pr_param <- ssgseaParam(pr_adj, hallmark_sets, normalize = TRUE)
pr_gsva <- gsva(pr_param, verbose = FALSE)

cat(sprintf("  TC ssGSEA: %d pathways x %d samples\n", nrow(tc_gsva), ncol(tc_gsva)))
cat(sprintf("  PR ssGSEA: %d pathways x %d samples\n", nrow(pr_gsva), ncol(pr_gsva)))

# ---- Test pathway scores between subtypes ----
test_pathway_subtypes <- function(gsva_mat, subtype_vec, omics_name) {
  results <- data.frame()
  cs1 <- names(subtype_vec)[subtype_vec == "CS1"]
  cs2 <- names(subtype_vec)[subtype_vec == "CS2"]
  cs1 <- intersect(cs1, colnames(gsva_mat))
  cs2 <- intersect(cs2, colnames(gsva_mat))
  
  for (pw in rownames(gsva_mat)) {
    v1 <- gsva_mat[pw, cs1]
    v2 <- gsva_mat[pw, cs2]
    tryCatch({
      wt <- wilcox.test(v1, v2)
      results <- rbind(results, data.frame(
        pathway = pw,
        omics = omics_name,
        mean_CS1 = mean(v1, na.rm = TRUE),
        mean_CS2 = mean(v2, na.rm = TRUE),
        diff_CS1_CS2 = mean(v1, na.rm = TRUE) - mean(v2, na.rm = TRUE),
        pvalue = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(results$pvalue), ]
  return(results)
}

pw_tc <- test_pathway_subtypes(tc_gsva, subtype_vec, "transcriptomics")
pw_pr <- test_pathway_subtypes(pr_gsva, subtype_vec, "proteomics")

cat(sprintf("  TC pathways sig (P<0.05): %d / %d\n", sum(pw_tc$pvalue < 0.05), nrow(pw_tc)))
cat(sprintf("  PR pathways sig (P<0.05): %d / %d\n", sum(pw_pr$pvalue < 0.05), nrow(pw_pr)))

# Show top pathways
cat("\n  Top subtype-differential Hallmark pathways (TC):\n")
for (i in 1:min(10, nrow(pw_tc))) {
  cat(sprintf("    %s: CS1=%.3f, CS2=%.3f, P=%.4f\n",
              gsub("HALLMARK_", "", pw_tc$pathway[i]),
              pw_tc$mean_CS1[i], pw_tc$mean_CS2[i], pw_tc$pvalue[i]))
}

pw_all <- rbind(pw_tc, pw_pr)
write.csv(pw_all, file.path(OUT_DIR, "subtype_pathway_scores.csv"), row.names = FALSE)


# ---- Pathway heatmap by subtype ----
# Select top variable pathways across subtypes
top_pw_tc <- head(pw_tc$pathway, 20)
top_pw_pr <- head(pw_pr$pathway, 20)

# Combined pathway heatmap (TC)
if (length(top_pw_tc) >= 5) {
  mat_pw <- tc_gsva[top_pw_tc, ]
  
  # Clean pathway names
  rownames(mat_pw) <- gsub("HALLMARK_", "", rownames(mat_pw))
  rownames(mat_pw) <- gsub("_", " ", rownames(mat_pw))
  
  # Order by subtype
  samp_order <- c(
    sort(names(subtype_vec)[subtype_vec == "CS1"]),
    sort(names(subtype_vec)[subtype_vec == "CS2"])
  )
  samp_order <- intersect(samp_order, colnames(mat_pw))
  mat_pw <- mat_pw[, samp_order]
  
  anno_col <- data.frame(Subtype = subtype_vec[samp_order])
  rownames(anno_col) <- samp_order
  
  pdf(file.path(FIG_DIR, "pathway_heatmap_TC.pdf"), width = 10, height = 8)
  pheatmap(mat_pw,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           annotation_col = anno_col,
           annotation_colors = list(Subtype = c("CS1" = "#3C5488", "CS2" = "#E64B35")),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 8,
           main = "Hallmark Pathway Scores by Subtype (Transcriptomics)",
           scale = "row")
  dev.off()
  cat("  Pathway heatmap (TC) saved.\n")
}

# PR pathway heatmap
if (length(top_pw_pr) >= 5) {
  mat_pw_pr <- pr_gsva[top_pw_pr, ]
  rownames(mat_pw_pr) <- gsub("HALLMARK_", "", rownames(mat_pw_pr))
  rownames(mat_pw_pr) <- gsub("_", " ", rownames(mat_pw_pr))
  samp_order_pr <- intersect(samp_order, colnames(mat_pw_pr))
  mat_pw_pr <- mat_pw_pr[, samp_order_pr]
  
  anno_col_pr <- data.frame(Subtype = subtype_vec[samp_order_pr])
  rownames(anno_col_pr) <- samp_order_pr
  
  pdf(file.path(FIG_DIR, "pathway_heatmap_PR.pdf"), width = 10, height = 8)
  pheatmap(mat_pw_pr,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           annotation_col = anno_col_pr,
           annotation_colors = list(Subtype = c("CS1" = "#3C5488", "CS2" = "#E64B35")),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 8,
           main = "Hallmark Pathway Scores by Subtype (Proteomics)",
           scale = "row")
  dev.off()
  cat("  Pathway heatmap (PR) saved.\n")
}


# ============================================================================
# 2. IMMUNE MICROENVIRONMENT BY SUBTYPE
# ============================================================================
cat("\n>>> 2. Immune microenvironment by subtype\n")

# Load ssGSEA immune scores from Phase 4
ssgsea_tc_file <- file.path(P4_DIR, "ssGSEA_transcriptomics.csv")
ssgsea_pr_file <- file.path(P4_DIR, "ssGSEA_proteomics.csv")

immune_tc <- read.csv(ssgsea_tc_file, check.names = FALSE, row.names = 1)
immune_pr <- read.csv(ssgsea_pr_file, check.names = FALSE, row.names = 1)

cat(sprintf("  Immune TC scores: %d cell types x %d samples\n", nrow(immune_tc), ncol(immune_tc)))
cat(sprintf("  Immune PR scores: %d cell types x %d samples\n", nrow(immune_pr), ncol(immune_pr)))

# Select Adjacent samples
adj_immune_tc <- immune_tc[, intersect(adj_samples, colnames(immune_tc)), drop = FALSE]
adj_immune_pr <- immune_pr[, intersect(adj_samples, colnames(immune_pr)), drop = FALSE]

# Test immune cell infiltration between subtypes
test_immune_subtypes <- function(immune_mat, subtype_vec, omics_name) {
  results <- data.frame()
  cs1 <- intersect(names(subtype_vec)[subtype_vec == "CS1"], colnames(immune_mat))
  cs2 <- intersect(names(subtype_vec)[subtype_vec == "CS2"], colnames(immune_mat))
  
  for (cell in rownames(immune_mat)) {
    v1 <- as.numeric(immune_mat[cell, cs1])
    v2 <- as.numeric(immune_mat[cell, cs2])
    tryCatch({
      wt <- wilcox.test(v1, v2)
      results <- rbind(results, data.frame(
        cell_type = cell,
        omics = omics_name,
        mean_CS1 = mean(v1, na.rm = TRUE),
        mean_CS2 = mean(v2, na.rm = TRUE),
        diff = mean(v1, na.rm = TRUE) - mean(v2, na.rm = TRUE),
        pvalue = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(results$pvalue), ]
  return(results)
}

imm_tc_test <- test_immune_subtypes(adj_immune_tc, subtype_vec, "transcriptomics")
imm_pr_test <- test_immune_subtypes(adj_immune_pr, subtype_vec, "proteomics")

cat(sprintf("  TC immune cells sig (P<0.05): %d / %d\n",
            sum(imm_tc_test$pvalue < 0.05), nrow(imm_tc_test)))
cat(sprintf("  PR immune cells sig (P<0.05): %d / %d\n",
            sum(imm_pr_test$pvalue < 0.05), nrow(imm_pr_test)))

# Show top immune cells
cat("\n  Top subtype-differential immune cells (TC):\n")
for (i in 1:min(8, nrow(imm_tc_test))) {
  cat(sprintf("    %s: CS1=%.3f, CS2=%.3f, P=%.4f %s\n",
              imm_tc_test$cell_type[i],
              imm_tc_test$mean_CS1[i], imm_tc_test$mean_CS2[i],
              imm_tc_test$pvalue[i],
              ifelse(imm_tc_test$pvalue[i] < 0.05, "*", "")))
}

imm_all <- rbind(imm_tc_test, imm_pr_test)
write.csv(imm_all, file.path(OUT_DIR, "subtype_immune_scores.csv"), row.names = FALSE)

# ---- Immune heatmap by subtype ----
# Combined TC + PR immune cell heatmap
immune_combined <- rbind(adj_immune_tc, adj_immune_pr)
shared_samps <- intersect(colnames(adj_immune_tc), colnames(adj_immune_pr))

if (length(shared_samps) >= 4) {
  # Use TC immune scores (more comprehensive)
  imm_plot <- as.matrix(adj_immune_tc[, shared_samps])
  
  samp_order <- c(
    sort(intersect(shared_samps, names(subtype_vec)[subtype_vec == "CS1"])),
    sort(intersect(shared_samps, names(subtype_vec)[subtype_vec == "CS2"]))
  )
  imm_plot <- imm_plot[, samp_order]
  
  anno_col <- data.frame(Subtype = subtype_vec[samp_order])
  rownames(anno_col) <- samp_order
  
  pdf(file.path(FIG_DIR, "immune_heatmap_subtype.pdf"), width = 10, height = 8)
  pheatmap(imm_plot,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           annotation_col = anno_col,
           annotation_colors = list(Subtype = c("CS1" = "#3C5488", "CS2" = "#E64B35")),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 8,
           main = "Immune Cell Infiltration by Subtype (TC ssGSEA)",
           scale = "row")
  dev.off()
  cat("  Immune heatmap saved.\n")
}

# ---- Immune boxplots (top cells) ----
top_immune_cells <- head(imm_tc_test$cell_type, 8)
imm_box_data <- data.frame()
for (cell in top_immune_cells) {
  vals <- as.numeric(adj_immune_tc[cell, ])
  samps <- colnames(adj_immune_tc)
  st <- subtype_vec[samps]
  imm_box_data <- rbind(imm_box_data, data.frame(
    cell_type = cell, score = vals, subtype = st, stringsAsFactors = FALSE
  ))
}

p_imm <- ggplot(imm_box_data, aes(x = subtype, y = score, fill = subtype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("CS1" = "#3C5488", "CS2" = "#E64B35")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  labs(title = "Immune Cell Infiltration by Molecular Subtype",
       x = "Subtype", y = "ssGSEA Score") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(size = 7))

ggsave(file.path(FIG_DIR, "immune_boxplots_subtype.pdf"), p_imm, width = 12, height = 6)
cat("  Immune boxplots saved.\n")


# ============================================================================
# 3. METABOLIC LANDSCAPE BY SUBTYPE
# ============================================================================
cat("\n>>> 3. Metabolic landscape by subtype\n")

# Load metabolomics annotation
met_anno_file <- file.path(PROC_DIR, "metabolomics_annotation.csv")
met_anno <- read.csv(met_anno_file, check.names = FALSE)
cat(sprintf("  Metabolite annotations: %d\n", nrow(met_anno)))

# Test metabolites between subtypes
cs1_met <- intersect(names(subtype_vec)[subtype_vec == "CS1"], colnames(met_adj))
cs2_met <- intersect(names(subtype_vec)[subtype_vec == "CS2"], colnames(met_adj))

met_results <- data.frame()
for (m in rownames(met_adj)) {
  v1 <- met_adj[m, cs1_met]
  v2 <- met_adj[m, cs2_met]
  tryCatch({
    wt <- wilcox.test(v1, v2)
    fc <- mean(v1, na.rm = TRUE) - mean(v2, na.rm = TRUE)
    met_results <- rbind(met_results, data.frame(
      compound = m,
      mean_CS1 = mean(v1, na.rm = TRUE),
      mean_CS2 = mean(v2, na.rm = TRUE),
      log2FC_CS1vsCS2 = fc,
      pvalue = wt$p.value,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}
met_results$padj <- p.adjust(met_results$pvalue, method = "BH")
met_results <- met_results[order(met_results$pvalue), ]

cat(sprintf("  Metabolites sig (P<0.05): %d / %d\n",
            sum(met_results$pvalue < 0.05), nrow(met_results)))

# Annotate metabolites
id_col <- intersect(c("Compound", "compound_id", "CompoundID"), colnames(met_anno))
if (length(id_col) > 0) {
  met_results <- merge(met_results, met_anno, by.x = "compound", by.y = id_col[1], all.x = TRUE)
}

write.csv(met_results, file.path(OUT_DIR, "subtype_metabolites.csv"), row.names = FALSE)

# Show top metabolites
cat("\n  Top subtype-differential metabolites:\n")
for (i in 1:min(10, nrow(met_results))) {
  name <- met_results$compound[i]
  # Try to get annotation name
  name_cols <- intersect(c("Metabolites", "Name", "metabolite_name"), colnames(met_results))
  if (length(name_cols) > 0 && !is.na(met_results[[name_cols[1]]][i])) {
    name <- paste0(met_results$compound[i], " (", met_results[[name_cols[1]]][i], ")")
  }
  cat(sprintf("    %s: FC=%.2f, P=%.4f\n", name, met_results$log2FC_CS1vsCS2[i], met_results$pvalue[i]))
}

# ---- Metabolite class enrichment by subtype ----
class_col <- intersect(c("ClassI", "Class", "class"), colnames(met_results))
if (length(class_col) > 0) {
  sig_mets <- met_results[met_results$pvalue < 0.05, ]
  if (nrow(sig_mets) > 0) {
    class_counts <- table(sig_mets[[class_col[1]]])
    class_counts <- sort(class_counts, decreasing = TRUE)
    
    cat("\n  Subtype-differential metabolite classes:\n")
    for (cl in names(class_counts)) {
      cat(sprintf("    %s: %d\n", cl, class_counts[cl]))
    }
    
    # Barplot of metabolite classes
    class_df <- as.data.frame(class_counts)
    colnames(class_df) <- c("Class", "Count")
    class_df <- class_df[class_df$Count > 0, ]
    
    if (nrow(class_df) >= 2) {
      p_class <- ggplot(class_df, aes(x = reorder(Class, Count), y = Count, fill = Class)) +
        geom_bar(stat = "identity", alpha = 0.8) +
        coord_flip() +
        labs(title = "Subtype-Differential Metabolite Classes",
             x = "", y = "Number of metabolites (P<0.05)") +
        theme_bw(base_size = 10) +
        theme(legend.position = "none")
      ggsave(file.path(FIG_DIR, "metabolite_class_barplot.pdf"), p_class, width = 8, height = 5)
      cat("  Metabolite class barplot saved.\n")
    }
  }
}

# ---- Top metabolite volcano-like plot ----
if (nrow(met_results) > 0) {
  met_results$neglog10P <- -log10(met_results$pvalue)
  met_results$sig <- ifelse(met_results$pvalue < 0.05 & abs(met_results$log2FC_CS1vsCS2) > 0.5,
                             "Significant", "NS")
  
  # Add labels for top hits
  top_label <- head(met_results$compound, 15)
  met_results$label <- ifelse(met_results$compound %in% top_label, met_results$compound, "")
  
  # Try to use annotation names
  name_cols <- intersect(c("Metabolites", "Name"), colnames(met_results))
  if (length(name_cols) > 0) {
    idx <- met_results$compound %in% top_label & !is.na(met_results[[name_cols[1]]])
    met_results$label[idx] <- met_results[[name_cols[1]]][idx]
  }
  
  p_vol <- ggplot(met_results, aes(x = log2FC_CS1vsCS2, y = neglog10P, color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
    scale_color_manual(values = c("NS" = "grey60", "Significant" = "#E64B35")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = "Metabolites: CS1 vs CS2",
         x = "log2(FC) CS1 vs CS2", y = "-log10(P-value)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, "metabolite_volcano_subtype.pdf"), p_vol, width = 8, height = 6)
  cat("  Metabolite volcano plot saved.\n")
}


# ============================================================================
# 4. KEY MOLECULAR MARKERS PER SUBTYPE
# ============================================================================
cat("\n>>> 4. Key molecular markers per subtype\n")

# Load subtype-specific differential features from Phase 6
diff_k2 <- read.csv(file.path(P6_DIR, "subtype_diff_features_K2.csv"), stringsAsFactors = FALSE)

# Top markers for each subtype (by absolute FC and significance)
markers_all <- data.frame()
for (s in c("CS1", "CS2")) {
  for (om in c("TC", "PR", "MET")) {
    sub_diff <- diff_k2[diff_k2$subtype == s & diff_k2$omics == om, ]
    sub_diff <- sub_diff[sub_diff$pvalue < 0.05, ]
    sub_diff <- sub_diff[order(-abs(sub_diff$log2FC)), ]
    top_markers <- head(sub_diff, 20)
    markers_all <- rbind(markers_all, top_markers)
  }
}

cat(sprintf("  Total top markers: %d\n", nrow(markers_all)))
write.csv(markers_all, file.path(OUT_DIR, "subtype_top_markers.csv"), row.names = FALSE)

# Show top markers per subtype
for (s in c("CS1", "CS2")) {
  cat(sprintf("\n  %s top markers:\n", s))
  for (om in c("TC", "PR")) {
    sub_m <- markers_all[markers_all$subtype == s & markers_all$omics == om, ]
    if (nrow(sub_m) > 0) {
      top5 <- head(sub_m$feature, 5)
      cat(sprintf("    %s: %s\n", om, paste(top5, collapse = ", ")))
    }
  }
}


# ============================================================================
# 5. COMPREHENSIVE ONCOPLOT-STYLE HEATMAP
# ============================================================================
cat("\n>>> 5. Comprehensive subtype summary heatmap\n")

# Select top features from each category for integrated visualization
# Top 15 genes, 15 proteins, top 5 pathways, top 5 immune cells

# ---- Prepare gene layer ----
tc_markers_cs1 <- diff_k2[diff_k2$subtype == "CS1" & diff_k2$omics == "TC" & diff_k2$pvalue < 0.05, ]
tc_markers_cs1 <- tc_markers_cs1[order(-abs(tc_markers_cs1$log2FC)), ]
tc_markers_cs2 <- diff_k2[diff_k2$subtype == "CS2" & diff_k2$omics == "TC" & diff_k2$pvalue < 0.05, ]
tc_markers_cs2 <- tc_markers_cs2[order(-abs(tc_markers_cs2$log2FC)), ]

top_genes <- unique(c(head(tc_markers_cs1$feature, 8), head(tc_markers_cs2$feature, 7)))
top_genes <- top_genes[top_genes %in% rownames(tc_adj)]

# ---- Prepare protein layer ----
pr_markers_cs1 <- diff_k2[diff_k2$subtype == "CS1" & diff_k2$omics == "PR" & diff_k2$pvalue < 0.05, ]
pr_markers_cs1 <- pr_markers_cs1[order(-abs(pr_markers_cs1$log2FC)), ]
pr_markers_cs2 <- diff_k2[diff_k2$subtype == "CS2" & diff_k2$omics == "PR" & diff_k2$pvalue < 0.05, ]
pr_markers_cs2 <- pr_markers_cs2[order(-abs(pr_markers_cs2$log2FC)), ]

top_proteins <- unique(c(head(pr_markers_cs1$feature, 8), head(pr_markers_cs2$feature, 7)))
top_proteins <- top_proteins[top_proteins %in% rownames(pr_adj)]

# ---- Build comprehensive heatmap ----
samp_order <- c(
  sort(names(subtype_vec)[subtype_vec == "CS1"]),
  sort(names(subtype_vec)[subtype_vec == "CS2"])
)

# Gene expression layer
if (length(top_genes) >= 5) {
  gene_mat <- t(scale(t(tc_adj[top_genes, samp_order])))
  gene_mat[gene_mat > 3] <- 3
  gene_mat[gene_mat < -3] <- -3
}

# Protein expression layer
if (length(top_proteins) >= 5) {
  prot_mat <- t(scale(t(pr_adj[top_proteins, samp_order])))
  prot_mat[prot_mat > 3] <- 3
  prot_mat[prot_mat < -3] <- -3
}

# Pathway scores layer
top_pathways <- head(pw_tc$pathway, 10)
pw_mat <- tc_gsva[top_pathways, samp_order]
pw_mat_scaled <- t(scale(t(pw_mat)))
pw_mat_scaled[pw_mat_scaled > 3] <- 3
pw_mat_scaled[pw_mat_scaled < -3] <- -3
rownames(pw_mat_scaled) <- gsub("HALLMARK_", "", rownames(pw_mat_scaled))
rownames(pw_mat_scaled) <- gsub("_", " ", rownames(pw_mat_scaled))

# Immune cell scores layer
top_cells <- head(imm_tc_test$cell_type, 8)
imm_mat <- as.matrix(adj_immune_tc[top_cells, samp_order])
imm_mat_scaled <- t(scale(t(imm_mat)))
imm_mat_scaled[imm_mat_scaled > 3] <- 3
imm_mat_scaled[imm_mat_scaled < -3] <- -3

# ---- Draw with ComplexHeatmap ----
col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

# Top annotation
ha_top <- HeatmapAnnotation(
  Subtype = subtype_vec[samp_order],
  col = list(Subtype = c("CS1" = "#3C5488", "CS2" = "#E64B35")),
  annotation_name_side = "left",
  show_legend = TRUE
)

# Build heatmap list
ht_list <- NULL

if (length(top_genes) >= 5) {
  ht_list <- Heatmap(gene_mat,
                      name = "z-score",
                      col = col_fun,
                      top_annotation = ha_top,
                      row_title = "Transcriptomics",
                      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      cluster_columns = FALSE,
                      cluster_rows = TRUE,
                      show_column_names = TRUE,
                      row_names_gp = gpar(fontsize = 7),
                      column_names_gp = gpar(fontsize = 8),
                      height = unit(length(top_genes) * 4, "mm"))
}

if (length(top_proteins) >= 5) {
  ht_prot <- Heatmap(prot_mat,
                      name = "z-score (PR)",
                      col = col_fun,
                      row_title = "Proteomics",
                      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      cluster_columns = FALSE,
                      cluster_rows = TRUE,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = 7),
                      height = unit(length(top_proteins) * 4, "mm"),
                      show_heatmap_legend = FALSE)
  if (is.null(ht_list)) ht_list <- ht_prot else ht_list <- ht_list %v% ht_prot
}

ht_pw <- Heatmap(pw_mat_scaled,
                  name = "z-score (PW)",
                  col = col_fun,
                  row_title = "Hallmark Pathways",
                  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                  cluster_columns = FALSE,
                  cluster_rows = TRUE,
                  show_column_names = FALSE,
                  row_names_gp = gpar(fontsize = 7),
                  height = unit(nrow(pw_mat_scaled) * 4, "mm"),
                  show_heatmap_legend = FALSE)

ht_imm <- Heatmap(imm_mat_scaled,
                   name = "z-score (IMM)",
                   col = col_fun,
                   row_title = "Immune Cells",
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   cluster_columns = FALSE,
                   cluster_rows = TRUE,
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 7),
                   height = unit(nrow(imm_mat_scaled) * 4, "mm"),
                   show_heatmap_legend = FALSE)

if (!is.null(ht_list)) {
  ht_list <- ht_list %v% ht_pw %v% ht_imm
} else {
  ht_list <- ht_pw %v% ht_imm
}

# Calculate total height
total_rows <- length(top_genes) + length(top_proteins) + nrow(pw_mat_scaled) + nrow(imm_mat_scaled)
fig_height <- max(10, total_rows * 0.3 + 3)

pdf(file.path(FIG_DIR, "comprehensive_subtype_heatmap.pdf"), width = 10, height = fig_height)
draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right",
     column_title = "Multi-dimensional Characterization of HAE Molecular Subtypes")
dev.off()
cat("  Comprehensive heatmap saved.\n")


# ============================================================================
# 6. SUBTYPE-SPECIFIC BIOLOGICAL THEMES
# ============================================================================
cat("\n>>> 6. Subtype-specific biological themes\n")

# Summarize CS1 vs CS2 pathway differences
cs1_up_pathways <- pw_tc[pw_tc$pvalue < 0.1 & pw_tc$diff_CS1_CS2 > 0, "pathway"]
cs2_up_pathways <- pw_tc[pw_tc$pvalue < 0.1 & pw_tc$diff_CS1_CS2 < 0, "pathway"]

cat(sprintf("  CS1-enriched pathways (P<0.1): %d\n", length(cs1_up_pathways)))
if (length(cs1_up_pathways) > 0) {
  cat(sprintf("    %s\n", paste(gsub("HALLMARK_", "", cs1_up_pathways), collapse = ", ")))
}

cat(sprintf("  CS2-enriched pathways (P<0.1): %d\n", length(cs2_up_pathways)))
if (length(cs2_up_pathways) > 0) {
  cat(sprintf("    %s\n", paste(gsub("HALLMARK_", "", cs2_up_pathways), collapse = ", ")))
}

# Summarize CS1 vs CS2 immune differences
cs1_up_immune <- imm_tc_test[imm_tc_test$pvalue < 0.1 & imm_tc_test$diff > 0, "cell_type"]
cs2_up_immune <- imm_tc_test[imm_tc_test$pvalue < 0.1 & imm_tc_test$diff < 0, "cell_type"]

cat(sprintf("\n  CS1-enriched immune cells (P<0.1): %d\n", length(cs1_up_immune)))
if (length(cs1_up_immune) > 0) cat(sprintf("    %s\n", paste(cs1_up_immune, collapse = ", ")))
cat(sprintf("  CS2-enriched immune cells (P<0.1): %d\n", length(cs2_up_immune)))
if (length(cs2_up_immune) > 0) cat(sprintf("    %s\n", paste(cs2_up_immune, collapse = ", ")))

# Create theme summary
theme_summary <- data.frame(
  subtype = character(),
  dimension = character(),
  theme = character(),
  direction = character(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(pw_tc)) {
  if (pw_tc$pvalue[i] < 0.1) {
    theme_summary <- rbind(theme_summary, data.frame(
      subtype = ifelse(pw_tc$diff_CS1_CS2[i] > 0, "CS1", "CS2"),
      dimension = "Pathway",
      theme = gsub("HALLMARK_", "", pw_tc$pathway[i]),
      direction = ifelse(pw_tc$diff_CS1_CS2[i] > 0, "Up in CS1", "Up in CS2"),
      pvalue = pw_tc$pvalue[i]
    ))
  }
}

for (i in 1:nrow(imm_tc_test)) {
  if (imm_tc_test$pvalue[i] < 0.1) {
    theme_summary <- rbind(theme_summary, data.frame(
      subtype = ifelse(imm_tc_test$diff[i] > 0, "CS1", "CS2"),
      dimension = "Immune",
      theme = imm_tc_test$cell_type[i],
      direction = ifelse(imm_tc_test$diff[i] > 0, "Up in CS1", "Up in CS2"),
      pvalue = imm_tc_test$pvalue[i]
    ))
  }
}

write.csv(theme_summary, file.path(OUT_DIR, "subtype_biological_themes.csv"), row.names = FALSE)
cat(sprintf("\n  Total biological themes identified: %d\n", nrow(theme_summary)))


# ============================================================================
# 7. DRUGGABLE TARGET ANALYSIS
# ============================================================================
cat("\n>>> 7. Druggable target analysis\n")

# Known druggable gene families
druggable_families <- list(
  kinases = c("EGFR", "ERBB2", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
              "PIK3CA", "PIK3CB", "AKT1", "AKT2", "MTOR", "CDK4", "CDK6", "JAK1", "JAK2",
              "SRC", "ABL1", "FLT3", "KIT", "PDGFRA", "PDGFRB", "FGFR1", "FGFR2",
              "MET", "ALK", "ROS1", "RET", "NTRK1", "VEGFR2", "BTK", "SYK"),
  immune_checkpoints = c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
                          "BTLA", "VSIR", "CD47", "SIRPA", "IDO1", "IDO2",
                          "TNFRSF4", "TNFRSF9", "TNFRSF18", "CD28", "ICOS"),
  metabolic_enzymes = c("FASN", "ACACA", "ACLY", "HMGCR", "SCD", "LDHA", "LDHB",
                         "PKM", "HK1", "HK2", "GLS", "IDH1", "IDH2",
                         "NAMPT", "NNMT", "ARG1", "ARG2"),
  epigenetic_regulators = c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2",
                             "EZH2", "KDM1A", "HDAC1", "HDAC2", "HDAC3",
                             "BRD4", "KAT2A", "KAT2B"),
  ECM_targets = c("MMP2", "MMP9", "MMP14", "ADAM17", "LOX", "LOXL2",
                   "COL1A1", "COL1A2", "COL3A1", "FN1", "TNC", "POSTN")
)

druggable_results <- data.frame()
for (family in names(druggable_families)) {
  genes <- druggable_families[[family]]
  
  # Check in transcriptomics
  tc_present <- intersect(genes, rownames(tc_adj))
  for (g in tc_present) {
    cs1_vals <- tc_adj[g, names(subtype_vec)[subtype_vec == "CS1"]]
    cs2_vals <- tc_adj[g, names(subtype_vec)[subtype_vec == "CS2"]]
    tryCatch({
      wt <- wilcox.test(cs1_vals, cs2_vals)
      druggable_results <- rbind(druggable_results, data.frame(
        gene = g, family = family, omics = "TC",
        mean_CS1 = mean(cs1_vals), mean_CS2 = mean(cs2_vals),
        log2FC = mean(cs1_vals) - mean(cs2_vals),
        pvalue = wt$p.value, stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  
  # Check in proteomics
  pr_present <- intersect(genes, rownames(pr_adj))
  for (g in pr_present) {
    cs1_vals <- pr_adj[g, names(subtype_vec)[subtype_vec == "CS1"]]
    cs2_vals <- pr_adj[g, names(subtype_vec)[subtype_vec == "CS2"]]
    tryCatch({
      wt <- wilcox.test(cs1_vals, cs2_vals)
      druggable_results <- rbind(druggable_results, data.frame(
        gene = g, family = family, omics = "PR",
        mean_CS1 = mean(cs1_vals), mean_CS2 = mean(cs2_vals),
        log2FC = mean(cs1_vals) - mean(cs2_vals),
        pvalue = wt$p.value, stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
}

druggable_results$padj <- p.adjust(druggable_results$pvalue, method = "BH")
druggable_results <- druggable_results[order(druggable_results$pvalue), ]
write.csv(druggable_results, file.path(OUT_DIR, "druggable_targets_by_subtype.csv"), row.names = FALSE)

sig_drugs <- druggable_results[druggable_results$pvalue < 0.05, ]
cat(sprintf("  Druggable targets tested: %d, significant (P<0.05): %d\n",
            nrow(druggable_results), nrow(sig_drugs)))

if (nrow(sig_drugs) > 0) {
  cat("  Significant druggable targets:\n")
  for (i in 1:min(15, nrow(sig_drugs))) {
    cat(sprintf("    %s (%s, %s): CS1=%.2f, CS2=%.2f, FC=%.2f, P=%.4f\n",
                sig_drugs$gene[i], sig_drugs$family[i], sig_drugs$omics[i],
                sig_drugs$mean_CS1[i], sig_drugs$mean_CS2[i],
                sig_drugs$log2FC[i], sig_drugs$pvalue[i]))
  }
}

# ---- Druggable target heatmap ----
if (nrow(sig_drugs) > 0) {
  # Select unique genes
  drug_genes <- unique(sig_drugs$gene)
  drug_tc <- intersect(drug_genes, rownames(tc_adj))
  
  if (length(drug_tc) >= 3) {
    drug_mat <- t(scale(t(tc_adj[drug_tc, samp_order, drop = FALSE])))
    drug_mat[drug_mat > 3] <- 3
    drug_mat[drug_mat < -3] <- -3
    
    # Family annotation
    gene_family <- sig_drugs$family[match(drug_tc, sig_drugs$gene)]
    anno_row <- data.frame(Family = gene_family)
    rownames(anno_row) <- drug_tc
    
    family_colors <- c(
      "kinases" = "#E64B35", "immune_checkpoints" = "#4DBBD5",
      "metabolic_enzymes" = "#00A087", "epigenetic_regulators" = "#F39B7F",
      "ECM_targets" = "#8491B4"
    )
    
    pdf(file.path(FIG_DIR, "druggable_targets_heatmap.pdf"), width = 10, height = max(5, length(drug_tc) * 0.4 + 2))
    pheatmap(drug_mat,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             annotation_col = data.frame(Subtype = subtype_vec[samp_order], row.names = samp_order),
             annotation_row = anno_row,
             annotation_colors = list(
               Subtype = c("CS1" = "#3C5488", "CS2" = "#E64B35"),
               Family = family_colors[unique(gene_family)]
             ),
             cluster_cols = FALSE,
             cluster_rows = TRUE,
             show_colnames = TRUE,
             fontsize_row = 8,
             fontsize_col = 8,
             main = "Druggable Targets by Molecular Subtype")
    dev.off()
    cat("  Druggable targets heatmap saved.\n")
  }
}


# ============================================================================
# 8. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 7 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Subtypes: CS1 (n=%d), CS2 (n=%d)\n",
            sum(subtype_vec == "CS1"), sum(subtype_vec == "CS2")))

cat(sprintf("\n1. Pathway profiling:\n"))
cat(sprintf("   TC Hallmark: %d sig (P<0.05)\n", sum(pw_tc$pvalue < 0.05)))
cat(sprintf("   PR Hallmark: %d sig (P<0.05)\n", sum(pw_pr$pvalue < 0.05)))

cat(sprintf("\n2. Immune microenvironment:\n"))
cat(sprintf("   TC immune cells sig: %d / %d\n", sum(imm_tc_test$pvalue < 0.05), nrow(imm_tc_test)))
cat(sprintf("   PR immune cells sig: %d / %d\n", sum(imm_pr_test$pvalue < 0.05), nrow(imm_pr_test)))

cat(sprintf("\n3. Metabolic landscape:\n"))
cat(sprintf("   Metabolites sig: %d / %d\n", sum(met_results$pvalue < 0.05), nrow(met_results)))

cat(sprintf("\n4. Druggable targets:\n"))
cat(sprintf("   Tested: %d, sig: %d\n", nrow(druggable_results), nrow(sig_drugs)))

cat(sprintf("\n5. Biological themes (P<0.1): %d\n", nrow(theme_summary)))

cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 7 COMPLETE.\n")
