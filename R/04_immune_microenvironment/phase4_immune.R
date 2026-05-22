#!/usr/bin/env Rscript
# ============================================================================
# Phase 4: Immune Microenvironment Analysis
# HAE Multi-omics Integration Study
# ============================================================================
# Strategy:
#   1. ssGSEA-based immune cell deconvolution (18 immune/stromal signatures)
#   2. MCP-counter deconvolution (proteomics compatible)
#   3. Immune checkpoint & exhaustion marker profiling
#   4. Immune pathway scoring across omics
#   5. Immune-metabolic crosstalk analysis
# ============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(limma)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
OUT_DIR  <- file.path(PROJECT, "analysis/results/phase4_immune")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 4: Immune Microenvironment\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA & CONVERT TO GENE SYMBOLS
# ============================================================================
cat(">>> 0. Loading data & converting IDs to gene symbols\n")

# Differential results (contain ID -> gene_name mappings)
degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# --- Transcriptomics: convert Ensembl ID -> gene symbol ---
tc_vst_raw <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                       check.names = FALSE, row.names = 1)
# Build mapping from DEGs (gene_id -> gene_name)
id2sym_tc <- setNames(degs$gene_name, degs$gene_id)
id2sym_tc <- id2sym_tc[!is.na(id2sym_tc) & id2sym_tc != ""]

# Map row names
new_names <- id2sym_tc[rownames(tc_vst_raw)]
keep <- !is.na(new_names) & new_names != ""
tc_vst <- tc_vst_raw[keep, ]
rownames(tc_vst) <- make.unique(new_names[keep])
cat(sprintf("  Transcriptomics: %d -> %d genes (symbol-mapped) x %d samples\n",
            nrow(tc_vst_raw), nrow(tc_vst), ncol(tc_vst)))

# --- Proteomics: convert Ensembl Protein ID -> gene symbol ---
pr_log2_raw <- read.csv(file.path(PROC_DIR, "proteomics_log2_norm.csv"),
                        check.names = FALSE, row.names = 1)
# Build mapping from DEPs (Protein -> gene_name)
id2sym_pr <- setNames(deps$gene_name, deps$Protein)
id2sym_pr <- id2sym_pr[!is.na(id2sym_pr) & id2sym_pr != "" & id2sym_pr != "_"]

new_names_pr <- id2sym_pr[rownames(pr_log2_raw)]
keep_pr <- !is.na(new_names_pr) & new_names_pr != ""
pr_log2 <- pr_log2_raw[keep_pr, ]
rownames(pr_log2) <- make.unique(new_names_pr[keep_pr])
cat(sprintf("  Proteomics: %d -> %d proteins (symbol-mapped) x %d samples\n",
            nrow(pr_log2_raw), nrow(pr_log2), ncol(pr_log2)))

# Sample correspondence
sample_corr <- read.csv(file.path(PROC_DIR, "sample_correspondence.csv"), check.names = FALSE)


# ============================================================================
# 1. IMMUNE CELL SIGNATURE GENE SETS (Charoentong et al. 2017)
# ============================================================================
cat("\n>>> 1. Preparing immune cell signatures\n")

# Comprehensive immune cell signatures for ssGSEA
# Based on Charoentong et al. (2017) and Bindea et al. (2013)
immune_signatures <- list(
  # Innate immunity
  Macrophages_M1 = c("NOS2", "TNF", "IL1B", "IL6", "IL12A", "IL12B", "CD80", "CD86",
                      "IRF5", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IDO1"),
  Macrophages_M2 = c("CD163", "MRC1", "MSR1", "CD68", "IL10", "TGFB1", "CCL17",
                      "CCL22", "ARG1", "VEGFA", "MMP9", "FN1"),
  NK_cells = c("NCR1", "KLRD1", "KLRC1", "KLRB1", "CD160", "KIR2DL1", "KIR2DL3",
                "KIR3DL1", "NKG7", "GNLY", "PRF1", "GZMB"),
  Dendritic_cells = c("CD1A", "CD1C", "CD1E", "FCER1A", "CLEC10A", "CD209",
                       "HLA-DQA1", "HLA-DQB1", "ITGAX", "BATF3", "IRF8"),
  Neutrophils = c("CEACAM8", "FPR1", "SIGLEC5", "CSF3R", "FCAR", "FCGR3B",
                   "CXCR1", "CXCR2", "S100A8", "S100A9", "S100A12"),
  Mast_cells = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "FCER1A", "HDC",
                  "KIT", "ENPP3", "GATA2"),
  Eosinophils = c("CCR3", "SIGLEC8", "PRG2", "EPX", "RNASE2", "RNASE3",
                   "IL5RA", "PTGDR2"),

  # Adaptive immunity - T cells
  CD8_T_cells = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "TBX21",
                   "EOMES", "NKG7", "KLRK1"),
  CD4_T_cells = c("CD4", "IL7R", "TCF7", "LEF1", "SELL", "CCR7", "CD27",
                   "CD28", "ICOS", "BCL6"),
  Th1_cells = c("TBX21", "IFNG", "TNF", "IL2", "STAT4", "IL12RB1", "IL12RB2",
                 "CCR5", "CXCR3"),
  Th2_cells = c("GATA3", "IL4", "IL5", "IL13", "IL4R", "CCR4", "CCR3",
                 "STAT6", "IL10"),
  Th17_cells = c("RORC", "IL17A", "IL17F", "IL22", "IL23R", "CCR6", "IL6R",
                  "STAT3", "BATF"),
  Treg_cells = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF18", "IKZF2",
                  "LRRC32", "IL10", "TGFB1"),
  Gamma_delta_T = c("TRGC1", "TRGC2", "TRDC", "TRDV1", "TRDV2", "TRDV3"),

  # Adaptive immunity - B cells
  B_cells = c("CD19", "CD79A", "CD79B", "MS4A1", "CD22", "BLNK", "PAX5",
               "BLK", "CR2", "TNFRSF13C"),
  Plasma_cells = c("SDC1", "XBP1", "PRDM1", "IRF4", "IGHG1", "IGHG2",
                    "IGHA1", "IGKC", "JCHAIN", "MZB1"),

  # Stromal
  Fibroblasts = c("FAP", "PDPN", "COL1A1", "COL1A2", "COL3A1", "ACTA2",
                   "FN1", "VIM", "TAGLN", "MYH11"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "KDR", "FLT1", "ENG", "PLVAP",
                   "ESAM", "CLDN5", "ERG")
)

cat(sprintf("  Defined %d immune cell signatures\n", length(immune_signatures)))


# ============================================================================
# 2. ssGSEA IMMUNE DECONVOLUTION — TRANSCRIPTOMICS
# ============================================================================
cat("\n>>> 2. ssGSEA Immune Deconvolution (Transcriptomics)\n")

# Convert gene names to match signature genes
tc_mat <- as.matrix(tc_vst)

# Run ssGSEA via GSVA
cat("  Running ssGSEA on transcriptomics...\n")
param_tc <- ssgseaParam(tc_mat, immune_signatures, normalize = TRUE)
ssgsea_tc <- gsva(param_tc, verbose = FALSE)
cat(sprintf("  ssGSEA matrix: %d cell types x %d samples\n",
            nrow(ssgsea_tc), ncol(ssgsea_tc)))

write.csv(ssgsea_tc, file.path(OUT_DIR, "ssGSEA_transcriptomics.csv"))


# ============================================================================
# 3. ssGSEA IMMUNE DECONVOLUTION — PROTEOMICS
# ============================================================================
cat("\n>>> 3. ssGSEA Immune Deconvolution (Proteomics)\n")

pr_mat <- as.matrix(pr_log2)
cat("  Running ssGSEA on proteomics...\n")
param_pr <- ssgseaParam(pr_mat, immune_signatures, normalize = TRUE)
ssgsea_pr <- gsva(param_pr, verbose = FALSE)
cat(sprintf("  ssGSEA matrix: %d cell types x %d samples\n",
            nrow(ssgsea_pr), ncol(ssgsea_pr)))

write.csv(ssgsea_pr, file.path(OUT_DIR, "ssGSEA_proteomics.csv"))


# ============================================================================
# 4. DIFFERENTIAL IMMUNE CELL ABUNDANCE (Adjacent vs Normal)
# ============================================================================
cat("\n>>> 4. Differential Immune Cell Abundance\n")

test_immune_diff <- function(ssgsea_mat, sample_info, label) {
  results <- data.frame()
  
  # Extract patient IDs for paired testing
  extract_id <- function(x) as.integer(gsub("Normal|Adjacent", "", x))
  
  for (ct in rownames(ssgsea_mat)) {
    scores <- ssgsea_mat[ct, ]
    
    # Match sample groups
    normal_samples <- sample_info$sample_id[sample_info$group == "Normal"]
    adjacent_samples <- sample_info$sample_id[sample_info$group == "Adjacent"]
    
    normal_scores <- scores[names(scores) %in% normal_samples]
    adjacent_scores <- scores[names(scores) %in% adjacent_samples]
    
    if (length(normal_scores) < 3 || length(adjacent_scores) < 3) next
    
    # Match pairs by patient ID for paired Wilcoxon test
    nor_ids <- extract_id(names(normal_scores))
    adj_ids <- extract_id(names(adjacent_scores))
    paired_ids <- intersect(nor_ids, adj_ids)
    
    if (length(paired_ids) >= 3) {
      # Use paired Wilcoxon (correct for matched patient samples)
      nor_paired <- normal_scores[paste0("Normal", paired_ids)]
      adj_paired <- adjacent_scores[paste0("Adjacent", paired_ids)]
      wt <- wilcox.test(adj_paired, nor_paired, paired = TRUE)
    } else {
      # Fallback to unpaired if pairing not possible
      wt <- wilcox.test(adjacent_scores, normal_scores, paired = FALSE)
    }
    
    results <- rbind(results, data.frame(
      cell_type = ct,
      mean_Normal = mean(normal_scores),
      mean_Adjacent = mean(adjacent_scores),
      diff = mean(adjacent_scores) - mean(normal_scores),
      pvalue = wt$p.value,
      stringsAsFactors = FALSE
    ))
  }
  
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results$direction <- ifelse(results$diff > 0, "Enriched in Adjacent", "Enriched in Normal")
  results <- results[order(results$pvalue), ]
  
  cat(sprintf("  %s: %d cell types, %d significant (P<0.05), %d FDR<0.1\n",
              label, nrow(results), sum(results$pvalue < 0.05),
              sum(results$padj < 0.1, na.rm = TRUE)))
  return(results)
}

# Build sample info from column names
# Columns are: Normal1..Normal14, Adjacent1..Adjacent14
tc_samples <- colnames(tc_vst)
tc_si <- data.frame(
  sample_id = tc_samples,
  group = ifelse(grepl("^Normal", tc_samples), "Normal", "Adjacent"),
  stringsAsFactors = FALSE
)

pr_samples <- colnames(pr_log2)
pr_si <- data.frame(
  sample_id = pr_samples,
  group = ifelse(grepl("^Normal", pr_samples), "Normal", "Adjacent"),
  stringsAsFactors = FALSE
)

immune_diff_tc <- test_immune_diff(ssgsea_tc, tc_si, "Transcriptomics")
immune_diff_pr <- test_immune_diff(ssgsea_pr, pr_si, "Proteomics")

write.csv(immune_diff_tc, file.path(OUT_DIR, "immune_diff_transcriptomics.csv"), row.names = FALSE)
write.csv(immune_diff_pr, file.path(OUT_DIR, "immune_diff_proteomics.csv"), row.names = FALSE)


# ============================================================================
# 5. VISUALIZATION: Immune Heatmap & Boxplots
# ============================================================================
cat("\n>>> 5. Immune Visualization\n")

# 5a. Heatmap of ssGSEA scores (transcriptomics)
col_ann_tc <- data.frame(
  Group = tc_si$group[match(colnames(ssgsea_tc), tc_si$sample_id)],
  row.names = colnames(ssgsea_tc)
)
ann_colors <- list(Group = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35"))

ssgsea_scaled <- t(scale(t(ssgsea_tc)))

pdf(file.path(FIG_DIR, "heatmap_ssGSEA_transcriptomics.pdf"),
    width = 10, height = 8)
pheatmap(
  ssgsea_scaled,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  annotation_col = col_ann_tc,
  annotation_colors = ann_colors,
  cluster_rows = TRUE, cluster_cols = TRUE,
  fontsize_row = 9,
  show_colnames = FALSE,
  main = "Immune Cell ssGSEA Scores (Transcriptomics)"
)
dev.off()
cat("  Transcriptomics ssGSEA heatmap saved.\n")

# 5b. Heatmap of ssGSEA scores (proteomics)
col_ann_pr <- data.frame(
  Group = pr_si$group[match(colnames(ssgsea_pr), pr_si$sample_id)],
  row.names = colnames(ssgsea_pr)
)

ssgsea_pr_scaled <- t(scale(t(ssgsea_pr)))

pdf(file.path(FIG_DIR, "heatmap_ssGSEA_proteomics.pdf"),
    width = 10, height = 8)
pheatmap(
  ssgsea_pr_scaled,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  annotation_col = col_ann_pr,
  annotation_colors = ann_colors,
  cluster_rows = TRUE, cluster_cols = TRUE,
  fontsize_row = 9,
  show_colnames = FALSE,
  main = "Immune Cell ssGSEA Scores (Proteomics)"
)
dev.off()
cat("  Proteomics ssGSEA heatmap saved.\n")

# 5c. Dotplot: cross-omics immune cell comparison
merged_diff <- merge(
  immune_diff_tc[, c("cell_type", "diff", "pvalue")],
  immune_diff_pr[, c("cell_type", "diff", "pvalue")],
  by = "cell_type", suffixes = c("_TC", "_PR")
)

merged_diff$cell_label <- gsub("_", " ", merged_diff$cell_type)
merged_diff$sig_TC <- merged_diff$pvalue_TC < 0.05
merged_diff$sig_PR <- merged_diff$pvalue_PR < 0.05
merged_diff$concordant <- sign(merged_diff$diff_TC) == sign(merged_diff$diff_PR)

p_immune <- ggplot(merged_diff, aes(diff_TC, diff_PR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey70") +
  geom_point(aes(color = concordant, shape = sig_TC | sig_PR), size = 3.5, alpha = 0.8) +
  geom_text_repel(aes(label = cell_label), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("TRUE" = "#00A087", "FALSE" = "#E64B35")) +
  labs(title = "Cross-omics Immune Cell Changes (Adjacent vs Normal)",
       x = "ssGSEA Difference (Transcriptomics)",
       y = "ssGSEA Difference (Proteomics)",
       color = "Concordant", shape = "P<0.05 in\neither omics") +
  theme_bw(base_size = 12) +
  coord_equal()
ggsave(file.path(FIG_DIR, "immune_cross_omics_scatter.pdf"), p_immune,
       width = 9, height = 8)
cat("  Cross-omics immune scatter saved.\n")

# 5d. Barplot: immune cell change summary
bar_data <- merged_diff %>%
  select(cell_label, diff_TC, diff_PR) %>%
  pivot_longer(cols = c(diff_TC, diff_PR), names_to = "omics", values_to = "diff") %>%
  mutate(omics = ifelse(omics == "diff_TC", "Transcriptomics", "Proteomics"))

p_bar <- ggplot(bar_data, aes(diff, reorder(cell_label, diff), fill = omics)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Transcriptomics" = "#3C5488", "Proteomics" = "#E64B35")) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30") +
  labs(title = "Immune Cell Enrichment (Adjacent vs Normal)",
       x = "ssGSEA Score Difference", y = NULL, fill = "Omics") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))
ggsave(file.path(FIG_DIR, "immune_cell_barplot.pdf"), p_bar, width = 9, height = 8)
cat("  Immune cell barplot saved.\n")


# ============================================================================
# 6. IMMUNE CHECKPOINT & EXHAUSTION MARKERS
# ============================================================================
cat("\n>>> 6. Immune Checkpoint & Exhaustion Markers\n")

checkpoint_genes <- c(
  # Inhibitory checkpoints
  "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
  "VSIR", "BTLA", "ADORA2A",
  # Co-stimulatory
  "CD28", "CD80", "CD86", "ICOS", "TNFRSF9", "TNFRSF4", "CD40LG", "CD40",
  # Exhaustion markers
  "TOX", "ENTPD1", "LAYN", "CXCL13",
  # HLA
  "HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-DRA", "HLA-DRB1"
)

extract_marker_expr <- function(de_results, genes, label) {
  matched <- de_results[de_results$gene_name %in% genes, ]
  matched$omics <- label
  matched$sig <- matched$P.Value < 0.05 & abs(matched$logFC) > 0.585
  return(matched[, c("gene_name", "logFC", "t", "P.Value", "adj.P.Val", "omics", "sig")])
}

ckpt_tc <- extract_marker_expr(degs, checkpoint_genes, "Transcriptomics")
ckpt_pr <- extract_marker_expr(deps, checkpoint_genes, "Proteomics")
ckpt_all <- rbind(ckpt_tc, ckpt_pr)

write.csv(ckpt_all, file.path(OUT_DIR, "checkpoint_marker_expression.csv"), row.names = FALSE)

cat(sprintf("  Checkpoint markers detected: TC=%d, PR=%d\n", nrow(ckpt_tc), nrow(ckpt_pr)))
cat(sprintf("  Significant (P<0.05, |FC|>1.5): TC=%d, PR=%d\n",
            sum(ckpt_tc$sig), sum(ckpt_pr$sig)))

# Heatmap: checkpoint markers across omics
if (nrow(ckpt_tc) > 0 && nrow(ckpt_pr) > 0) {
  tc_fc <- setNames(ckpt_tc$logFC, ckpt_tc$gene_name)
  pr_fc <- setNames(ckpt_pr$logFC, ckpt_pr$gene_name)
  shared <- intersect(names(tc_fc), names(pr_fc))
  
  if (length(shared) > 1) {
    mat <- cbind(mRNA = tc_fc[shared], Protein = pr_fc[shared])
    rownames(mat) <- shared
    
    # Significance annotation
    tc_sig_map <- setNames(ckpt_tc$sig, ckpt_tc$gene_name)
    pr_sig_map <- setNames(ckpt_pr$sig, ckpt_pr$gene_name)
    
    ann_row <- data.frame(
      Category = ifelse(shared %in% c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3",
                                        "HAVCR2","TIGIT","VSIR","BTLA","ADORA2A"),
                          "Inhibitory",
                   ifelse(shared %in% c("HLA-A","HLA-B","HLA-C","B2M","HLA-DRA","HLA-DRB1"),
                          "HLA",
                   ifelse(shared %in% c("TOX","ENTPD1","LAYN","CXCL13"),
                          "Exhaustion", "Co-stimulatory"))),
      row.names = shared
    )
    
    pdf(file.path(FIG_DIR, "heatmap_checkpoint_markers.pdf"),
        width = 7, height = max(5, length(shared) * 0.35 + 2))
    pheatmap(
      mat,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      cluster_cols = FALSE,
      annotation_row = ann_row,
      fontsize_row = 9,
      main = "Immune Checkpoint & Exhaustion Markers (logFC)",
      breaks = seq(-max(abs(mat), na.rm=TRUE), max(abs(mat), na.rm=TRUE), length.out = 101)
    )
    dev.off()
    cat("  Checkpoint marker heatmap saved.\n")
  }
}


# ============================================================================
# 7. IMMUNE PATHWAY SCORING (Hallmark immune pathways)
# ============================================================================
cat("\n>>> 7. Immune Pathway Scoring\n")

# Read Phase 2 GSEA results for immune-related pathways
gsea_tc_hm <- read.csv(file.path(PROJECT, "analysis/results/phase2_enrichment/GSEA_transcriptomics_Hallmark.csv"),
                       check.names = FALSE)
gsea_pr_hm <- read.csv(file.path(PROJECT, "analysis/results/phase2_enrichment/GSEA_proteomics_Hallmark.csv"),
                       check.names = FALSE)

# Immune-related Hallmark pathways
immune_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_COAGULATION",
  "HALLMARK_KRAS_SIGNALING_UP"
)

immune_gsea_tc <- gsea_tc_hm[gsea_tc_hm$pathway %in% immune_pathways, ]
immune_gsea_pr <- gsea_pr_hm[gsea_pr_hm$pathway %in% immune_pathways, ]

if (nrow(immune_gsea_tc) > 0 && nrow(immune_gsea_pr) > 0) {
  immune_merged <- merge(
    immune_gsea_tc[, c("pathway", "NES", "padj")],
    immune_gsea_pr[, c("pathway", "NES", "padj")],
    by = "pathway", suffixes = c("_TC", "_PR")
  )
  immune_merged$pathway_label <- gsub("^HALLMARK_", "", immune_merged$pathway)
  immune_merged$pathway_label <- gsub("_", " ", immune_merged$pathway_label)
  
  write.csv(immune_merged, file.path(OUT_DIR, "immune_hallmark_NES.csv"), row.names = FALSE)
  
  # Paired barplot of immune pathway NES
  im_long <- immune_merged %>%
    select(pathway_label, NES_TC, NES_PR) %>%
    pivot_longer(cols = c(NES_TC, NES_PR), names_to = "omics", values_to = "NES") %>%
    mutate(omics = ifelse(omics == "NES_TC", "Transcriptomics", "Proteomics"))
  
  p_immune_pw <- ggplot(im_long, aes(NES, reorder(pathway_label, NES), fill = omics)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
    scale_fill_manual(values = c("Transcriptomics" = "#3C5488", "Proteomics" = "#E64B35")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = "Immune-Related Hallmark Pathway NES",
         x = "Normalized Enrichment Score", y = NULL, fill = "Omics") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file.path(FIG_DIR, "immune_hallmark_NES_barplot.pdf"), p_immune_pw,
         width = 9, height = 6)
  cat("  Immune Hallmark NES barplot saved.\n")
}


# ============================================================================
# 8. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 4 SUMMARY\n")
cat("========================================\n")

cat(sprintf("ssGSEA deconvolution: %d cell types analyzed\n", length(immune_signatures)))
cat("\nDifferential immune cells (P<0.05):\n")
cat(sprintf("  Transcriptomics: %d / %d\n",
            sum(immune_diff_tc$pvalue < 0.05), nrow(immune_diff_tc)))
cat(sprintf("  Proteomics: %d / %d\n",
            sum(immune_diff_pr$pvalue < 0.05), nrow(immune_diff_pr)))

cat(sprintf("\nCheckpoint markers detected across both omics: %d\n",
            length(intersect(ckpt_tc$gene_name, ckpt_pr$gene_name))))

cat(sprintf("\nCross-omics concordant immune cells: %d / %d\n",
            sum(merged_diff$concordant), nrow(merged_diff)))

cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 4 COMPLETE.\n")
