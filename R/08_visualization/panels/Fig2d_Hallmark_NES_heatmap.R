#!/usr/bin/env Rscript
# =============================================================================
# Fig2d_Hallmark_NES_heatmap.R -- Cross-omics Hallmark NES Heatmap
# Figure 2d | 90 x 85 mm | ComplexHeatmap
# =============================================================================

slot <- get_panel_slot("Fig2d_Hallmark_NES_heatmap.pdf")
W <- slot$w; H <- slot$h

hall <- read.csv(file.path(RES, "phase2_enrichment/cross_omics_convergent_Hallmark.csv"),
                 stringsAsFactors = FALSE)
hall$pathway_clean <- clean_pathway(hall$pathway)

nes_mat <- matrix(NA, nrow = nrow(hall), ncol = 2)
rownames(nes_mat) <- hall$pathway_clean

tc_nes <- find_col(hall, c("NES_TC", "NES_transcriptomics", "NES.TC"))
pr_nes <- find_col(hall, c("NES_PR", "NES_proteomics", "NES.PR"))
if (!is.null(tc_nes) && !is.null(pr_nes)) {
  nes_mat[, 1] <- hall[[tc_nes]]
  nes_mat[, 2] <- hall[[pr_nes]]
} else {
  nes_cols <- grep("NES", colnames(hall), value = TRUE)
  if (length(nes_cols) >= 2) {
    nes_mat[, 1] <- hall[[nes_cols[1]]]
    nes_mat[, 2] <- hall[[nes_cols[2]]]
  }
}
colnames(nes_mat) <- c("Transcriptomics", "Proteomics")

# Filter by significance if padj columns exist
padj_cols <- grep("padj|P_adj|FDR|q.val", colnames(hall), value = TRUE)
if (length(padj_cols) >= 1) {
  sig_idx <- apply(hall[, padj_cols, drop = FALSE], 1, function(x) any(x < 0.25, na.rm = TRUE))
  if (sum(sig_idx) > 5) nes_mat <- nes_mat[sig_idx, , drop = FALSE]
}
keep <- complete.cases(nes_mat)
nes_mat <- nes_mat[keep, , drop = FALSE]
nes_mat <- nes_mat[order(rowMeans(nes_mat)), , drop = FALSE]

s <- sp("Fig2d_Hallmark_NES_heatmap.pdf", W, H)
cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
ht <- Heatmap(nes_mat, name = "NES", col = col_nes,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, row_names_side = "left",
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_names_gp = gp_row_names(7), column_names_gp = gp_col_names(9),
  column_names_rot = 0, border = TRUE, rect_gp = gp_border(),
  width = unit(40, "mm"), heatmap_legend_param = std_legend_param())
draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()
cat("  -> Fig2d_Hallmark_NES_heatmap.pdf\n")
