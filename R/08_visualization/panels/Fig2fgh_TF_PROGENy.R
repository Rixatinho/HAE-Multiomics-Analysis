#!/usr/bin/env Rscript
# =============================================================================
# Fig2fgh_TF_PROGENy.R -- DoRothEA TF Activity & PROGENy Pathway Activity
# =============================================================================
# Figure 2f: DoRothEA TF activity heatmap (59 x 75 mm)
# Figure 2g: PROGENy pathway activity heatmap (59 x 75 mm)
# Figure 2h: PROGENy pathway activity boxplot (59 x 75 mm)
# =============================================================================
# Usage:
#   /Users/rishat/miniforge3/envs/multiomics/bin/Rscript Fig2fgh_TF_PROGENy.R
# =============================================================================

cat("=== Figure 2f-h: DoRothEA TF & PROGENy Pathway Activity ===\n")

# =============================================================================
# Load Libraries
# =============================================================================
suppressPackageStartupMessages({
  library(decoupleR)
  library(OmnipathR)
  library(progeny)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(grid)
})

# =============================================================================
# Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "analysis/results")
PROC <- file.path(BASE, "analysis/data/processed")
OUT  <- file.path(BASE, "图片/Figure_2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants (NC/J Hep standard)
# =============================================================================
FONT_FAMILY <- "Helvetica"
MM_PER_INCH <- 25.4
W_3COL <- 59  # 3-column width
H_STD  <- 75  # standard height

# Register fonts
pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# Color palette
COL_NORMAL   <- "#4DBBD5"
COL_ADJACENT <- "#E64B35"

# =============================================================================
# Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 9, face = "bold", hjust = 0),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.key.size = unit(2, "mm"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# Helper functions
# =============================================================================
gp_row <- function(sz = 6) gpar(fontsize = sz, fontfamily = FONT_FAMILY)
gp_col <- function(sz = 6) gpar(fontsize = sz, fontfamily = FONT_FAMILY, fontface = "bold")

col_zscore <- colorRamp2(c(-2, 0, 2), c("#3C5488", "#FFFFFF", "#E64B35"))

# =============================================================================
# Load Data
# =============================================================================
cat(">>> Loading data...\n")

# Load DEGs for gene mapping
degs <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"), 
                 check.names = FALSE, stringsAsFactors = FALSE)

# Load VST expression matrix
tc_vst_raw <- read.csv(file.path(PROC, "transcriptomics_vst_paired.csv"),
                       check.names = FALSE, row.names = 1)

# Map ENSG to gene symbol
id2sym <- setNames(degs$gene_name, degs$gene_id)
id2sym <- id2sym[!is.na(id2sym) & id2sym != ""]
tc_sym <- id2sym[rownames(tc_vst_raw)]
keep <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst_raw[keep, ])
rownames(tc_mat) <- tc_sym[keep]

# Define groups
group_vec <- ifelse(grepl("^Normal", colnames(tc_mat)), "Normal", "Adjacent")
names(group_vec) <- colnames(tc_mat)

adj_samps <- names(group_vec)[group_vec == "Adjacent"]
nor_samps <- names(group_vec)[group_vec == "Normal"]

cat(sprintf("  TC matrix: %d genes x %d samples\n", nrow(tc_mat), ncol(tc_mat)))

# =============================================================================
# 1. TF Activity (DoRothEA via decoupleR)
# =============================================================================
cat("\n>>> 1. TF activity inference (DoRothEA)...\n")

# Get DoRothEA regulons
dorothea_net <- tryCatch({
  get_dorothea(organism = "human", levels = c("A", "B", "C"))
}, error = function(e) {
  cat("  DoRothEA error:", e$message, "\n")
  cat("  Trying collectri...\n")
  get_collectri(organism = "human")
})

cat(sprintf("  Regulon: %d TFs\n", length(unique(dorothea_net$source))))

# Run ULM
tf_activity <- run_ulm(mat = tc_mat, net = dorothea_net, .source = "source",
                       .target = "target", .mor = "mor", minsize = 5)

# Pivot to matrix
tf_mat <- tf_activity %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  as.data.frame()
rownames(tf_mat) <- tf_mat$source
tf_mat$source <- NULL
tf_mat <- as.matrix(tf_mat)

cat(sprintf("  TF matrix: %d TFs x %d samples\n", nrow(tf_mat), ncol(tf_mat)))

# Test TF differences: Adjacent vs Normal
tf_diff <- data.frame()
for (tf in rownames(tf_mat)) {
  tryCatch({
    wt <- wilcox.test(tf_mat[tf, adj_samps], tf_mat[tf, nor_samps])
    fc <- mean(tf_mat[tf, adj_samps]) - mean(tf_mat[tf, nor_samps])
    tf_diff <- rbind(tf_diff, data.frame(
      TF = tf, diff = fc, pvalue = wt$p.value, stringsAsFactors = FALSE))
  }, error = function(e) NULL)
}
tf_diff$padj <- p.adjust(tf_diff$pvalue, method = "BH")
tf_diff <- tf_diff[order(tf_diff$pvalue), ]

# =============================================================================
# Figure 2f: Top TF Activity Heatmap
# =============================================================================
cat("\n>>> Figure 2f: TF activity heatmap...\n")

top_tfs <- head(tf_diff$TF, 20)
tf_plot <- tf_mat[top_tfs, , drop = FALSE]

# Scale rows
tf_plot_z <- t(scale(t(tf_plot)))
tf_plot_z[tf_plot_z > 2] <- 2
tf_plot_z[tf_plot_z < -2] <- -2

# Sample order
samp_order <- c(sort(nor_samps), sort(adj_samps))
samp_order <- intersect(samp_order, colnames(tf_plot_z))
tf_plot_z <- tf_plot_z[, samp_order]

# Annotation
ha <- HeatmapAnnotation(
  Group = group_vec[samp_order],
  col = list(Group = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT)),
  annotation_name_gp = gpar(fontsize = 6, fontfamily = FONT_FAMILY),
  simple_anno_size = unit(2, "mm"),
  show_legend = TRUE
)

# Save PDF
pdf_path <- file.path(OUT, "Fig2f_DoRothEA_TF.pdf")
cairo_pdf(pdf_path, width = W_3COL / MM_PER_INCH, height = H_STD / MM_PER_INCH, family = FONT_FAMILY)

ht <- Heatmap(tf_plot_z, name = "Z-score",
              col = col_zscore,
              top_annotation = ha,
              cluster_rows = TRUE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = FALSE,
              row_names_side = "right",
              row_names_gp = gp_row(5),
              column_title = "DoRothEA TF Activity",
              column_title_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY, fontface = "bold"),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6, fontfamily = FONT_FAMILY),
                labels_gp = gpar(fontsize = 5, fontfamily = FONT_FAMILY),
                legend_height = unit(15, "mm"), grid_width = unit(2, "mm")
              ),
              width = unit(35, "mm"))

draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
cat(sprintf("  Saved: %s\n", pdf_path))

# =============================================================================
# 2. PROGENy Pathway Activity
# =============================================================================
cat("\n>>> 2. PROGENy pathway activity...\n")

pw_scores <- progeny(tc_mat, scale = FALSE, organism = "Human", top = 500,
                     perm = 1000, verbose = FALSE)
pw_mat <- t(pw_scores)
rownames(pw_mat) <- gsub("\\.", "-", rownames(pw_mat))

# Remove zero-variance pathways
pw_var <- apply(pw_mat, 1, var, na.rm = TRUE)
pw_mat <- pw_mat[pw_var > 1e-10, , drop = FALSE]

cat(sprintf("  Pathway matrix: %d pathways x %d samples\n", nrow(pw_mat), ncol(pw_mat)))

# =============================================================================
# Figure 2g: PROGENy Heatmap
# =============================================================================
cat("\n>>> Figure 2g: PROGENy heatmap...\n")

pw_plot_z <- t(scale(t(pw_mat)))
pw_plot_z[pw_plot_z > 2] <- 2
pw_plot_z[pw_plot_z < -2] <- -2
pw_plot_z <- pw_plot_z[, samp_order]

pdf_path <- file.path(OUT, "Fig2g_PROGENy_heatmap.pdf")
cairo_pdf(pdf_path, width = W_3COL / MM_PER_INCH, height = H_STD / MM_PER_INCH, family = FONT_FAMILY)

ht <- Heatmap(pw_plot_z, name = "Z-score",
              col = col_zscore,
              top_annotation = ha,
              cluster_rows = TRUE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = FALSE,
              row_names_side = "right",
              row_names_gp = gp_row(6),
              column_title = "PROGENy Pathway Activity",
              column_title_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY, fontface = "bold"),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6, fontfamily = FONT_FAMILY),
                labels_gp = gpar(fontsize = 5, fontfamily = FONT_FAMILY),
                legend_height = unit(15, "mm"), grid_width = unit(2, "mm")
              ),
              width = unit(35, "mm"))

draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
cat(sprintf("  Saved: %s\n", pdf_path))

# =============================================================================
# Figure 2h: PROGENy Boxplot
# =============================================================================
cat("\n>>> Figure 2h: PROGENy boxplot...\n")

# Build long-format data
pw_long <- data.frame()
for (pw in rownames(pw_mat)) {
  pw_long <- rbind(pw_long, data.frame(
    pathway = pw,
    score = c(pw_mat[pw, nor_samps], pw_mat[pw, adj_samps]),
    group = c(rep("Normal", length(nor_samps)), rep("Adjacent", length(adj_samps))),
    stringsAsFactors = FALSE
  ))
}

# Calculate p-values for annotation
pw_pvals <- sapply(unique(pw_long$pathway), function(p) {
  sub <- pw_long[pw_long$pathway == p, ]
  tryCatch({
    wilcox.test(score ~ group, data = sub)$p.value
  }, error = function(e) NA)
})
pw_sig <- names(pw_pvals)[pw_pvals < 0.05]

p <- ggplot(pw_long, aes(x = reorder(pathway, score, FUN = median), y = score, fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT)) +
  coord_flip() +
  labs(title = "PROGENy Activity", x = "", y = "Activity Score", fill = "Group") +
  theme_nc +
  theme(
    legend.position = "right",
    legend.key.size = unit(2, "mm"),
    axis.text.y = element_text(size = 6, face = ifelse(levels(reorder(pw_long$pathway, pw_long$score, FUN = median)) %in% pw_sig, "bold", "plain"))
  )

pdf_path <- file.path(OUT, "Fig2h_PROGENy_boxplot.pdf")
ggsave(pdf_path, p, width = W_3COL / MM_PER_INCH, height = H_STD / MM_PER_INCH, 
       device = cairo_pdf, family = FONT_FAMILY)
cat(sprintf("  Saved: %s\n", pdf_path))

# =============================================================================
# Summary
# =============================================================================
cat("\n=== Figure 2f-h Complete ===\n")
cat("Output files:\n")
cat(sprintf("  - Fig2f_DoRothEA_TF.pdf (%d TFs)\n", length(top_tfs)))
cat(sprintf("  - Fig2g_PROGENy_heatmap.pdf (%d pathways)\n", nrow(pw_mat)))
cat(sprintf("  - Fig2h_PROGENy_boxplot.pdf (%d significant)\n", length(pw_sig)))
