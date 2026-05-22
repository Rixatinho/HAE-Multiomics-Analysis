#!/usr/bin/env Rscript
# =============================================================================
# Figure_7_standalone.R
# Title: Machine Learning Biomarker Discovery & External Validation
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Nature Communications - Cairo PDF | Arial 8-12pt (8pt floor, gradient hierarchy)
# =============================================================================
# Figure 7: Machine learning reveals candidate multi-omics biomarker panels
# 5 panels:
#   a - Top single-feature classifier per omics (faceted lollipop, AUC + 95% CI)
#   b - Cross-cohort effect-size bubble matrix (top 10 genes x 6 cohorts + pooled, 130 samples)
#   c - External validation AUC barplot (5 independent cohorts)
#   d - ML ensemble classification (CS1 vs CS2 LOOCV, RF/LASSO/SVM/XGB/Stack)
#   e - Feature overlap Venn diagram (MOFA hubs / RF / Differential -> FGG consensus)
# =============================================================================
# Usage:
#   /Users/rishat/miniforge3/envs/multiomics/bin/Rscript Figure_7_standalone.R
# =============================================================================

cat("=== Figure 7: Machine Learning Biomarker Discovery & Validation ===\n")
cat("  Loading libraries...\n")


# =============================================================================
# SECTION 1: Library Imports + Arial Font Registration
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(ggsci)
  library(metafor)
})

# --- Register Arial in R's PostScript/PDF font databases ---
# cairo_pdf natively embeds ArialMT TrueType from the system, but R's grid
# package uses the PostScript font database for text-width calculation. Without
# this registration, grid emits harmless but noisy "font family 'Arial' not
# found in PostScript font database" warnings.  Mapping Arial -> Helvetica
# metrics is safe because the two fonts are metrically identical.
pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
cat("  Arial registered in PDF/PostScript font databases\n")


# =============================================================================
# SECTION 2: Project Paths (EDIT THESE IF RUNNING ON A DIFFERENT MACHINE)
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/main/Figure_7")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants (J Hep Standard)
# =============================================================================
FONT_FAMILY <- "Arial"      # cairo_pdf maps to Arial
FONT_GRID   <- "Arial"      # unified for ComplexHeatmap
MM_PER_INCH <- 25.4

# J Hep page dimensions (mm)
W_SINGLE  <- 89    # single column
W_DOUBLE  <- 183   # double column
W_HALF    <- 89    # single column = 89mm
H_STD     <- 85    # standard height
H_TALL    <- 100   # tall panel
H_MAX     <- 240

ASSEMBLY_DPI <- 600   # pixel-exact rendering constant

mm2in <- function(mm) mm / MM_PER_INCH
in2mm <- function(inch) inch * MM_PER_INCH
# Typography hierarchy (minimum 8pt rule)
FS_TITLE      <- 10    # plot.title, column_title (bold)
FS_SUBTITLE   <- 9     # plot.subtitle
FS_AXIS_TITLE <- 9     # axis.title (bold)
FS_AXIS_TEXT  <- 8     # axis.text, strip.text (minimum)
FS_LEGEND_T   <- 8     # legend.title (bold)
FS_LEGEND_L   <- 8     # legend.text
FS_ANNO       <- 8     # annotation name
FS_ROW_NAME   <- 8     # heatmap row/column names
FS_CELL       <- 8     # heatmap cell text
FS_TAG        <- 12    # panel tag (A, B, C...) — aligned with Figure 6
FS_GEOM_TEXT  <- 2.82  # geom_text size (= 8pt in mm)

# =============================================================================
# SECTION 4: Color Palette (ggsci JCO / Journal of Clinical Oncology)
# =============================================================================
PAL_CAT <- pal_jco("default")(10)

# Semantic colors
COL_UP       <- "#CD534CFF"    # up-regulated (NPG red)
COL_DOWN     <- "#0073C2FF"    # down-regulated (NPG cyan)
COL_NS       <- "#868686FF"    # not significant
COL_NA       <- "#F0F0F0"    # NA / background

# Omics layer colors
COL_TC       <- "#0073C2FF"    # transcriptomics (navy)
COL_PR       <- "#CD534CFF"    # proteomics (red)
COL_MT       <- "#EFC000FF"    # metabolomics (teal)

# Group colors
COL_NORMAL   <- "#7AA6DCFF"    # normal tissue (cyan)
COL_ADJACENT <- "#CD534CFF"    # adjacent/disease tissue (red)

# Subtype colors
COL_CS1      <- "#CD534CFF"    # consensus subtype 1 (red) / Adjacent enriched
COL_CS2      <- "#0073C2FF"    # consensus subtype 2 (navy) / Normal enriched

# =============================================================================
# SECTION 5: Color Scale Functions (circlize::colorRamp2)
# =============================================================================
col_div     <- colorRamp2(c(-2, 0, 2),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_nes     <- colorRamp2(c(-3, 0, 3),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_cor     <- colorRamp2(c(-1, 0, 1),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_zscore  <- colorRamp2(c(-2.5, 0, 2.5), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_immune  <- colorRamp2(c(-2, 0, 2),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 6: ggplot2 Theme (theme_bw base, J Hep style)
# =============================================================================
# J Hep requirement: ALL text elements must explicitly set family = "Arial"
# Font sizes: axis.text 8pt, axis.title 9pt, plot.title 10pt bold
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text               = element_text(family = FONT_FAMILY, size = 8),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.line          = element_line(linewidth = 0.4, color = "black"),
    axis.ticks         = element_line(linewidth = 0.4, color = "black"),
    axis.ticks.length  = unit(1.5, "mm"),
    axis.text          = element_text(family = FONT_FAMILY, size = 8, color = "black"),
    axis.title         = element_text(family = FONT_FAMILY, size = 9, face = "bold", color = "black"),
    plot.title         = element_text(family = FONT_FAMILY, size = 10, face = "bold", hjust = 0.5),
    plot.title.position = "panel",
    legend.text        = element_text(family = FONT_FAMILY, size = 8),
    legend.title       = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    legend.key.size    = unit(3.5, "mm"),
    legend.background  = element_blank(),
    legend.box.background = element_blank(),
    legend.spacing.y   = unit(1, "mm"),
    strip.text         = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    strip.background   = element_rect(fill = "grey96", linewidth = 0.3, color = "grey80"),
    plot.margin        = margin(5, 5, 5, 5, "mm")
  )
theme_set(theme_nc)
# Alias for compatibility
theme_pub <- theme_nc


# =============================================================================
# SECTION 7: ComplexHeatmap gpar Factory Functions
# =============================================================================
gp_row_names <- function(size = 8, italic = FALSE) {
  gpar(fontsize = size, fontfamily = FONT_GRID,
       fontface = if (italic) "italic" else "plain")
}

gp_col_names <- function(size = 8, bold = TRUE) {
  gpar(fontsize = size, fontfamily = FONT_GRID,
       fontface = if (bold) "bold" else "plain")
}

gp_legend_title <- function(size = 8) {
  gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
}

gp_legend_labels <- function(size = 8) {
  gpar(fontsize = size, fontfamily = FONT_GRID)
}

gp_cell_text <- function(size = 8) {
  gpar(fontsize = size, fontfamily = FONT_GRID)
}

gp_anno_name <- function(size = 8) {
  gpar(fontsize = size, fontfamily = FONT_GRID)
}

gp_border <- function() {
  gpar(col = "white", lwd = 0.3)
}

gp_row_title <- function(size = 8) {
  gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
}

std_legend_param <- function() {
  list(
    title_gp      = gp_legend_title(),
    labels_gp     = gp_legend_labels(),
    legend_height = unit(20, "mm"),
    grid_width    = unit(3, "mm")
  )
}

# =============================================================================
# SECTION 8: Helper Functions
# =============================================================================

# --- Save helper ---
sp <- function(filename, w_mm, h_mm, out_dir = NULL) {
  if (is.null(out_dir)) out_dir <- OUT
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  list(path = file.path(out_dir, filename),
       width = target_w_px / ASSEMBLY_DPI, height = target_h_px / ASSEMBLY_DPI)
}

save_pdf <- function(filename, w_mm, h_mm, expr, out_dir = NULL) {
  s <- sp(filename, w_mm, h_mm, out_dir)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  force(expr)
  dev.off()
  cat("  PDF saved:", s$path, "\n")
  invisible(s$path)
}

# --- CSV reader ---
read_csv_safe <- function(filepath, ...) {
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE, ...)
  if (ncol(df) > 1) {
    first_col <- df[[1]]
    if (is.character(first_col) || is.factor(first_col)) {
      rownames(df) <- make.unique(as.character(first_col))
      df[[1]] <- NULL
    }
  }
  df
}

# --- Matrix scaling ---
clamp_matrix <- function(mat, lim = 2) {
  mat[mat > lim] <- lim; mat[mat < -lim] <- -lim; mat
}

scale_rows <- function(mat, lim = 2) {
  mat <- t(scale(t(as.matrix(mat)))); mat[is.na(mat)] <- 0; clamp_matrix(mat, lim)
}

# --- Column finder ---
find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else NULL
}



# =============================================================================
# SECTION 9: Data Loading
# =============================================================================
cat("  Loading biomarker and ML data...\n")

cv_summary <- read.csv(file.path(RES, "phase8_biomarker/CV_fold_AUC_summary.csv"),
                       stringsAsFactors = FALSE)
roc_ci     <- read.csv(file.path(RES, "phase8_biomarker/ROC_confidence_intervals.csv"),
                       stringsAsFactors = FALSE)
forest_df  <- read.csv(file.path(RES, "enhancement14_meta_analysis/forest_plot_data.csv"),
                       stringsAsFactors = FALSE)
val_auc    <- read.csv(file.path(RES, "enhancement14_meta_analysis/signature_validation_AUC.csv"),
                       stringsAsFactors = FALSE)
ml_comp    <- read.csv(file.path(RES, "enhancement17_ml_ensemble/model_comparison.csv"),
                       stringsAsFactors = FALSE)
rf_feat    <- read.csv(file.path(RES, "enhancement17_ml_ensemble/rf_top20_features.csv"),
                       stringsAsFactors = FALSE)
lasso_feat <- read.csv(file.path(RES, "enhancement17_ml_ensemble/lasso_top_features.csv"),
                       stringsAsFactors = FALSE)
mofa_wt    <- read.csv(file.path(RES, "phase5_integration/MOFA2_top_weights_per_factor.csv"),
                       stringsAsFactors = FALSE)
deg_sig    <- tryCatch(read.csv(file.path(RES, "phase1_diff/DEGs_significant.csv"),
                                stringsAsFactors = FALSE), error = function(e) NULL)

cat("  Data loaded.\n")



# --- Panel object registry (for optimized assembly) ---
# Strategy A: ggplot objects kept directly (best quality - vector)
# Strategy B: ComplexHeatmap captured via grid.grabExpr (single rasterization)
# Strategy C: External/pre-rendered PDF read at 600 DPI (fallback)
obj_A <- obj_B <- obj_C <- obj_D <- obj_E <- NULL
# =============================================================================
# SECTION 10: Render Figure 7 -- Machine Learning Biomarker Discovery
# =============================================================================
cat("\n--- Rendering Figure 7: Machine Learning Biomarker Discovery ---\n")

# -----------------------------------------------------------------------------
# Fig 7c: Top single-feature ROC per omics (horizontal lollipop with 95% CI)  [moved from old 7a]
# -----------------------------------------------------------------------------
cat("  Fig7c: Top single-feature classifiers\n")
tryCatch({
  omics_map <- c("TC" = "Transcriptomics", "PR" = "Proteomics", "MET" = "Metabolomics")
  omics_cols <- c("Transcriptomics" = COL_TC, "Proteomics" = COL_PR, "Metabolomics" = COL_MT)

  # Top 3 features per omics by AUC
  roc_top <- roc_ci %>%
    dplyr::mutate(omics_full = unname(omics_map[omics])) %>%
    dplyr::filter(!is.na(omics_full)) %>%
    dplyr::group_by(omics_full) %>%
    dplyr::arrange(dplyr::desc(AUC), .by_group = TRUE) %>%
    dplyr::slice_head(n = 3) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      omics_full = factor(omics_full,
                          levels = c("Transcriptomics", "Proteomics", "Metabolomics"))
    ) %>%
    dplyr::arrange(omics_full, AUC) %>%
    dplyr::mutate(feature = factor(feature, levels = feature))

  p <- ggplot(roc_top, aes(x = AUC, y = feature, color = omics_full)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_segment(aes(x = CI_lower, xend = CI_upper, yend = feature), linewidth = 0.5) +
    geom_point(size = 2.2) +
    facet_grid(omics_full ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(values = omics_cols, guide = "none") +
    scale_x_continuous(limits = c(0.4, 1.05), breaks = seq(0.4, 1, 0.2)) +
    labs(title = "Top single-feature classifiers",
         x = "AUC (95% CI)", y = NULL) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text.y.left = element_text(angle = 0, hjust = 1,
                                            face = "bold", size = 8),
          panel.spacing.y = unit(1, "mm"),
          axis.text.y = element_text(size = 8, family = FONT_FAMILY))
  save_pdf("Fig7c_top_features.pdf", 91, 80, print(p))
  obj_A <<- p
}, error = function(e) cat("    ERROR Fig7b:", e$message, "\n"))


# -----------------------------------------------------------------------------
# Fig 7b: Meta-analysis forest plot (6 cohorts)
# -----------------------------------------------------------------------------
cat("  Fig7a: Cross-cohort effect-size overview (bubble heatmap)\n")
tryCatch({
  # Top 10 genes by # datasets covered
  gene_count <- table(forest_df$gene)
  top_genes <- names(sort(gene_count, decreasing = TRUE))[1:min(10, length(gene_count))]
  fdf <- forest_df[forest_df$gene %in% top_genes &
                     !is.na(forest_df$logFC) & !is.na(forest_df$SE) &
                     forest_df$SE > 0, ]

  # Random-effects meta-analysis per gene (REML)
  meta_summary <- do.call(rbind, lapply(top_genes, function(g) {
    sub <- fdf[fdf$gene == g, ]
    if (nrow(sub) < 2) return(NULL)
    res <- tryCatch(metafor::rma(yi = sub$logFC, sei = sub$SE,
                                  method = "REML", control = list(maxiter = 200)),
                    error = function(e) NULL,
                    warning = function(w) NULL)
    if (is.null(res)) return(NULL)
    data.frame(gene = g,
               pooled  = as.numeric(res$beta),
               I2      = res$I2,
               pooled_p = res$pval,
               n_studies = nrow(sub),
               stringsAsFactors = FALSE)
  }))

  # Caps for visual scale
  CAP_FC <- 5    # |logFC| display cap
  CAP_LP <- 8    # -log10(P) size cap

  # Per-study cells
  study_cells <- data.frame(
    gene       = fdf$gene,
    dataset    = fdf$dataset,
    logFC_cap  = pmax(pmin(fdf$logFC, CAP_FC), -CAP_FC),
    nlp        = pmin(-log10(pmax(fdf$pvalue, 1e-15)), CAP_LP),
    is_pooled  = FALSE,
    stringsAsFactors = FALSE
  )

  # Pooled cells (right-most "Pooled" column)
  pooled_cells <- data.frame(
    gene       = meta_summary$gene,
    dataset    = "Pooled",
    logFC_cap  = pmax(pmin(meta_summary$pooled, CAP_FC), -CAP_FC),
    nlp        = pmin(-log10(pmax(meta_summary$pooled_p, 1e-15)), CAP_LP),
    is_pooled  = TRUE,
    stringsAsFactors = FALSE
  )

  cells <- rbind(study_cells, pooled_cells)

  ds_order  <- c("HAE_our", "GSE184297", "GSE124362", "GSE278225",
                 "GSE24376_10984", "GSE24376_10985", "Pooled")
  ds_order <- c("HAE_our", "GSE184297", "GSE124362", "GSE278225",
                "GSE24376_10984", "GSE24376_10985", "Pooled", "I2")
  cells$dataset <- factor(cells$dataset, levels = ds_order)
  # Order genes by absolute pooled effect (most informative on top)
  ord <- meta_summary$gene[order(-abs(meta_summary$pooled))]
  cells$gene <- factor(cells$gene, levels = rev(ord))

  # I2 annotation strip (dedicated x-axis column, drawn inside panel)
  i2_lab <- data.frame(
    gene    = factor(meta_summary$gene, levels = rev(ord)),
    dataset = factor("I2", levels = ds_order),
    label   = sprintf("I\u00b2=%.0f%%", meta_summary$I2)
  )

  p <- ggplot(cells, aes(x = dataset, y = gene)) +
    geom_vline(xintercept = 6.5, color = "grey50",
               linetype = "dashed", linewidth = 0.3) +
    geom_point(aes(color = logFC_cap, size = nlp,
                   shape = is_pooled)) +
    scale_x_discrete(drop = FALSE,
                     labels = c("HAE_our" = "HAE_our",
                                "GSE184297" = "GSE184297",
                                "GSE124362" = "GSE124362",
                                "GSE278225" = "GSE278225",
                                "GSE24376_10984" = "GSE24376_10984",
                                "GSE24376_10985" = "GSE24376_10985",
                                "Pooled" = "Pooled",
                                "I2" = "")) +
    scale_color_gradient2(low = COL_DOWN, mid = "white", high = COL_UP,
                          midpoint = 0, limits = c(-CAP_FC, CAP_FC),
                          name = "log2 FC", oob = scales::squish,
                          breaks = c(-4, -2, 0, 2, 4),
                          guide = guide_colourbar(raster = FALSE, nbin = 50,
                                                  barwidth = unit(3, "mm"),
                                                  barheight = unit(20, "mm"))) +
    scale_size_continuous(range = c(1, 5), limits = c(0, CAP_LP),
                          name = "-log10(P)",
                          breaks = c(2, 4, 6, 8)) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18),
                       guide = "none") +
    geom_text(data = i2_lab,
              aes(x = dataset, y = gene, label = label),
              hjust = 0.5, size = FS_GEOM_TEXT,
              family = FONT_FAMILY, color = "grey30") +
    labs(title = "Cross-cohort effect-size overview: top 10 genes \u00d7 6 cohorts + meta-pool",
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8, face = "italic"),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
          plot.margin = margin(5, 5, 5, 5, "mm"),
          legend.position = "right",
          legend.box = "vertical",
          legend.key.size = unit(3, "mm"))
  save_pdf("Fig7a_overview_bubble.pdf", 183, 85, print(p))
  obj_B <<- p
}, error = function(e) cat("    ERROR Fig7a:", e$message, "\n"))

# -----------------------------------------------------------------------------
# Fig 7b: External validation AUC  [moved from old 7c] barplot
# -----------------------------------------------------------------------------
cat("  Fig7b: External validation AUC\n")
tryCatch({
  val_auc$dataset <- factor(val_auc$dataset, levels = val_auc$dataset)
  val_auc$sig <- val_auc$wilcox_P < 0.05
  val_auc$sig_label <- ifelse(val_auc$wilcox_P < 0.001, "***",
                       ifelse(val_auc$wilcox_P < 0.01,  "**",
                       ifelse(val_auc$wilcox_P < 0.05,  "*", "ns")))

  p <- ggplot(val_auc, aes(x = dataset, y = AUC, fill = sig)) +
    geom_col(width = 0.65, alpha = 0.9, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = AUC_CI_lower, ymax = pmin(AUC_CI_upper, 1.05)),
                  width = 0.2, linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", AUC), y = 0.05),
              family = FONT_FAMILY, size = FS_GEOM_TEXT,
              color = "white", fontface = "bold") +
    geom_text(aes(label = sig_label,
                  y = pmin(AUC_CI_upper, 1.05) + 0.04),
              family = FONT_FAMILY, size = FS_GEOM_TEXT, fontface = "bold") +
    scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = COL_PR),
                      labels = c("FALSE" = "ns (P\u22650.05)",
                                 "TRUE"  = "validated (P<0.05)"),
                      name = NULL) +
    scale_y_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.2)) +
    labs(title = "External validation across 5 independent cohorts",
         x = NULL, y = "AUC (signature score)") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          legend.position = "top", legend.key.size = unit(3, "mm"))
  save_pdf("Fig7b_validation_AUC.pdf", 92, 80, print(p))
  obj_C <<- p
}, error = function(e) cat("    ERROR Fig7b:", e$message, "\n"))


# -----------------------------------------------------------------------------
# Fig 7d: ML ensemble classification performance
# -----------------------------------------------------------------------------
cat("  Fig7d: ML ensemble performance\n")
tryCatch({
  # Sort models by AUC descending; highlight best (Stacking Ensemble)
  ml_order <- ml_comp$Model[order(-ml_comp$AUC)]
  ml_long <- reshape2::melt(ml_comp, id.vars = c("Model", "CV_Scheme"),
                            variable.name = "Metric", value.name = "Value")
  ml_long$Model  <- factor(ml_long$Model,  levels = ml_order)
  ml_long$Metric <- factor(ml_long$Metric, levels = c("Accuracy", "AUC", "F1"))

  p <- ggplot(ml_long, aes(x = Model, y = Value, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.2f", Value)),
              position = position_dodge(width = 0.7), vjust = -0.4,
              family = FONT_FAMILY, size = FS_GEOM_TEXT) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = c("Accuracy" = COL_TC, "AUC" = COL_PR, "F1" = COL_MT)) +
    scale_y_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.2)) +
    labs(title = "Subtype CS1 vs CS2 (LOOCV)",
         x = NULL, y = "Score", fill = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          plot.title  = element_text(size = FS_TITLE, face = "bold", hjust = 0.5),
          legend.position = "top", legend.key.size = unit(3, "mm"))
  save_pdf("Fig7d_model_comparison.pdf", 91, 80, print(p))
  obj_D <<- p
}, error = function(e) cat("    ERROR Fig7d:", e$message, "\n"))


# -----------------------------------------------------------------------------
# Fig 7e: Feature overlap Venn diagram (MOFA / RF / DE)
# Correct feature sets for consensus analysis:
#   MOFA: cross-factor hub features (appearing in 2+ MOFA factors)
#   RF:   all features from condition-level RF model (gene-level unique)
#   DE:   DEPs with nominal P < 0.05 (gene-level unique)
# Result: FGG is the sole feature appearing in all three sets
# -----------------------------------------------------------------------------
cat("  Fig7e: Feature Venn diagram\n")
tryCatch({
  # 1. MOFA set: cross-factor hub features (features in 2+ MOFA factors)
  hub_file <- file.path(RES, "enhancement22_dl_fusion/cross_factor_hub_features.csv")
  hub <- read.csv(hub_file, stringsAsFactors = FALSE)
  mofa_set <- hub$feature

  # 2. RF set: all features from condition-level RF model (gene-level unique)
  rf_cond_file <- file.path(RES, "enhancement22_dl_fusion/rf_importance_condition.csv")
  rf_cond <- read.csv(rf_cond_file, stringsAsFactors = FALSE)
  rf_set <- unique(gsub("_TC$|_PR$|_MET$", "", rf_cond$feature))

  # 3. Diff set: DEPs with nominal P < 0.05 (gene-level unique)
  dep_file <- file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv")
  dep_full <- read.csv(dep_file, stringsAsFactors = FALSE)
  diff_set <- unique(dep_full$gene_name[dep_full$P.Value < 0.05])

  cat("    MOFA hubs:", length(mofa_set), "| RF condition:", length(rf_set),
      "| DEPs P<0.05:", length(diff_set), "\n")

  # Verify FGG consensus
  all_three <- Reduce(intersect, list(mofa_set, rf_set, diff_set))
  cat("    Consensus (all 3):", paste(all_three, collapse = ", "), "\n")

  # Build Venn using VennDiagram
  suppressPackageStartupMessages(library(VennDiagram))
  venn_list <- list(
    "MOFA2\nhubs" = mofa_set,
    "RF\nimportance" = rf_set,
    "Differential\nexpression" = diff_set
  )

  s <- sp("Fig7e_feature_venn.pdf", 92, 80)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  grid::grid.newpage()
  vp <- venn.diagram(
    venn_list, filename = NULL,
    fill = c(COL_TC, COL_PR, COL_MT), alpha = 0.35,
    cat.cex = 0.75, cex = 0.7,
    cat.fontfamily = FONT_FAMILY,
    fontfamily = FONT_FAMILY,
    cat.dist = c(0.06, 0.06, 0.04),
    main = "Multi-method feature consensus",
    main.fontfamily = FONT_FAMILY,
    main.cex = 0.85,
    sub = paste0("Center: ", paste(all_three, collapse = ", ")),
    sub.fontfamily = FONT_FAMILY,
    sub.cex = 0.75
  )
  grid::grid.draw(vp)
  dev.off()
  obj_E <<- vp  # Store Venn grob for vector assembly
  cat("    -> Fig7e_feature_venn.pdf\n")
}, error = function(e) cat("    ERROR Fig7e:", e$message, "\n"))

# =============================================================================
# SECTION 11: Composite Figure Assembly (Pure Vector Grid Viewport — AI-editable)
# =============================================================================
cat("\n--- Assembling composite Figure_7 (VECTOR, 183x245mm) ---\n")

tryCatch({
  DPI <- 600
  W_TOTAL <- 183; H_TOTAL <- 245

  # Row 1 (85mm): A(183) full-width  [overview bubble]
  # Row 2 (80mm): B(92) + C(91) = 183  [validation AUC | top single-feature]
  # Row 3 (80mm): D(91) + E(92) = 183  [unchanged]
  # New layout: A full-width top, B+C middle row, D+E bottom row
  W_B <- 92; W_C <- W_TOTAL - W_B  # middle row split (B left = validation AUC, C right = top features)
  W_D <- 91; W_E <- W_TOTAL - W_D  # bottom row split unchanged
  H1 <- 85; H2 <- 80; H3 <- H_TOTAL - H1 - H2  # row heights: A=85, B/C=80, D/E=80

  # --- Render function: places all panels + labels into current device ---
  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))

    # Row 1 (top): panel a full-width (overview bubble) — y = H2 + H3 from bottom
    y_row1 <- H2 + H3
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row1, "mm"),
      width = unit(W_TOTAL, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_B, newpage = FALSE)  # obj_B holds the overview-bubble plot, now panel A
    grid::popViewport()

    # Row 2 (middle): panel b (validation AUC, was obj_C) | panel c (top features, was obj_A)
    y_row2 <- H3

    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_B, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_C, newpage = FALSE)  # obj_C holds the validation-AUC plot, now panel B
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x = unit(W_B, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_C, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_A, newpage = FALSE)  # obj_A holds the top-feature plot, now panel C
    grid::popViewport()

    # Row 3 (bottom): panels d & e — y = 0
    # Panel d
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W_D, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_D, newpage = FALSE)
    grid::popViewport()

    # Panel e (VennDiagram grob — use grid.draw)
    grid::pushViewport(grid::viewport(
      x = unit(W_D, "mm"), y = unit(0, "mm"),
      width = unit(W_E, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")
    ))
    grid::grid.draw(obj_E)
    grid::popViewport()

    # --- Panel labels (bold, FS_TAG = 12pt) ---
    # Layout: a=top-left (full-width row1), b=row2-left, c=row2-right, d=row3-left, e=row3-right
    label_data <- data.frame(
      text = c("a", "b", "c", "d", "e"),
      x_mm = c(1, 1, W_B + 1, 1, W_D + 1),
      y_mm = c(H_TOTAL - 1, y_row2 + H2 - 1, y_row2 + H2 - 1, H3 - 1, H3 - 1),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(label_data))) {
      grid::grid.text(
        label = label_data$text[i],
        x = unit(label_data$x_mm[i], "mm"),
        y = unit(label_data$y_mm[i], "mm"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY)
      )
    }

    grid::popViewport()  # pop full-page viewport
  }

  # --- Save vector PDF (AI-editable) ---
  cairo_pdf(file.path(OUT, "Figure_7.pdf"),
            width = W_TOTAL / 25.4, height = H_TOTAL / 25.4, family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> Figure_7.pdf (VECTOR, AI-editable)\n")

  # --- Save PNG (600 DPI for review) ---
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "Figure_7.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_7.png\n")

  # --- Save TIFF (600 DPI for journal submission) ---
  grDevices::tiff(file.path(OUT, "Figure_7.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_7.tiff\n")

  cat("  Figure_7 DONE (183x245mm, VECTOR PDF + 600DPI PNG/TIFF)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

# =============================================================================
cat("\n=== Figure 7 rendering complete ===\n")
