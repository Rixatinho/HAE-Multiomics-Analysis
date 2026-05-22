#!/usr/bin/env Rscript
# =============================================================================
# Figure_5_standalone.R
# Complete, Self-Contained Code for Figure 5
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Journal of Hepatology (Elsevier) - TIFF 600dpi + PDF (optimized assembly)
# =============================================================================
# Figure 5: Cell Communication and Master Regulators
# 7 panels:
#   (a) BayesPrism cell type composition stacked barplot
#   (b) CellChat pathway communication strength heatmap
#   (c) Ligand-receptor pair differential heatmap
#   (d) Ligand-receptor bubble plot
#   (e) TCR/BCR diversity boxplot
#   (f) GRN master regulator network
#   (g) TF centrality ranking bar plot
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript Figure_5_standalone.R
# =============================================================================

cat("=== Figure 5: Cell Communication and Master Regulators ===\n")
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
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggrepel)
  library(RColorBrewer)
  library(ggpubr)
  library(ggsci)
  library(magick)
})

# --- Register Arial in R's PostScript/PDF font databases ---
pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
cat("  Arial registered in PDF/PostScript font databases\n")


# =============================================================================
# SECTION 2: Project Paths (EDIT THESE IF RUNNING ON A DIFFERENT MACHINE)
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/main/Figure_5")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants (J Hepatol Standard)
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4

W_SINGLE  <- 89
W_DOUBLE  <- 183
W_HALF    <- 89
H_STD     <- 85
H_TALL    <- 100
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
FS_TAG        <- 16    # panel tag (A, B, C...) — gradient rule: 16pt > title 10pt > body 8pt
FS_GEOM_TEXT  <- 2.82  # geom_text size (= 8pt in mm)

# =============================================================================
# SECTION 4: Color Palette (ggsci JCO / Journal of Clinical Oncology)
# =============================================================================
PAL_CAT <- pal_jco("default")(10)

COL_UP       <- "#CD534CFF"
COL_DOWN     <- "#0073C2FF"
COL_NS       <- "#868686FF"
COL_NA       <- "#F0F0F0"

COL_TC       <- "#0073C2FF"
COL_PR       <- "#CD534CFF"
COL_MT       <- "#EFC000FF"

COL_NORMAL   <- "#7AA6DCFF"
COL_ADJACENT <- "#CD534CFF"

COL_CS1      <- "#CD534CFF"
COL_CS2      <- "#0073C2FF"

# =============================================================================
# SECTION 5: Color Scale Functions (circlize::colorRamp2)
# =============================================================================
col_div     <- colorRamp2(c(-2, 0, 2),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_nes     <- colorRamp2(c(-3, 0, 3),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_cor     <- colorRamp2(c(-1, 0, 1),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_zscore  <- colorRamp2(c(-2.5, 0, 2.5), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 6: ggplot2 Theme (theme_bw base, J Hepatol style)
# =============================================================================
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
    legend.spacing.x   = unit(2, "mm"),
    legend.margin      = margin(2, 2, 2, 2),
    strip.text         = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    strip.background   = element_rect(fill = "grey96", linewidth = 0.5, color = "grey80"),
    plot.margin        = margin(4, 4, 2, 4, "mm")
  )
theme_set(theme_nc)
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

# --- Save helper (pixel-exact: ensures PDF rasterized at 600 DPI = exact slot pixels) ---
sp <- function(filename, w_mm, h_mm, out_dir = NULL) {
  if (is.null(out_dir)) out_dir <- OUT
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  list(path = file.path(out_dir, filename),
       width = target_w_px / ASSEMBLY_DPI, height = target_h_px / ASSEMBLY_DPI)
}

save_pdf <- function(filename, w, h, expr) {
  fpath <- file.path(OUT, filename)
  target_w_px <- round(w * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h * ASSEMBLY_DPI / 25.4)
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  force(expr)
  dev.off()
  cat(sprintf("  Panel saved: %s (%d×%d px)\n", basename(fpath), target_w_px, target_h_px))
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

# --- Group helpers ---
make_group_vector <- function(sample_names) {
  ifelse(grepl("Normal|^N\\d", sample_names), "Normal", "Adjacent")
}

# --- Clean cell type names ---
clean_celltype <- function(x) {
  x <- gsub("_", " ", x)
  x
}

# --- Clean clinical variable names ---
clean_varname <- function(x) {
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(x)
  x <- gsub(" Pct$", " (%)", x)
  x <- gsub("^Wbc$", "WBC", x); x <- gsub("^Alt$", "ALT", x)
  x <- gsub("^Ast$", "AST", x); x <- gsub("^Ggt$", "GGT", x)
  x
}


# =============================================================================
# SECTION 9: Data Loading
# =============================================================================
DECONV_DIR <- file.path(RES, "enhancement16_deconvolution")
COMM_DIR   <- file.path(RES, "enhancement11_cell_communication")
TCR_DIR    <- file.path(RES, "enhancement21_tcr_bcr")
GRN_DIR    <- file.path(RES, "enhancement15_grn")

cat("  Data loaded inline in render logic below.\n")

# --- Panel object registry (for optimized assembly) ---
# Strategy A: ggplot objects kept directly (best quality - vector)
# Strategy B: ComplexHeatmap captured via grid.grabExpr (single rasterization)
# Strategy C: External/pre-rendered PDF read at 600 DPI (fallback)
obj_A <- obj_B <- obj_C <- obj_D <- obj_E <- obj_F <- obj_G <- NULL


# =============================================================================
# SECTION 10: Render Figure 5 -- Cell Communication & Master Regulators
# =============================================================================
cat("\n--- Rendering Figure 5: Cell Communication & Master Regulators ---\n")
cat("\n=== FIGURE 5 ===\n")


# ---------------------------------------------------------------------------
# Fig 5a: BayesPrism cell type composition stacked barplot
# ---------------------------------------------------------------------------
cat("  Fig5a\n")
tryCatch({
  bp_theta <- read.csv(file.path(DECONV_DIR, "bayesprism_cell_proportions.csv"),
                       stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

  props_long <- data.frame(
    sample = rep(rownames(bp_theta), ncol(bp_theta)),
    cell_type = rep(colnames(bp_theta), each = nrow(bp_theta)),
    proportion = as.vector(as.matrix(bp_theta)),
    stringsAsFactors = FALSE
  )
  props_long$condition <- ifelse(grepl("Normal", props_long$sample), "Normal", "Adjacent")

  n_types <- length(unique(props_long$cell_type))
  type_colors <- colorRampPalette(brewer.pal(min(12, n_types), "Set3"))(n_types)
  names(type_colors) <- unique(props_long$cell_type)

  props_long$cell_type_clean <- clean_celltype(props_long$cell_type)

  # --- Compute per-condition mean proportions + sort cell types by |delta| ---
  props_mean <- props_long %>%
    group_by(condition, cell_type_clean) %>%
    summarise(mean_prop = mean(proportion), .groups = "drop")

  props_wide <- props_mean %>%
    pivot_wider(names_from = condition, values_from = mean_prop, values_fill = 0) %>%
    mutate(delta = Adjacent - Normal)

  # Sort: cell types with largest absolute difference at top
  ct_ord_a <- props_wide %>%
    arrange(desc(abs(delta))) %>%
    pull(cell_type_clean)
  props_wide$cell_type_clean <- factor(props_wide$cell_type_clean, levels = ct_ord_a)
  props_long$ct_f             <- factor(props_long$cell_type_clean, levels = ct_ord_a)

  # --- Wilcoxon test per cell type → y-axis label colouring ---
  pvals_a <- sapply(ct_ord_a, function(ct) {
    a <- props_long$proportion[props_long$cell_type_clean == ct & props_long$condition == "Adjacent"]
    n <- props_long$proportion[props_long$cell_type_clean == ct & props_long$condition == "Normal"]
    if (length(a) < 2 || length(n) < 2) return(1)
    suppressWarnings(wilcox.test(a, n)$p.value)
  })
  delta_a <- props_wide$delta[match(ct_ord_a, as.character(props_wide$cell_type_clean))]
  # Diagnostic output — shows raw p-values and deltas
  cat("    [Panel A] Wilcoxon raw p-values and deltas per cell type:\n")
  for (i in seq_along(ct_ord_a)) {
    cat(sprintf("      %s: p=%.4f  delta=%.4f\n", ct_ord_a[i], pvals_a[i], delta_a[i]))
  }
  # Use raw p < 0.05 (no BH correction — small n, visualisation only)
  y_col_a <- ifelse(pvals_a < 0.05 & delta_a > 0, COL_ADJACENT,
             ifelse(pvals_a < 0.05 & delta_a < 0, COL_NORMAL, "grey30"))
  cat(sprintf("    [Panel A] coloured: red=%d blue=%d grey=%d\n",
              sum(pvals_a < 0.05 & delta_a > 0),
              sum(pvals_a < 0.05 & delta_a < 0),
              sum(pvals_a >= 0.05)))

  # --- Dumbbell (Cleveland) dot plot: group means + individual sample jitter ---
  p4a <- ggplot(props_wide, aes(y = cell_type_clean)) +
    # Background: individual sample jitter (no labels — only group colour matters)
    geom_jitter(data = props_long,
                aes(x = proportion, y = ct_f, color = condition),
                height = 0.18, size = 0.55, alpha = 0.38, show.legend = FALSE) +
    # Connecting segment between group means
    geom_segment(aes(x = Normal, xend = Adjacent,
                     y = cell_type_clean, yend = cell_type_clean),
                 color = "grey62", linewidth = 0.55) +
    # Group mean dots
    geom_point(aes(x = Normal,    color = "Normal"),   size = 2.8, alpha = 0.95) +
    geom_point(aes(x = Adjacent,  color = "Adjacent"), size = 2.8, alpha = 0.95) +
    scale_color_manual(
      values = c("Adjacent" = COL_ADJACENT, "Normal" = COL_NORMAL),
      name   = "Group") +
    scale_x_continuous(
      labels = function(x) paste0(round(x * 100, 1), "%"),
      expand = expansion(mult = c(0.02, 0.10))) +
    labs(x = "Cell Type Proportion", y = NULL,
         title = "BayesPrism Cell Type Composition") +
    theme(axis.text.y  = element_text(size = 7.5, color = y_col_a),
          axis.text.x  = element_text(size = 8),
          axis.title.x = element_text(size = 8.5, face = "bold"),
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 8, face = "bold"),
          legend.key.size  = unit(3.5, "mm"),
          legend.position  = "right",
          plot.margin  = margin(3, 4, 1, 4, "mm"))
  save_pdf("Fig5a_stacked_barplot.pdf", 183, 55, print(p4a))
  obj_A <- p4a  # Strategy A: keep ggplot object for direct assembly
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# ---------------------------------------------------------------------------
# Fig 5b: CellChat pathway-level differential signalling bar plot
#   Replaces the previous 5x5 sender-receiver communication-strength heatmap,
#   which did not match ¶34 (specific pathway scores) nor ¶161-(b) legend
#   ("Bar plot displays pathway score differences").
# ---------------------------------------------------------------------------
cat("  Fig5b\n")
tryCatch({
  pw <- read.csv(file.path(RES, "enhancement16_deconvolution",
                            "pathway_communication_summary.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
  pw <- pw[order(-pw$diff), ]
  pw$pathway <- factor(pw$pathway, levels = rev(pw$pathway))
  pw$direction <- ifelse(pw$diff > 0, "Up in Adjacent", "Down in Adjacent")
  pw$lab_x <- ifelse(pw$diff >= 0, pw$diff + 0.02, pw$diff - 0.02)
  pw$lab_h <- ifelse(pw$diff >= 0, 0, 1)
  cat(sprintf("    Panel B: %d pathways; UP=%d DOWN=%d; diff range [%.3f, %.3f]\n",
              nrow(pw), sum(pw$diff > 0), sum(pw$diff < 0),
              min(pw$diff), max(pw$diff)))
  cat(sprintf("    Top suppressed: %s (%.3f) | Top activated: %s (%.3f)\n",
              as.character(pw$pathway[which.min(pw$diff)]), min(pw$diff),
              as.character(pw$pathway[which.max(pw$diff)]), max(pw$diff)))

  p4b <- ggplot(pw, aes(y = pathway, x = diff, fill = direction)) +
    geom_col(width = 0.72, alpha = 0.92) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_text(aes(x = lab_x, label = sprintf("%+.2f", diff),
                  hjust = lab_h),
              size = 2.0, family = FONT_FAMILY, color = "grey25") +
    scale_fill_manual(values = c("Up in Adjacent"   = COL_ADJACENT,
                                  "Down in Adjacent" = COL_NORMAL),
                       guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    labs(title = "CellChat Pathways",
         x = "\u0394 Comm. Score (Adj \u2212 Nor)", y = NULL) +
    theme(plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
          axis.text.y  = element_text(size = 7),
          axis.text.x  = element_text(size = 7.5),
          axis.title.x = element_text(size = 8, face = "bold"),
          plot.margin  = margin(2, 3, 1, 2, "mm"))
  save_pdf("Fig5b_cell_communication.pdf", 53, 76, print(p4b))
  obj_B <- p4b   # Strategy A: pure ggplot, vector assembly
  grob_B <- NULL # disable old grob path
  cat("    -> saved (pathway-level diff bar plot)\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# ---------------------------------------------------------------------------
# Fig 5c: Ligand-receptor pair differential heatmap (per-pair centred)
# ---------------------------------------------------------------------------
cat("  Fig5c\n")
tryCatch({
  lr_diff <- read.csv(file.path(COMM_DIR, "lr_pair_differential.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  lr_diff <- lr_diff[order(lr_diff$pvalue), ]

  # Top 25 most significant L-R pairs
  top_lr <- head(lr_diff, 25)
  lr_mat <- as.matrix(top_lr[, c("mean_Nor", "mean_Adj")])
  colnames(lr_mat) <- c("Normal", "Adjacent")
  rownames(lr_mat) <- paste0(top_lr$ligand, " -> ", top_lr$receptor)

  # Per-pair centering: guarantees directional color consistency within each block.
  # For each L-R pair: dev = value − midpoint(Adjacent, Normal)
  #   → "UP in Adjacent" pair: Adjacent dev > 0 (red), Normal dev < 0 (blue), no exceptions
  #   → "UP in Normal"   pair: Normal dev > 0 (red), Adjacent dev < 0 (blue), no exceptions
  #   → Color saturation ∝ |Adjacent − Normal|: larger delta = more vivid contrast
  lr_mat_t  <- t(lr_mat)
  lr_mat_t  <- lr_mat_t[c("Adjacent", "Normal"), , drop = FALSE]
  pair_mid  <- colMeans(lr_mat_t)                        # per-pair midpoint
  lr_dev_t  <- sweep(lr_mat_t, 2, pair_mid, "-")         # deviation from midpoint
  dev_scale <- max(abs(lr_dev_t), na.rm = TRUE)          # global scale to [-1, +1]
  lr_dev_t  <- lr_dev_t / dev_scale
  lr_dev_t[is.na(lr_dev_t)] <- 0

  # ---- Data-driven column ordering for maximum visual contrast ----
  # delta = Adjacent dev − Normal dev = (Adjacent − Normal) / dev_scale  (same sign as raw delta)
  delta_z_c <- lr_dev_t["Adjacent", ] - lr_dev_t["Normal", ]
  # Three-key sort: direction block → category (same category together) → delta descending
  dir_vec_c <- ifelse(delta_z_c > 0, "UP in Adjacent", "UP in Normal")
  col_ord_c <- order(
    match(dir_vec_c, c("UP in Adjacent", "UP in Normal")),  # 1st: direction block
    top_lr$category,                                         # 2nd: category (alphabetical)
    -delta_z_c                                               # 3rd: delta descending within category
  )
  lr_dev_t  <- lr_dev_t[, col_ord_c, drop = FALSE]
  top_lr    <- top_lr[col_ord_c, ]
  delta_z_c <- delta_z_c[col_ord_c]

  # Split columns into two visually separated direction blocks
  dir_split_c <- factor(
    ifelse(delta_z_c > 0, "UP in Adjacent", "UP in Normal"),
    levels = c("UP in Adjacent", "UP in Normal"))

  # Fixed colour range [-1, +1] (full saturation = max raw delta in dataset)
  col_div_c <- colorRamp2(c(-1, 0, 1), c(COL_DOWN, "#FFFFFF", COL_UP))
  cat(sprintf("    C: delta [%.3f, %.3f]; dev_scale=%.4f; UP-Adj=%d UP-Nor=%d\n",
              min(delta_z_c), max(delta_z_c), dev_scale,
              sum(delta_z_c > 0), sum(delta_z_c <= 0)))

  # Column annotation: category + direction (built AFTER reordering top_lr)
  cat_col <- PAL_CAT[seq_len(length(unique(top_lr$category)))]
  names(cat_col) <- unique(top_lr$category)

  col_anno_c <- HeatmapAnnotation(
    Category  = top_lr$category,
    col = list(Category = cat_col),
    annotation_name_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
    annotation_name_side = "left",
    simple_anno_size = unit(2.5, "mm"),
    annotation_legend_param = list(
      Category = list(
        title_gp = gpar(fontsize = 7, fontfamily = FONT_GRID, fontface = "bold"),
        labels_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
        grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"),
        nrow = 4))
  )

  s <- sp("Fig5c_LR_heatmap.pdf", 130, 76)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  ht4c <- Heatmap(lr_dev_t, name = "Norm. Delta", col = col_div_c,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    row_names_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 6.5, fontfamily = FONT_GRID),
    column_names_rot = 45,
    column_names_side = "bottom",
    top_annotation = col_anno_c,
    rect_gp = gpar(col = "white", lwd = 0.3),
    height = unit(8, "mm"),
    column_split = dir_split_c,
    cluster_column_slices = FALSE,
    column_title = c("UP in Adjacent", "UP in Normal"),
    column_title_gp = gpar(fontsize = 7.5, fontface = "bold",
                            fontfamily = FONT_GRID,
                            col = c(COL_ADJACENT, COL_NORMAL)),
    column_gap = unit(2, "mm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 7, fontfamily = FONT_GRID, fontface = "bold"),
      labels_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
      legend_width = unit(20, "mm"), grid_height = unit(2.5, "mm"),
      direction = "horizontal"))
  draw(ht4c, padding = unit(c(5, 3, 8, 3), "mm"),
       column_title = "L-R Pair Differential Profile (top 25 by P)",
       column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_GRID),
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom",
       merge_legends = TRUE); dev.off()
  # Strategy B: capture ComplexHeatmap as grob for direct assembly
  grob_C <- grid::grid.grabExpr({
    draw(ht4c, padding = unit(c(5, 3, 8, 3), "mm"),
         column_title = "L-R Pair Differential Profile (top 25 by P)",
         column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_GRID),
         heatmap_legend_side = "bottom",
         annotation_legend_side = "bottom",
         merge_legends = TRUE)
  }, width = s$width, height = s$height)
  obj_C <- patchwork::wrap_elements(grob_C)
  cat("    -> saved\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# ---------------------------------------------------------------------------
# Fig 5d: Ligand-receptor lollipop plot (top 25 by |diff|)
# ---------------------------------------------------------------------------
cat("  Fig5d\n")
tryCatch({
  lr_diff <- read.csv(file.path(COMM_DIR, "lr_pair_differential.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  lr_diff <- lr_diff[order(-abs(lr_diff$diff)), ]
  top25 <- head(lr_diff, 25)

  # ---- Lollipop chart: directly visualise delta-score (Adjacent - Normal) ----
  # Sort rows by signed differential (descending): UP-in-Adjacent at top
  top25$diff_score <- top25$mean_Adj - top25$mean_Nor
  top25 <- top25[order(top25$diff_score, decreasing = TRUE), ]
  top25$LR_pair <- paste0(top25$ligand, " -> ", top25$receptor)
  top25$LR_pair <- factor(top25$LR_pair, levels = rev(top25$LR_pair))

  # Y-axis label colours: red = UP in Adjacent, blue = UP in Normal
  # factor levels are bottom-to-top = rev(sorted diff_score order)
  y_colors <- rev(ifelse(top25$diff_score > 0, COL_ADJACENT, COL_NORMAL))

  # Category colour palette
  cats     <- sort(unique(top25$category))
  cat_cols <- setNames(PAL_CAT[seq_along(cats)], cats)

  p4d <- ggplot(top25, aes(y = LR_pair, x = diff_score)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.4, color = "grey45") +
    geom_segment(aes(x = 0, xend = diff_score,
                     y = LR_pair, yend = LR_pair),
                 color = "grey78", linewidth = 0.5) +
    geom_point(aes(color = category, size = abs(diff_score)), alpha = 0.92) +
    scale_color_manual(values = cat_cols, name = "Category") +
    scale_size_continuous(range = c(1.2, 5),
                          name = "|\u0394 Score|",
                          guide = guide_legend(override.aes = list(alpha = 1))) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    labs(x = "\u0394 Comm. Score (Adjacent \u2212 Normal)",
         y = NULL, title = "Top L-R Pairs by \u0394 Score (top 25)") +
    theme(plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
          axis.text.y  = element_text(size = 7.5, color = y_colors, face = "plain"),
          axis.text.x  = element_text(size = 8),
          axis.title.x = element_text(size = 8.5, face = "bold"),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 7.5, face = "bold"),
          legend.key.size  = unit(3, "mm"),
          legend.spacing.y = unit(0.5, "mm"),
          plot.margin  = margin(3, 4, 1, 2, "mm"))
  save_pdf("Fig5d_LR_bubble.pdf", 90, 76, print(p4d))
  obj_D <- p4d
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# ---------------------------------------------------------------------------
# Fig 5e: TCR/BCR immune repertoire expression heatmap
# Enhancement: rows sorted by Adjacent-Normal delta-z for maximum visual contrast;
#   columns sorted within group by mean expression; tight adaptive colour range;
#   right annotation with delta-z barplot showing direction and magnitude.
# ---------------------------------------------------------------------------
cat("  Fig5e\n")
tryCatch({
  ig_tr <- read.csv(file.path(TCR_DIR, "ig_tr_family_expression.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)

  expr_mat <- as.matrix(ig_tr[, c("IGH", "IGK", "IGL", "TRB", "TRA", "TRG", "TRD")])
  rownames(expr_mat) <- ig_tr$sample
  expr_mat <- t(expr_mat)   # genes as rows, samples as columns

  # Row-wise z-score
  expr_z <- t(scale(t(expr_mat)))
  expr_z[is.na(expr_z)] <- 0
  expr_z <- clamp_matrix(expr_z, 2)

  # ---- STEP 1: Group membership (before any reordering) ----
  is_adj <- !grepl("Normal", colnames(expr_z))   # TRUE = Adjacent

  # ---- STEP 2: Per-gene delta-z (Adjacent mean - Normal mean) → row order ----
  delta_z <- rowMeans(expr_z[,  is_adj, drop = FALSE]) -
             rowMeans(expr_z[, !is_adj, drop = FALSE])
  row_ord <- order(delta_z, decreasing = TRUE)   # most UP-in-Adjacent at top

  # ---- STEP 3: Sort columns within each group by decreasing mean expression ----
  adj_srt <- which( is_adj)[order(-colMeans(expr_z[,  is_adj, drop = FALSE]))]
  nor_srt <- which(!is_adj)[order(-colMeans(expr_z[, !is_adj, drop = FALSE]))]
  col_ord <- c(adj_srt, nor_srt)

  # ---- STEP 4: Apply orderings ----
  expr_z_f  <- expr_z[row_ord, col_ord]
  dz_f      <- delta_z[row_ord]

  # ---- STEP 5: Tight adaptive colour range (2nd-98th pct) for max contrast ----
  z_lim     <- max(abs(quantile(expr_z_f, c(0.02, 0.98))), 1.5)
  col_div_e <- colorRamp2(c(-z_lim, 0, z_lim),
                           c(COL_DOWN, "#FFFFFF", COL_UP))
  cat(sprintf("    delta-z range: [%.2f, %.2f]; colour lim: +-%.2f\n",
              min(dz_f), max(dz_f), z_lim))

  # ---- Column split (condition) + thin colour bar (legend suppressed — titles are clear) ----
  cond_vec <- factor(
    ifelse(grepl("Normal", colnames(expr_z_f)), "Normal", "Adjacent"),
    levels = c("Adjacent", "Normal"))

  col_anno <- HeatmapAnnotation(
    Condition = cond_vec,
    col = list(Condition = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT)),
    annotation_name_gp  = gpar(fontsize = 7, fontfamily = FONT_GRID),
    annotation_name_side = "right",
    show_legend = FALSE,
    simple_anno_size = unit(2.5, "mm"))

  # ---- Row annotation: Adj-Nor delta-z barplot only ----
  row_anno <- rowAnnotation(
    "Adj-Nor" = anno_barplot(
      dz_f, baseline = 0,
      gp   = gpar(fill = ifelse(dz_f > 0, COL_ADJACENT, COL_NORMAL), col = NA),
      width = unit(14, "mm"),
      axis_param = list(gp = gpar(fontsize = 5.5, fontfamily = FONT_GRID))
    ),
    annotation_name_gp  = gpar(fontsize = 6.5, fontfamily = FONT_GRID),
    annotation_name_rot = 0
  )

  s <- sp("Fig5e_ig_tr_expression_heatmap.pdf", 93, 76)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  ht4e <- Heatmap(expr_z_f, name = "Z-score", col = col_div_e,
    cluster_rows    = FALSE,   # sorted by delta-z: UP-in-Adj at top
    cluster_columns = FALSE,   # sorted within group by mean expression
    show_row_dend = FALSE, show_column_dend = FALSE,
    row_names_gp    = gp_row_names(7, italic = TRUE),
    row_names_side  = "left",
    show_column_names = FALSE,
    top_annotation   = col_anno,
    right_annotation = row_anno,
    column_split    = cond_vec,
    cluster_column_slices = FALSE,
    column_title    = c("Adjacent", "Normal"),
    column_title_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    column_gap      = unit(2, "mm"),
    heatmap_legend_param = list(
      title_gp  = gpar(fontsize = 7, fontfamily = FONT_GRID, fontface = "bold"),
      labels_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
      legend_width = unit(20, "mm"), grid_height = unit(2.5, "mm"),
      direction = "horizontal"),
    row_names_max_width = unit(50, "mm"))
  draw(ht4e, padding = unit(c(1, 4, 6, 3), "mm"),
       column_title = "TCR/BCR Immune Repertoire Expression",
       column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_GRID),
       heatmap_legend_side = "bottom")
  dev.off()
  # Strategy B: capture ComplexHeatmap as grob for direct assembly
  grob_E <- grid::grid.grabExpr({
    draw(ht4e, padding = unit(c(1, 4, 6, 3), "mm"),
         column_title = "TCR/BCR Immune Repertoire Expression",
         column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_GRID),
         heatmap_legend_side = "bottom")
  }, width = s$width, height = s$height)
  obj_E <- patchwork::wrap_elements(grob_E)
  cat("    -> saved\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# ---------------------------------------------------------------------------
# Fig 5f: GRN master regulator network
# ---------------------------------------------------------------------------
cat("  Fig5f\n")
tryCatch({
  links <- read.csv(file.path(GRN_DIR, "grn_significant_links.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  nodes <- read.csv(file.path(GRN_DIR, "grn_node_metrics.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)

  # Top 20 master regulators by PageRank
  top_tfs <- nodes[nodes$IsTF == TRUE, ]
  top_tfs <- head(top_tfs[order(-top_tfs$PageRank), ], 20)

  # Filter links to include only edges from top TFs
  net_links <- links[links$TF %in% top_tfs$Gene, ]

  # Limit targets to top-weighted edges per TF (max 5 targets each)
  net_links <- net_links %>%
    group_by(TF) %>%
    slice_max(order_by = Weight, n = 5) %>%
    ungroup()

  # Build igraph
  g <- graph_from_data_frame(net_links[, c("TF", "Target", "Weight")],
                             directed = TRUE, vertices = NULL)

  # Node attributes — enhance TF visual distinction
  V(g)$is_tf <- V(g)$name %in% top_tfs$Gene
  V(g)$size <- ifelse(V(g)$is_tf, 5, 1.2)
  V(g)$color <- ifelse(V(g)$is_tf, COL_UP, COL_NS)
  V(g)$label <- ifelse(V(g)$is_tf, V(g)$name, "")

  tg <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(node_type = ifelse(is_tf, "TF", "Target"))

  set.seed(42)  # reproducible layout
  p4f <- ggraph(tg, layout = "fr") +
    geom_edge_arc(aes(alpha = Weight), strength = 0.1,
                  color = "grey70", linewidth = 0.2,
                  arrow = arrow(length = unit(0.8, "mm"), type = "closed")) +
    geom_node_point(aes(size = size, color = node_type),
                    alpha = 0.85, stroke = 0.3) +
    ggrepel::geom_text_repel(
      aes(x = x, y = y, label = label),
      size = 2.3, family = FONT_FAMILY, fontface = "bold.italic",
      color = "black",
      bg.color = "white", bg.r = 0.12,
      box.padding = 0.25, point.padding = 0.15,
      min.segment.length = 0, segment.size = 0.2, segment.color = "grey50",
      max.overlaps = Inf, force = 2, seed = 42
    ) +
    scale_color_manual(values = c("TF" = COL_UP, "Target" = COL_NS), name = "Node Type") +
    scale_size_identity() +
    scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
    labs(title = "GRN: Top 20 Master Regulators") +
    coord_cartesian(clip = "off") +
    theme_void(base_family = FONT_FAMILY) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                     family = FONT_FAMILY),
          legend.text = element_text(size = 8, family = FONT_FAMILY),
          legend.title = element_text(size = 8, face = "bold", family = FONT_FAMILY),
          plot.margin = margin(3, 4, 1, 4, "mm"))
  save_pdf("Fig5f_master_regulator_network.pdf", 78, 70, print(p4f))
  obj_F <- p4f  # Strategy A: keep ggraph/ggplot object for direct assembly
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# ---------------------------------------------------------------------------
# Fig 5g: TF centrality ranking bar plot
# ---------------------------------------------------------------------------
cat("  Fig5g\n")
tryCatch({
  nodes <- read.csv(file.path(GRN_DIR, "grn_node_metrics.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)

  tf_nodes <- nodes[nodes$IsTF == TRUE, ]
  tf_nodes <- head(tf_nodes[order(-tf_nodes$PageRank), ], 20)
  tf_nodes$Gene <- factor(tf_nodes$Gene, levels = rev(tf_nodes$Gene))

  p4g <- ggplot(tf_nodes, aes(x = Gene, y = PageRank)) +
    geom_col(aes(fill = TotalDegree), width = 0.75) +
    geom_text(aes(label = TotalDegree), hjust = -0.2, size = 2.5,
              family = FONT_FAMILY) +
    scale_fill_gradient(low = PAL_CAT[5], high = COL_UP, name = "Total Degree") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_flip() +
    labs(title = "TF Centrality Ranking (PageRank)",
         x = NULL, y = "PageRank Score") +
    theme(axis.text.y = element_text(face = "italic", size = 8),
          axis.text.x = element_text(size = 8),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, face = "bold"),
          plot.margin = margin(3, 4, 1, 4, "mm"))
  save_pdf("Fig5g_tf_centrality_ranking.pdf", 105, 70, print(p4g))
  obj_G <- p4g  # Strategy A: keep ggplot object for direct assembly
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================

# =============================================================================
# SECTION 11: Composite Figure Assembly (magick zero-distortion, 183x277mm)
# =============================================================================
cat("\n--- Assembling composite Figure_5 (magick zero-distortion, 183x277mm) ---\n")

tryCatch({
  W_TOTAL <- 183   # mm (J Hepatol double-column)
  H_TOTAL <- 277   # mm (55+76+76+70: H4 extended for F/G label spacing)
  DPI     <- 600
  px_per_mm <- DPI / 25.4

  # --- Row heights (must sum to 256) ---
  H1 <- 55   # Row 1: A (BayesPrism dumbbell plot, increased for cell-type row count)
  H2 <- 76   # Row 2: B + C (+10mm for Panel C legend clearance)
  H3 <- 76   # Row 3: D + E (+6mm for Panel D row annotation spacing)
  H4 <- H_TOTAL - H1 - H2 - H3   # Row 4: F + G (GRN + TF centrality) = 70

  # --- Column widths per row (each row must sum to 183) ---
  W_A <- W_TOTAL                              # Row 1: A full
  W_B <- 53; W_C <- W_TOTAL - W_B             # Row 2: 53+130=183
  W_D <- 90; W_E <- W_TOTAL - W_D             # Row 3: 90+93=183
  W_F <- 78; W_G <- W_TOTAL - W_F             # Row 4: 78+105=183

  # --- Pixel dimensions ---
  px_W <- round(W_TOTAL * px_per_mm)
  px_H <- round(H_TOTAL * px_per_mm)

  px_H1 <- round(H1 * px_per_mm)
  px_H2 <- round(H2 * px_per_mm)
  px_H3 <- round(H3 * px_per_mm)
  px_H4 <- px_H - px_H1 - px_H2 - px_H3   # remainder avoids rounding gaps

  px_WA <- px_W
  px_WB <- round(W_B * px_per_mm)
  px_WC <- px_W - px_WB
  px_WD <- round(W_D * px_per_mm)
  px_WE <- px_W - px_WD
  px_WF <- round(W_F * px_per_mm)
  px_WG <- px_W - px_WF

  cat(sprintf("  Layout: H1=%dmm, H2=%dmm, H3=%dmm, H4=%dmm\n", H1, H2, H3, H4))
  cat(sprintf("  Row2 W: B=%dmm, C=%dmm\n", W_B, W_C))
  cat(sprintf("  Row3 W: D=%dmm, E=%dmm\n", W_D, W_E))
  cat(sprintf("  Row4 W: F=%dmm, G=%dmm\n", W_F, W_G))
  cat(sprintf("  Canvas: %d x %d px @ %d DPI\n", px_W, px_H, DPI))

  # --- Read each panel PDF (rendered at exact target dimensions) ---
  read_panel <- function(fn) {
    img <- magick::image_read_pdf(file.path(OUT, fn), density = DPI)
    info <- magick::image_info(img)
    cat(sprintf("    %s: %d x %d px\n", fn, info$width, info$height))
    img
  }

  img_A <- read_panel("Fig5a_stacked_barplot.pdf")
  img_B <- read_panel("Fig5b_cell_communication.pdf")
  img_C <- read_panel("Fig5c_LR_heatmap.pdf")
  img_D <- read_panel("Fig5d_LR_bubble.pdf")
  img_E <- read_panel("Fig5e_ig_tr_expression_heatmap.pdf")
  img_F <- read_panel("Fig5f_master_regulator_network.pdf")
  img_G <- read_panel("Fig5g_tf_centrality_ranking.pdf")

  # --- Build white canvas and composite panels at exact pixel offsets ---
  canvas <- magick::image_blank(px_W, px_H, color = "white")

  y1 <- 0
  y2 <- px_H1
  y3 <- px_H1 + px_H2
  y4 <- px_H1 + px_H2 + px_H3

  # Row 1: A full width
  canvas <- magick::image_composite(canvas, img_A, offset = sprintf("+%d+%d", 0, y1))
  # Row 2: B + C
  canvas <- magick::image_composite(canvas, img_B, offset = sprintf("+%d+%d", 0, y2))
  canvas <- magick::image_composite(canvas, img_C, offset = sprintf("+%d+%d", px_WB, y2))
  # Row 3: D + E
  canvas <- magick::image_composite(canvas, img_D, offset = sprintf("+%d+%d", 0, y3))
  canvas <- magick::image_composite(canvas, img_E, offset = sprintf("+%d+%d", px_WD, y3))
  # Row 4: F + G
  canvas <- magick::image_composite(canvas, img_F, offset = sprintf("+%d+%d", 0, y4))
  canvas <- magick::image_composite(canvas, img_G, offset = sprintf("+%d+%d", px_WF, y4))

  info <- magick::image_info(canvas)
  cat(sprintf("  Composite: %d x %d px (expected %d x %d)\n",
              info$width, info$height, px_W, px_H))

  # --- Vector panel label overlay (FS_TAG=16 pt, bold) ---
  raster_img <- grDevices::as.raster(canvas)

  render_with_labels <- function() {
    grid::grid.raster(raster_img, interpolate = FALSE,
                      width = unit(W_TOTAL, "mm"),
                      height = unit(H_TOTAL, "mm"))
    # Label positions (mm from top-left)
    label_data <- data.frame(
      text = c("a", "b", "c", "d", "e", "f", "g"),
      x_mm = c(2,                        # a: row 1 left
               2,                        # b: row 2 col 1 (B)
               W_B + 2,                  # c: row 2 col 2 (C)
               2,                        # d: row 3 col 1 (D)
               W_D + 2,                  # e: row 3 col 2 (E)
               2,                        # f: row 4 col 1 (F)
               W_F + 2),                 # g: row 4 col 2 (G)
      y_mm = c(2,                        # a: top of page
               H1 + 2,                   # b: top of row 2
               H1 + 2,                   # c: top of row 2
               H1 + H2 + 2,              # d: top of row 3
               H1 + H2 + 2,              # e: top of row 3
               H1 + H2 + H3 + 2,         # f: top of row 4
               H1 + H2 + H3 + 2),        # g: top of row 4
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(label_data))) {
      grid::grid.text(
        label = label_data$text[i],
        x = grid::unit(label_data$x_mm[i], "mm"),
        y = grid::unit(H_TOTAL - label_data$y_mm[i], "mm"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = FS_TAG, fontface = "bold",
                        fontfamily = FONT_FAMILY)
      )
    }
  }

  # --- Output PDF ---
  cairo_pdf(file.path(OUT, "Figure_5.pdf"),
            width = W_TOTAL/25.4, height = H_TOTAL/25.4, family = FONT_FAMILY)
  render_with_labels()
  dev.off()
  cat("  -> Figure_5.pdf\n")

  # --- Output PNG (600 DPI, review) ---
  grDevices::png(file.path(OUT, "Figure_5.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_with_labels()
  dev.off()
  cat("  -> Figure_5.png\n")

  # --- Output TIFF (600 DPI, submission) ---
  grDevices::tiff(file.path(OUT, "Figure_5.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_with_labels()
  dev.off()
  cat("  -> Figure_5.tiff\n")

  cat("  Composite Figure_5 DONE (183x277mm, 600 DPI, magick zero-distortion)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

# =============================================================================
# SECTION 12: Pure Vector PDF Assembly for Adobe Illustrator
# Strategy: cairo_pdf + grid::viewport — each panel rendered directly into its
# target viewport; all text, lines and colour blocks remain as editable vectors.
# Grobs (B/C/E) were captured at their exact target dimensions, so grid.draw()
# fills the matching viewport without rescaling artefacts.
# =============================================================================
cat("\n--- Assembling Figure_5_AI.pdf (pure vector, Illustrator-editable) ---\n")

tryCatch({
  # Layout constants — must mirror SECTION 11
  W_TOTAL <- 183; H_TOTAL <- 277
  H1 <- 55; H2 <- 76; H3 <- 76; H4 <- H_TOTAL - H1 - H2 - H3   # = 70
  W_B <- 53;  W_C <- W_TOTAL - W_B    # 53 + 130 = 183
  W_D <- 90;  W_E <- W_TOTAL - W_D    # 90 +  93 = 183
  W_F <- 78;  W_G <- W_TOTAL - W_F    # 78 + 105 = 183

  # Helper: viewport anchored at (x, y) bottom-left corner, in mm
  vp_mm <- function(x, y, w, h)
    grid::viewport(x      = grid::unit(x, "mm"),
                   y      = grid::unit(y, "mm"),
                   width  = grid::unit(w, "mm"),
                   height = grid::unit(h, "mm"),
                   just   = c("left", "bottom"))

  # y from page bottom (origin = bottom-left)
  vp_A <- vp_mm(0,         H2+H3+H4, W_TOTAL, H1)   # row 1
  vp_B <- vp_mm(0,         H3+H4,    W_B,     H2)   # row 2 left
  vp_C <- vp_mm(W_B,       H3+H4,    W_C,     H2)   # row 2 right
  vp_D <- vp_mm(0,         H4,       W_D,     H3)   # row 3 left
  vp_E <- vp_mm(W_D,       H4,       W_E,     H3)   # row 3 right
  vp_F <- vp_mm(0,         0,        W_F,     H4)   # row 4 left
  vp_G <- vp_mm(W_F,       0,        W_G,     H4)   # row 4 right

  cairo_pdf(file.path(OUT, "Figure_5_AI.pdf"),
            width  = W_TOTAL / 25.4,
            height = H_TOTAL / 25.4,
            family = FONT_FAMILY)
  grid::grid.newpage()

  # --- ggplot panels (print scales automatically to viewport) ---
  if (!is.null(obj_A)) { print(obj_A, vp = vp_A); cat("  [AI] a\n") }
  if (!is.null(obj_B)) { print(obj_B, vp = vp_B); cat("  [AI] b\n") }
  if (!is.null(obj_D)) { print(obj_D, vp = vp_D); cat("  [AI] d\n") }
  if (!is.null(obj_F)) { print(obj_F, vp = vp_F); cat("  [AI] f\n") }
  if (!is.null(obj_G)) { print(obj_G, vp = vp_G); cat("  [AI] g\n") }

  # --- ComplexHeatmap grobs (captured at matching dimensions → no rescaling) ---
  draw_grob <- function(grob, vp, label) {
    if (exists(deparse(substitute(grob))) && !is.null(grob)) {
      grid::pushViewport(vp)
      grid::grid.draw(grob)
      grid::popViewport()
      cat(sprintf("  [AI] %s\n", label))
    }
  }
  if (exists("grob_B") && !is.null(grob_B)) {
    grid::pushViewport(vp_B); grid::grid.draw(grob_B); grid::popViewport()
    cat("  [AI] b\n")
  }
  if (exists("grob_C") && !is.null(grob_C)) {
    grid::pushViewport(vp_C); grid::grid.draw(grob_C); grid::popViewport()
    cat("  [AI] c\n")
  }
  if (exists("grob_E") && !is.null(grob_E)) {
    grid::pushViewport(vp_E); grid::grid.draw(grob_E); grid::popViewport()
    cat("  [AI] e\n")
  }

  # --- Vector panel labels a–g (y_top = distance from page top) ---
  lbl <- data.frame(
    lab  = c("a", "b",    "c",    "d",       "e",       "f",          "g"),
    x    = c(2,   2,      W_B+2,  2,         W_D+2,     2,            W_F+2),
    ytop = c(2,   H1+2,   H1+2,   H1+H2+2,   H1+H2+2,   H1+H2+H3+2,   H1+H2+H3+2),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(lbl))) {
    grid::grid.text(
      label = lbl$lab[i],
      x     = grid::unit(lbl$x[i],              "mm"),
      y     = grid::unit(H_TOTAL - lbl$ytop[i], "mm"),
      just  = c("left", "top"),
      gp    = grid::gpar(fontsize   = FS_TAG,
                         fontface   = "bold",
                         fontfamily = FONT_FAMILY)
    )
  }

  dev.off()
  cat("  -> Figure_5_AI.pdf (pure vector, Adobe Illustrator editable)\n")

}, error = function(e) cat("  ERROR AI assembly:", e$message, "\n"))

# =============================================================================
cat("\n=== Figure 5 rendering complete ===\n")