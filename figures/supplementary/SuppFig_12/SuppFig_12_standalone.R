#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_12_standalone.R
# Supplementary Figure 12: Extended Biomarker Validation & Paired-Difference Approach
# HAE Multi-omics Study | Target: EBioMedicine (Lancet family)
# =============================================================================
# Panels (8):
#   A = Cross-dataset direction consistency heatmap (top consensus DEGs)
#   B = Concordance threshold curve (proportion of genes at each agreement level)
#   C = Meta-analysis cohort characteristics bubble plot
#   D = Signature validation AUC across cohorts (lollipop + CI)
#   E = Paired-difference LOPO-CV top 15 features (lollipop + CI)
#   F = Paired delta heatmap (top 20 discriminative features)
#   G = Feature-level AUC distribution by omics layer (boxplot + jitter)
#   H = Bootstrap LASSO stability selection (lollipop)
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_12_standalone.R
# =============================================================================

cat("=== Supplementary Figure 12: Extended Biomarker Validation ===\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(ggsci)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "02_analysis/results")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_12")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM <- 25.4

FS_GEOM <- 2.82   # geom_text size (≈8pt)
FS_TAG  <- 12     # panel label font size (bold, top-left)

COL_UP   <- "#CD534CFF"; COL_DOWN <- "#0073C2FF"; COL_NS <- "#868686FF"
COL_TC   <- "#0073C2FF"; COL_PR   <- "#CD534CFF"; COL_MT <- "#EFC000FF"
PAL_CAT  <- pal_jco("default")(10)

# =============================================================================
# Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text             = element_text(family = FONT_FAMILY, size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks       = element_line(linewidth = 0.3, color = "black"),
    axis.text        = element_text(size = 8, color = "black", family = FONT_FAMILY),
    axis.title       = element_text(size = 9, face = "bold", family = FONT_FAMILY),
    plot.title       = element_text(size = 10, face = "bold", hjust = 0,
                                    family = FONT_FAMILY),
    legend.text      = element_text(size = 8, family = FONT_FAMILY),
    legend.title     = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    legend.key.size  = unit(3, "mm"),
    legend.background = element_blank(),
    strip.text       = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    plot.margin      = margin(2, 3, 2, 3, "mm")
  )
theme_set(theme_nc)

gp_rn <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID)
gp_cn <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
ht_opt$message <- FALSE

save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fp <- file.path(OUT, filename)
  cairo_pdf(fp, width = w_mm / MM, height = h_mm / MM, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat("  ->", filename, sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

mm2in <- function(mm) mm / MM

# =============================================================================
# PANEL A: Cross-dataset direction consistency heatmap
# =============================================================================
cat("\n--- Panel A: Direction consistency heatmap ---\n")
ht_A <- NULL
tryCatch({
  meta <- read.csv(file.path(RES, "enhancement14_meta_analysis/meta_analysis_consensus_DEGs.csv"),
                   stringsAsFactors = FALSE)
  lfc_cols <- grep("^logFC_", colnames(meta), value = TRUE)
  datasets <- gsub("^logFC_", "", lfc_cols)

  # Select top 18 genes by n_valid then |mean_logFC| (reduced to avoid label overlap
  # within composite-figure row height of 62mm)
  meta$n_valid <- rowSums(!is.na(meta[, lfc_cols]) & meta[, lfc_cols] != "NA")
  meta_top <- head(meta[order(-meta$n_valid, -abs(as.numeric(meta$mean_logFC))), ], 18)

  mat <- as.matrix(meta_top[, lfc_cols])
  mat[mat == "NA"] <- NA
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- meta_top$gene
  colnames(mat) <- datasets
  mat[mat >  3] <-  3
  mat[mat < -3] <- -3

  col_fun <- colorRamp2(c(-3, 0, 3), c(COL_DOWN, "white", COL_UP))

  save_panel_pdf("Supp12a_direction_heatmap.pdf", 120, 110, {
    ht <- Heatmap(mat,
                  name          = "logFC",
                  col           = col_fun,
                  cluster_rows    = TRUE,
                  cluster_columns = FALSE,
                  show_row_dend   = FALSE,
                  na_col          = "#F0F0F0",
                  row_names_gp    = gp_rn(7),
                  column_names_gp = gp_cn(7),
                  column_names_rot = 45,
                  column_names_side = "bottom",
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    v <- mat[i, j]
                    if (!is.na(v))
                      grid.text(sprintf("%.1f", v), x, y,
                                gp = gpar(fontsize = 6, fontfamily = FONT_GRID,
                                          col = ifelse(abs(v) > 1.8, "white", "black")))
                  },
                  column_title    = "Cross-dataset expression direction",
                  column_title_gp = gpar(fontsize = 10, fontfamily = FONT_GRID,
                                         fontface = "bold"),
                  heatmap_legend_param = list(
                    title_gp  = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
                    labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID)
                  ))
    draw(ht)
    ht_A <<- ht
  })
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL B: Concordance threshold curve (replaces bar chart)
# =============================================================================
cat("\n--- Panel B: Concordance threshold curve ---\n")
p_B <- NULL
tryCatch({
  meta <- read.csv(file.path(RES, "enhancement14_meta_analysis/meta_analysis_consensus_DEGs.csv"),
                   stringsAsFactors = FALSE)
  agree      <- as.numeric(meta$direction_agree)
  thresholds <- c(0.67, 0.75, 0.80, 0.83, 1.00)
  props      <- sapply(thresholds, function(th) mean(agree >= th, na.rm = TRUE))
  df_line <- data.frame(
    Threshold  = thresholds * 100,
    Proportion = props * 100,
    Label      = sprintf("%.1f%%", props * 100)
  )

  p_B <- ggplot(df_line, aes(x = Threshold, y = Proportion)) +
    geom_area(fill = PAL_CAT[1], alpha = 0.15) +
    geom_line(color = PAL_CAT[1], linewidth = 0.8) +
    geom_point(color = PAL_CAT[1], fill = "white", shape = 21,
               size = 3, stroke = 1.2) +
    geom_text(aes(label = Label), vjust = -0.6, size = FS_GEOM,
              family = FONT_FAMILY, color = "black") +
    scale_x_continuous(breaks = thresholds * 100,
                       labels = paste0(thresholds * 100, "%")) +
    scale_y_continuous(limits = c(0, 115),
                       breaks = seq(0, 100, 25),
                       labels = paste0(seq(0, 100, 25), "%")) +
    labs(title  = "Directional concordance across cohorts",
         x      = "Agreement threshold",
         y      = "Proportion of genes") +
    theme(axis.text.x = element_text(size = 8))
  save_panel_pdf("Supp12b_concordance.pdf", 89, 75, print(p_B))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL C: Dataset characteristics bubble
# =============================================================================
cat("\n--- Panel C: Dataset characteristics ---\n")
p_C <- NULL
tryCatch({
  ds <- read.csv(file.path(RES, "enhancement14_meta_analysis/effective_sample_summary.csv"),
                 stringsAsFactors = FALSE)
  ds$species_label <- ifelse(ds$species == "human", "Human", "Mouse")

  p_C <- ggplot(ds, aes(x = reorder(dataset, -n_total), y = n_total,
                         size = n_total, color = species_label)) +
    geom_point(alpha = 0.75) +
    scale_size_continuous(range = c(4, 14), name = "N samples") +
    scale_color_manual(values = c("Human" = COL_TC, "Mouse" = COL_MT),
                       name = "Species") +
    geom_text(aes(label = n_total), size = FS_GEOM, color = "black",
              family = FONT_FAMILY) +
    labs(title = "Meta-analysis cohort overview (130 samples)",
         x = NULL, y = "Sample size") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  save_panel_pdf("Supp12c_dataset_bubble.pdf", 89, 80, print(p_C))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL D: Signature validation AUC – lollipop + CI
#          Below-chance cohorts (AUC < 0.5) shown as open circles
# =============================================================================
cat("\n--- Panel D: Signature validation AUC (lollipop) ---\n")
p_D <- NULL
tryCatch({
  val <- read.csv(file.path(RES, "enhancement14_meta_analysis/signature_validation_AUC.csv"),
                  stringsAsFactors = FALSE)
  val$dataset    <- factor(val$dataset, levels = val$dataset[order(val$AUC)])
  val$species_col <- ifelse(val$species == "human", "Human", "Mouse")
  val$above_chance <- val$AUC >= 0.5

  p_D <- ggplot(val, aes(y = dataset, x = AUC)) +
    # stem from AUC = 0.5 (reference, not 0)
    geom_segment(aes(x = 0.5, xend = AUC, yend = dataset,
                     color = species_col),
                 linewidth = 0.7) +
    # CI error bar
    geom_errorbarh(aes(xmin = AUC_CI_lower, xmax = pmin(AUC_CI_upper, 1),
                       color = species_col),
                   height = 0.2, linewidth = 0.4) +
    # filled dot = above chance; open dot = below chance
    geom_point(data = subset(val,  above_chance),
               aes(color = species_col), size = 3.2) +
    geom_point(data = subset(val, !above_chance),
               aes(color = species_col), shape = 21, fill = "white",
               size = 3.2, stroke = 1.1) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    scale_color_manual(values = c("Human" = COL_TC, "Mouse" = COL_MT),
                       name = "Species") +
    scale_x_continuous(limits = c(0, 1.08),
                       breaks  = seq(0, 1, 0.25)) +
    labs(title = "Signature validation across cohorts",
         x = "AUC (95% CI)", y = NULL) +
    theme(legend.position = "top")
  save_panel_pdf("Supp12d_validation_AUC.pdf", 89, 80, print(p_D))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL E: Paired-difference LOPO-CV top features – lollipop + CI
# =============================================================================
cat("\n--- Panel E: Paired LOPO-CV top features (lollipop) ---\n")
p_E <- NULL
tryCatch({
  lopo     <- read.csv(file.path(RES, "optimization_biomarker/single_feature_lopo_cv.csv"),
                       stringsAsFactors = FALSE)
  lopo$lopo_auc <- as.numeric(lopo$lopo_auc)
  lopo_top <- head(lopo[order(-lopo$lopo_auc), ], 15)
  lopo_top$feature <- factor(lopo_top$feature, levels = rev(lopo_top$feature))

  omics_cols  <- c("TC"  = COL_TC, "PR"  = COL_PR, "MET" = COL_MT)
  omics_labs  <- c("TC"  = "Transcriptomics",
                   "PR"  = "Proteomics",
                   "MET" = "Metabolomics")

  p_E <- ggplot(lopo_top, aes(x = lopo_auc, y = feature, color = omics)) +
    geom_segment(aes(x = 0.5, xend = lopo_auc, yend = feature),
                 color = "grey78", linewidth = 0.55) +
    geom_errorbarh(aes(xmin = as.numeric(ci_low), xmax = as.numeric(ci_high)),
                   height = 0.28, linewidth = 0.4) +
    geom_point(size = 3.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.35) +
    scale_color_manual(values = omics_cols, labels = omics_labs,
                       name = "Omics") +
    scale_x_continuous(limits = c(0.4, 1.05),
                       breaks = seq(0.4, 1.0, 0.1)) +
    labs(title = "Paired-difference LOPO-CV AUC",
         x = "AUC (95% CI)", y = NULL) +
    theme(legend.position = "top")
  save_panel_pdf("Supp12e_paired_LOPO.pdf", 89, 95, print(p_E))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL F: Paired delta heatmap (top 20 features)
# =============================================================================
cat("\n--- Panel F: Paired delta heatmap ---\n")
ht_F <- NULL
tryCatch({
  delta     <- read.csv(file.path(RES, "optimization_biomarker/paired_delta_wilcoxon_results.csv"),
                        stringsAsFactors = FALSE)
  delta_top <- head(delta[order(delta$pvalue), ], 20)
  # Re-order rows by omics layer (TC -> PR -> MET) while preserving p-value rank within
  delta_top$omics <- factor(delta_top$omics, levels = c("TC", "PR", "MET"))
  delta_top <- delta_top[order(delta_top$omics, delta_top$pvalue), ]

  mat_delta <- matrix(as.numeric(delta_top$effect_size), ncol = 1)
  rownames(mat_delta) <- delta_top$feature
  colnames(mat_delta) <- "Effect size"
  mat_delta[mat_delta >  3] <-  3
  mat_delta[mat_delta < -3] <- -3

  col_fun   <- colorRamp2(c(0, 1.5, 3), c("white", "#EFC000FF", COL_UP))
  omics_col <- c("TC" = COL_TC, "PR" = COL_PR, "MET" = COL_MT)

  save_panel_pdf("Supp12f_delta_heatmap.pdf", 70, 100, {
    ra <- rowAnnotation(
      Omics = as.character(delta_top$omics),
      col   = list(Omics = omics_col),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        title_gp  = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
        labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID)
      )
    )
    ht <- Heatmap(mat_delta,
                  name             = "Effect\nsize",
                  col              = col_fun,
                  cluster_rows     = FALSE,
                  cluster_columns  = FALSE,
                  row_names_gp     = gp_rn(8),
                  show_column_names = FALSE,
                  row_split        = as.character(delta_top$omics),
                  row_title        = NULL,
                  row_gap          = unit(1.5, "mm"),
                  left_annotation  = ra,
                  column_title     = "Top paired-difference features",
                  column_title_gp  = gpar(fontsize = 10, fontfamily = FONT_GRID,
                                          fontface = "bold"),
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    grid.text(sprintf("%.2f", mat_delta[i, j]), x, y,
                              gp = gpar(fontsize = 8, fontfamily = FONT_GRID))
                  },
                  heatmap_legend_param = list(
                    title_gp  = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
                    labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID)
                  ))
    draw(ht)
    ht_F <<- ht
  })
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL G: Feature-level AUC distribution by omics layer
# =============================================================================
cat("\n--- Panel G: AUC by omics layer ---\n")
p_G <- NULL
tryCatch({
  lopo      <- read.csv(file.path(RES, "optimization_biomarker/single_feature_lopo_cv.csv"),
                        stringsAsFactors = FALSE)
  lopo$lopo_auc <- as.numeric(lopo$lopo_auc)
  lopo$omics_label <- c("TC"  = "Transcriptomics",
                         "PR"  = "Proteomics",
                         "MET" = "Metabolomics")[lopo$omics]
  omics_cols_full <- c("Transcriptomics" = COL_TC,
                        "Proteomics"      = COL_PR,
                        "Metabolomics"    = COL_MT)

  p_G <- ggplot(lopo, aes(x = omics_label, y = lopo_auc, fill = omics_label)) +
    geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 0.8) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed",
               color = COL_UP, linewidth = 0.35) +
    annotate("text", x = 3.4, y = 0.815, label = "AUC = 0.8",
             size = FS_GEOM, color = COL_UP, hjust = 1, family = FONT_FAMILY) +
    scale_fill_manual(values = omics_cols_full) +
    scale_y_continuous(limits = c(0.55, 1.0),
                       breaks = seq(0.6, 1.0, 0.1)) +
    labs(title = "Feature-level AUC distribution by omics",
         x = NULL, y = "LOPO-CV AUC") +
    theme(legend.position = "none")
  save_panel_pdf("Supp12g_auc_by_omics.pdf", 89, 75, print(p_G))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL H: Stability selection frequency – lollipop
# =============================================================================
cat("\n--- Panel H: Stability selection (lollipop) ---\n")
p_H <- NULL
tryCatch({
  stab <- read.csv(file.path(RES, "optimization_biomarker/stability_selection_subtype.csv"),
                   stringsAsFactors = FALSE)
  stab$selection_frequency <- as.numeric(stab$selection_frequency)
  stab <- stab[order(-stab$selection_frequency), ]
  stab$feature <- factor(stab$feature, levels = rev(stab$feature))
  omics_cols <- c("TC" = COL_TC, "PR" = COL_PR, "MET" = COL_MT)
  omics_labs <- c("TC" = "Transcriptomics", "PR" = "Proteomics",
                  "MET" = "Metabolomics")

  p_H <- ggplot(stab, aes(x = selection_frequency, y = feature, color = omics)) +
    geom_segment(aes(x = 0, xend = selection_frequency, yend = feature),
                 color = "grey75", linewidth = 0.8) +
    geom_point(size = 4.5) +
    geom_text(aes(label = paste0(selection_frequency, "%")),
              hjust = -0.4, size = FS_GEOM, family = FONT_FAMILY, color = "black") +
    geom_vline(xintercept = 50, linetype = "dashed",
               color = COL_UP, linewidth = 0.45) +
    annotate("text", x = 52, y = 0.55, label = "50% threshold",
             size = FS_GEOM, color = COL_UP, hjust = 0, family = FONT_FAMILY) +
    scale_color_manual(values = omics_cols, labels = omics_labs,
                       name = "Omics") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(title = "Bootstrap LASSO stability selection",
         x = "Selection frequency (%)", y = NULL) +
    theme(legend.position = "top",
          axis.text.y = element_text(face = "bold"))
  save_panel_pdf("Supp12h_stability_freq.pdf", 89, 65, print(p_H))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# Composite Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite SuppFig_12 (vector, 183x245mm) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  H1 <- 62; H2 <- 62; H3 <- 62; H4 <- H_TOTAL - H1 - H2 - H3  # = 59mm
  W_L <- 91; W_R <- W_TOTAL - W_L                               # = 92mm

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width  = unit(W_TOTAL, "mm"),
      height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL),
      yscale = c(0, H_TOTAL)
    ))

    # Row 1: A (heatmap, left) | B (line+area, right)
    grid::pushViewport(grid::viewport(
      x=unit(0,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    if (!is.null(ht_A)) draw(ht_A, newpage=FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x=unit(W_L,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    if (!is.null(p_B)) print(p_B, newpage=FALSE)
    grid::popViewport()

    # Row 2: C (bubble, left) | D (lollipop AUC, right)
    grid::pushViewport(grid::viewport(
      x=unit(0,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    if (!is.null(p_C)) print(p_C, newpage=FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x=unit(W_L,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    if (!is.null(p_D)) print(p_D, newpage=FALSE)
    grid::popViewport()

    # Row 3: E (lollipop LOPO, left) | F (delta heatmap, right)
    grid::pushViewport(grid::viewport(
      x=unit(0,"mm"), y=unit(H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    if (!is.null(p_E)) print(p_E, newpage=FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x=unit(W_L,"mm"), y=unit(H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    if (!is.null(ht_F)) draw(ht_F, newpage=FALSE)
    grid::popViewport()

    # Row 4: G (boxplot, left) | H (lollipop stability, right)
    grid::pushViewport(grid::viewport(
      x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_L,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    if (!is.null(p_G)) print(p_G, newpage=FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x=unit(W_L,"mm"), y=unit(0,"mm"),
      width=unit(W_R,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    if (!is.null(p_H)) print(p_H, newpage=FALSE)
    grid::popViewport()

    # Panel labels a–h
    label_data <- data.frame(
      text  = c("A", "B", "C", "D", "E", "F", "G", "H"),
      x_mm  = c(2, W_L+2, 2, W_L+2, 2, W_L+2, 2, W_L+2),
      y_mm  = c(H2+H3+H4+H1-2, H2+H3+H4+H1-2,
                H3+H4+H2-2,    H3+H4+H2-2,
                H4+H3-2,       H4+H3-2,
                H4-2,          H4-2),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(label_data))) {
      grid::grid.text(
        label = label_data$text[i],
        x     = unit(label_data$x_mm[i], "mm"),
        y     = unit(label_data$y_mm[i], "mm"),
        just  = c("left", "top"),
        gp    = grid::gpar(fontsize = FS_TAG, fontface = "bold",
                           fontfamily = FONT_FAMILY)
      )
    }
    grid::popViewport()
  }

  # PDF (vector)
  cairo_pdf(file.path(OUT, "SuppFig_12.pdf"),
            width  = W_TOTAL / 25.4,
            height = H_TOTAL / 25.4,
            family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> SuppFig_12.pdf (vector)\n")

  # PNG (preview)
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "SuppFig_12.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> SuppFig_12.png\n")

  # TIFF (submission)
  grDevices::tiff(file.path(OUT, "SuppFig_12.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> SuppFig_12.tiff\n")

  cat("  SuppFig_12 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 12 rendering complete ===\n")
