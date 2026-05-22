#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_08_standalone.R  (v2 zero-distortion polish: 2026-05-21)
# Supplementary Figure 8: Regulatory architecture, ligand–receptor crosstalk,
#   and deconvolution validation supporting Fig. 5
# HAE Multi-omics Study | Nature Communications
# =============================================================================
# Polish v2 fixes:
#   - Panel A: cluster_rows=FALSE (no dendrogram); M0-M13 ordered numerically;
#              full-width 183mm × 55mm to eliminate row-label overlap.
#   - Panel B: bold.italic gene labels at 3.2pt with white halo for clarity.
#   - Panel D: 89×70mm (20 TFs at 3.5mm/row, axis 7pt) — no overlap.
#   - Panel F: 89×70mm (20 pathways at 3.5mm/row, axis 7pt) — no overlap.
#   - Panel H: full-width 183×30mm (9 cells at ~3.3mm/row) — eliminates
#              bottom-right whitespace; more bar runway for 0.97-0.99 values.
#   - Panels saved at EXACT viewport target dimensions → zero distortion.
# Panels (a–h):
#   a = WGCNA transcriptomic module-trait correlation heatmap (M0–M13 vs 15 traits)
#   b = WGCNA hub gene scatter (|kME| vs |GS|), top-10 hubs labelled
#   c = Condition-specific GRN edge comparison (Peri-lesional vs Normal vs Shared)
#   d = Top-20 TF regulatory degree change (Adjacent − Normal)
#   e = Ligand–receptor pair counts by pathway category (10 categories)
#   f = Pathway-level mean L–R communication score (enhanced vs reduced)
#   g = BayesPrism deconvolution cross-validation (Spearman vs ssGSEA)
#   h = Cross-method deconvolution sensitivity (mean concordance, 9 cell types)
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_08_standalone.R
# =============================================================================

cat("=== Supplementary Figure 8: Regulatory architecture & deconvolution validation ===\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggsci)
  library(ggpubr)
  library(ggrepel)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "02_analysis/results")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_08")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FS_TAG <- 12; FONT_GRID <- "Arial"; MM <- 25.4

# Final composite layout (target dimensions for each panel)
W_TOTAL <- 183; H_TOTAL <- 275
W_LEFT  <- 89;  W_RIGHT <- 94
HA <- 55; HBC <- 50; HDE <- 70; HFG <- 70; HH <- 30

COL_UP <- "#CD534CFF"; COL_DOWN <- "#0073C2FF"; COL_NS <- "#868686FF"
COL_NORMAL <- "#7AA6DCFF"; COL_ADJACENT <- "#CD534CFF"
PAL_CAT <- pal_jco("default")(10)

col_div <- colorRamp2(c(-2, 0, 2), c(COL_DOWN, "#FFFFFF", COL_UP))
col_cor <- colorRamp2(c(-1, 0, 1), c(COL_DOWN, "#FFFFFF", COL_UP))

# =============================================================================
# Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(2.5, "mm"),
    legend.background = element_blank(),
    plot.margin = margin(3, 3, 3, 3, "mm")
  )
theme_set(theme_nc)
ht_opt$message <- FALSE

# =============================================================================
# Helpers
# =============================================================================
save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  cairo_pdf(fpath, width = w_mm/MM, height = h_mm/MM, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat("  ->", basename(fpath), "(", w_mm, "x", h_mm, "mm)\n")
}

gp_rn <- function(sz = 7.5) gpar(fontsize = sz, fontfamily = FONT_GRID)
gp_cn <- function(sz = 7.5) gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
gp_lt <- function(sz = 8)   gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
gp_ll <- function(sz = 7.5) gpar(fontsize = sz, fontfamily = FONT_GRID)
std_lp <- function() list(title_gp = gp_lt(), labels_gp = gp_ll(),
                          legend_height = unit(20, "mm"), grid_width = unit(3, "mm"))

pretty_cat <- function(x) gsub("_", " ", x)

# =============================================================================
# Panel A: WGCNA module-trait correlation heatmap (no dendrogram, full-width)
# =============================================================================
cat("  Panel A: WGCNA module-trait\n")
ht_A <- NULL
tryCatch({
  cor_mat <- read.csv(file.path(RES, "enhancement2_wgcna/TC_module_trait_cor.csv"),
                      stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  pval_mat <- read.csv(file.path(RES, "enhancement2_wgcna/TC_module_trait_pval_BH.csv"),
                       stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  cor_m <- as.matrix(cor_mat)
  pval_m <- as.matrix(pval_mat)

  sig_stars <- ifelse(pval_m < 0.001, "***",
               ifelse(pval_m < 0.01,  "**",
               ifelse(pval_m < 0.05,  "*", "")))

  # Order modules numerically: M0, M1, ..., M13
  mod_order <- order(as.numeric(rownames(cor_m)))
  cor_m <- cor_m[mod_order, , drop = FALSE]
  pval_m <- pval_m[mod_order, , drop = FALSE]
  sig_stars <- sig_stars[mod_order, , drop = FALSE]
  rownames(cor_m) <- paste0("M", rownames(cor_m))

  # Save Panel A at EXACT composite target (183 × 55 mm)
  s <- file.path(OUT, "Supp08a_WGCNA_module_trait.pdf")
  cairo_pdf(s, width = W_TOTAL/MM, height = HA/MM, family = FONT_FAMILY)
  ht <- Heatmap(cor_m, name = "Correlation", col = col_cor,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    row_names_gp = gp_rn(7), row_names_side = "left",
    column_names_gp = gp_rn(7), column_names_side = "bottom",
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (nchar(sig_stars[i, j]) > 0)
        grid.text(sig_stars[i, j], x, y, gp = gpar(fontsize = 7, fontfamily = FONT_GRID))
    },
    column_title = "Module-trait correlation (WGCNA)",
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
    column_title_side = "top",
    heatmap_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll(),
                                 legend_height = unit(15, "mm"),
                                 grid_width = unit(2.5, "mm")),
    row_names_max_width = unit(15, "mm"),
    width = unit(W_TOTAL - 40, "mm"),
    height = unit(HA - 22, "mm"))
  draw(ht, padding = unit(c(2, 4, 2, 4), "mm")); dev.off()
  ht_A <<- ht
  cat("  -> Supp08a_WGCNA_module_trait.pdf ( 183 x 55 mm )\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Panel B: WGCNA hub gene scatter (bold.italic + halo for clarity)
# =============================================================================
cat("  Panel B: WGCNA hub scatter\n")
p_B <- NULL
tryCatch({
  hubs <- read.csv(file.path(RES, "enhancement2_wgcna/TC_hub_genes.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  hubs$is_hub <- hubs$hub_score >= quantile(hubs$hub_score, 0.9)
  hubs_ann <- hubs[!grepl("^ENSG", hubs$gene), ]
  top_genes <- head(hubs_ann[order(-hubs_ann$hub_score), ], 10)

  p <- ggplot(hubs, aes(x = abs_kME, y = abs_GS)) +
    geom_point(aes(color = module), alpha = 0.5, size = 1) +
    geom_point(data = top_genes, aes(x = abs_kME, y = abs_GS),
               color = COL_UP, size = 2.2, stroke = 0.5) +
    geom_text_repel(data = top_genes, aes(x = abs_kME, y = abs_GS, label = gene),
                    size = 3.2, family = FONT_FAMILY, fontface = "bold.italic",
                    bg.color = "white", bg.r = 0.18,
                    box.padding = 0.35, point.padding = 0.25,
                    max.overlaps = 20, segment.size = 0.25,
                    segment.color = "grey40", min.segment.length = 0) +
    scale_color_manual(values = setNames(PAL_CAT[seq_len(length(unique(hubs$module)))],
                                         unique(hubs$module))) +
    labs(x = "|kME| (Module Membership)", y = "|GS| (Gene Significance)",
         title = "WGCNA hub gene identification", color = "Module") +
    theme(legend.position = "right",
          legend.key.size = unit(2.5, "mm"))
  save_panel_pdf("Supp08b_WGCNA_hub_scatter.pdf", W_LEFT, HBC, print(p))
  p_B <<- p
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Panel C: GRN condition-specific edge barplot
# =============================================================================
cat("  Panel C: GRN edge comparison\n")
p_C <- NULL
tryCatch({
  comp <- read.csv(file.path(RES, "enhancement15_grn/grn_condition_comparison.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  adj_only <- comp$Value[comp$Metric == "Adjacent_Only"]
  nor_only <- comp$Value[comp$Metric == "Normal_Only"]
  shared   <- comp$Value[comp$Metric == "Shared"]

  df <- data.frame(
    Category = c("Peri-lesional\nonly", "Normal\nonly", "Shared"),
    Edges = c(adj_only, nor_only, shared),
    stringsAsFactors = FALSE)
  df$Category <- factor(df$Category, levels = df$Category)

  p <- ggplot(df, aes(x = Category, y = Edges, fill = Category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = format(Edges, big.mark = ",")),
              vjust = -0.4, size = 3, fontface = "bold", family = FONT_FAMILY) +
    scale_fill_manual(values = c(COL_ADJACENT, COL_NORMAL, COL_NS)) +
    labs(x = NULL, y = "Number of regulatory edges",
         title = "Condition-specific regulatory rewiring") +
    theme(legend.position = "none") +
    ylim(0, max(df$Edges) * 1.18)
  save_panel_pdf("Supp08c_GRN_edge_comparison.pdf", W_RIGHT, HBC, print(p))
  p_C <<- p
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Panel D: Top-20 TF differential activity (89×70mm, 7pt italic)
# =============================================================================
cat("  Panel D: Differential TF activity\n")
p_D <- NULL
tryCatch({
  tf_diff <- read.csv(file.path(RES, "enhancement15_grn/grn_tf_differential_activity.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  tf_diff <- tf_diff[order(-abs(tf_diff$Diff)), ]
  top20 <- head(tf_diff, 20)
  top20$TF <- factor(top20$TF, levels = rev(top20$TF))
  top20$direction <- ifelse(top20$Diff > 0, "Gained", "Lost")

  p <- ggplot(top20, aes(x = TF, y = Diff, fill = direction)) +
    geom_col(width = 0.78) +
    geom_text(aes(label = Diff),
              hjust = ifelse(top20$Diff > 0, -0.15, 1.15),
              size = 2.4, family = FONT_FAMILY) +
    scale_fill_manual(values = c("Gained" = COL_UP, "Lost" = COL_DOWN),
                      name = "Edge change") +
    coord_flip() +
    labs(x = NULL, y = "Degree difference (Peri-lesional − Normal)",
         title = "TF regulatory degree change") +
    theme(axis.text.y = element_text(face = "italic", size = 7),
          axis.text.x = element_text(size = 7),
          legend.position = c(0.85, 0.20),
          legend.background = element_rect(fill = "white", color = "grey80",
                                           linewidth = 0.2))
  save_panel_pdf("Supp08d_TF_degree_diff.pdf", W_LEFT, HDE, print(p))
  p_D <<- p
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Panel E: LR category barplot (94×70mm)
# =============================================================================
cat("  Panel E: LR category barplot\n")
p_E <- NULL
tryCatch({
  lr_cat <- read.csv(file.path(RES, "enhancement11_cell_communication/lr_category_summary.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  cat_col <- colnames(lr_cat)[1]
  n_col <- intersect(c("n_pairs", "count", "n"), colnames(lr_cat))
  if (length(n_col) == 0) n_col <- colnames(lr_cat)[2]
  lr_cat$category <- pretty_cat(lr_cat[[cat_col]])
  lr_cat$count <- as.numeric(lr_cat[[n_col[1]]])
  lr_cat <- lr_cat[order(-lr_cat$count), ]
  lr_cat$category <- factor(lr_cat$category, levels = rev(lr_cat$category))

  p <- ggplot(lr_cat, aes(x = category, y = count)) +
    geom_col(fill = PAL_CAT[1], width = 0.72) +
    geom_text(aes(label = count), hjust = -0.25,
              size = 2.6, family = FONT_FAMILY) +
    coord_flip() +
    labs(x = NULL, y = "Number of L–R pairs",
         title = "L–R pairs by pathway category") +
    theme(axis.text.y = element_text(size = 7.5),
          axis.text.x = element_text(size = 7.5)) +
    ylim(0, max(lr_cat$count) * 1.18)
  save_panel_pdf("Supp08e_LR_category.pdf", W_RIGHT, HDE, print(p))
  p_E <<- p
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Panel F: pathway communication scores (89×70mm, 7pt)
# =============================================================================
cat("  Panel F: Pathway communication scores\n")
p_F <- NULL
tryCatch({
  pw_comm <- read.csv(file.path(RES, "ligand_receptor_network/pathway_communication_summary.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  pw_comm$pathway <- rownames(pw_comm)
  pw_comm <- pw_comm[order(pw_comm$mean_LR_score), ]
  pw_comm$pathway <- factor(pw_comm$pathway, levels = pw_comm$pathway)
  pw_comm$direction <- ifelse(pw_comm$mean_LR_score > 0, "Enhanced", "Reduced")

  p <- ggplot(pw_comm, aes(x = pathway, y = mean_LR_score, fill = direction)) +
    geom_col(width = 0.78) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c("Enhanced" = COL_UP, "Reduced" = COL_DOWN)) +
    coord_flip() +
    labs(x = NULL, y = "Mean L–R communication score",
         title = "Pathway-level communication change", fill = "Direction") +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          legend.position = c(0.85, 0.20),
          legend.background = element_rect(fill = "white", color = "grey80",
                                           linewidth = 0.2))
  save_panel_pdf("Supp08f_pathway_communication.pdf", W_LEFT, HFG, print(p))
  p_F <<- p
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Panel G: Deconvolution validation (94×70mm)
# =============================================================================
cat("  Panel G: Deconvolution validation\n")
p_G <- NULL
tryCatch({
  val <- read.csv(file.path(RES, "enhancement16_deconvolution/bayesprism_validation_vs_ssGSEA.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  val$cell_label <- gsub("_", " ", val$bayesprism_type)
  val$sig <- ifelse(val$pvalue < 0.05, "Significant", "NS")

  p <- ggplot(val, aes(x = reorder(cell_label, correlation), y = correlation, fill = sig)) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", correlation)),
              hjust = ifelse(val$correlation > 0, -0.15, 1.15),
              size = 2.6, family = FONT_FAMILY) +
    scale_fill_manual(values = c("Significant" = PAL_CAT[1], "NS" = COL_NS)) +
    coord_flip() +
    labs(x = NULL, y = "Spearman correlation (BayesPrism vs ssGSEA)",
         title = "Deconvolution cross-validation", fill = NULL) +
    theme(axis.text.y = element_text(size = 7.5),
          axis.text.x = element_text(size = 7.5),
          legend.position = c(0.82, 0.20),
          legend.background = element_rect(fill = "white", color = "grey80",
                                           linewidth = 0.2)) +
    ylim(-0.2, 1.15)
  save_panel_pdf("Supp08g_deconv_validation.pdf", W_RIGHT, HFG, print(p))
  p_G <<- p
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Panel H: Cross-method concordance (FULL-WIDTH 183×30mm to remove whitespace)
# =============================================================================
cat("  Panel H: Sensitivity analysis (full-width)\n")
p_H <- NULL
tryCatch({
  sens <- read.csv(file.path(RES, "upgrade_deconv_sensitivity/sensitivity_analysis_full_table.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  sens$cell_label <- gsub("_", " ", sens$cell_type)
  sens <- sens[order(-sens$mean_concordance), ]
  sens$cell_label <- factor(sens$cell_label, levels = sens$cell_label)

  # Full-width vertical bars (cell types on x-axis) maximises horizontal real-estate
  p <- ggplot(sens, aes(x = cell_label, y = mean_concordance)) +
    geom_col(fill = PAL_CAT[2], width = 0.66) +
    geom_text(aes(label = sprintf("%.3f", mean_concordance)),
              vjust = -0.5, size = 2.6, fontface = "bold", family = FONT_FAMILY) +
    labs(x = NULL, y = "Concordance",
         title = "Deconvolution sensitivity (multi-method agreement across 9 cell types)") +
    theme(axis.text.x = element_text(size = 7.5, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 7.5),
          axis.title.y = element_text(size = 8, face = "bold",
                                      margin = margin(r = 1.5, unit = "mm")),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          plot.margin = margin(2, 4, 2, 6, "mm")) +
    coord_cartesian(ylim = c(0, 1.18), clip = "off")
  save_panel_pdf("Supp08h_sensitivity.pdf", W_TOTAL, HH, print(p))
  p_H <<- p
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================
# Composite Assembly — Pure-Vector Grid Viewport (zero distortion)
# Each Panel was rendered at EXACT viewport target dims → no scaling.
# =============================================================================
cat("\n--- Assembling composite SuppFig_08 (vector, 183x275mm, 600DPI) ---\n")

tryCatch({
  DPI <- 600
  # Y stack (top → bottom): A | BC | DE | FG | H
  yA  <- HBC + HDE + HFG + HH    # = 220
  yBC <- HDE + HFG + HH          # = 170
  yDE <- HFG + HH                # = 100
  yFG <- HH                      # = 30
  yH  <- 0

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    # Panel A (Row1 full-width heatmap, 183×55)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(yA,"mm"),
      width=unit(W_TOTAL,"mm"), height=unit(HA,"mm"), just=c("left","bottom")))
    draw(ht_A, padding = unit(c(2, 8, 2, 4), "mm"), newpage = FALSE)
    grid::popViewport()
    # Panel B (Row2 left, 89×50)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(yBC,"mm"),
      width=unit(W_LEFT,"mm"), height=unit(HBC,"mm"), just=c("left","bottom")))
    print(p_B, newpage=FALSE); grid::popViewport()
    # Panel C (Row2 right, 94×50)
    grid::pushViewport(grid::viewport(x=unit(W_LEFT,"mm"), y=unit(yBC,"mm"),
      width=unit(W_RIGHT,"mm"), height=unit(HBC,"mm"), just=c("left","bottom")))
    print(p_C, newpage=FALSE); grid::popViewport()
    # Panel D (Row3 left, 89×70)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(yDE,"mm"),
      width=unit(W_LEFT,"mm"), height=unit(HDE,"mm"), just=c("left","bottom")))
    print(p_D, newpage=FALSE); grid::popViewport()
    # Panel E (Row3 right, 94×70)
    grid::pushViewport(grid::viewport(x=unit(W_LEFT,"mm"), y=unit(yDE,"mm"),
      width=unit(W_RIGHT,"mm"), height=unit(HDE,"mm"), just=c("left","bottom")))
    print(p_E, newpage=FALSE); grid::popViewport()
    # Panel F (Row4 left, 89×70)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(yFG,"mm"),
      width=unit(W_LEFT,"mm"), height=unit(HFG,"mm"), just=c("left","bottom")))
    print(p_F, newpage=FALSE); grid::popViewport()
    # Panel G (Row4 right, 94×70)
    grid::pushViewport(grid::viewport(x=unit(W_LEFT,"mm"), y=unit(yFG,"mm"),
      width=unit(W_RIGHT,"mm"), height=unit(HFG,"mm"), just=c("left","bottom")))
    print(p_G, newpage=FALSE); grid::popViewport()
    # Panel H (Row5 FULL-WIDTH 183×30) — eliminates bottom-right whitespace
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(yH,"mm"),
      width=unit(W_TOTAL,"mm"), height=unit(HH,"mm"), just=c("left","bottom")))
    print(p_H, newpage=FALSE); grid::popViewport()

    # Panel labels — top-aligned within each row
    label_data <- data.frame(
      text = c("a","b","c","d","e","f","g","h"),
      x_mm = c(2,
               2,    W_LEFT+2,
               2,    W_LEFT+2,
               2,    W_LEFT+2,
               2),
      y_mm = c(H_TOTAL-2,
               yA-2,    yA-2,
               yBC-2,   yBC-2,
               yDE-2,   yDE-2,
               yFG-2),
      stringsAsFactors = FALSE)
    for (i in seq_len(nrow(label_data))) {
      grid::grid.text(label = label_data$text[i],
        x = unit(label_data$x_mm[i], "mm"),
        y = unit(label_data$y_mm[i], "mm"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY))
    }
    grid::popViewport()
  }

  # Vector PDF (AI-editable)
  cairo_pdf(file.path(OUT, "SuppFig_08.pdf"),
            width = W_TOTAL/25.4, height = H_TOTAL/25.4, family = FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_08.pdf (vector)\n")

  # PNG for review
  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT, "SuppFig_08.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_08.png\n")

  # TIFF for journal submission
  grDevices::tiff(file.path(OUT, "SuppFig_08.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_08.tiff\n")

  cat("  SuppFig_08 v2 DONE — zero-distortion, vector AI-editable.\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 8 v2 polish render complete ===\n")
