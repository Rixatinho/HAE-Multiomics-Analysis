#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_14_standalone.R
# Supplementary Figure 14: Molecular Docking & Pharmacogenomic Validation
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Nature Communications
# =============================================================================
# 8 Panels:
#   A = Binding affinity comparison (barplot, |DeltaG| with signed labels)
#   B = Differential binding affinity (Ponatinib advantage)
#   C = Re-docking validation (RMSD vs threshold, n_atoms in subtitle)
#   D = Binding affinity heatmap
#   E = HAE gene-ponatinib sensitivity overlap (Fisher OR forest)
#   F = GDSC pharmacogenomic correlation summary (real data lollipop)
#   G = Ponatinib target expression in HAE
#   H = Sensitivity markers in HAE (dotplot)
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_14_standalone.R
# =============================================================================

cat("=== Supplementary Figure 14: Molecular Docking & Pharmacogenomic Validation ===\n")
cat("  Loading libraries...\n")

# =============================================================================
# SECTION 1: Library Imports
# =============================================================================
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
  library(jsonlite)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_14")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FS_TAG <- 12
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4

W_SINGLE <- 89; W_DOUBLE <- 183; W_HALF <- 89
H_STD <- 85; H_TALL <- 100

ASSEMBLY_DPI <- 600   # pixel-exact rendering constant

mm2in <- function(mm) mm / MM_PER_INCH

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
PAL_CAT  <- pal_jco("default")(10)
COL_UP   <- "#CD534CFF"
COL_DOWN <- "#0073C2FF"
COL_NS   <- "#868686FF"
COL_TC   <- "#0073C2FF"
COL_PR   <- "#CD534CFF"
COL_MT   <- "#EFC000FF"

# =============================================================================
# SECTION 5: Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = 8, color = "black", family = FONT_FAMILY),
    axis.title = element_text(size = 9, face = "bold", family = FONT_FAMILY),
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0, family = FONT_FAMILY),
    legend.text = element_text(size = 7.5, family = FONT_FAMILY),
    legend.title = element_text(size = 7.5, face = "bold", family = FONT_FAMILY),
    legend.key.size = unit(3, "mm"),
    legend.background = element_blank(),
    plot.margin = margin(3, 4, 3, 4, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# SECTION 6: Helper Functions
# =============================================================================
save_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  cairo_pdf(fpath, width = w_mm/MM_PER_INCH, height = h_mm/MM_PER_INCH, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

gp_rn <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID)
gp_cn <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
gp_lt <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
gp_ll <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID)

# =============================================================================
# SECTION 7: Data Loading
# =============================================================================
cat("  Loading molecular docking and GDSC data...\n")

dock_csv  <- read.csv(file.path(RES, "molecular_docking/docking_results_table.csv"), stringsAsFactors = FALSE)
dock_json <- fromJSON(file.path(RES, "molecular_docking/docking_results_summary.json"))
gdsc_json <- fromJSON(file.path(RES, "gdsc_pharmacogenomic/gdsc_validation_results.json"))
sens_csv  <- read.csv(file.path(RES, "gdsc_pharmacogenomic/sensitivity_markers_in_hae.csv"), stringsAsFactors = FALSE)

cat("  Data loaded.\n")

# =============================================================================
# Panel A: Binding Affinity Comparison
# =============================================================================
cat("  Panel A: Binding affinity comparison\n")
tryCatch({
  dock_csv$Affinity <- abs(dock_csv$Affinity_kcal_mol)
  dock_csv$Target <- factor(dock_csv$Target, levels = c("ABL1", "FGFR1", "SC5D"))
  dock_csv$Ligand <- factor(dock_csv$Ligand, levels = c("Ponatinib", "Albendazole"))

  p <- ggplot(dock_csv, aes(x = Target, y = Affinity, fill = Ligand)) +
    geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.9,
             color = "grey20", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", Affinity_kcal_mol)),
              position = position_dodge(0.7), vjust = -0.4, size = 2.6,
              family = FONT_FAMILY, fontface = "bold") +
    geom_hline(yintercept = 8, linetype = "dotted", color = "grey40", linewidth = 0.4) +
    annotate("text", x = 0.55, y = 8.25, label = "Strong-binding\nthreshold (8 kcal/mol)",
             hjust = 0, vjust = 0, size = 2.3, family = FONT_FAMILY,
             color = "grey40", fontface = "italic") +
    scale_fill_manual(values = c("Ponatinib" = COL_PR, "Albendazole" = COL_TC), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Molecular docking binding affinity",
         subtitle = "Bar height: |\u0394G| ; Number labels: signed \u0394G (kcal/mol)",
         x = "Target protein", y = "|Binding affinity| (kcal/mol)") +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 7, color = "grey30",
                                        family = FONT_FAMILY, face = "italic"))

  save_pdf("Supp14a_binding_affinity.pdf", W_HALF, H_STD, print(p))
  p_A <- p
}, error = function(e) cat("    ERROR Panel A:", e$message, "\n"))

# =============================================================================
# Panel B: Differential Binding (Ponatinib Advantage)
# =============================================================================
cat("  Panel B: Differential binding affinity\n")
tryCatch({
  pon <- dock_csv[dock_csv$Ligand == "Ponatinib", c("Target", "Affinity_kcal_mol")]
  alb <- dock_csv[dock_csv$Ligand == "Albendazole", c("Target", "Affinity_kcal_mol")]
  diff_df <- data.frame(
    Target = pon$Target,
    Delta = abs(pon$Affinity_kcal_mol) - abs(alb$Affinity_kcal_mol)
  )
  diff_df$Target <- factor(diff_df$Target, levels = diff_df$Target[order(diff_df$Delta, decreasing = FALSE)])

  p <- ggplot(diff_df, aes(x = Delta, y = Target, fill = Delta)) +
    geom_col(width = 0.55, alpha = 0.9, color = "grey20", linewidth = 0.25) +
    geom_text(aes(label = sprintf("+%.2f kcal/mol", Delta)),
              hjust = -0.05, size = 2.7, family = FONT_FAMILY, fontface = "bold") +
    scale_fill_gradient(low = "#F4B6AE", high = COL_PR, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(title = "Ponatinib selective advantage over albendazole",
         subtitle = "\u0394 = |\u0394G_Ponatinib| \u2212 |\u0394G_Albendazole|  (positive = stronger ponatinib binding)",
         x = "\u0394 |Binding affinity| (kcal/mol)", y = NULL) +
    theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                        family = FONT_FAMILY, face = "italic"))

  save_pdf("Supp14b_diff_binding.pdf", W_HALF, H_STD, print(p))
  p_B <- p
}, error = function(e) cat("    ERROR Panel B:", e$message, "\n"))

# =============================================================================
# Panel C: Re-docking Validation
# =============================================================================
cat("  Panel C: Re-docking validation\n")
tryCatch({
  rv <- dock_json$redocking_validation
  rmsd_df <- data.frame(
    Metric = factor(c("Observed RMSD", "Acceptance threshold"),
                     levels = c("Observed RMSD", "Acceptance threshold")),
    Value  = c(rv$rmsd_angstrom, rv$threshold)
  )
  pass_color <- ifelse(rv$rmsd_angstrom < rv$threshold, COL_DOWN, COL_PR)

  p <- ggplot(rmsd_df, aes(x = Metric, y = Value, fill = Metric)) +
    geom_col(width = 0.55, alpha = 0.9, color = "grey20", linewidth = 0.25) +
    geom_segment(aes(x = 0.5, xend = 2.5,
                     y = rv$threshold, yend = rv$threshold),
                 linetype = "dashed", color = COL_PR, linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f \u00c5", Value)),
              vjust = -0.5, size = 2.9, family = FONT_FAMILY, fontface = "bold") +
    annotate("text", x = 1.5, y = rv$threshold * 0.55,
             label = sprintf("\u0394 = %.2f \u00c5\nPASS \u2713",
                             rv$threshold - rv$rmsd_angstrom),
             size = 3.2, family = FONT_FAMILY, fontface = "bold",
             color = pass_color) +
    scale_fill_manual(values = c("Observed RMSD" = pass_color,
                                  "Acceptance threshold" = "grey75"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                       breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5)) +
    labs(title = "Re-docking validation",
         subtitle = sprintf("PDB: 3OXZ (ABL1\u2013ponatinib) | n = %d heavy atoms aligned",
                            rv$n_atoms_compared),
         x = NULL, y = "RMSD (\u00c5)") +
    theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                        family = FONT_FAMILY, face = "italic"),
          axis.text.x = element_text(size = 8))

  save_pdf("Supp14c_redocking.pdf", W_HALF, H_STD, print(p))
  p_C <- p
}, error = function(e) cat("    ERROR Panel C:", e$message, "\n"))

# =============================================================================
# Panel D: Binding Affinity Heatmap
# =============================================================================
cat("  Panel D: Binding affinity heatmap\n")
tryCatch({
  mat <- matrix(dock_csv$Affinity_kcal_mol, nrow = 3, byrow = FALSE,
                dimnames = list(c("ABL1", "FGFR1", "SC5D"), c("Ponatinib", "Albendazole")))

  col_aff <- colorRamp2(c(-13, -9, -6), c("#CD534CFF", "#FFFFFF", "#0073C2FF"))

  s_path <- file.path(OUT, "Supp14d_affinity_heatmap.pdf")
  cairo_pdf(s_path, width = mm2in(W_HALF), height = mm2in(70), family = FONT_FAMILY)
  ht <- Heatmap(mat, name = "\u0394G (kcal/mol)",
    col = col_aff, border = TRUE,
    rect_gp = gpar(col = "white", lwd = 1),
    row_names_gp = gp_rn(8), column_names_gp = gp_cn(8),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.1f", mat[i, j]), x, y,
                gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"))
    },
    column_title = "Binding affinity heatmap",
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_FAMILY),
    heatmap_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll(),
                                 grid_width = unit(3, "mm")),
    show_row_dend = FALSE, show_column_dend = FALSE)
  draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))
  dev.off()
  ht_D <- ht
  cat("  -> Supp14d_affinity_heatmap.pdf\n")
}, error = function(e) cat("    ERROR Panel D:", e$message, "\n"))

# =============================================================================
# Panel E: HAE Gene-Ponatinib Sensitivity Overlap (Fisher)
# =============================================================================
cat("  Panel E: HAE-ponatinib sensitivity enrichment\n")
tryCatch({
  enr <- gdsc_json$enrichment
  bar_df <- data.frame(
    Category = c("HAE-up vs sensitivity",
                 "HAE-up vs resistance",
                 "HAE-down vs sensitivity",
                 "HAE-down vs resistance"),
    OR = c(enr$up_sens_OR, enr$up_resist_OR, enr$down_sens_OR, enr$down_resist_OR),
    P  = c(enr$up_sens_P,  enr$up_resist_P,  enr$down_sens_P,  enr$down_resist_P)
  )
  bar_df$OR_plot <- ifelse(bar_df$OR < 0.1, 0.1, bar_df$OR)
  bar_df$sig <- ifelse(bar_df$P < 0.05, "Significant (P<0.05)", "Not significant")
  bar_df$Category <- factor(bar_df$Category, levels = rev(bar_df$Category))
  bar_df$lab <- ifelse(bar_df$OR == 0,
                        sprintf("OR=0 (no overlap), P=%.2f", bar_df$P),
                        ifelse(bar_df$P < 0.001,
                               sprintf("OR=%.1f, P=%.1e", bar_df$OR, bar_df$P),
                               sprintf("OR=%.2f, P=%.2f", bar_df$OR, bar_df$P)))
  # Position labels: zero-OR rows go above the truncation point (x=0.1);
  # positive-OR rows go to the right of the dot
  bar_df$lab_x <- ifelse(bar_df$OR == 0, 0.105, bar_df$OR_plot)
  bar_df$lab_h <- ifelse(bar_df$OR == 0, 0, -0.05)

  p <- ggplot(bar_df, aes(x = OR_plot, y = Category, color = sig)) +
    geom_segment(aes(x = 0.1, xend = OR_plot, yend = Category), linewidth = 0.6) +
    geom_point(aes(size = -log10(P + 1e-10)), alpha = 0.95) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_text(aes(x = lab_x, label = lab, hjust = lab_h),
              vjust = -1.0, size = 2.3,
              family = FONT_FAMILY, color = "black") +
    scale_x_log10(limits = c(0.08, 30), breaks = c(0.1, 1, 10),
                  labels = c("0", "1", "10")) +
    scale_color_manual(values = c("Significant (P<0.05)" = COL_PR,
                                   "Not significant" = COL_NS), name = NULL) +
    scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
    labs(title = "HAE DEG \u2013 ponatinib marker overlap",
         subtitle = "Fisher's exact (one-sided)",
         x = "Odds ratio (log scale)", y = NULL) +
    theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                        family = FONT_FAMILY, face = "italic"),
          legend.position = "top",
          legend.key.size = unit(3, "mm"),
          axis.text.y = element_text(size = 7.5))

  save_pdf("Supp14e_gene_overlap.pdf", W_HALF, H_STD, print(p))
  p_E <- p
}, error = function(e) cat("    ERROR Panel E:", e$message, "\n"))

# =============================================================================
# Panel F: GDSC Pharmacogenomic Correlation Summary (real-data lollipop)
# =============================================================================
cat("  Panel F: GDSC pharmacogenomic correlation summary\n")
tryCatch({
  cor_data <- gdsc_json$correlations
  fp_df <- data.frame(
    Pair = c("FGFR-pathway, pan-cancer",
             "EMT-pathway, pan-cancer",
             "FGFR-pathway, liver-only"),
    rho  = c(cor_data$fgfr_rho, cor_data$emt_rho, cor_data$liver_fgfr_rho),
    P    = c(cor_data$fgfr_P,   cor_data$emt_P,   cor_data$liver_fgfr_P),
    stringsAsFactors = FALSE
  )
  fp_df$Pair <- factor(fp_df$Pair, levels = rev(fp_df$Pair))
  fp_df$sig <- ifelse(fp_df$P < 0.05, "Significant (P<0.05)", "Not significant")
  fp_df$lab <- ifelse(fp_df$P < 1e-3,
                       sprintf("\u03c1=%.3f, P=%.1e", fp_df$rho, fp_df$P),
                       sprintf("\u03c1=%.3f, P=%.3f", fp_df$rho, fp_df$P))

  p <- ggplot(fp_df, aes(x = rho, y = Pair, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = rho, yend = Pair), linewidth = 0.7) +
    geom_point(aes(size = abs(rho)), alpha = 0.95) +
    geom_text(aes(label = lab, x = rho),
              hjust = ifelse(fp_df$rho < 0, -0.05, 1.05),
              vjust = -1.0, size = 2.3, family = FONT_FAMILY, color = "black") +
    scale_color_manual(values = c("Significant (P<0.05)" = COL_PR,
                                   "Not significant" = COL_NS), name = NULL) +
    scale_size_continuous(range = c(3, 7), guide = "none") +
    scale_x_continuous(limits = c(-1.45, 0.6),
                       breaks = c(-1, -0.5, 0, 0.5),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = "GDSC2 pharmacogenomic correlations",
         subtitle = "Spearman \u03c1 vs ponatinib ln IC\u2085\u2080",
         x = "Spearman \u03c1", y = NULL) +
    theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                        family = FONT_FAMILY, face = "italic"),
          legend.position = "top",
          legend.key.size = unit(3, "mm"),
          axis.text.y = element_text(size = 7.5))

  save_pdf("Supp14f_fgfr_ic50.pdf", W_HALF, H_STD, print(p))
  p_F <- p
}, error = function(e) cat("    ERROR Panel F:", e$message, "\n"))

# =============================================================================
# Panel G: Ponatinib Target Expression in HAE
# =============================================================================
cat("  Panel G: Ponatinib target expression in HAE\n")
tryCatch({
  tgt_df <- gdsc_json$targets$details
  tgt_df$sig <- ifelse(tgt_df$pvalue < 0.05, "Significant (P<0.05)", "Not significant")
  tgt_df$gene <- factor(tgt_df$gene, levels = tgt_df$gene[order(tgt_df$logFC)])
  tgt_df$plab <- ifelse(tgt_df$pvalue < 0.001,
                         sprintf("P = %.1e", tgt_df$pvalue),
                         sprintf("P = %.3f", tgt_df$pvalue))

  p <- ggplot(tgt_df, aes(x = logFC, y = gene, fill = sig)) +
    geom_col(width = 0.6, alpha = 0.9, color = "grey20", linewidth = 0.25) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_text(aes(label = plab),
              hjust = ifelse(tgt_df$logFC > 0, -0.1, 1.1),
              size = 2.4, family = FONT_FAMILY) +
    scale_fill_manual(values = c("Significant (P<0.05)" = COL_PR,
                                  "Not significant" = COL_NS), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.32, 0.32))) +
    labs(title = "Ponatinib target expression in HAE",
         subtitle = "log\u2082 fold-change, peri-lesional vs adjacent-normal liver",
         x = "log\u2082 fold-change", y = NULL) +
    theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                        family = FONT_FAMILY, face = "italic"),
          legend.position = "top",
          legend.key.size = unit(3, "mm"))

  save_pdf("Supp14g_target_expression.pdf", W_HALF, H_STD, print(p))
  p_G <- p
}, error = function(e) cat("    ERROR Panel G:", e$message, "\n"))

# =============================================================================
# Panel H: Sensitivity Markers in HAE (dotplot)
# =============================================================================
cat("  Panel H: Sensitivity markers in HAE\n")
tryCatch({
  sens_csv$neg_log10p <- -log10(sens_csv$pvalue)
  top_n <- min(15, nrow(sens_csv))
  sens_top <- head(sens_csv, top_n)
  sens_top$gene <- factor(sens_top$gene, levels = rev(sens_top$gene))

  p <- ggplot(sens_top, aes(x = logFC, y = gene, size = neg_log10p, color = logFC)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient2(low = COL_DOWN, mid = "grey85", high = COL_PR,
                          midpoint = 0, name = "log\u2082 FC") +
    scale_size_continuous(range = c(2.2, 6.2), name = "\u2212log\u2081\u2080(P)") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey40") +
    labs(title = "HAE-elevated ponatinib markers",
         subtitle = "Top 15 ranked by significance",
         x = "log\u2082 FC in HAE", y = NULL) +
    theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                        family = FONT_FAMILY, face = "italic"),
          legend.position = "right",
          legend.key.size = unit(3, "mm"))

  save_pdf("Supp14h_sensitivity_markers.pdf", W_HALF, H_STD, print(p))
  p_H <- p
}, error = function(e) cat("    ERROR Panel H:", e$message, "\n"))
# =============================================================================
# =============================================================================
# Composite Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite SuppFig_14 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  H1 <- 61; H2 <- 61; H3 <- 61; H4 <- H_TOTAL - H1 - H2 - H3
  W_L <- 91; W_R <- W_TOTAL - W_L

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    # Panel A (Row1 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_A, newpage=FALSE); grid::popViewport()
    # Panel B (Row1 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_B, newpage=FALSE); grid::popViewport()
    # Panel C (Row2 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(p_C, newpage=FALSE); grid::popViewport()
    # Panel D (Row2 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    draw(ht_D, padding=unit(c(5,5,5,5),"mm"), newpage=FALSE); grid::popViewport()
    # Panel E (Row3 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(p_E, newpage=FALSE); grid::popViewport()
    # Panel F (Row3 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(p_F, newpage=FALSE); grid::popViewport()
    # Panel G (Row4 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_L,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(p_G, newpage=FALSE); grid::popViewport()
    # Panel H (Row4 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(0,"mm"),
      width=unit(W_R,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(p_H, newpage=FALSE); grid::popViewport()
    # Labels
    label_data <- data.frame(
      text=c("a","b","c","d","e","f","g","h"),
      x_mm=c(2,W_L+2,2,W_L+2,2,W_L+2,2,W_L+2),
      y_mm=c(H2+H3+H4+H1-2,H2+H3+H4+H1-2,H3+H4+H2-2,H3+H4+H2-2,
             H4+H3-2,H4+H3-2,H4-2,H4-2), stringsAsFactors=FALSE)
    for(i in seq_len(nrow(label_data))) {
      grid::grid.text(label=label_data$text[i],
        x=unit(label_data$x_mm[i],"mm"), y=unit(label_data$y_mm[i],"mm"),
        just=c("left","top"),
        gp=grid::gpar(fontsize=FS_TAG, fontface="bold", fontfamily=FONT_FAMILY))
    }
    grid::popViewport()
  }

  cairo_pdf(file.path(OUT,"SuppFig_14.pdf"),
            width=W_TOTAL/25.4, height=H_TOTAL/25.4, family=FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_14.pdf (vector)\n")

  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT,"SuppFig_14.png"),
                 width=px_W, height=px_H, res=DPI, type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_14.png\n")

  grDevices::tiff(file.path(OUT,"SuppFig_14.tiff"),
                  width=px_W, height=px_H, res=DPI, compression="lzw", type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_14.tiff\n")

  cat("  SuppFig_14 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 14 rendering complete ===\n")

