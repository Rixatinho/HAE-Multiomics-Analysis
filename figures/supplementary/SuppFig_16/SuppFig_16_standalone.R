#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_16_standalone.R
# Supplementary Figure 16: Independent Single-cell Validation (HRA000553)
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: EBioMedicine (Lancet family)
# =============================================================================
# 4 Panels:
#   A = UMAP of 83,921 CD45+ cells coloured by lineage annotation (16 types)
#   B = Lineage composition across 11 samples (tile heatmap, sqrt scale)
#   C = Paired PL vs AN dot plots for 6 lineages with largest |dz|
#   D = Paired PL vs AN module-score dot plots (6 immune function modules)
# Data: 02_analysis/results/enhancement49_singlecell_validation/
# Usage: /Users/rishat/miniforge3/envs/multiomics/bin/Rscript SuppFig_16_standalone.R
# =============================================================================

cat("=== Supplementary Figure 16: Single-cell Validation ===\n")
cat("  Loading libraries...\n")

# =============================================================================
# SECTION 1: Library Imports
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(ggsci)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results/enhancement49_singlecell_validation")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_16")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FS_TAG <- 12
FONT_GRID <- "Arial"
MM_PER_INCH <- 25.4
DPI <- 600

# Double-column layout, two rows
W_L <- 88; W_R <- 95
H1  <- 78; H2  <- 78
W_TOTAL <- W_L + W_R
H_TOTAL <- H1 + H2

mm2in <- function(mm) mm / MM_PER_INCH

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
PAL_NPG   <- pal_npg("nrc")(10)
PAL_D3    <- pal_d3("category20")(20)
PAL_CAT16 <- c(PAL_NPG, PAL_D3[seq(2, 20, 2)])[1:16]
PAL_PT    <- PAL_NPG[c(1, 3, 4)]  # patients P02 P03 P04
COL_TISSUE <- c("AN" = "#0073C2FF", "PL" = "#CD534CFF")

# =============================================================================
# SECTION 5: Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = 7, color = "black", family = FONT_FAMILY),
    axis.title = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    plot.title = element_text(size = 9, face = "bold", hjust = 0, family = FONT_FAMILY),
    legend.text = element_text(size = 6, family = FONT_FAMILY),
    legend.title = element_text(size = 7, face = "bold", family = FONT_FAMILY),
    legend.key.size = unit(2.2, "mm"),
    legend.background = element_blank(),
    strip.text = element_text(size = 6.5, face = "bold", family = FONT_FAMILY),
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.3),
    plot.margin = margin(1.5, 3, 1.5, 3, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# SECTION 6: Helper Functions
# =============================================================================
save_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  cairo_pdf(fpath, width = w_mm / MM_PER_INCH, height = h_mm / MM_PER_INCH,
            family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

# =============================================================================
# SECTION 7: Data Loading
# =============================================================================
cat("  Loading single-cell validation outputs...\n")

meta <- read.csv(file.path(RES, "cell_metadata_umap.csv"),
                 stringsAsFactors = FALSE)
comp <- read.csv(file.path(RES, "composition_per_sample.csv"),
                 stringsAsFactors = FALSE)
mods <- read.csv(file.path(RES, "module_scores_per_sample.csv"),
                 stringsAsFactors = FALSE)
comp_test <- read.csv(file.path(RES, "composition_PL_vs_AN_paired_wilcoxon.csv"),
                      stringsAsFactors = FALSE)
mod_test  <- read.csv(file.path(RES, "module_PL_vs_AN_paired_wilcoxon.csv"),
                      stringsAsFactors = FALSE)

TISSUE_ORD <- c("PB", "PL", "AN")
TISSUE_LBL  <- c(PB = "Peripheral\nblood", PL = "Peri-lesional", AN = "Adjacent")
meta$tissue  <- factor(meta$tissue, levels = TISSUE_ORD)
comp$tissue  <- factor(comp$tissue, levels = TISSUE_ORD)
mods$tissue  <- factor(mods$tissue, levels = TISSUE_ORD)

cat(sprintf("  meta: %d cells | comp: %d rows | mods: %d rows\n",
            nrow(meta), nrow(comp), nrow(mods)))

# =============================================================================
# SECTION 8: Panel A - UMAP by lineage
# =============================================================================
cat("  Panel A: UMAP...\n")

ct_levels <- names(sort(table(meta$cell_type), decreasing = TRUE))
meta$cell_type <- factor(meta$cell_type, levels = rev(ct_levels))
names(PAL_CAT16) <- ct_levels

p_A <- ggplot(meta, aes(x = umap1, y = umap2, colour = cell_type)) +
  geom_point(size = 0.15, alpha = 0.65, stroke = 0) +
  scale_colour_manual(values = PAL_CAT16, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.2, alpha = 1),
                               ncol = 2)) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "83,921 CD45+ cells, 16 lineages") +
  theme(legend.position = "bottom",
        axis.text = element_blank(), axis.ticks = element_blank())

save_pdf("Supp16a_umap.pdf", W_L, H1, print(p_A))

# =============================================================================
# SECTION 9: Panel B - Composition heatmap (tile)
# =============================================================================
cat("  Panel B: composition heatmap...\n")

# Complete the grid (absent combos = 0 cells)
comp_grid <- tidyr::expand_grid(
  patient = sort(unique(comp$patient)),
  tissue  = TISSUE_ORD,
  majority_voting = ct_levels) %>%
  dplyr::left_join(comp, by = c("patient", "tissue", "majority_voting")) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0L),
                frac = tidyr::replace_na(frac, 0))

comp_grid$sample <- paste0(comp_grid$patient, "_", comp_grid$tissue)
samp_order <- comp_grid %>%
  dplyr::distinct(patient, tissue, sample) %>%
  dplyr::arrange(patient, tissue) %>%
  dplyr::pull(sample)
comp_grid$sample <- factor(comp_grid$sample, levels = samp_order)
comp_grid$majority_voting <- factor(comp_grid$majority_voting, levels = rev(ct_levels))
comp_grid$frac_sqrt <- sqrt(comp_grid$frac)

p_B <- ggplot(comp_grid,
              aes(x = sample, y = majority_voting, fill = frac_sqrt)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "grey97", high = "#0073C2FF",
                      name = expression(sqrt(fraction))) +
  labs(x = NULL, y = NULL,
       title = "Lineage composition per sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 6),
        legend.position = "right")

save_pdf("Supp16b_composition_heatmap.pdf", W_R, H1, print(p_B))

# =============================================================================
# SECTION 10: Panel C - Paired dot plots, 6 lineages with largest |dz|
# =============================================================================
cat("  Panel C: paired lineage dot plots...\n")

top_ct <- comp_test %>%
  dplyr::arrange(dplyr::desc(abs(cohens_dz))) %>%
  dplyr::slice(1:6) %>%
  dplyr::pull(cell_type)

comp_pair <- comp %>%
  dplyr::filter(tissue %in% c("AN", "PL"), majority_voting %in% top_ct) %>%
  dplyr::mutate(tissue = factor(tissue, levels = c("AN", "PL")))

dz_lookup <- setNames(comp_test$cohens_dz, comp_test$cell_type)
comp_pair$facet_lbl <- paste0(comp_pair$majority_voting,
                              "\ndz = ", sprintf("%.2f", dz_lookup[comp_pair$majority_voting]))
facet_ord <- top_ct[order(-abs(dz_lookup[top_ct]))]
comp_pair$facet_lbl <- factor(comp_pair$facet_lbl,
  levels = paste0(facet_ord, "\ndz = ", sprintf("%.2f", dz_lookup[facet_ord])))

p_C <- ggplot(comp_pair,
              aes(x = tissue, y = frac, group = patient, colour = patient)) +
  geom_line(colour = "grey55", linewidth = 0.3) +
  geom_point(size = 1.6) +
  facet_wrap(~facet_lbl, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = setNames(PAL_PT, c("P02", "P03", "P04")),
                      name = "Patient") +
  scale_x_discrete(labels = c("AN" = "AN", "PL" = "PL")) +
  labs(x = NULL, y = "Fraction of CD45+ cells",
       title = "Lineage abundance, PL vs AN (paired)") +
  theme(legend.position = "bottom")

save_pdf("Supp16c_lineage_paired.pdf", W_L, H2, print(p_C))

# =============================================================================
# SECTION 11: Panel D - Paired module-score dot plots
# =============================================================================
cat("  Panel D: module-score dot plots...\n")

mod_labs <- c(score_exhaustion = "Exhaustion",
              score_IFNg_response = "IFN-\u03b3 response",
              score_inflammatory_response = "Inflammatory\nresponse",
              score_complement = "Complement",
              score_TCR_genes = "TCR module",
              score_BCR_genes = "BCR module")

mods_pair <- mods %>%
  dplyr::filter(tissue %in% c("AN", "PL")) %>%
  dplyr::mutate(tissue = factor(tissue, levels = c("AN", "PL"))) %>%
  tidyr::pivot_longer(dplyr::starts_with("score_"),
                      names_to = "module", values_to = "score")

dz_mod <- setNames(mod_test$cohens_dz, mod_test$module)
mods_pair$facet_lbl <- paste0(mod_labs[mods_pair$module],
                              "\ndz = ", sprintf("%.2f", dz_mod[mods_pair$module]))
mod_ord <- names(mod_labs)[order(-abs(dz_mod[names(mod_labs)]))]
mods_pair$facet_lbl <- factor(mods_pair$facet_lbl,
  levels = paste0(mod_labs[mod_ord], "\ndz = ", sprintf("%.2f", dz_mod[mod_ord])))

p_D <- ggplot(mods_pair,
              aes(x = tissue, y = score, group = patient, colour = patient)) +
  geom_line(colour = "grey55", linewidth = 0.3) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70",
             linetype = "dashed") +
  facet_wrap(~facet_lbl, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = setNames(PAL_PT, c("P02", "P03", "P04")),
                      name = "Patient") +
  labs(x = NULL, y = "Module score (per-sample mean)",
       title = "Immune module scores, PL vs AN (paired)") +
  theme(legend.position = "bottom")

save_pdf("Supp16d_module_paired.pdf", W_R, H2, print(p_D))

# =============================================================================
# SECTION 12: Composite Assembly (vector PDF + PNG + TIFF)
# =============================================================================
cat("  Assembling composite...\n")

render_vector_composite <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = 0, y = 0, width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
    just = c("left", "bottom"), gp = grid::gpar(fontsize = 8)))
  # Panel A (Row1 left)
  grid::pushViewport(grid::viewport(x = unit(0, "mm"), y = unit(H2, "mm"),
    width = unit(W_L, "mm"), height = unit(H1, "mm"), just = c("left", "bottom")))
  print(p_A, newpage = FALSE); grid::popViewport()
  # Panel B (Row1 right)
  grid::pushViewport(grid::viewport(x = unit(W_L, "mm"), y = unit(H2, "mm"),
    width = unit(W_R, "mm"), height = unit(H1, "mm"), just = c("left", "bottom")))
  print(p_B, newpage = FALSE); grid::popViewport()
  # Panel C (Row2 left)
  grid::pushViewport(grid::viewport(x = unit(0, "mm"), y = unit(0, "mm"),
    width = unit(W_L, "mm"), height = unit(H2, "mm"), just = c("left", "bottom")))
  print(p_C, newpage = FALSE); grid::popViewport()
  # Panel D (Row2 right)
  grid::pushViewport(grid::viewport(x = unit(W_L, "mm"), y = unit(0, "mm"),
    width = unit(W_R, "mm"), height = unit(H2, "mm"), just = c("left", "bottom")))
  print(p_D, newpage = FALSE); grid::popViewport()
  # Labels
  label_data <- data.frame(
    text = c("A", "B", "C", "D"),
    x_mm = c(2, W_L + 2, 2, W_L + 2),
    y_mm = c(H_TOTAL - 2, H_TOTAL - 2, H2 - 2, H2 - 2),
    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(label_data))) {
    grid::grid.text(label = label_data$text[i],
      x = unit(label_data$x_mm[i], "mm"), y = unit(label_data$y_mm[i], "mm"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = FS_TAG, fontface = "bold",
                      fontfamily = FONT_FAMILY))
  }
  grid::popViewport()
}

tryCatch({
  cairo_pdf(file.path(OUT, "SuppFig_16.pdf"),
            width = W_TOTAL / 25.4, height = H_TOTAL / 25.4,
            family = FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_16.pdf (vector)\n")

  px_W <- round(W_TOTAL * DPI / 25.4); px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "SuppFig_16.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_16.png\n")

  grDevices::tiff(file.path(OUT, "SuppFig_16.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_16.tiff\n")

  cat("  SuppFig_16 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 16 rendering complete ===\n")
