#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_11_standalone.R
# Supplementary Figure 11: Extended Molecular Subtyping and Network Validation
# HAE Multi-omics Study | Target: EBioMedicine (Lancet family)
# =============================================================================
# Panels (8):
#   A = Internal cluster validation indices (CH/DB/Dunn/Silhouette, all vote K=2)
#   B = Silhouette per-sample barplot K=2
#   C = Clinical association (lollipop, -log10P)
#   D = Protein-metabolite correlation distribution
#   E = Hub proteins lollipop by metabolite connectivity
#   F = Druggable target category summary
#   G = CS1 vs CS2 cell type dumbbell plot
#   H = Multi-method subtype concordance
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_11_standalone.R
# =============================================================================

cat("=== Supplementary Figure 11: Extended Subtyping & Network ===\n")

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
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_11")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM <- 25.4

FS_GEOM <- 2.82; FS_TAG <- 12
COL_UP <- "#CD534CFF"; COL_DOWN <- "#0073C2FF"
COL_CS1 <- "#CD534CFF"; COL_CS2 <- "#0073C2FF"
COL_TC <- "#0073C2FF"; COL_PR <- "#CD534CFF"; COL_MT <- "#EFC000FF"
PAL_CAT <- pal_jco("default")(10)

# =============================================================================
# Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text = element_text(family = FONT_FAMILY, size = 8),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = 8, color = "black", family = FONT_FAMILY),
    axis.title = element_text(size = 9, face = "bold", family = FONT_FAMILY),
    plot.title = element_text(size = 10, face = "bold", hjust = 0, family = FONT_FAMILY),
    legend.text = element_text(size = 8, family = FONT_FAMILY),
    legend.title = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    legend.key.size = unit(3, "mm"),
    legend.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    plot.margin = margin(2, 3, 2, 3, "mm")
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
# PANEL A: Internal cluster validation indices (all 4 vote K=2)
# FIX: strip.text >= 8pt (was 7pt — NC violation)
# =============================================================================
cat("\n--- Panel A: Internal validation indices ---\n")
tryCatch({
  cvi <- read.csv(file.path(RES, "subtype_robustness_v2/cluster_validation_indices.csv"),
                  stringsAsFactors = FALSE)

  cvi_long <- cvi %>%
    dplyr::select(K, CH, DB, Dunn, Silhouette) %>%
    pivot_longer(-K, names_to = "Index", values_to = "Value")
  cvi_long$Index <- factor(cvi_long$Index,
    levels = c("CH", "DB", "Dunn", "Silhouette"),
    labels = c("Calinski-Harabasz\n(higher = better)",
               "Davies-Bouldin\n(lower = better)",
               "Dunn Index\n(higher = better)",
               "Silhouette\n(higher = better)"))

  opt_points <- cvi_long %>%
    group_by(Index) %>%
    filter(
      (grepl("lower", as.character(Index)) & Value == min(Value)) |
      (!grepl("lower", as.character(Index)) & Value == max(Value))
    ) %>%
    slice(1) %>%
    ungroup()

  p <- ggplot(cvi_long, aes(x = factor(K), y = Value, group = 1)) +
    geom_line(linewidth = 0.7, color = "grey40") +
    geom_point(size = 2, color = "grey40") +
    geom_point(data = opt_points, size = 3.5, color = COL_UP, shape = 18) +
    facet_wrap(~Index, scales = "free_y", nrow = 2) +
    labs(title = "Internal validation (optimal K vote: 4/4 = K=2)",
         x = "K", y = "Index value") +
    # FIX: strip.text must be >= 8pt (NC compliance)
    theme(strip.text = element_text(size = 8, face = "bold", family = FONT_FAMILY),
          plot.title = element_text(size = 9, face = "bold"))
  save_panel_pdf("Supp11a_validation_indices.pdf", 89, 75, print(p))
  p_A <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL B: Silhouette per-sample K=2  (canonical plot, keep as bar)
# =============================================================================
cat("\n--- Panel B: Silhouette per-sample ---\n")
tryCatch({
  sil <- read.csv(file.path(RES, "phase6_subtyping/silhouette_per_sample_K2.csv"),
                  stringsAsFactors = FALSE)
  sil <- sil[order(sil$cluster, -sil$sil_width), ]
  sil$sample <- factor(sil$sample, levels = sil$sample)
  sil$cluster <- factor(sil$cluster)

  p <- ggplot(sil, aes(x = sample, y = sil_width, fill = cluster)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("1" = COL_CS1, "2" = COL_CS2),
                      labels = c("CS1", "CS2"), name = "Subtype") +
    geom_hline(yintercept = mean(sil$sil_width), linetype = "dashed",
               linewidth = 0.4, color = "grey40") +
    annotate("text", x = nrow(sil) * 0.8, y = mean(sil$sil_width) + 0.03,
             label = sprintf("Mean = %.3f", mean(sil$sil_width)),
             size = FS_GEOM, family = FONT_FAMILY) +
    labs(title = "Silhouette widths (K=2)", x = "Sample", y = "Silhouette width") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  save_panel_pdf("Supp11b_silhouette_K2.pdf", 89, 75, print(p))
  p_B <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL C: Clinical association — LOLLIPOP chart
# FIX1: bar -> lollipop (eliminate 3rd consecutive horizontal bar)
# FIX2: plain "-log10(P)" avoids subscript font-size violation (was 6.3pt)
# =============================================================================
cat("\n--- Panel C: Clinical association (lollipop) ---\n")
tryCatch({
  clin <- read.csv(file.path(RES, "phase6_subtyping/clinical_association_K2.csv"),
                   stringsAsFactors = FALSE)
  clin$neg_log10p <- -log10(clin$pvalue + 1e-10)
  clin <- clin[order(clin$pvalue), ]

  var_labels <- c(
    WBC = "WBC", total_bilirubin = "Total bilirubin", ABZ_treatment = "ABZ treatment",
    ethnicity = "Ethnicity", surgery_history = "Surgery history",
    disease_duration = "Disease duration", ALT = "ALT", albumin = "Albumin",
    lesion_location = "Lesion location", GGT = "GGT",
    neutrophil_pct = "Neutrophil %", lymphocyte_pct = "Lymphocyte %",
    age = "Age", endemic_area = "Endemic area", cholesterol = "Cholesterol",
    lesion_lodation = "Lesion location"
  )
  clin$var_label <- ifelse(clin$variable %in% names(var_labels),
                           var_labels[clin$variable], gsub("_", " ", clin$variable))
  clin$var_label <- factor(clin$var_label, levels = rev(head(clin$var_label, 15)))
  clin$sig <- ifelse(clin$pvalue < 0.05, "Significant", "Not significant")
  clin15 <- head(clin, 15)

  thresh <- -log10(0.05)

  p <- ggplot(clin15, aes(x = neg_log10p, y = var_label, color = sig)) +
    # Lollipop stem
    geom_segment(aes(x = 0, xend = neg_log10p, yend = var_label),
                 linewidth = 0.6) +
    # Lollipop head
    geom_point(size = 3) +
    scale_color_manual(values = c("Significant" = COL_UP, "Not significant" = PAL_CAT[1]),
                       guide = "none") +
    geom_vline(xintercept = thresh, linetype = "dashed",
               linewidth = 0.3, color = "grey40") +
    # FIX: plain string avoids subscript size violation
    annotate("text", x = thresh + 0.05, y = 0.55,
             label = "P = 0.05", size = FS_GEOM, family = FONT_FAMILY,
             hjust = 0, color = "grey40") +
    labs(title = "Clinical variable association",
         x = "-log10(P)", y = NULL) +
    theme(axis.text.y = element_text(size = 8))
  save_panel_pdf("Supp11c_clinical_assoc.pdf", 89, 80, print(p))
  p_C <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL D: Protein-metabolite correlation distribution (histogram, keep)
# =============================================================================
cat("\n--- Panel D: Protein-metabolite correlations ---\n")
tryCatch({
  pm <- read.csv(file.path(RES, "protein_metabolite_network/protein_metabolite_correlations.csv"),
                 stringsAsFactors = FALSE)
  p <- ggplot(pm, aes(x = rho)) +
    geom_histogram(binwidth = 0.05, fill = PAL_CAT[2], color = "white", linewidth = 0.2) +
    geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed",
               linewidth = 0.4, color = COL_UP) +
    annotate("text", x = 0.63, y = max(table(cut(pm$rho, breaks = seq(-1, 1, 0.05)))) * 0.85,
             label = sprintf("n = %d\nedges (FDR < 0.05)", nrow(pm)),
             size = FS_GEOM, family = FONT_FAMILY, hjust = 0) +
    labs(title = "Protein-metabolite correlations",
         x = "Spearman rho", y = "Count")
  save_panel_pdf("Supp11d_pm_corr_dist.pdf", 89, 75, print(p))
  p_D <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL E: Hub proteins — LOLLIPOP chart (replaces bar chart)
# More elegant ranking visualization for top-journal standard
# =============================================================================
cat("\n--- Panel E: Hub proteins (lollipop) ---\n")
tryCatch({
  pm <- read.csv(file.path(RES, "protein_metabolite_network/protein_metabolite_correlations.csv"),
                 stringsAsFactors = FALSE)
  hub <- pm %>% count(protein, name = "degree") %>% arrange(desc(degree)) %>% head(20)
  hub$protein <- factor(hub$protein, levels = rev(hub$protein))

  p <- ggplot(hub, aes(x = degree, y = protein, color = degree)) +
    # Lollipop stem (light grey)
    geom_segment(aes(x = 0, xend = degree, yend = protein),
                 linewidth = 0.5, color = "grey75") +
    # Lollipop head with gradient
    geom_point(size = 3.5) +
    scale_color_gradient(low = "#B8D4E3", high = "#1B7837", guide = "none") +
    geom_text(aes(label = degree, color = NULL), hjust = -0.5,
              size = FS_GEOM, family = FONT_FAMILY, color = "grey30") +
    labs(title = "Top hub proteins (by metabolite degree)",
         x = "Degree", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)))
  save_panel_pdf("Supp11e_hub_proteins.pdf", 89, 85, print(p))
  p_E <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL F: Druggable target categories (bar, 5 categories, distinct colors — keep)
# =============================================================================
cat("\n--- Panel F: Druggable target categories ---\n")
tryCatch({
  drug <- read.csv(file.path(RES, "phase7_characterization/druggable_targets_by_subtype.csv"),
                   stringsAsFactors = FALSE)
  fam_summary <- drug %>% count(family, name = "count") %>% arrange(desc(count))

  fam_labels <- c(
    kinases = "Kinases", metabolic_enzymes = "Metabolic enzymes",
    ECM_targets = "ECM targets", immune_checkpoints = "Immune checkpoints",
    epigenetic_regulators = "Epigenetic regulators"
  )
  fam_summary$fam_label <- ifelse(fam_summary$family %in% names(fam_labels),
                                  fam_labels[fam_summary$family],
                                  gsub("_", " ", fam_summary$family))
  fam_summary$fam_label <- factor(fam_summary$fam_label,
                                  levels = rev(fam_summary$fam_label))

  p <- ggplot(fam_summary, aes(x = count, y = fam_label, fill = fam_label)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = setNames(PAL_CAT[1:nrow(fam_summary)],
                                        levels(fam_summary$fam_label)),
                      guide = "none") +
    geom_text(aes(label = count), hjust = -0.2, size = FS_GEOM, family = FONT_FAMILY) +
    labs(title = "Druggable targets by category",
         x = "Number of targets", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)))
  save_panel_pdf("Supp11f_drug_category.pdf", 89, 75, print(p))
  p_F <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL G: CS1 vs CS2 — DUMBBELL PLOT (replaces grouped bar chart)
# Connects CS1 and CS2 scores per cell type with a segment + colored dots.
# Sorted by absolute difference; significance asterisks at right margin.
# =============================================================================
cat("\n--- Panel G: Cell type dumbbell (CS1 vs CS2) ---\n")
tryCatch({
  imm <- read.csv(file.path(RES, "phase7_characterization/subtype_immune_scores.csv"),
                  stringsAsFactors = FALSE)
  imm$cell_type <- gsub("_", " ", imm$cell_type)
  imm <- imm[order(imm$pvalue), ]
  imm <- imm[!duplicated(imm$cell_type), ]
  # Sort by absolute difference (most divergent at top)
  imm <- imm[order(abs(imm$diff), decreasing = TRUE), ]
  imm_top <- head(imm, 12)

  imm_top$sig_label <- ifelse(imm_top$padj < 0.001, "***",
                       ifelse(imm_top$padj < 0.01,  "**",
                       ifelse(imm_top$padj < 0.05,  "*", "")))

  ct_levels <- rev(imm_top$cell_type)
  imm_top$cell_type <- factor(imm_top$cell_type, levels = ct_levels)

  # Long format for dots
  imm_long <- imm_top %>%
    dplyr::select(cell_type, mean_CS1, mean_CS2) %>%
    pivot_longer(-cell_type, names_to = "subtype", values_to = "score")
  imm_long$subtype <- gsub("mean_", "", imm_long$subtype)
  imm_long$cell_type <- factor(imm_long$cell_type, levels = ct_levels)

  # Significance at right margin
  x_max <- max(c(imm_top$mean_CS1, imm_top$mean_CS2), na.rm = TRUE)
  sig_df <- imm_top %>%
    filter(sig_label != "") %>%
    mutate(x_pos = x_max * 1.08)

  p <- ggplot() +
    # Connecting segment between CS1 and CS2
    geom_segment(data = imm_top,
                 aes(x = mean_CS1, xend = mean_CS2,
                     y = cell_type, yend = cell_type),
                 linewidth = 0.9, color = "grey65") +
    # Zero reference line
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.3, color = "grey50") +
    # Dots for CS1 and CS2
    geom_point(data = imm_long,
               aes(x = score, y = cell_type, color = subtype),
               size = 3.2) +
    scale_color_manual(values = c(CS1 = COL_CS1, CS2 = COL_CS2), name = "Subtype") +
    labs(title = "Cell type enrichment by subtype",
         x = "Mean enrichment score", y = NULL) +
    theme(legend.position = "top",
          legend.margin = margin(0, 0, 0, 0))
  if (nrow(sig_df) > 0) {
    p <- p + geom_text(data = sig_df,
                       aes(x = x_pos, y = cell_type, label = sig_label),
                       inherit.aes = FALSE, size = 3, family = FONT_FAMILY,
                       vjust = 0.75, color = "grey30")
  }
  save_panel_pdf("Supp11g_celltype_subtype.pdf", 89, 85, print(p))
  p_G <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL H: Multi-method subtype concordance (tile, keep)
# =============================================================================
cat("\n--- Panel H: Multi-method concordance ---\n")
tryCatch({
  mm <- read.csv(file.path(RES, "enhancement22_dl_fusion/multi_method_comparison.csv"),
                 stringsAsFactors = FALSE)
  mm <- mm[order(mm$ConsensusCluster, mm$sample), ]
  sample_order <- mm$sample

  mm_long <- mm %>%
    pivot_longer(-sample, names_to = "method", values_to = "assignment")
  mm_long$method <- gsub("_cluster|_prediction", "", mm_long$method)
  mm_long$method <- gsub("ConsensusCluster", "Consensus", mm_long$method)
  mm_long$assignment <- gsub("CS|MOFA_C|SNF_C|RF_CS", "", mm_long$assignment)
  mm_long$assignment <- as.integer(mm_long$assignment)
  mm_long$sample <- factor(mm_long$sample, levels = sample_order)
  mm_long$method <- factor(mm_long$method,
    levels = c("Consensus", "MOFA", "SNF", "RF"))

  p <- ggplot(mm_long, aes(x = method, y = sample, fill = factor(assignment))) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("1" = COL_CS1, "2" = COL_CS2),
                      labels = c("CS1/C1", "CS2/C2"), name = "Assignment") +
    labs(title = "Multi-method subtype agreement",
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  save_panel_pdf("Supp11h_method_concordance.pdf", 89, 80, print(p))
  p_H <- p
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# Composite Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite SuppFig_11 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  H1 <- 61; H2 <- 61; H3 <- 61; H4 <- H_TOTAL - H1 - H2 - H3
  W_L <- 91; W_R <- W_TOTAL - W_L

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_A, newpage=FALSE); grid::popViewport()
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_B, newpage=FALSE); grid::popViewport()
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(p_C, newpage=FALSE); grid::popViewport()
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(p_D, newpage=FALSE); grid::popViewport()
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(p_E, newpage=FALSE); grid::popViewport()
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(p_F, newpage=FALSE); grid::popViewport()
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_L,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(p_G, newpage=FALSE); grid::popViewport()
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(0,"mm"),
      width=unit(W_R,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(p_H, newpage=FALSE); grid::popViewport()
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

  cairo_pdf(file.path(OUT,"SuppFig_11.pdf"),
            width=W_TOTAL/25.4, height=H_TOTAL/25.4, family=FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_11.pdf (vector)\n")

  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT,"SuppFig_11.png"),
                 width=px_W, height=px_H, res=DPI, type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_11.png\n")

  grDevices::tiff(file.path(OUT,"SuppFig_11.tiff"),
                  width=px_W, height=px_H, res=DPI, compression="lzw", type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_11.tiff\n")

  cat("  SuppFig_11 DONE\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 11 rendering complete ===\n")

# Clean up stale files
stale <- list.files(OUT, pattern = "^Supp10[a-h]_|^Supp11a_K_metrics", full.names = TRUE)
if (length(stale) > 0) { file.remove(stale); cat("  Cleaned", length(stale), "stale files\n") }
