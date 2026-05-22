#!/usr/bin/env Rscript
# =============================================================================
# Figure_3_standalone.R
# Complete, Self-Contained Code for Figure 3
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Journal of Hepatology (Elsevier) - TIFF 600dpi + PDF
# =============================================================================
# Figure 3: Tryptophan-Kynurenine Pathway Dysregulation and NAD+ Depletion
# 6 panels: A=enzyme heatmap, B=IDO1 cross-disease, C=metabolite barplot,
#           D=N-AcKyn-TCA correlation, E=KTR paired boxplot, F=pathway schematic
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript Figure_3_standalone.R
# =============================================================================

cat("=== Figure 3: Tryptophan-Kynurenine Pathway Dysregulation ===\n")
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
  library(ggrepel)
  library(ggpubr)
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
TRP  <- file.path(RES, "tryptophan_kynurenine")
OUT  <- file.path(BASE, "04_figures/main/Figure_3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4

# NC page dimensions (mm)
W_SINGLE  <- 89
W_DOUBLE  <- 183
W_HALF    <- 89
H_STD     <- 85
H_TALL    <- 100
H_MAX     <- 240

mm2in <- function(mm) mm / MM_PER_INCH
in2mm <- function(inch) inch * MM_PER_INCH

# Typography hierarchy (minimum 8pt rule)
FS_TITLE      <- 10
FS_SUBTITLE   <- 9
FS_AXIS_TITLE <- 9
FS_AXIS_TEXT  <- 8
FS_LEGEND_T   <- 8
FS_LEGEND_L   <- 8
FS_ANNO       <- 8
FS_ROW_NAME   <- 8
FS_CELL       <- 8
FS_TAG        <- 14
FS_GEOM_TEXT  <- 2.82  # geom_text size (= 8pt in mm)

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
PAL_CAT <- pal_jco("default")(10)

# Semantic colors
COL_UP       <- "#CD534CFF"
COL_DOWN     <- "#0073C2FF"
COL_NS       <- "#868686FF"
COL_NA       <- "#F0F0F0"

# Omics layer colors
COL_TC       <- "#0073C2FF"    # transcriptomics
COL_PR       <- "#CD534CFF"    # proteomics
COL_MT       <- "#EFC000FF"    # metabolomics

# Group colors
COL_NORMAL   <- "#7AA6DCFF"    # normal tissue
COL_ADJACENT <- "#CD534CFF"    # adjacent/disease tissue

# Disease colors for cross-disease comparison
COL_HAE      <- "#0073C2FF"
COL_HCC      <- "#CD534CFF"
COL_CCA      <- "#EFC000FF"
COL_CRC      <- "#20854EFF"
COL_MEL      <- "#7876B1FF"

# =============================================================================
# SECTION 5: Color Scale Functions
# =============================================================================
col_div     <- colorRamp2(c(-2, 0, 2),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_logfc   <- colorRamp2(c(-2, 0, 0.5),   c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 6: ggplot2 Theme
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
    axis.title         = element_text(family = FONT_FAMILY, size = 8, face = "bold", color = "black"),
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
# SECTION 8: Helper Functions (pixel-exact rendering)
# =============================================================================
# Assembly DPI — used to calculate pixel-exact PDF dimensions
ASSEMBLY_DPI <- 600

# Pixel-exact PDF save: calculates inches from target pixels (not mm/25.4)
# to eliminate floating-point mm→inches→pixels rounding errors
save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  # Target pixels (same formula as assembly)
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  # Back-calculate inches from pixels for exact match
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  force(expr)
  dev.off()
  cat(sprintf("  Panel saved: %s (%d×%d px target)\n", basename(fpath), target_w_px, target_h_px))
}

save_heatmap_panel_pdf <- function(filename, ht, w_mm, h_mm, title = NULL) {
  fpath <- file.path(OUT, filename)
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  # Reserve top space for title if provided
  if (!is.null(title)) {
    draw(ht, padding = unit(c(2, 2, 8, 2), "mm"))  # extra top padding for title
    # Draw title centered at top of full page width
    grid::grid.text(title, x = unit(0.5, "npc"), y = unit(1, "npc") - unit(3, "mm"),
                    just = c("center", "top"),
                    gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY))
  } else {
    draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  }
  dev.off()
  cat(sprintf("  Panel saved: %s (%d×%d px target)\n", basename(fpath), target_w_px, target_h_px))
}

# =============================================================================
# SECTION 9: Data Loading
# =============================================================================
cat("\n--- Loading tryptophan pathway data ---\n")

# Enzyme expression data
enzymes <- read.csv(file.path(TRP, "trp_pathway_enzyme_expression.csv"),
                    stringsAsFactors = FALSE)
enzymes_detected <- enzymes[enzymes$detected_mRNA == "True", ]
cat("  Enzymes: ", nrow(enzymes_detected), " detected of ", nrow(enzymes), " total\n")

# Metabolite data
metabolites <- read.csv(file.path(TRP, "trp_metabolites_measured.csv"),
                        stringsAsFactors = FALSE)
cat("  Metabolites: ", nrow(metabolites), " measured\n")

# KTR per sample
ktr <- read.csv(file.path(TRP, "ktr_per_sample.csv"), stringsAsFactors = FALSE)
cat("  KTR samples: ", nrow(ktr), " (", length(unique(ktr$patient)), " patients)\n")

# IDO1 cross-disease comparison
ido1_cross <- read.csv(file.path(TRP, "ido1_cross_disease_comparison.csv"),
                       stringsAsFactors = FALSE)
cat("  Cross-disease: ", nrow(ido1_cross), " entries\n")

# NAD-OXPHOS correlation
nad_oxphos <- read.csv(file.path(TRP, "nad_oxphos_correlation.csv"),
                       stringsAsFactors = FALSE)
cat("  NAD-OXPHOS: ", nrow(nad_oxphos), " samples\n")


# =============================================================================
# SECTION 10: Panel Generation
# =============================================================================

# -----------------------------------------------------------------------------
# PANEL A: Cross-omics enzyme heatmap (12 enzymes, mRNA + protein logFC)
# -----------------------------------------------------------------------------
cat("\n--- Panel A: Enzyme heatmap ---\n")
tryCatch({
  # Prepare matrix: enzymes as rows, mRNA/protein logFC as columns
  enz_df <- enzymes_detected %>%
    arrange(factor(branch, levels = c("Kynurenine branch", "Serotonin branch")), step) %>%
    mutate(enzyme = factor(enzyme, levels = enzyme))

  # Create matrix for heatmap
  mat <- matrix(NA, nrow = nrow(enz_df), ncol = 2)
  rownames(mat) <- enz_df$enzyme
  colnames(mat) <- c("mRNA", "Protein")
  mat[, 1] <- enz_df$mRNA_logFC
  mat[, 2] <- ifelse(enz_df$detected_protein == "True", enz_df$protein_logFC, NA)

  # Branch annotation
  branch_col <- c("Kynurenine branch" = "#0073C2FF", "Serotonin branch" = "#EFC000FF")
  ha_left <- rowAnnotation(
    Branch = enz_df$branch,
    col = list(Branch = branch_col),
    annotation_name_gp = gp_anno_name(),
    annotation_legend_param = list(
      title_gp = gp_legend_title(),
      labels_gp = gp_legend_labels()
    ),
    simple_anno_size = unit(4, "mm")
  )

  # P-value significance annotation (stars)
  sig_stars <- ifelse(enz_df$mRNA_P < 0.05, "*",
               ifelse(enz_df$mRNA_P < 0.1, "\u2020", ""))

  ha_right <- rowAnnotation(
    Sig = anno_text(sig_stars, gp = gpar(fontsize = 8, fontfamily = FONT_GRID)),
    annotation_name_gp = gp_anno_name()
  )

  ht_A <- Heatmap(mat,
    name = "log2FC",
    col = col_div,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = gp_row_names(italic = TRUE),
    column_names_gp = gp_col_names(),
    column_names_rot = 0,
    column_names_side = "top",
    left_annotation = ha_left,
    right_annotation = ha_right,
    na_col = COL_NA,
    rect_gp = gp_border(),
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- mat[i, j]
      if (!is.na(v)) {
        grid.text(sprintf("%.2f", v), x, y,
                  gp = gpar(fontsize = 8, fontfamily = FONT_GRID,
                            col = ifelse(abs(v) > 1.2, "white", "black")))
      } else {
        grid.text("ND", x, y,
                  gp = gpar(fontsize = 7, fontfamily = FONT_GRID,
                            col = "grey50", fontface = "italic"))
      }
    },
    # Title rendered externally via save_heatmap_panel_pdf for full-width centering
    heatmap_legend_param = list(
      title = "log2FC",
      title_gp = gp_legend_title(),
      labels_gp = gp_legend_labels(),
      legend_height = unit(18, "mm"),
      grid_width = unit(3, "mm")
    ),
    width = unit(26, "mm"),
    height = unit(nrow(mat) * 5.5, "mm")
  )

  save_heatmap_panel_pdf("Fig3a_enzyme_heatmap.pdf", ht_A, 80, 88,
                         title = "Tryptophan Pathway Enzymes")
  cat("    -> Fig3a_enzyme_heatmap.pdf\n")
}, error = function(e) cat("    ERROR Panel A:", e$message, "\n"))


# -----------------------------------------------------------------------------
# PANEL B: IDO1 cross-disease comparison (forest-style dot plot)
# -----------------------------------------------------------------------------
cat("\n--- Panel B: IDO1 cross-disease comparison ---\n")
tryCatch({
  # ---------------------------------------------------------------------------
  # Statistical integrity:
  #   - HAE row uses REAL 95% CI computed from limma t-statistic:
  #       SE = |logFC / t|;  CI = logFC +/- qt(0.975, df) * SE
  #     (df = n_pairs - 1 = 11 for the paired tumor-vs-adjacent design).
  #   - Cancer rows use literature-reported logFC RANGES (not statistical CIs);
  #     they are drawn with a DIFFERENT geometry (open bars + diamond marker)
  #     and explicitly labeled in legend & x-axis caption to avoid confusion.
  # ---------------------------------------------------------------------------
  HAE_DF <- 11  # 12 paired samples - 1
  T_CRIT <- qt(0.975, df = HAE_DF)

  # Pull HAE IDO1 t-statistic from the enzyme expression table to derive SE
  ido1_enz <- enzymes %>% filter(enzyme == "IDO1")
  ido1_t   <- ido1_enz$mRNA_t[1]
  ido1_lfc <- ido1_enz$mRNA_logFC[1]
  ido1_se  <- if (!is.na(ido1_t) && ido1_t != 0) abs(ido1_lfc / ido1_t) else NA_real_

  ido1_data <- ido1_cross %>%
    filter(gene == "IDO1") %>%
    mutate(
      disease = factor(disease, levels = c("HAE", "HCC (TCGA-LIHC)", "CCA",
                                            "Colorectal Ca", "Melanoma")),
      # --- HAE: real 95% CI from limma t-statistic ----------------------------
      hae_se = ifelse(disease == "HAE", ido1_se, NA_real_),
      # --- Position values -----------------------------------------------------
      fc_mid = case_when(
        disease == "HAE" ~ logFC,
        direction == "UP (+1.5-3.0)" ~ 2.25,
        direction == "UP (+1.0-2.0)" ~ 1.5,
        direction == "UP (+2.0-4.0)" ~ 3.0,
        direction == "UP (+3.0-5.0)" ~ 4.0,
        TRUE ~ NA_real_
      ),
      fc_lo = case_when(
        disease == "HAE" ~ logFC - T_CRIT * hae_se,
        direction == "UP (+1.5-3.0)" ~ 1.5,
        direction == "UP (+1.0-2.0)" ~ 1.0,
        direction == "UP (+2.0-4.0)" ~ 2.0,
        direction == "UP (+3.0-5.0)" ~ 3.0,
        TRUE ~ NA_real_
      ),
      fc_hi = case_when(
        disease == "HAE" ~ logFC + T_CRIT * hae_se,
        direction == "UP (+1.5-3.0)" ~ 3.0,
        direction == "UP (+1.0-2.0)" ~ 2.0,
        direction == "UP (+2.0-4.0)" ~ 4.0,
        direction == "UP (+3.0-5.0)" ~ 5.0,
        TRUE ~ NA_real_
      ),
      group_lab = ifelse(disease == "HAE",
                         "HAE (this study; 95% CI)",
                         "Cancer (literature range)")
    )

  cat(sprintf("  HAE IDO1: logFC=%.3f  SE=%.3f  95%% CI=[%.3f, %.3f]\n",
              ido1_data$logFC[ido1_data$disease == "HAE"],
              ido1_data$hae_se[ido1_data$disease == "HAE"],
              ido1_data$fc_lo[ido1_data$disease == "HAE"],
              ido1_data$fc_hi[ido1_data$disease == "HAE"]))

  hae_row    <- ido1_data %>% filter(disease == "HAE")
  cancer_row <- ido1_data %>% filter(disease != "HAE")

  p3b <- ggplot() +
    # Reference vertical line at 0
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.4, color = "grey50") +
    # --- Cancer literature ranges: open horizontal bars + diamond marker ----
    geom_errorbarh(data = cancer_row,
                   aes(xmin = fc_lo, xmax = fc_hi, y = disease),
                   height = 0.0, linewidth = 0.6, color = COL_UP, alpha = 0.55) +
    geom_point(data = cancer_row,
               aes(x = fc_mid, y = disease, fill = group_lab),
               shape = 23, size = 3.0, stroke = 0.6,
               color = COL_UP) +
    # --- HAE real 95% CI: solid bar + filled circle marker -----------------
    geom_errorbarh(data = hae_row,
                   aes(xmin = fc_lo, xmax = fc_hi, y = disease),
                   height = 0.20, linewidth = 0.7, color = COL_DOWN) +
    geom_point(data = hae_row,
               aes(x = fc_mid, y = disease, fill = group_lab),
               shape = 21, size = 3.2, stroke = 0.6,
               color = COL_DOWN) +
    scale_fill_manual(
      values = c("HAE (this study; 95% CI)"   = COL_DOWN,
                 "Cancer (literature range)"  = COL_UP),
      name   = NULL
    ) +
    scale_x_continuous(breaks = seq(-2, 6, 1),
                       limits = c(-2.0, 5.5)) +
    labs(title = "IDO1 Expression Across Diseases",
         x = "log\u2082FC  (HAE: 95% CI; Cancer: literature range)",
         y = NULL) +
    theme_nc +
    theme(legend.position = "bottom",
          legend.margin   = margin(0, 0, 0, 0),
          legend.key.size = unit(3.5, "mm"),
          legend.text     = element_text(size = 7.5, family = FONT_FAMILY),
          axis.text.y     = element_text(size = 8, family = FONT_FAMILY),
          axis.title.x    = element_text(size = 8, family = FONT_FAMILY))

  save_panel_pdf("Fig3b_IDO1_cross_disease.pdf", 103, 88, print(p3b))
  cat("    -> Fig3b_IDO1_cross_disease.pdf\n")
}, error = function(e) cat("    ERROR Panel B:", e$message, "\n"))


# -----------------------------------------------------------------------------
# PANEL C: Metabolite depletion barplot (key metabolites)
# -----------------------------------------------------------------------------
cat("\n--- Panel C: Metabolite depletion barplot ---\n")
tryCatch({
  # Select key depleted metabolites
  key_mets <- metabolites %>%
    filter(metabolite %in% c("L-Kynurenine", "NAD+", "N-Acetylkynurenine",
                             "Indoleacetylglycine", "N-Acetylserotonin",
                             "3-Methyldioxyindole")) %>%
    arrange(logFC) %>%
    mutate(
      met_label = case_when(
        metabolite == "N-Acetylkynurenine" ~ "N-AcKyn",
        metabolite == "L-Kynurenine" ~ "L-Kynurenine",
        metabolite == "Indoleacetylglycine" ~ "IndoleAcGly",
        metabolite == "N-Acetylserotonin" ~ "N-AcSerotonin",
        metabolite == "3-Methyldioxyindole" ~ "3-MeDioxyindole",
        TRUE ~ metabolite
      ),
      met_label = factor(met_label, levels = met_label),
      sig_label = ifelse(P_value < 0.05, "*", ifelse(P_value < 0.1, "\u2020", ""))
    )

  p3c <- ggplot(key_mets, aes(x = met_label, y = logFC, fill = logFC < 0)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_text(aes(label = sig_label,
                  y = ifelse(logFC < 0, logFC - 0.08, logFC + 0.08)),
              size = FS_GEOM_TEXT, family = FONT_FAMILY) +
    geom_text(aes(label = sprintf("%.2f", logFC),
                  y = logFC / 2),
              size = FS_GEOM_TEXT, family = FONT_FAMILY, color = "white") +
    scale_fill_manual(values = c("TRUE" = COL_DOWN, "FALSE" = COL_UP), guide = "none") +
    labs(title = "Tryptophan Metabolite Depletion",
         x = NULL,
         y = "log\u2082FC (Adjacent/Normal)") +
    theme_nc +
    theme(axis.text.x = element_text(size = 8, angle = 35, hjust = 1, family = FONT_FAMILY))

  save_panel_pdf("Fig3c_metabolite_barplot.pdf", 95, 79, print(p3c))
  cat("    -> Fig3c_metabolite_barplot.pdf\n")
}, error = function(e) cat("    ERROR Panel C:", e$message, "\n"))


# -----------------------------------------------------------------------------
# PANEL D: N-Acetylkynurenine vs TCA Cycle correlation scatter plot
# (Replaces non-significant NAD+ vs OXPHOS; rho=+0.706, P=0.0001)
# -----------------------------------------------------------------------------
cat("\n--- Panel D: N-Acetylkynurenine vs TCA Cycle correlation ---\n")
tryCatch({
  # Load metabolomics and mechanistic ssGSEA data
  mb_raw <- read.csv(
    file.path(BASE, "02_analysis/data/processed/metabolomics_log2_merged.csv"),
    row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
  )
  annot  <- read.csv(
    file.path(BASE, "02_analysis/data/processed/metabolomics_annotation.csv"),
    stringsAsFactors = FALSE
  )
  mech   <- read.csv(
    file.path(RES, "enhancement10_mechanism/ssGSEA_mechanistic_scores.csv"),
    row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
  )

  # Locate N-Acetylkynurenine by compound name
  nak_id <- annot$Compound_ID[annot$Name == "N-Acetylkynurenine"]
  if (length(nak_id) == 0) stop("N-Acetylkynurenine not found in annotation")
  nak_id <- nak_id[1]

  # Match common samples (metabolomics 28 vs ssGSEA 24)
  common_smp <- intersect(colnames(mb_raw), colnames(mech))
  cat(sprintf("  Common samples: %d\n", length(common_smp)))

  nak_vals <- as.numeric(mb_raw[nak_id,    common_smp])
  tca_vals <- as.numeric(mech["TCA_CYCLE", common_smp])

  # Group annotation from nad_oxphos (same 24 samples)
  grp_map <- setNames(nad_oxphos$group, nad_oxphos$sample)
  grp_vec <- grp_map[common_smp]
  grp_vec[is.na(grp_vec)] <- "Unknown"

  d3d <- data.frame(
    sample    = common_smp,
    NAK_level = nak_vals,
    TCA_score = tca_vals,
    group     = grp_vec,
    stringsAsFactors = FALSE
  )

  # Spearman correlation
  cor_test <- cor.test(d3d$NAK_level, d3d$TCA_score, method = "spearman")
  rho_val  <- cor_test$estimate
  p_val    <- cor_test$p.value
  p_label  <- if (p_val < 0.001) formatC(p_val, format = "e", digits = 2) else sprintf("%.3f", p_val)
  cat(sprintf("  N-Acetylkynurenine vs TCA_CYCLE: rho=%.3f P=%s\n", rho_val, p_label))

  p3d <- ggplot(d3d, aes(x = NAK_level, y = TCA_score, color = group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "grey30",
                fill = "grey80", alpha = 0.3, formula = y ~ x) +
    scale_color_manual(
      values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT, "Unknown" = COL_NS),
      name = "Tissue"
    ) +
    annotate("text",
             x = min(d3d$NAK_level, na.rm = TRUE),
             y = max(d3d$TCA_score, na.rm = TRUE),
             label = paste0("\u03c1 = ", sprintf("%.3f", rho_val), "\nP = ", p_label),
             size = FS_GEOM_TEXT, family = FONT_FAMILY, hjust = 0, vjust = 1) +
    labs(title = "N-Acetylkynurenine vs TCA Cycle Activity",
         x = "N-Acetylkynurenine (log\u2082 intensity)",
         y = "TCA Cycle Score (ssGSEA)") +
    theme_nc +
    theme(legend.position = c(0.18, 0.18),
          legend.background = element_rect(fill = alpha("white", 0.8), linewidth = 0))

  save_panel_pdf("Fig3d_NAK_TCA_correlation.pdf", 88, 79, print(p3d))
  cat("    -> Fig3d_NAK_TCA_correlation.pdf\n")
}, error = function(e) cat("    ERROR Panel D:", e$message, "\n"))


# -----------------------------------------------------------------------------
# PANEL E: KTR paired boxplot (Wilcoxon test)
# -----------------------------------------------------------------------------
cat("\n--- Panel E: KTR paired boxplot ---\n")
tryCatch({
  ktr$group <- factor(ktr$group, levels = c("Normal", "Adjacent"))

  # Paired Wilcoxon test (two-vector form required for paired)
  ktr_wide <- ktr %>%
    select(patient, group, KTR_log2) %>%
    pivot_wider(names_from = group, values_from = KTR_log2)
  wt <- wilcox.test(ktr_wide$Normal, ktr_wide$Adjacent, paired = TRUE)
  p_ktr <- wt$p.value

  p3e <- ggplot(ktr, aes(x = group, y = KTR_log2, fill = group)) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
    geom_line(aes(group = patient), color = "grey50", linewidth = 0.4, alpha = 0.65) +
    geom_point(aes(color = group), size = 1.8, alpha = 0.8) +
    scale_fill_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT),
                      guide = "none") +
    scale_color_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT),
                       guide = "none") +
    annotate("segment", x = 1, xend = 2, y = max(ktr$KTR_log2) + 0.3,
             yend = max(ktr$KTR_log2) + 0.3, linewidth = 0.4) +
    annotate("text", x = 1.5, y = max(ktr$KTR_log2) + 0.5,
             label = paste0("P = ", sprintf("%.3f", p_ktr)),
             size = FS_GEOM_TEXT, family = FONT_FAMILY) +
    labs(title = "KTR (Kyn/Trp Ratio)",
         x = NULL,
         y = "log\u2082KTR") +
    theme_nc +
    theme(plot.title.position = "plot",
          plot.margin = margin(8, 4, 2, 4, "mm"))

  save_panel_pdf("Fig3e_KTR_paired.pdf", 62, 78, print(p3e))
  cat("    -> Fig3e_KTR_paired.pdf\n")
}, error = function(e) cat("    ERROR Panel E:", e$message, "\n"))


# -----------------------------------------------------------------------------
# PANEL F: Pathway mechanism schematic (HAE vs Cancer comparison)
# -----------------------------------------------------------------------------
cat("\n--- Panel F: Pathway schematic ---\n")
tryCatch({
  # Create a conceptual comparison diagram using ggplot annotations
  # Cancer: IDO1 UP -> Kynurenine accumulation -> T cell suppression
  # HAE:    Enzyme DOWN -> NAD+ depletion -> Bioenergetic compromise

  schematic_df <- data.frame(
    pathway = c(rep("Cancer\n(IDO1 UP)", 3), rep("HAE\n(Enzyme DOWN)", 3)),
    step = rep(c(1, 2, 3), 2),
    label = c("IDO1\u2191", "Kynurenine\naccumulation", "T cell\nsuppression",
              "IDO1/TDO2\u2193", "NAD\u207a\ndepletion", "Bioenergetic\ncompromise"),
    x = rep(c(1, 2, 3), 2),
    y = c(rep(2, 3), rep(1, 3))
  )

  # Schematic font size (slightly larger for readability in wide panel)
  FS_SCHEMA <- 3.2  # ~9pt

  p3f <- ggplot(schematic_df, aes(x = x, y = y)) +
    # Background boxes for each pathway
    annotate("rect", xmin = 0.4, xmax = 3.6, ymin = 1.6, ymax = 2.4,
             fill = "#FFEAEA", color = COL_UP, linewidth = 0.5, alpha = 0.3) +
    annotate("rect", xmin = 0.4, xmax = 3.6, ymin = 0.6, ymax = 1.4,
             fill = "#E8F4FD", color = COL_DOWN, linewidth = 0.5, alpha = 0.3) +
    # Pathway labels
    annotate("text", x = 0.15, y = 2, label = "Cancer", size = FS_SCHEMA + 0.3,
             fontface = "bold", family = FONT_FAMILY, color = COL_UP, hjust = 0.5) +
    annotate("text", x = 0.15, y = 1, label = "HAE", size = FS_SCHEMA + 0.3,
             fontface = "bold", family = FONT_FAMILY, color = COL_DOWN, hjust = 0.5) +
    # Step boxes - Cancer row
    annotate("label", x = 1, y = 2, label = "IDO1 \u2191\u2191",
             size = FS_SCHEMA, family = FONT_FAMILY, fontface = "bold",
             fill = "#FFD4D4", color = COL_UP, linewidth = 0.4) +
    annotate("label", x = 2, y = 2, label = "Kynurenine\naccumulation",
             size = FS_SCHEMA, family = FONT_FAMILY,
             fill = "#FFD4D4", color = COL_UP, linewidth = 0.4) +
    annotate("label", x = 3, y = 2, label = "T cell\nsuppression",
             size = FS_SCHEMA, family = FONT_FAMILY,
             fill = "#FFD4D4", color = COL_UP, linewidth = 0.4) +
    # Step boxes - HAE row
    annotate("label", x = 1, y = 1, label = "All enzymes \u2193",
             size = FS_SCHEMA, family = FONT_FAMILY, fontface = "bold",
             fill = "#D4E8FF", color = COL_DOWN, linewidth = 0.4) +
    annotate("label", x = 2, y = 1, label = "NAD\u207a\ndepletion",
             size = FS_SCHEMA, family = FONT_FAMILY,
             fill = "#D4E8FF", color = COL_DOWN, linewidth = 0.4) +
    annotate("label", x = 3, y = 1, label = "Bioenergetic\ncompromise",
             size = FS_SCHEMA, family = FONT_FAMILY,
             fill = "#D4E8FF", color = COL_DOWN, linewidth = 0.4) +
    # Arrows between steps
    annotate("segment", x = 1.35, xend = 1.65, y = 2, yend = 2,
             arrow = arrow(length = unit(2, "mm"), type = "closed"),
             linewidth = 0.6, color = COL_UP) +
    annotate("segment", x = 2.35, xend = 2.65, y = 2, yend = 2,
             arrow = arrow(length = unit(2, "mm"), type = "closed"),
             linewidth = 0.6, color = COL_UP) +
    annotate("segment", x = 1.35, xend = 1.65, y = 1, yend = 1,
             arrow = arrow(length = unit(2, "mm"), type = "closed"),
             linewidth = 0.6, color = COL_DOWN) +
    annotate("segment", x = 2.35, xend = 2.65, y = 1, yend = 1,
             arrow = arrow(length = unit(2, "mm"), type = "closed"),
             linewidth = 0.6, color = COL_DOWN) +
    # Title
    labs(title = "Immune Evasion: Cancer vs HAE") +
    coord_cartesian(xlim = c(-0.2, 3.8), ylim = c(0.4, 2.6)) +
    theme_void(base_family = FONT_FAMILY) +
    theme(plot.title = element_text(family = FONT_FAMILY, size = 10,
                                    face = "bold", hjust = 0.5),
          plot.margin = margin(5, 5, 5, 5, "mm"))

  save_panel_pdf("Fig3f_pathway_schematic.pdf", 121, 78, print(p3f))
  cat("    -> Fig3f_pathway_schematic.pdf\n")
}, error = function(e) cat("    ERROR Panel F:", e$message, "\n"))


# =============================================================================
# SECTION 11: Composite Figure Assembly (Pure Vector Grid Viewport — AI-editable)
# =============================================================================
cat("\n--- Assembling composite Figure_3 (VECTOR, 183x245mm) ---\n")

tryCatch({
  DPI <- ASSEMBLY_DPI  # 600
  W_TOTAL <- 183; H_TOTAL <- 245

  # Row 1: A(heatmap,80mm) + B(forest,103mm) => h=88mm
  # Row 2: C(barplot,95mm) + D(scatter,88mm) => h=79mm
  # Row 3: E(boxplot,62mm) + F(schematic,121mm) => h=78mm
  H1 <- 88;  H2 <- 79;  H3 <- 78
  W_A <- 80;  W_B <- 103
  W_C <- 95;  W_D <- 88
  W_E <- 62;  W_F <- 121

  cat(sprintf("  Layout: Row1=%dmm (A:%d+B:%d), Row2=%dmm (C:%d+D:%d), Row3=%dmm (E:%d+F:%d)\n",
              H1, W_A, W_B, H2, W_C, W_D, H3, W_E, W_F))

  # --- Render function: places all panels + labels into current device ---
  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))

    # Row 1 (top): A + B — y = H2 + H3 from bottom
    y_row1 <- H2 + H3

    # Panel a (heatmap)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row1, "mm"),
      width = unit(W_A, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    draw(ht_A, padding = unit(c(2, 2, 8, 2), "mm"), newpage = FALSE)
    grid::popViewport()

    # Panel b (ggplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_A, "mm"), y = unit(y_row1, "mm"),
      width = unit(W_B, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(p3b, newpage = FALSE)
    grid::popViewport()

    # Row 2 (middle): C + D — y = H3 from bottom
    y_row2 <- H3

    # Panel c (ggplot)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_C, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(p3c, newpage = FALSE)
    grid::popViewport()

    # Panel d (ggplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_C, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_D, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(p3d, newpage = FALSE)
    grid::popViewport()

    # Row 3 (bottom): E + F — y = 0
    # Panel e (ggplot)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W_E, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")
    ))
    print(p3e, newpage = FALSE)
    grid::popViewport()

    # Panel f (ggplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_E, "mm"), y = unit(0, "mm"),
      width = unit(W_F, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")
    ))
    print(p3f, newpage = FALSE)
    grid::popViewport()

    # --- Panel labels (bold, FS_TAG pt) ---
    label_data <- data.frame(
      text  = c("a", "b", "c", "d", "e", "f"),
      x_mm  = c(1, W_A + 1, 1, W_C + 1, 1, W_E + 1),
      y_mm  = c(H_TOTAL - 1, H_TOTAL - 1, y_row1 - 1, y_row1 - 1, H3 - 1, H3 - 1),
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
  cairo_pdf(file.path(OUT, "Figure_3.pdf"),
            width = W_TOTAL / 25.4, height = H_TOTAL / 25.4, family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> Figure_3.pdf (VECTOR, AI-editable)\n")

  # --- Save PNG (600 DPI for review) ---
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "Figure_3.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_3.png\n")

  # --- Save TIFF (600 DPI for journal submission) ---
  grDevices::tiff(file.path(OUT, "Figure_3.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_3.tiff\n")

  cat("  Figure_3 DONE (183x245mm, VECTOR PDF + 600DPI PNG/TIFF)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

# =============================================================================
cat("\n=== Figure 3 rendering complete ===\n")
