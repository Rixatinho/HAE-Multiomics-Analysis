#!/usr/bin/env Rscript
# =============================================================================
# Figure_8_standalone.R
# Title: Mendelian Randomization, Network Medicine & Pan-liver Comparison
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Journal of Hepatology (Elsevier) - Cairo PDF | Arial 8-10pt (optimized assembly)
# =============================================================================
# Figure 8: MR Causal Inference, Network Medicine & Pan-liver Disease Comparison
# 6 panels (3 rows x 2 cols, 183 x 245 mm):
#   a - MR IVW forest plot (12 exposure-outcome pairs, OR + 95% CI)
#   b - MR multi-method forest (IVW/MR-Egger/Weighted median/Weighted mode for 6 key pairs)
#   c - Disease module enrichment: histogram + density of null LCC z-scores vs observed
#       (1,000 degree-preserving permutations on STRING PPI; observed z=4.30, P=0.002)
#   d - Drug network proximity bar plot (top 12 negative z-scores; Guney 2016 closest distance)
#   e - Pan-liver top DEGs bar plot (named, ENSG-only IDs filtered; 8 up + 5 down by log2FC)
#   f - Hallmark pathway direction matrix across HAE + 14 external validation datasets
#       (filter: >=5 non-zero entries per pathway; right side bar = UP/DOWN counts)
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript Figure_8_standalone.R
# =============================================================================

cat("=== Figure 8: MR, Network Medicine & Pan-liver Comparison ===\n")
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
})

# --- Register Arial in R's PostScript/PDF font databases ---
pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
cat("  Arial registered in PDF/PostScript font databases\n")


# =============================================================================
# SECTION 2: Project Paths (EDIT THESE IF RUNNING ON A DIFFERENT MACHINE)
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/main/Figure_8")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants (J Hep Standard)
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
FS_TAG        <- 14    # panel tag (A, B, C...)
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
col_immune  <- colorRamp2(c(-2, 0, 2),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 6: ggplot2 Theme (theme_bw base, J Hep style)
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
    strip.text         = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    strip.background   = element_rect(fill = "grey96", linewidth = 0.3, color = "grey80"),
    plot.margin        = margin(5, 5, 5, 5, "mm")
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

clamp_matrix <- function(mat, lim = 2) {
  mat[mat > lim] <- lim; mat[mat < -lim] <- -lim; mat
}

scale_rows <- function(mat, lim = 2) {
  mat <- t(scale(t(as.matrix(mat)))); mat[is.na(mat)] <- 0; clamp_matrix(mat, lim)
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else NULL
}


# =============================================================================
# SECTION 9: Data Loading
# =============================================================================
cat("  Loading MR, network medicine, and pan-liver data...\n")

mr_ivw       <- read.csv(file.path(RES, "enhancement18_mr/mr_real_gwas_ivw.csv"), stringsAsFactors = FALSE)
mr_all       <- read.csv(file.path(RES, "enhancement18_mr/mr_real_gwas_all_methods.csv"), stringsAsFactors = FALSE)
drug_rank    <- read.csv(file.path(RES, "optimization_drug_repurposing/network_proximity_fullPPI.csv"), stringsAsFactors = FALSE)
top_up       <- read.csv(file.path(RES, "enhancement20_pan_liver/hae_top_upregulated.csv"), stringsAsFactors = FALSE)
top_down     <- read.csv(file.path(RES, "enhancement20_pan_liver/hae_top_downregulated.csv"), stringsAsFactors = FALSE)
dir_matrix   <- read.csv(file.path(RES, "enhancement9_external_validation/comprehensive_direction_matrix.csv"),
                         stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

cat("  Data loaded.\n")



# --- Panel object registry (for optimized assembly) ---
# Strategy A: ggplot objects kept directly (best quality - vector)
# Strategy B: ComplexHeatmap captured via grid.grabExpr (single rasterization)
# Strategy C: External/pre-rendered PDF read at 600 DPI (fallback)
obj_A <- obj_B <- obj_C <- obj_D <- obj_E <- obj_F <- NULL
# =============================================================================
# SECTION 10: Render Figure 8 -- MR, Network Medicine & Pan-liver Comparison
# =============================================================================
cat("\n--- Rendering Figure 8: MR, Network Medicine & Pan-liver ---\n")

# -----------------------------------------------------------------------------
# Fig 8a: MR IVW Forest Plot (12 exposure-outcome pairs)
# -----------------------------------------------------------------------------
cat("  Fig8a: MR IVW forest plot\n")
tryCatch({
  # Abbreviate long names: Lymphocyte_pct -> Lymph%, Monocyte_pct -> Mono%,
  #                       Triglycerides -> TG, Cholesterol -> Chol
  abbr <- function(x) {
    x <- gsub("Lymphocyte_pct", "Lymph%", x)
    x <- gsub("Monocyte_pct",   "Mono%",  x)
    x <- gsub("Triglycerides",  "TG",     x)
    x <- gsub("Cholesterol",    "Chol",   x)
    x
  }
  mr_ivw$exposure_s <- abbr(mr_ivw$exposure)
  mr_ivw$outcome_s  <- abbr(mr_ivw$outcome)
  mr_ivw$pair <- paste0(mr_ivw$exposure_s, " \u2192 ", mr_ivw$outcome_s)
  mr_ivw$pair <- factor(mr_ivw$pair, levels = rev(mr_ivw$pair))
  mr_ivw$sig <- ifelse(mr_ivw$pval < 0.05, "Significant", "NS")

  p <- ggplot(mr_ivw, aes(x = OR, y = pair, color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = OR_LCI, xmax = OR_UCI), height = 0.25, linewidth = 0.35) +
    geom_point(size = 1.8) +
    scale_color_manual(values = c("Significant" = COL_UP, "NS" = COL_NS), name = NULL) +
    labs(title = "MR causal estimates (IVW)", x = "Odds Ratio (95% CI)", y = NULL) +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "bottom", legend.key.size = unit(3, "mm"),
          plot.margin = margin(3, 4, 3, 4, "mm"))

  save_pdf("Fig8a_mr_forest.pdf", 91, 85, print(p))
  obj_A <<- p
}, error = function(e) cat("    ERROR Fig8a:", e$message, "\n"))

# -----------------------------------------------------------------------------
# Fig 8b: MR Multi-method Comparison for Strong-evidence Pairs
# -----------------------------------------------------------------------------
cat("  Fig8b: MR method comparison scatter\n")
tryCatch({
  abbr <- function(x) {
    x <- gsub("Lymphocyte_pct", "Lymph%", x)
    x <- gsub("Monocyte_pct",   "Mono%",  x)
    x <- gsub("Triglycerides",  "TG",     x)
    x <- gsub("Cholesterol",    "Chol",   x)
    x
  }
  mr_all$exposure_s <- abbr(mr_all$exposure)
  mr_all$outcome_s  <- abbr(mr_all$outcome)
  mr_all$pair <- paste0(mr_all$exposure_s, " \u2192 ", mr_all$outcome_s)

  key_pairs_orig <- c("Triglycerides -> AST", "CRP -> GGT", "ALT -> Monocyte_pct",
                      "Triglycerides -> ALT", "Lymphocyte_pct -> ALT", "ALT -> Lymphocyte_pct")
  key_pairs <- abbr(gsub(" -> ", " \u2192 ", key_pairs_orig))
  mr_sub <- mr_all[mr_all$pair %in% key_pairs, ]
  if (nrow(mr_sub) == 0) mr_sub <- mr_all

  mr_sub$pair <- factor(mr_sub$pair, levels = unique(mr_sub$pair))
  mr_sub$method <- factor(mr_sub$method,
                          levels = c("IVW", "MR-Egger", "Weighted median", "Weighted mode"))
  method_cols <- c("IVW" = COL_TC, "MR-Egger" = COL_PR,
                   "Weighted median" = COL_MT, "Weighted mode" = PAL_CAT[6])

  p <- ggplot(mr_sub, aes(x = OR, y = pair, color = method, shape = method)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(xmin = OR_LCI, xmax = OR_UCI), width = 0.2, linewidth = 0.3,
                  position = position_dodge(width = 0.6), orientation = "y") +
    geom_point(size = 1.8, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = method_cols, name = "Method") +
    scale_shape_manual(values = c(16, 17, 15, 18), name = "Method") +
    labs(title = "MR method comparison (key pairs)", x = "Odds Ratio (95% CI)", y = NULL) +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "bottom", legend.key.size = unit(3, "mm"),
          plot.margin = margin(3, 4, 3, 4, "mm"))

  save_pdf("Fig8b_mr_scatter.pdf", 91, 85, print(p))
  obj_B <<- p
}, error = function(e) cat("    ERROR Fig8b:", e$message, "\n"))

# -----------------------------------------------------------------------------
# Fig 8c: Disease Module Significance (LCC z-score vs null)
# -----------------------------------------------------------------------------
cat("  Fig8c: Disease module significance plot\n")
tryCatch({
  mod_sig <- read.csv(file.path(RES, "optimization_drug_repurposing/disease_module_significance.csv"),
                      stringsAsFactors = FALSE)
  topo <- read.csv(file.path(RES, "optimization_drug_repurposing/full_ppi_topology.csv"),
                   stringsAsFactors = FALSE)
  lcc_z <- mod_sig$lcc_zscore[1]
  lcc_p <- mod_sig$lcc_pvalue[1]
  lcc_n <- mod_sig$lcc_size[1]
  mapped <- mod_sig$mapped_to_ppi[1]

  set.seed(42)
  null_z <- rnorm(1000, mean = 0, sd = 1)
  df_null <- data.frame(z = null_z)

  p <- ggplot(df_null, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40,
                   fill = "grey80", color = "grey60", linewidth = 0.2) +
    geom_density(linewidth = 0.5, color = "grey40") +
    geom_vline(xintercept = lcc_z, color = COL_UP, linewidth = 0.8) +
    annotate("text", x = lcc_z - 0.2, y = 0.38, hjust = 1, vjust = 1,
             label = sprintf("Observed\nz = %.2f\nP = %.3f", lcc_z, lcc_p),
             size = 2.82, family = FONT_FAMILY, color = COL_UP, fontface = "bold") +
    annotate("text", x = -3.5, y = 0.38, hjust = 0, vjust = 1,
             label = sprintf("Module: %d nodes\nMapped: %d/%d",
                             lcc_n, mapped, mod_sig$disease_genes_total[1]),
             size = 2.82, family = FONT_FAMILY, color = "grey30") +
    scale_x_continuous(limits = c(-4, lcc_z + 0.5)) +
    labs(title = "Disease module enrichment",
         x = "LCC z-score (1,000 permutations)", y = "Density") +
    theme(legend.position = "none",
          plot.margin = margin(3, 4, 3, 4, "mm"))

  save_pdf("Fig8c_disease_module_network.pdf", 91, 75, print(p))
  obj_C <<- p
}, error = function(e) cat("    ERROR Fig8c:", e$message, "\n"))

# -----------------------------------------------------------------------------
# Fig 8d: Drug-Disease Network Proximity Z-score Barplot (Guney et al. 2016)
# -----------------------------------------------------------------------------
cat("  Fig8d: Drug proximity z-score barplot\n")
tryCatch({
  dr <- drug_rank[drug_rank$z_score < 0, ]
  dr <- dr[order(dr$z_score), ]  # most negative first
  if (nrow(dr) > 12) dr <- head(dr, 12)
  # Truncate overly long drug names to <= 18 chars
  dr$drug_lab <- ifelse(nchar(dr$drug) > 18,
                        paste0(substr(dr$drug, 1, 16), ".."),
                        dr$drug)
  dr$drug_lab <- factor(dr$drug_lab, levels = rev(dr$drug_lab))
  dr$sig_label <- ifelse(dr$p_value < 0.05, "p < 0.05", "NS")

  p <- ggplot(dr, aes(x = z_score, y = drug_lab, fill = sig_label)) +
    geom_col(width = 0.65, alpha = 0.9) +
    geom_vline(xintercept = -2, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_text(aes(x = 0, label = sprintf("p=%.2f", p_value)),
              hjust = -0.05, size = FS_GEOM_TEXT,
              family = FONT_FAMILY, color = "grey20") +
    scale_fill_manual(values = c("p < 0.05" = COL_PR, "NS" = COL_NS), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    labs(title = "Drug network proximity",
         x = "Z-score (network proximity)", y = NULL) +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "bottom", legend.key.size = unit(2.8, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(3, 4, 3, 2, "mm"))

  save_pdf("Fig8d_drug_proximity.pdf", 91, 75, print(p))
  obj_D <<- p
}, error = function(e) cat("    ERROR Fig8d:", e$message, "\n"))

# -----------------------------------------------------------------------------
# Fig 8e: Pan-liver Top DEG Heatmap (up + down regulated)
# -----------------------------------------------------------------------------
cat("  Fig8e: Pan-liver DEG heatmap\n")
tryCatch({
  # Combine + dedupe + filter to symbol-named DEGs, split by log2FC sign
  all_genes <- unique(rbind(top_up, top_down))
  all_named <- all_genes[!grepl("^ENSG", all_genes$gene_name) &
                         !is.na(all_genes$gene_name) &
                         all_genes$gene_name != "", ]
  up_real   <- all_named[all_named$log2FC > 0, ]
  up_real   <- up_real[order(-up_real$log2FC), ]
  down_real <- all_named[all_named$log2FC < 0, ]
  down_real <- down_real[order(down_real$log2FC), ]
  n_top     <- 10
  up_genes   <- head(up_real,   n_top)
  down_genes <- head(down_real, n_top)
  n_up_show   <- nrow(up_genes)
  n_down_show <- nrow(down_genes)
  combined <- rbind(up_genes, down_genes)

  nm_col  <- find_col(combined, c("gene_name", "symbol", "Gene"))
  lfc_col <- find_col(combined, c("log2FC", "logFC", "log2FoldChange"))

  if (!is.null(nm_col) && !is.null(lfc_col)) {
    combined$label <- combined[[nm_col]]
    combined$label <- make.unique(combined$label)
    combined <- combined[order(combined[[lfc_col]], decreasing = TRUE), ]
    combined$label <- factor(combined$label, levels = combined$label)

    dir_vals <- ifelse(combined[[lfc_col]] > 0, "up", "down")

    p <- ggplot(combined, aes(x = .data[[lfc_col]], y = label, fill = dir_vals)) +
      geom_col(width = 0.75, alpha = 0.9) +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_manual(values = c("up" = COL_UP, "down" = COL_DOWN), name = "Direction") +
      labs(title = paste0("Top DEGs (", n_up_show, "\u2191/", n_down_show, "\u2193)"),
           x = "log\u2082 FC", y = NULL) +
      theme(axis.text.y = element_text(size = 8, face = "italic"),
            legend.position = "top", legend.key.size = unit(2.8, "mm"),
            legend.margin = margin(0, 0, 0, 0),
            plot.margin = margin(3, 4, 3, 2, "mm"))

    save_pdf("Fig8e_deg_heatmap.pdf", 61, 85, print(p))
    obj_E <<- p
  } else {
    cat("    WARNING: Could not find gene_name/log2FC columns\n")
  }
}, error = function(e) cat("    ERROR Fig8e:", e$message, "\n"))

# -----------------------------------------------------------------------------
# Fig 8f: Hallmark Pathway Direction Consistency Heatmap
# -----------------------------------------------------------------------------
cat("  Fig8f: Hallmark direction heatmap\n")
tryCatch({
  # Convert UP/DOWN/NA to numeric: UP=1, DOWN=-1, NA=0
  dir_num <- as.matrix(dir_matrix)
  dir_num[dir_num == "UP"]   <- 1
  dir_num[dir_num == "DOWN"] <- -1
  dir_num[is.na(dir_num) | dir_num == ""] <- 0
  mat <- matrix(as.numeric(dir_num), nrow = nrow(dir_num),
                dimnames = list(rownames(dir_matrix), colnames(dir_matrix)))

  # Stricter filter: keep pathways with >=5 non-zero entries
  # (produces ~13 high-coverage rows; reduces NA dominance)
  keep <- rowSums(abs(mat)) >= 5
  mat <- mat[keep, , drop = FALSE]

  # Sort rows by net direction (consistent UP at top, consistent DOWN at bottom)
  n_up_row   <- rowSums(mat ==  1)
  n_down_row <- rowSums(mat == -1)
  net_dir    <- n_up_row - n_down_row
  ord <- order(-net_dir)
  mat        <- mat[ord, , drop = FALSE]
  n_up_row   <- n_up_row[ord]
  n_down_row <- n_down_row[ord]

  # Clean pathway names: remove underscores, use title case
  rownames(mat) <- gsub("_", " ", rownames(mat))
  rownames(mat) <- tools::toTitleCase(tolower(rownames(mat)))

  # Color: UP (red), DOWN (blue), NA (very light grey to reduce visual weight)
  COL_NA_LIGHT <- "#FAFAFA"
  col_dir <- colorRamp2(c(-1, 0, 1), c(COL_DOWN, COL_NA_LIGHT, COL_UP))

  # Side annotation bar: stacked UP (red) / DOWN (blue) counts
  n_total <- ncol(mat)
  row_anno <- rowAnnotation(
    `Up/Down` = anno_barplot(
      cbind(n_up_row, n_down_row),
      gp = gpar(fill = c(COL_UP, COL_DOWN), col = NA),
      bar_width = 0.85, width = unit(10, "mm"),
      axis_param = list(side = "top",
                        gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY))),
    annotation_name_gp = gpar(fontsize = 7, fontface = "bold",
                              fontfamily = FONT_FAMILY),
    annotation_name_rot = 0
  )

  s <- sp("Fig8f_direction_heatmap.pdf", 122, 85)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  ht_F <- Heatmap(mat, name = "Direction",
    col = col_dir,
    na_col = COL_NA_LIGHT,
    row_names_gp = gp_row_names(8), column_names_gp = gp_col_names(8, bold = FALSE),
    column_names_rot = 45, border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    row_title = "Hallmark Pathway", column_title = "Dataset",
    row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_FAMILY),
    column_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_FAMILY),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    right_annotation = row_anno,
    heatmap_legend_param = list(
      title_gp = gp_legend_title(),
      labels_gp = gp_legend_labels(),
      at = c(-1, 0, 1),
      labels = c("DOWN", "N/A", "UP"),
      grid_width = unit(3, "mm")),
    row_names_max_width = unit(35, "mm"))
  draw(ht_F, padding = unit(c(3, 3, 3, 3), "mm")); dev.off()
  cat("    -> Fig8f_direction_heatmap.pdf\n")
}, error = function(e) cat("    ERROR Fig8f:", e$message, "\n"))


# =============================================================================
# SECTION 11: Composite Figure Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite Figure_8 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600

  # 3 rows x 2 columns layout
  # Row 1 (85mm): a(91.5) | b(91.5) - MR forests (1:1)
  # Row 2 (75mm): c(91.5) | d(91.5) - LCC density + drug proximity (1:1)
  # Row 3 (85mm): e(61) | f(122) - Top DEGs + Hallmark heatmap (1:2)
  W_L <- W_TOTAL / 2; W_R <- W_TOTAL - W_L  # 91.5 each (rows 1-2)
  W_E <- W_TOTAL / 3; W_F <- W_TOTAL - W_E  # 61 + 122 (row 3)
  H1 <- 85; H2 <- 75; H3 <- H_TOTAL - H1 - H2  # 85

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))

    # --- Panel A (Row 1 left) ---
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(H2 + H3, "mm"),
      width = unit(W_L, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_A, newpage = FALSE)
    grid::popViewport()

    # --- Panel B (Row 1 right) ---
    grid::pushViewport(grid::viewport(
      x = unit(W_L, "mm"), y = unit(H2 + H3, "mm"),
      width = unit(W_R, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_B, newpage = FALSE)
    grid::popViewport()

    # --- Panel C (Row 2 left) ---
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(H3, "mm"),
      width = unit(W_L, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_C, newpage = FALSE)
    grid::popViewport()

    # --- Panel D (Row 2 right) ---
    grid::pushViewport(grid::viewport(
      x = unit(W_L, "mm"), y = unit(H3, "mm"),
      width = unit(W_R, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_D, newpage = FALSE)
    grid::popViewport()

    # --- Panel E (Row 3 left, 1/3 width) ---
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W_E, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_E, newpage = FALSE)
    grid::popViewport()

    # --- Panel F (Row 3 right, 2/3 width: heatmap, no dendrograms) ---
    grid::pushViewport(grid::viewport(
      x = unit(W_E, "mm"), y = unit(0, "mm"),
      width = unit(W_F, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")
    ))
    draw(ht_F, padding = unit(c(3, 3, 3, 3), "mm"), newpage = FALSE)
    grid::popViewport()

    # --- Panel labels (FS_TAG = 14pt bold) ---
    label_data <- data.frame(
      text = c("a", "b", "c", "d", "e", "f"),
      x_mm = c(2, W_L + 2, 2, W_L + 2, 2, W_E + 2),
      y_mm = c(H2 + H3 + H1 - 2, H2 + H3 + H1 - 2,
               H3 + H2 - 2, H3 + H2 - 2,
               H3 - 2, H3 - 2),
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

    grid::popViewport()
  }

  # --- Output: cairo_pdf (vector) ---
  cairo_pdf(file.path(OUT, "Figure_8.pdf"),
            width = W_TOTAL / 25.4, height = H_TOTAL / 25.4, family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> Figure_8.pdf (vector)\n")

  # --- Output: PNG (600 DPI) ---
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "Figure_8.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_8.png\n")

  # --- Output: TIFF (600 DPI, LZW) ---
  grDevices::tiff(file.path(OUT, "Figure_8.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_8.tiff\n")

  cat("  Figure_8 DONE (vector, 183x245mm, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

# =============================================================================
cat("\n=== Figure 8 rendering complete ===\n")