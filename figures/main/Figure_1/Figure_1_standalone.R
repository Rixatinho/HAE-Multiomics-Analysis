#!/usr/bin/env Rscript
# =============================================================================
# Figure_1_standalone.R  (v2 — Plan B Restructure)
# Complete, Self-Contained Code for Figure 1
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: EBioMedicine (Lancet family) — Pure Vector PDF + PNG + TIFF
# =============================================================================
# Figure 1 panels (6 total):
#   a: Workflow diagram (AI-generated, not coded here)
#   b: Multi-omics differential overview (diverging bar chart)
#   c: mRNA-Protein discordance scatter (log2FC concordance)
#   d: Top integrated features heatmap (DEGs + DEPs, Z-scored)
#   e: Cross-omics Hallmark pathway heatmap (enhanced, with concordance)
#   f: Metabolic suppression convergence (lollipop chart)
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript Figure_1_standalone.R
# =============================================================================

cat("=== Figure 1: Multi-omics Molecular Landscape of HAE (v2) ===\n")
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
  library(ggrepel)
})

# --- Register Arial ---
pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
ht_opt$TITLE_PADDING <- unit(c(4, 4), "pt")
ht_opt$message <- FALSE

# =============================================================================
# SECTION 2: Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/main/Figure_1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
DPI <- 600
W_TOTAL <- 183  # EBioMedicine double-column width (mm)
H_TOTAL <- 222  # Total height for b-f composite

# Content-adaptive row heights (mm)
H_R1 <- 72   # Row 1: overview + scatter
H_R2 <- 76   # Row 2: pathway heatmap (full-width)
H_R3 <- 74   # Row 3: integrated heatmap + metabolic bar

# Content-adaptive column widths (mm)
W_B <- 85; W_C <- W_TOTAL - W_B  # Row 1: 85 + 98
W_E <- 100; W_F <- W_TOTAL - W_E  # Row 3: 100 + 83

mm2in <- function(mm) mm / 25.4

# Typography (EBioMedicine minimum 8pt)
FS_TITLE    <- 10
FS_SUBTITLE <- 9
FS_AXIS_T   <- 9
FS_AXIS     <- 8
FS_LEGEND_T <- 8
FS_LEGEND   <- 8
FS_TAG      <- 14
FS_GEOM     <- 2.85  # geom_text size (8.1pt)

# =============================================================================
# SECTION 4: Color Palette (ggsci JCO — consistent with Figures 2-6)
# =============================================================================
# Primary palette (JCO: CD534C = red, 0073C2 = navy)
COL_UP       <- "#CD534CFF"    # Warm red (upregulated)
COL_DOWN     <- "#0073C2FF"    # Steel blue (downregulated)
COL_NS       <- "#868686FF"    # Darker grey (not significant)
COL_CONCORD  <- "#55A868"    # Sage green (concordant)
COL_DISCORD  <- "#DD8452"    # Amber (discordant)

# Omics layer accent colors
COL_TC <- "#0073C2FF"    # Transcriptomics (navy)
COL_PR <- "#CD534CFF"    # Proteomics (red)
COL_MT <- "#EFC000FF"    # Metabolomics (gold)

# Group colors
COL_NORMAL   <- "#7AA6DCFF"
COL_ADJACENT <- "#CD534CFF"

# Diverging scale
COL_DIVERGE <- c("#0073C2FF", "#FFFFFF", "#CD534CFF")

# =============================================================================
# SECTION 5: Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text               = element_text(family = FONT_FAMILY, size = 8),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.line          = element_blank(),
    axis.ticks         = element_line(linewidth = 0.35, color = "black"),
    axis.ticks.length  = unit(1.2, "mm"),
    axis.text          = element_text(family = FONT_FAMILY, size = 8, color = "black"),
    axis.title         = element_text(family = FONT_FAMILY, size = 9, face = "bold"),
    plot.title         = element_text(family = FONT_FAMILY, size = 10, face = "bold", hjust = 0),
    plot.title.position = "plot",
    legend.text        = element_text(family = FONT_FAMILY, size = 8),
    legend.title       = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    legend.key.size    = unit(3.5, "mm"),
    legend.background  = element_blank(),
    legend.box.background = element_blank(),
    strip.text         = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    strip.background   = element_rect(fill = "grey96", linewidth = 0.4),
    plot.margin        = margin(4, 4, 4, 4, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# SECTION 6: Helper functions
# =============================================================================
save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  w_px <- round(w_mm * DPI / 25.4)
  h_px <- round(h_mm * DPI / 25.4)
  cairo_pdf(fpath, width = w_px / DPI, height = h_px / DPI, family = FONT_FAMILY)
  force(expr)
  dev.off()
  cat(sprintf("  [OK] %s (%d x %d px)\n", basename(fpath), w_px, h_px))
}

gp_row <- function(sz = 8, italic = FALSE) {
  gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = if (italic) "italic" else "plain")
}
gp_col <- function(sz = 8, bold = TRUE) {
  gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = if (bold) "bold" else "plain")
}
gp_title <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
gp_legend_t <- function() gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold")
gp_legend_l <- function() gpar(fontsize = 8, fontfamily = FONT_GRID)

# --- ENSG -> Gene Symbol mapping ---
build_gene_map <- function() {
  map <- data.frame(id = character(), symbol = character(), stringsAsFactors = FALSE)
  tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    db <- AnnotationDbi::select(org.Hs.eg.db,
      keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
      columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
    db <- db[!is.na(db$SYMBOL) & !duplicated(db$ENSEMBL), ]
    map <- data.frame(id = db$ENSEMBL, symbol = db$SYMBOL, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  # Supplement from DEGs/DEPs files
  for (f in c("phase1_diff/DEGs_Adjacent_vs_Normal.csv", "phase1_diff/DEPs_Adjacent_vs_Normal.csv")) {
    fp <- file.path(RES, f)
    if (file.exists(fp)) {
      d <- read.csv(fp, stringsAsFactors = FALSE)
      id_col <- intersect(c("gene_id", "Protein"), colnames(d))[1]
      nm_col <- intersect(c("gene_name", "symbol"), colnames(d))[1]
      if (!is.na(id_col) && !is.na(nm_col)) {
        sub <- data.frame(id = gsub("\\.[0-9]+$", "", d[[id_col]]),
                          symbol = d[[nm_col]], stringsAsFactors = FALSE)
        sub <- sub[!is.na(sub$symbol) & sub$symbol != "" & sub$symbol != "_" &
                   !grepl("^ENS", sub$symbol), ]
        map <- rbind(map, sub)
      }
    }
  }
  map[!duplicated(map$id), ]
}

# =============================================================================
# SECTION 7: Data Loading
# =============================================================================
cat("  Loading data...\n")

tc_mat <- read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"),
                   row.names = 1, check.names = FALSE)
pr_mat <- read.csv(file.path(DATA, "proteomics_log2_norm.csv"),
                   row.names = 1, check.names = FALSE)

degs <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
deps <- read.csv(file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
dems <- read.csv(file.path(RES, "phase1_diff/DEMs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)

gsea_tc <- read.csv(file.path(RES, "phase2_enrichment/GSEA_transcriptomics_Hallmark.csv"), stringsAsFactors = FALSE)
gsea_pr <- read.csv(file.path(RES, "phase2_enrichment/GSEA_proteomics_Hallmark.csv"), stringsAsFactors = FALSE)
convergent <- read.csv(file.path(RES, "phase2_enrichment/cross_omics_convergent_Hallmark.csv"), stringsAsFactors = FALSE)

gene_map <- build_gene_map()
cat(sprintf("  Data loaded: DEGs=%d, DEPs=%d, DEMs=%d, Gene map=%d\n",
            nrow(degs), nrow(deps), nrow(dems), nrow(gene_map)))

# =============================================================================
# SECTION 8: Panel B — Multi-omics Differential Overview
# =============================================================================
cat("\n--- Panel B: Multi-omics Differential Overview ---\n")

tryCatch({
  # Calculate sig counts using unified Methods threshold:
  # P < 0.05 & |log2FC| > 0.585 (i.e. fold change > 1.5x)
  LFC_CUT <- log2(1.5)  # = 0.5849625
  tc_up   <- sum(degs$P.Value < 0.05 & degs$logFC >  LFC_CUT)
  tc_down <- sum(degs$P.Value < 0.05 & degs$logFC < -LFC_CUT)
  pr_up   <- sum(deps$P.Value < 0.05 & deps$logFC >  LFC_CUT)
  pr_down <- sum(deps$P.Value < 0.05 & deps$logFC < -LFC_CUT)
  mt_up   <- sum(dems$P.Value < 0.05 & dems$logFC >  LFC_CUT)
  mt_down <- sum(dems$P.Value < 0.05 & dems$logFC < -LFC_CUT)
  cat(sprintf("  Panel B counts (P<0.05 & |log2FC|>0.585):\n"))
  cat(sprintf("    TC: up=%d, down=%d (total=%d)\n", tc_up, tc_down, tc_up + tc_down))
  cat(sprintf("    PR: up=%d, down=%d (total=%d)\n", pr_up, pr_down, pr_up + pr_down))
  cat(sprintf("    MT: up=%d, down=%d (total=%d)\n", mt_up, mt_down, mt_up + mt_down))

  df_b <- data.frame(
    Layer = factor(rep(c("Transcriptomics", "Proteomics", "Metabolomics"), each = 2),
                   levels = c("Transcriptomics", "Proteomics", "Metabolomics")),
    Direction = rep(c("Up-regulated", "Down-regulated"), 3),
    Count = c(tc_up, -tc_down, pr_up, -pr_down, mt_up, -mt_down),
    Label = c(tc_up, tc_down, pr_up, pr_down, mt_up, mt_down),
    stringsAsFactors = FALSE
  )

  # Total features per layer (integrated into x-axis labels)
  totals <- data.frame(
    Layer = factor(c("Transcriptomics", "Proteomics", "Metabolomics"),
                   levels = c("Transcriptomics", "Proteomics", "Metabolomics")),
    Total = c(nrow(degs), nrow(deps), nrow(dems)),
    stringsAsFactors = FALSE
  )
  # Build x-axis labels with totals
  layer_labels <- setNames(
    paste0(levels(totals$Layer), "\n(n = ", format(totals$Total, big.mark = ","), ")"),
    levels(totals$Layer)
  )

  max_val <- max(abs(df_b$Count)) * 1.25

  p_b <- ggplot(df_b, aes(x = Layer, y = Count, fill = Direction)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    geom_text(aes(label = Label,
                  y = ifelse(Count > 0, Count + max_val * 0.05, Count - max_val * 0.05)),
              size = FS_GEOM, family = FONT_FAMILY, fontface = "bold") +
    scale_fill_manual(values = c("Up-regulated" = COL_UP, "Down-regulated" = COL_DOWN),
                      name = NULL) +
    scale_x_discrete(labels = layer_labels) +
    scale_y_continuous(limits = c(-max_val, max_val),
                       labels = function(x) abs(x),
                       breaks = pretty(c(-max_val, max_val), n = 6)) +
    labs(title = "Differential features by omics layer",
         x = NULL,
         y = expression(paste("No. features (", italic(P), " < 0.05, |", log[2], "FC| > 0.585)"))) +
    theme_nc +
    theme(
      legend.position = "top",
      legend.justification = "right",
      legend.direction = "horizontal",
      legend.key.size = unit(3, "mm"),
      legend.margin = margin(0, 0, -2, 0, "mm"),
      axis.text.x = element_text(size = 8, face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey92")
    )

  save_panel_pdf("Fig1b_diff_overview.pdf", W_B, H_R1, print(p_b))
}, error = function(e) cat("  ERROR Panel B:", e$message, "\n"))

# =============================================================================
# SECTION 9: Panel C — mRNA-Protein Discordance Scatter
# =============================================================================
cat("\n--- Panel C: mRNA-Protein Discordance Scatter ---\n")

tryCatch({
  # Match genes between TC and PR by gene_name
  degs_c <- degs
  degs_c$symbol <- gene_map$symbol[match(gsub("\\.[0-9]+$", "", degs_c$gene_id), gene_map$id)]
  degs_c$symbol <- ifelse(is.na(degs_c$symbol), degs_c$gene_name, degs_c$symbol)
  degs_c <- degs_c[!is.na(degs_c$symbol) & degs_c$symbol != "" & !grepl("^ENS", degs_c$symbol), ]
  degs_c <- degs_c[!duplicated(degs_c$symbol), ]

  deps_c <- deps
  deps_c$symbol <- deps_c$gene_name
  deps_c <- deps_c[!is.na(deps_c$symbol) & deps_c$symbol != "" & deps_c$symbol != "_", ]
  deps_c <- deps_c[!duplicated(deps_c$symbol), ]

  common_genes <- intersect(degs_c$symbol, deps_c$symbol)
  cat(sprintf("  Common genes for scatter: %d\n", length(common_genes)))

  scatter_df <- data.frame(
    gene = common_genes,
    tc_logFC = degs_c$logFC[match(common_genes, degs_c$symbol)],
    pr_logFC = deps_c$logFC[match(common_genes, deps_c$symbol)],
    stringsAsFactors = FALSE
  )

  # Classify concordance quadrants
  scatter_df$quadrant <- with(scatter_df, case_when(
    tc_logFC > 0 & pr_logFC > 0 ~ "Concordant Up",
    tc_logFC < 0 & pr_logFC < 0 ~ "Concordant Down",
    TRUE ~ "Discordant"
  ))

  # Calculate Pearson correlation
  cor_val <- cor(scatter_df$tc_logFC, scatter_df$pr_logFC, use = "complete.obs")
  n_discord <- sum(scatter_df$quadrant == "Discordant")
  pct_discord <- round(100 * n_discord / nrow(scatter_df), 1)

  # Identify top discordant genes for labeling
  scatter_df$discord_score <- abs(scatter_df$tc_logFC) + abs(scatter_df$pr_logFC)
  top_discord <- scatter_df %>%
    filter(quadrant == "Discordant") %>%
    arrange(desc(discord_score)) %>%
    head(8)
  # Also label a few top concordant
  top_concord <- scatter_df %>%
    filter(quadrant != "Discordant") %>%
    arrange(desc(discord_score)) %>%
    head(4)
  label_genes <- rbind(top_discord, top_concord)

  scatter_df$show_label <- ifelse(scatter_df$gene %in% label_genes$gene, scatter_df$gene, "")

  p_c <- ggplot(scatter_df, aes(x = tc_logFC, y = pr_logFC)) +
    # Quadrant shading
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = COL_DISCORD, alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = COL_DISCORD, alpha = 0.04) +
    # Points
    geom_point(aes(color = quadrant), size = 0.3, alpha = 0.4) +
    # Reference lines
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    # Trend line
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6,
                color = "grey30", linetype = "solid") +
    # Labels
    geom_label_repel(
      data = scatter_df[scatter_df$show_label != "", ],
      aes(label = show_label), size = 2.2, family = FONT_FAMILY,
      fontface = "italic", max.overlaps = 20, force = 3,
      segment.size = 0.2, segment.color = "grey50",
      label.size = 0.12, label.padding = unit(0.1, "lines"),
      fill = alpha("white", 0.85), seed = 42
    ) +
    scale_color_manual(
      values = c("Concordant Up" = COL_UP, "Concordant Down" = COL_DOWN, "Discordant" = COL_DISCORD),
      name = NULL
    ) +
    # Annotation: correlation + discordance %
    annotate("text", x = Inf, y = -Inf,
             label = sprintf("r = %.2f\nDiscordant: %s%%\n(%d / %d genes)",
                             cor_val, pct_discord, n_discord, nrow(scatter_df)),
             hjust = 1.1, vjust = -0.3, size = 2.5, family = FONT_FAMILY, color = "grey30") +
    labs(title = "mRNA-Protein concordance",
         x = expression(log[2]~FC~"(Transcriptomics)"),
         y = expression(log[2]~FC~"(Proteomics)")) +
    coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-3, 3)) +
    theme_nc +
    theme(
      legend.position = c(0.82, 0.92),
      legend.direction = "vertical",
      legend.key.size = unit(2.5, "mm"),
      legend.spacing.y = unit(0.5, "mm")
    )

  save_panel_pdf("Fig1c_discordance_scatter.pdf", W_C, H_R1, print(p_c))
}, error = function(e) cat("  ERROR Panel C:", e$message, "\n"))

# =============================================================================
# SECTION 10: Panel E — Cross-omics Hallmark Pathway Heatmap (Row 3 left)
# =============================================================================
cat("\n--- Panel E: Cross-omics Hallmark Pathway Heatmap ---\n")

tryCatch({
  # Use convergent data (both TC and PR significant)
  conv <- convergent
  conv$pathway_clean <- gsub("^HALLMARK_", "", conv$pathway)
  conv$pathway_clean <- gsub("_", " ", conv$pathway_clean)
  conv$pathway_clean <- tools::toTitleCase(tolower(conv$pathway_clean))
  # Fix biological nomenclature (capitalize known abbreviations)
  conv$pathway_clean <- gsub("\\bIl6\\b", "IL-6", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bIl2\\b", "IL-2", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bJak\\b", "JAK", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bStat3\\b", "STAT3", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bStat5\\b", "STAT5", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bTnfa\\b", "TNF-a", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bNfkb\\b", "NF-kB", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bMtorc1\\b", "mTORC1", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bDna\\b", "DNA", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bUv\\b", "UV", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bE2f\\b", "E2F", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bG2m\\b", "G2M", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bKras\\b", "KRAS", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bMyc\\b", "MYC", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bP53\\b", "p53", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bPi3k\\b", "PI3K", conv$pathway_clean)
  conv$pathway_clean <- gsub("\\bAkt\\b", "AKT", conv$pathway_clean)

  # Select top pathways: all with both TC padj < 0.05
  conv_sig <- conv[conv$TC_padj < 0.05 & conv$PR_padj < 0.05, ]
  # Sort by maximum absolute NES
  conv_sig$max_nes <- pmax(abs(conv_sig$TC_NES), abs(conv_sig$PR_NES))
  conv_sig <- conv_sig[order(conv_sig$max_nes, decreasing = TRUE), ]
  # Take top 15
  if (nrow(conv_sig) > 15) conv_sig <- head(conv_sig, 15)

  # Concordance annotation
  concord_vec <- ifelse(conv_sig$concordant, "Concordant", "Discordant")

  # --- Sort rows logically: by Regulation, then by mean NES within each group ---
  conv_sig$mean_nes <- rowMeans(conv_sig[, c("TC_NES", "PR_NES")])
  conv_sig$reg_group <- concord_vec
  # Concordant first (both same direction), ordered by mean NES (most positive first)
  # Then Discordant, ordered by mean NES
  conc_rows <- conv_sig[conv_sig$reg_group == "Concordant", ]
  disc_rows <- conv_sig[conv_sig$reg_group == "Discordant", ]
  conc_rows <- conc_rows[order(conc_rows$mean_nes, decreasing = TRUE), ]
  disc_rows <- disc_rows[order(disc_rows$mean_nes, decreasing = TRUE), ]
  conv_sig <- rbind(conc_rows, disc_rows)
  concord_vec <- conv_sig$reg_group

  mat_d <- as.matrix(conv_sig[, c("TC_NES", "PR_NES")])
  rownames(mat_d) <- conv_sig$pathway_clean
  colnames(mat_d) <- c("mRNA", "Protein")

  # Regulation annotation (right side)
  reg_factor <- factor(concord_vec, levels = c("Concordant", "Discordant"))

  ha_right <- rowAnnotation(
    Regulation = reg_factor,
    col = list(Regulation = c("Concordant" = COL_CONCORD, "Discordant" = COL_DISCORD)),
    simple_anno_size = unit(3.5, "mm"),
    show_legend = TRUE,
    annotation_name_gp = gp_title(8),
    annotation_legend_param = list(
      Regulation = list(
        title = "Regulation",
        title_gp = gp_legend_t(),
        labels_gp = gp_legend_l()
      )
    )
  )

  col_fn_d <- colorRamp2(c(-3, 0, 3), COL_DIVERGE)

  ht_d <- Heatmap(mat_d, name = "NES", col = col_fn_d,
    row_names_gp = gp_row(8),
    column_names_gp = gp_col(9, bold = TRUE),
    column_names_rot = 0,
    column_names_side = "top",
    column_names_centered = TRUE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_split = reg_factor, cluster_row_slices = FALSE,
    row_gap = unit(2, "mm"),
    row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_GRID),
    show_row_dend = FALSE,
    right_annotation = ha_right,
    column_title = "Hallmark Enrichment",
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
    rect_gp = gpar(col = "white", lwd = 0.5),
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- mat_d[i, j]
      grid.text(sprintf("%.1f", v), x, y,
                gp = gpar(fontsize = 7.5, fontfamily = FONT_GRID,
                          col = ifelse(abs(v) > 2.0, "white", "black")))
    },
    heatmap_legend_param = list(
      title_gp = gp_legend_t(),
      labels_gp = gp_legend_l(),
      legend_height = unit(22, "mm"),
      grid_width = unit(3.5, "mm")
    ),
    width = unit(35, "mm"),
    row_names_max_width = unit(38, "mm")
  )

  save_panel_pdf("Fig1e_hallmark_crossomics.pdf", W_E, H_R3, {
    draw(ht_d, padding = unit(c(6, 4, 2, 10), "mm"), merge_legend = TRUE)
  })
}, error = function(e) cat("  ERROR Panel E:", e$message, "\n"))

# =============================================================================
# SECTION 11: Panel D — Multi-omics Convergence Dumbbell Plot (Full-width)
# =============================================================================
cat("\n--- Panel D: Multi-omics Convergence Dumbbell Plot ---\n")

tryCatch({
  # --- Get gene symbols for both omics ---
  degs_sym <- degs$gene_name
  if (is.null(degs_sym)) degs_sym <- gene_map$symbol[match(gsub("\\.[0-9]+$", "", degs$gene_id), gene_map$id)]
  pr_name_col <- intersect(c("gene_name", "symbol"), colnames(deps))[1]
  deps_sym <- deps[[pr_name_col]]

  # --- Build logFC and P-value named vectors ---
  deg_fc  <- setNames(degs$logFC, degs_sym)
  deg_p   <- setNames(degs$P.Value, degs_sym)
  dep_fc  <- setNames(deps$logFC, deps_sym)
  dep_p   <- setNames(deps$P.Value, deps_sym)

  # --- Find DEG-significant genes (manuscript definition: P<0.05 & |log2FC|>0.585)
  #     present in proteomics, with concordant direction at the protein level ---
  shared <- intersect(
    degs_sym[!is.na(degs_sym) & degs_sym != "" & !grepl("^ENS", degs_sym)],
    deps_sym[!is.na(deps_sym) & deps_sym != "" & !grepl("^ENS", deps_sym)]
  )
  LFC_CUT <- log2(1.5)
  deg_sig <- shared[deg_p[shared] < 0.05 & abs(deg_fc[shared]) > LFC_CUT]
  deg_sig <- deg_sig[!is.na(deg_sig)]
  # Keep only genes whose proteomic logFC is in the same direction (convergent)
  concordant <- deg_sig[sign(deg_fc[deg_sig]) == sign(dep_fc[deg_sig])]
  concordant <- concordant[!is.na(concordant)]
  cat(sprintf("  DEG-significant in shared: %d\n", length(deg_sig)))
  cat(sprintf("  Concordant (mRNA+protein same direction): %d\n", length(concordant)))

  # --- Balanced selection: 10 Up + 10 Down by combined significance ---
  conc_up <- concordant[deg_fc[concordant] > 0]
  conc_dn <- concordant[deg_fc[concordant] < 0]

  score_fn <- function(genes) {
    -log10(pmin(deg_p[genes], 1, na.rm = TRUE)) - log10(pmin(dep_p[genes], 1, na.rm = TRUE))
  }
  n_per <- 10
  top_up <- names(sort(score_fn(conc_up), decreasing = TRUE))[1:min(n_per, length(conc_up))]
  top_dn <- names(sort(score_fn(conc_dn), decreasing = TRUE))[1:min(n_per, length(conc_dn))]
  selected <- c(top_up, top_dn)
  selected <- selected[!is.na(selected)]
  cat(sprintf("  Selected: %d (Up=%d, Down=%d)\n", length(selected), length(top_up), length(top_dn)))

  # --- Build data frame for dumbbell plot ---
  df_dumb <- data.frame(
    gene = rep(selected, 2),
    logFC = c(deg_fc[selected], dep_fc[selected]),
    Omics = rep(c("mRNA", "Protein"), each = length(selected)),
    Direction = ifelse(rep(deg_fc[selected] > 0, 2), "Up", "Down"),
    stringsAsFactors = FALSE
  )

  # Order genes: Up group by mRNA logFC descending, Down group by mRNA logFC ascending
  up_order <- top_up[order(deg_fc[top_up], decreasing = TRUE)]
  dn_order <- top_dn[order(deg_fc[top_dn], decreasing = FALSE)]
  gene_order <- c(dn_order, up_order)  # Down at bottom, Up at top
  df_dumb$gene <- factor(df_dumb$gene, levels = gene_order)

  # Segment data for connecting lines
  df_seg <- data.frame(
    gene = selected,
    x_mRNA = deg_fc[selected],
    x_Protein = dep_fc[selected],
    Direction = ifelse(deg_fc[selected] > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )
  df_seg$gene <- factor(df_seg$gene, levels = gene_order)

  # --- Plot ---
  p_d <- ggplot() +
    # Vertical zero line
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
    # Connecting segments
    geom_segment(data = df_seg,
                 aes(x = x_mRNA, xend = x_Protein, y = gene, yend = gene, color = Direction),
                 linewidth = 0.6, alpha = 0.6) +
    # Points with shape mapped to Omics
    geom_point(data = df_dumb,
               aes(x = logFC, y = gene, color = Direction, shape = Omics),
               size = 2.8) +
    scale_color_manual(values = c("Up" = COL_UP, "Down" = COL_DOWN), name = "Direction") +
    scale_shape_manual(values = c("mRNA" = 16, "Protein" = 17), name = "Omics") +
    scale_x_continuous(name = expression(log[2]~"Fold Change"), expand = expansion(mult = 0.05)) +
    labs(y = NULL, title = "Multi-omics Convergent Dysregulation") +
    facet_grid(Direction ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_bw(base_size = 8, base_family = FONT_GRID) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8, face = "italic", family = FONT_GRID),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.major.x = element_line(linewidth = 0.15, color = "grey92"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text.y.left = element_text(angle = 0, size = 9, face = "bold", family = FONT_GRID),
      strip.placement = "outside",
      legend.position = "right",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7.5),
      plot.margin = margin(4, 6, 4, 4, "mm")
    ) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

  save_panel_pdf("Fig1d_top_features.pdf", W_TOTAL, H_R2, print(p_d))
}, error = function(e) cat("  ERROR Panel D:", e$message, "\n"))

# =============================================================================
# SECTION 12: Panel F — Metabolic Suppression Convergence
# =============================================================================
cat("\n--- Panel F: Metabolic Suppression Convergence ---\n")

tryCatch({
  # Select metabolic pathways from Hallmark (concordantly downregulated)
  metabolic_keywords <- c("XENOBIOTIC", "BILE_ACID", "FATTY_ACID", "HEME",
                          "PEROXISOME", "ADIPOGENESIS", "CHOLESTEROL",
                          "OXIDATIVE_PHOSPHORYLATION", "SPERMATOGENESIS")
  metab_conv <- convergent[grepl(paste(metabolic_keywords, collapse = "|"), convergent$pathway), ]
  # Also add from GSEA TC (metabolic pathways that are suppressed)
  metab_all <- gsea_tc[grepl(paste(metabolic_keywords, collapse = "|"), gsea_tc$pathway) &
                       gsea_tc$padj < 0.1, ]

  # Merge TC and PR NES for metabolic pathways
  df_f <- data.frame(
    pathway = metab_conv$pathway,
    TC_NES = metab_conv$TC_NES,
    PR_NES = metab_conv$PR_NES,
    concordant = metab_conv$concordant,
    stringsAsFactors = FALSE
  )
  # Clean names
  df_f$pathway_clean <- gsub("^HALLMARK_", "", df_f$pathway)
  df_f$pathway_clean <- gsub("_", " ", df_f$pathway_clean)
  df_f$pathway_clean <- tools::toTitleCase(tolower(df_f$pathway_clean))
  # Sort by TC NES
  df_f <- df_f[order(df_f$TC_NES), ]
  df_f$pathway_clean <- factor(df_f$pathway_clean, levels = df_f$pathway_clean)

  # Reshape for dumbbell plot
  df_long <- df_f %>%
    dplyr::select(pathway_clean, TC_NES, PR_NES) %>%
    pivot_longer(cols = c(TC_NES, PR_NES), names_to = "Layer", values_to = "NES") %>%
    mutate(Layer = ifelse(Layer == "TC_NES", "Transcriptomics", "Proteomics"))

  p_f <- ggplot(df_long, aes(x = NES, y = pathway_clean)) +
    # Dumbbell connecting segments (TC to PR for each pathway)
    geom_segment(data = df_f,
                 aes(x = TC_NES, xend = PR_NES, y = pathway_clean, yend = pathway_clean),
                 color = "grey60", linewidth = 0.8, inherit.aes = FALSE) +
    # Stem from zero to closer point
    geom_segment(data = df_f,
                 aes(x = 0, xend = pmax(TC_NES, PR_NES),
                     y = pathway_clean, yend = pathway_clean),
                 color = "grey88", linewidth = 0.3, linetype = "dotted", inherit.aes = FALSE) +
    # Points for each layer (larger, with border)
    geom_point(aes(color = Layer, shape = Layer), size = 3.2, stroke = 0.4) +
    # Suppression zone shading
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = COL_DOWN, alpha = 0.03) +
    scale_color_manual(values = c("Transcriptomics" = COL_TC, "Proteomics" = COL_PR),
                       name = "Omics Layer") +
    scale_shape_manual(values = c("Transcriptomics" = 16, "Proteomics" = 17),
                       name = "Omics Layer") +
    scale_x_continuous(breaks = seq(-3, 0, by = 1), limits = c(-3.5, 0),
                       expand = expansion(mult = c(0.02, 0.03))) +
    labs(title = "Metabolic pathway suppression",
         x = "Normalized Enrichment Score",
         y = NULL) +
    theme_nc +
    theme(
      legend.position = c(0.78, 0.18),
      legend.direction = "vertical",
      legend.key.size = unit(3.5, "mm"),
      legend.background = element_rect(fill = alpha("white", 0.9), linewidth = 0.3, color = "grey80"),
      panel.grid.major.x = element_line(linewidth = 0.15, color = "grey90"),
      panel.grid.major.y = element_line(linewidth = 0.1, color = "grey95"),
      axis.text.y = element_text(size = 8)
    )

  save_panel_pdf("Fig1f_metabolic_suppression.pdf", W_F, H_R3, print(p_f))
}, error = function(e) cat("  ERROR Panel F:", e$message, "\n"))

# =============================================================================
# SECTION 13: Composite Assembly (Pure Vector — Grid Viewports)
# =============================================================================
cat("\n--- Assembling composite Figure_1 b-f ---\n")

tryCatch({
  render_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))

    # --- Row 1: b (left) + c (right), y from H_R2+H_R3 to top ---
    y_r1 <- H_R2 + H_R3
    # Panel b
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_r1, "mm"),
      width = unit(W_B, "mm"), height = unit(H_R1, "mm"),
      just = c("left", "bottom")
    ))
    print(p_b, newpage = FALSE)
    grid::popViewport()
    # Panel c
    grid::pushViewport(grid::viewport(
      x = unit(W_B, "mm"), y = unit(y_r1, "mm"),
      width = unit(W_C, "mm"), height = unit(H_R1, "mm"),
      just = c("left", "bottom")
    ))
    print(p_c, newpage = FALSE)
    grid::popViewport()

    # --- Row 2: d (full-width) — Top features heatmap, y from H_R3 to H_R3+H_R2 ---
    y_r2 <- H_R3
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_r2, "mm"),
      width = unit(W_TOTAL, "mm"), height = unit(H_R2, "mm"),
      just = c("left", "bottom")
    ))
    print(p_d, newpage = FALSE)
    grid::popViewport()

    # --- Row 3: e (left, Hallmark) + f (right), y from 0 to H_R3 ---
    # Panel e
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W_E, "mm"), height = unit(H_R3, "mm"),
      just = c("left", "bottom")
    ))
    draw(ht_d, padding = unit(c(2, 4, 2, 4), "mm"), merge_legend = TRUE, newpage = FALSE)
    grid::popViewport()
    # Panel f
    grid::pushViewport(grid::viewport(
      x = unit(W_E, "mm"), y = unit(0, "mm"),
      width = unit(W_F, "mm"), height = unit(H_R3, "mm"),
      just = c("left", "bottom")
    ))
    print(p_f, newpage = FALSE)
    grid::popViewport()

    # --- Panel labels ---
    lbl <- data.frame(
      t = c("B", "C", "D", "E", "F"),
      x = c(1, W_B + 1, 1, 1, W_E + 1),
      y = c(H_TOTAL - 1, H_TOTAL - 1, y_r1 - 1, H_R3 - 1, H_R3 - 1),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(lbl))) {
      grid::grid.text(lbl$t[i], x = unit(lbl$x[i], "mm"), y = unit(lbl$y[i], "mm"),
                      just = c("left", "top"),
                      gp = gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY))
    }

    grid::popViewport()
  }

  # PDF
  w_in <- mm2in(W_TOTAL); h_in <- mm2in(H_TOTAL)
  cairo_pdf(file.path(OUT, "Figure_1_bf.pdf"), width = w_in, height = h_in, family = FONT_FAMILY)
  render_composite()
  dev.off()
  cat("  [OK] Figure_1_bf.pdf\n")

  # PNG
  png(file.path(OUT, "Figure_1_bf.png"), width = round(W_TOTAL * DPI / 25.4),
      height = round(H_TOTAL * DPI / 25.4), res = DPI, type = "cairo")
  render_composite()
  dev.off()
  cat("  [OK] Figure_1_bf.png\n")

  # TIFF
  tiff(file.path(OUT, "Figure_1_bf.tiff"), width = round(W_TOTAL * DPI / 25.4),
       height = round(H_TOTAL * DPI / 25.4), res = DPI, compression = "lzw", type = "cairo")
  render_composite()
  dev.off()
  cat("  [OK] Figure_1_bf.tiff\n")

}, error = function(e) cat("  ERROR Assembly:", e$message, "\n"))

cat("\n=== Figure 1 (v2) complete ===\n")
