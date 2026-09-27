#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_03_standalone.R  (v2: pixel-exact + magick canvas composite)
# Supplementary Figure 3: Pathway Enrichment + Reactome Validation
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: EBioMedicine (Lancet family) - vector PDF + PNG/TIFF (600 DPI) composite
# =============================================================================
# Final 8-panel layout (4 rows x 2 cols, 183x245mm canvas):
#   Row1 (60mm):  a = TC Hallmark GSEA dotplot       | b = PR Hallmark GSEA dotplot
#   Row2 (60mm):  c = GO-BP TC NES heatmap (top 15)  | d = GO-BP PR NES bar (top 12)
#   Row3 (60mm):  e = PR KEGG dotplot (top 12)       | f = Reactome category bar (top 12)
#   Row4 (65mm):  g = Immune-metabolic Spearman heat | h = Cross-cohort validation bar
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_03_standalone.R
# =============================================================================

cat("=== Supplementary Figure 3: Pathway Enrichment (v2) ===\n")
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
  library(ggsci)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_03")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants (NC: minimum 8pt)
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4
mm2in <- function(mm) mm / MM_PER_INCH

ASSEMBLY_DPI <- 600

# Final composite canvas
W_TOTAL <- 183
H_TOTAL <- 245

# Per-panel target dimensions (mm)
W_LEFT  <- 91;  W_RIGHT <- 92
H_ROW1  <- 60;  H_ROW2  <- 60;  H_ROW3 <- 60;  H_ROW4 <- 65

# Font sizes
FS_TITLE      <- 9
FS_AXIS_TITLE <- 8
FS_AXIS_TEXT  <- 7
FS_LEGEND_T   <- 7
FS_LEGEND_L   <- 7
FS_ROW_NAME   <- 7
FS_TAG        <- 12

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
COL_UP   <- "#CD534CFF"
COL_DOWN <- "#0073C2FF"
COL_VAL  <- "#EFC000FF"

col_nes <- colorRamp2(c(-3, 0, 3), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_cor <- colorRamp2(c(-0.6, 0, 0.6), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 5: Theme  (compact for 91x60mm slots)
# =============================================================================
theme_nc <- theme_bw(base_size = 7, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.4, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text  = element_text(size = FS_AXIS_TEXT, color = "black"),
    axis.title = element_text(size = FS_AXIS_TITLE, face = "bold"),
    plot.title = element_text(size = FS_TITLE, face = "bold", hjust = 0),
    legend.text  = element_text(size = FS_LEGEND_L),
    legend.title = element_text(size = FS_LEGEND_T, face = "bold"),
    legend.key.size = unit(2.2, "mm"),
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0, "mm"),
    legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm")
  )
theme_set(theme_nc)

ht_opt$message <- FALSE

# =============================================================================
# SECTION 6: Helper Functions
# =============================================================================
gp_anno_name <- function(size = FS_AXIS_TEXT) gpar(fontsize = size, fontfamily = FONT_GRID)
gp_legend_title  <- function(size = FS_LEGEND_T) gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
gp_legend_labels <- function(size = FS_LEGEND_L) gpar(fontsize = size, fontfamily = FONT_GRID)

std_legend_param <- function() {
  list(title_gp = gp_legend_title(), labels_gp = gp_legend_labels(),
       legend_height = unit(15, "mm"), grid_width = unit(2.5, "mm"))
}

# pixel-exact PDF saver for ggplot panels
save_panel_pdf <- function(filename, w_mm, h_mm, plot_obj) {
  fpath <- file.path(OUT, filename)
  cairo_pdf(fpath, width = mm2in(w_mm), height = mm2in(h_mm), family = FONT_FAMILY)
  print(plot_obj)
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

# pixel-exact PDF saver for ComplexHeatmap panels (no padding by default)
save_heatmap_panel_pdf <- function(filename, w_mm, h_mm, ht_obj, padding = c(2, 2, 2, 2)) {
  fpath <- file.path(OUT, filename)
  cairo_pdf(fpath, width = mm2in(w_mm), height = mm2in(h_mm), family = FONT_FAMILY)
  draw(ht_obj, padding = unit(padding, "mm"))
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else NULL
}

clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_|^REACTOME_|^KEGG_|^GOBP_|^GO_|^WP_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# truncate long pathway names with ellipsis
trim_label <- function(x, max_len = 38) {
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 1), "\u2026"), x)
}

# --- Compact GSEA dotplot builder (sized for 91x60mm slot) ---
make_gsea_dotplot <- function(gsea_file, title, n_top = 12) {
  df <- read.csv(gsea_file, stringsAsFactors = FALSE)
  name_col <- find_col(df, c("Description", "pathway_name", "pathway", "term", "ID"))
  p_col    <- find_col(df, c("pval", "pvalue", "p.adjust", "padj", "qvalue", "P.Value"))
  nes_col  <- find_col(df, c("NES", "enrichmentScore"))
  size_col <- find_col(df, c("setSize", "Count", "count", "size"))
  if (is.null(name_col) || is.null(p_col)) return(NULL)

  df$pathway_name <- trim_label(clean_pathway(df[[name_col]]), 38)
  df$neg_log10p <- -log10(df[[p_col]] + 1e-300)
  df$direction  <- if (!is.null(nes_col)) ifelse(df[[nes_col]] > 0, "Up", "Down") else "Up"
  df <- head(df[order(df$neg_log10p, decreasing = TRUE), ], n_top)

  p <- ggplot(df, aes(x = neg_log10p, y = reorder(pathway_name, neg_log10p)))
  if (!is.null(size_col)) {
    p <- p + geom_point(aes(size = .data[[size_col]], fill = direction),
                        shape = 21, color = "black", stroke = 0.25, alpha = 0.85) +
      scale_size_continuous(range = c(1.5, 4), name = "Count")
  } else {
    p <- p + geom_point(aes(fill = direction), size = 2.5,
                        shape = 21, color = "black", stroke = 0.25, alpha = 0.85)
  }
  p + scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), name = "Dir.") +
    labs(title = title, x = "\u2013log\u2081\u2080(P)", y = NULL) +
    guides(fill = guide_legend(override.aes = list(size = 2.5))) +
    theme_nc +
    theme(axis.text.y = element_text(size = FS_AXIS_TEXT, hjust = 1),
          legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(0.5, "mm"))
}

# =============================================================================
# SECTION 7: Data Paths
# =============================================================================
cat("  Locating GSEA results...\n")

gsea_tc_hallmark <- file.path(RES, "phase2_enrichment/GSEA_transcriptomics_Hallmark.csv")
gsea_pr_hallmark <- file.path(RES, "phase2_enrichment/GSEA_proteomics_Hallmark.csv")
gsea_tc_gobp     <- file.path(RES, "phase2_enrichment/GSEA_transcriptomics_GOBP.csv")
gsea_pr_gobp     <- file.path(RES, "phase2_enrichment/GSEA_proteomics_GOBP.csv")
gsea_pr_kegg     <- file.path(RES, "phase2_enrichment/GSEA_proteomics_KEGG.csv")
reactome_tc_file <- file.path(RES, "phase2_enrichment/GSEA_transcriptomics_Reactome.csv")
im_corr_file     <- file.path(RES, "enhancement10_mechanism/immune_metabolic_correlation.csv")
meta_file        <- file.path(RES, "enhancement14_meta_analysis/meta_analysis_consensus_DEGs.csv")

for (f in c(gsea_tc_hallmark, gsea_pr_hallmark, gsea_tc_gobp, gsea_pr_gobp,
            gsea_pr_kegg, reactome_tc_file, im_corr_file, meta_file)) {
  cat("  ", ifelse(file.exists(f), "[OK]", "[MISSING]"), basename(f), "\n")
}

# =============================================================================
# SECTION 8: Panel Renderers (pixel-exact + object registry for vector composite)
# =============================================================================
cat("\n--- Rendering panels (pixel-exact) ---\n")

# Object registry: every panel object (ggplot or Heatmap) is stored here so we
# can re-draw the composite as PURE VECTOR (no rasterisation).
PANELS <- list()

# --- Panel a: TC Hallmark GSEA dotplot (91x60mm) ---
tryCatch({
  p <- make_gsea_dotplot(gsea_tc_hallmark, "Transcriptomics: Hallmark GSEA", n_top = 12)
  if (!is.null(p)) {
    save_panel_pdf("Supp03a_TC_Hallmark_GSEA.pdf", W_LEFT, H_ROW1, p)
    PANELS$a <- p
  }
  cat("  [OK] Panel a\n")
}, error = function(e) cat("  ERROR Panel a:", e$message, "\n"))

# --- Panel b: PR Hallmark GSEA dotplot (92x60mm) ---
tryCatch({
  p <- make_gsea_dotplot(gsea_pr_hallmark, "Proteomics: Hallmark GSEA", n_top = 12)
  if (!is.null(p)) {
    save_panel_pdf("Supp03b_PR_Hallmark_GSEA.pdf", W_RIGHT, H_ROW1, p)
    PANELS$b <- p
  }
  cat("  [OK] Panel b\n")
}, error = function(e) cat("  ERROR Panel b:", e$message, "\n"))

# --- Panel c: GO-BP TC NES heatmap (91x60mm) ---
tryCatch({
  df <- read.csv(gsea_tc_gobp, stringsAsFactors = FALSE)
  name_col <- find_col(df, c("Description", "pathway_name", "pathway", "term", "ID"))
  nes_col  <- find_col(df, c("NES", "enrichmentScore"))
  p_col    <- find_col(df, c("pval", "pvalue", "p.adjust", "padj"))

  if (!is.null(name_col) && !is.null(nes_col) && !is.null(p_col)) {
    df_sig <- df[df[[p_col]] < 0.05, ]
    df_sig <- head(df_sig[order(abs(df_sig[[nes_col]]), decreasing = TRUE), ], 15)
    df_sig$pw <- trim_label(clean_pathway(df_sig[[name_col]]), 38)

    mat <- matrix(df_sig[[nes_col]], ncol = 1)
    rownames(mat) <- df_sig$pw
    colnames(mat) <- "TC"
    mat[mat >  3] <-  3
    mat[mat < -3] <- -3

    ht <- Heatmap(mat, name = "NES", col = col_nes,
      row_names_gp = gpar(fontsize = FS_ROW_NAME, fontfamily = FONT_GRID),
      show_column_names = TRUE,
      column_names_gp = gpar(fontsize = FS_AXIS_TEXT, fontface = "bold", fontfamily = FONT_GRID),
      column_names_rot = 0, column_names_centered = TRUE,
      cluster_rows = TRUE, cluster_columns = FALSE,
      show_row_dend = FALSE,
      border = TRUE, rect_gp = gpar(col = "white", lwd = 0.3),
      column_title = "GO-BP Enrichment (TC)",
      column_title_gp = gpar(fontsize = FS_TITLE, fontface = "bold", fontfamily = FONT_GRID),
      width = unit(8, "mm"),
      heatmap_legend_param = std_legend_param(),
      row_names_max_width = unit(55, "mm"))

    save_heatmap_panel_pdf("Supp03c_GOBP_TC_heatmap.pdf", W_LEFT, H_ROW2, ht,
                            padding = c(1, 1, 1, 1))
    PANELS$c <- ht
    cat("  [OK] Panel c\n")
  }
}, error = function(e) cat("  ERROR Panel c:", e$message, "\n"))

# --- Panel d: GO-BP PR enrichment barplot (92x60mm) ---
tryCatch({
  df <- read.csv(gsea_pr_gobp, stringsAsFactors = FALSE)
  name_col <- find_col(df, c("Description", "pathway_name", "pathway", "term", "ID"))
  nes_col  <- find_col(df, c("NES", "enrichmentScore"))
  p_col    <- find_col(df, c("pval", "pvalue", "p.adjust", "padj"))

  if (!is.null(name_col) && !is.null(nes_col) && !is.null(p_col)) {
    df_sig <- df[df[[p_col]] < 0.05, ]
    df_sig <- head(df_sig[order(abs(df_sig[[nes_col]]), decreasing = TRUE), ], 12)
    df_sig$pw <- trim_label(clean_pathway(df_sig[[name_col]]), 40)
    df_sig$direction <- ifelse(df_sig[[nes_col]] > 0, "Up", "Down")

    p <- ggplot(df_sig, aes(x = .data[[nes_col]], y = reorder(pw, .data[[nes_col]]),
                             fill = direction)) +
      geom_col(color = "black", linewidth = 0.2, alpha = 0.85, width = 0.7) +
      scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), name = "Dir.") +
      geom_vline(xintercept = 0, linewidth = 0.4) +
      labs(title = "GO-BP Enrichment (PR)", x = "NES", y = NULL) +
      theme_nc +
      theme(axis.text.y = element_text(size = FS_AXIS_TEXT, hjust = 1),
            legend.position = "right")
    save_panel_pdf("Supp03d_GOBP_PR_barplot.pdf", W_RIGHT, H_ROW2, p)
    PANELS$d <- p
    cat("  [OK] Panel d\n")
  }
}, error = function(e) cat("  ERROR Panel d:", e$message, "\n"))

# --- Panel e: PR KEGG enrichment dotplot (91x60mm) ---
tryCatch({
  p <- make_gsea_dotplot(gsea_pr_kegg, "Proteomics: KEGG", n_top = 12)
  if (!is.null(p)) {
    save_panel_pdf("Supp03e_PR_KEGG.pdf", W_LEFT, H_ROW3, p)
    PANELS$e <- p
  }
  cat("  [OK] Panel e\n")
}, error = function(e) cat("  ERROR Panel e:", e$message, "\n"))

# --- Panel f: Reactome category mean-NES barplot (92x60mm) ---
tryCatch({
  reactome_tc <- read.csv(reactome_tc_file, stringsAsFactors = FALSE)
  reactome_tc_sig <- reactome_tc[reactome_tc$padj < 0.05, ]

  extract_category <- function(pathway) {
    pw <- gsub("REACTOME_", "", pathway)
    parts <- str_split(pw, "_")[[1]]
    paste(parts[1:min(2, length(parts))], collapse = "_")
  }
  reactome_tc_sig$category <- sapply(reactome_tc_sig$pathway, extract_category)

  category_summary <- reactome_tc_sig %>%
    group_by(category) %>%
    summarise(n_pathways = n(), mean_NES = mean(NES), min_padj = min(padj),
              .groups = "drop") %>%
    arrange(desc(abs(mean_NES)))

  top_cats <- head(category_summary, 12)
  top_cats$direction <- ifelse(top_cats$mean_NES > 0, "Up", "Down")
  top_cats$category_disp <- trim_label(clean_pathway(top_cats$category), 38)

  p <- ggplot(top_cats, aes(x = reorder(category_disp, abs(mean_NES)),
                             y = mean_NES, fill = direction)) +
    geom_col(alpha = 0.85, width = 0.7, color = "black", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), name = "Dir.") +
    labs(title = "Reactome Categories (TC)",
         x = NULL, y = "Mean NES") +
    theme_nc +
    theme(axis.text.y = element_text(size = FS_AXIS_TEXT),
          legend.position = "right")
  save_panel_pdf("Supp03f_Reactome_categories.pdf", W_RIGHT, H_ROW3, p)
  PANELS$f <- p
  cat("  [OK] Panel f\n")
}, error = function(e) cat("  ERROR Panel f:", e$message, "\n"))

# --- Panel g: Immune-metabolic coupling heatmap (91x65mm) ---
tryCatch({
  if (!file.exists(im_corr_file)) {
    cat("  SKIP Panel g: immune_metabolic_correlation.csv not found\n")
  } else {
    immune_metab_corr <- read.csv(im_corr_file, row.names = 1, check.names = FALSE)
    mat <- as.matrix(immune_metab_corr)
    mat[mat >  0.6] <-  0.6
    mat[mat < -0.6] <- -0.6

    if (nrow(mat) > 15) {
      rv <- apply(mat, 1, function(x) var(x, na.rm = TRUE))
      mat <- mat[order(rv, decreasing = TRUE)[1:15], , drop = FALSE]
    }
    if (ncol(mat) > 12) {
      cv <- apply(mat, 2, function(x) var(x, na.rm = TRUE))
      mat <- mat[, order(cv, decreasing = TRUE)[1:12], drop = FALSE]
    }

    rownames(mat) <- trim_label(rownames(mat), 28)
    colnames(mat) <- trim_label(colnames(mat), 22)

    ht <- Heatmap(mat, name = "Spearman r", col = col_cor,
      cluster_rows = TRUE, cluster_columns = TRUE,
      show_row_dend = FALSE, show_column_dend = FALSE,
      show_row_names = TRUE, show_column_names = TRUE,
      row_names_gp    = gpar(fontsize = FS_ROW_NAME, fontfamily = FONT_GRID),
      column_names_gp = gpar(fontsize = FS_AXIS_TEXT, fontfamily = FONT_GRID),
      column_names_rot = 45,
      column_title = "Immune-Metabolic Coupling",
      column_title_gp = gpar(fontsize = FS_TITLE, fontface = "bold", fontfamily = FONT_GRID),
      heatmap_legend_param = std_legend_param(),
      row_names_max_width    = unit(32, "mm"),
      column_names_max_height = unit(20, "mm"))
    save_heatmap_panel_pdf("Supp03g_immune_metabolic.pdf", W_LEFT, H_ROW4, ht,
                            padding = c(1, 1, 1, 1))
    PANELS$g <- ht
    cat("  [OK] Panel g\n")
  }
}, error = function(e) cat("  ERROR Panel g:", e$message, "\n"))

# --- Panel h: Cross-cohort validation rate barplot (92x65mm) ---
tryCatch({
  if (!file.exists(meta_file)) {
    cat("  SKIP Panel h: meta_analysis_consensus_DEGs.csv not found\n")
  } else {
    reactome_tc <- read.csv(reactome_tc_file, stringsAsFactors = FALSE)
    reactome_tc_sig <- reactome_tc[reactome_tc$padj < 0.05, ]
    top_pathways <- head(reactome_tc_sig[order(reactome_tc_sig$padj), ], 30)

    parse_leading_edge <- function(le_string) {
      if (is.na(le_string) || le_string == "") return(character(0))
      str_split(le_string, ";")[[1]]
    }
    pathway_genes <- lapply(seq_len(nrow(top_pathways)), function(i) {
      parse_leading_edge(top_pathways$leadingEdge[i])
    })
    names(pathway_genes) <- top_pathways$pathway

    meta_results <- read.csv(meta_file, stringsAsFactors = FALSE)
    leading_genes <- unique(unlist(pathway_genes))
    gene_col <- find_col(meta_results, c("gene", "Gene", "symbol", "gene_symbol"))
    padj_col <- find_col(meta_results, c("fisher_padj", "padj", "p.adjust", "FDR"))

    if (!is.null(gene_col) && !is.null(padj_col)) {
      consensus_validated <- meta_results[meta_results[[gene_col]] %in% leading_genes &
                                           meta_results[[padj_col]] < 0.05, ]
      pathway_validation <- sapply(names(pathway_genes), function(pw) {
        genes <- pathway_genes[[pw]]
        if (length(genes) == 0) return(0)
        sum(genes %in% consensus_validated[[gene_col]]) / length(genes)
      })

      val_df <- data.frame(
        pathway = names(pathway_validation),
        validation_rate = pathway_validation,
        stringsAsFactors = FALSE) %>%
        arrange(desc(validation_rate))

      top_val <- head(val_df, 12)
      top_val$pathway_short <- trim_label(clean_pathway(top_val$pathway), 40)

      p <- ggplot(top_val, aes(x = reorder(pathway_short, validation_rate),
                                y = validation_rate * 100)) +
        geom_col(fill = COL_VAL, alpha = 0.9, width = 0.7,
                 color = "black", linewidth = 0.2) +
        geom_hline(yintercept = mean(pathway_validation) * 100,
                   linetype = "dashed", color = "grey40", linewidth = 0.4) +
        coord_flip() +
        labs(title = "Cross-Cohort Validation",
             x = NULL, y = "Validation Rate (%)") +
        theme_nc +
        theme(axis.text.y = element_text(size = FS_AXIS_TEXT))
      save_panel_pdf("Supp03h_pathway_validation.pdf", W_RIGHT, H_ROW4, p)
      PANELS$h <- p
      cat("  [OK] Panel h\n")
    }
  }
}, error = function(e) cat("  ERROR Panel h:", e$message, "\n"))

# =============================================================================
# SECTION 9: Pure-Vector Composite (viewport-based, AI-editable)
# =============================================================================
cat("\n--- Assembling SuppFig_03 composite (pure vector) ---\n")

# Heatmap padding (mm) per cell to keep title/legend inside its slot
HT_PADDING <- c(2.5, 1.5, 2.5, 1.5)  # bottom, left, top, right

draw_panel_in_cell <- function(obj, row, col) {
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  if (inherits(obj, "Heatmap") || inherits(obj, "HeatmapList")) {
    draw(obj, newpage = FALSE, padding = unit(HT_PADDING, "mm"))
  } else if (inherits(obj, "ggplot")) {
    print(obj, newpage = FALSE)
  }
  popViewport()
}

render_vector_composite <- function() {
  grid.newpage()
  pushViewport(viewport(
    width  = unit(W_TOTAL, "mm"),
    height = unit(H_TOTAL, "mm"),
    layout = grid.layout(
      nrow = 4, ncol = 2,
      heights = unit(c(H_ROW1, H_ROW2, H_ROW3, H_ROW4), "mm"),
      widths  = unit(c(W_LEFT, W_RIGHT), "mm")
    )
  ))

  # Place all 8 panels
  if (!is.null(PANELS$a)) draw_panel_in_cell(PANELS$a, 1, 1)
  if (!is.null(PANELS$b)) draw_panel_in_cell(PANELS$b, 1, 2)
  if (!is.null(PANELS$c)) draw_panel_in_cell(PANELS$c, 2, 1)
  if (!is.null(PANELS$d)) draw_panel_in_cell(PANELS$d, 2, 2)
  if (!is.null(PANELS$e)) draw_panel_in_cell(PANELS$e, 3, 1)
  if (!is.null(PANELS$f)) draw_panel_in_cell(PANELS$f, 3, 2)
  if (!is.null(PANELS$g)) draw_panel_in_cell(PANELS$g, 4, 1)
  if (!is.null(PANELS$h)) draw_panel_in_cell(PANELS$h, 4, 2)

  # Panel labels (a-h) overlaid as vector text
  label_positions <- list(
    list("A", 1, 1), list("B", 1, 2),
    list("C", 2, 1), list("D", 2, 2),
    list("E", 3, 1), list("F", 3, 2),
    list("G", 4, 1), list("H", 4, 2)
  )
  for (lp in label_positions) {
    pushViewport(viewport(layout.pos.row = lp[[2]], layout.pos.col = lp[[3]]))
    grid.text(lp[[1]],
              x = unit(1.2, "mm"),
              y = unit(1, "npc") - unit(1.2, "mm"),
              just = c("left", "top"),
              gp = gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY))
    popViewport()
  }
  popViewport()
}

# PDF (pure vector, AI-editable)
cairo_pdf(file.path(OUT, "SuppFig_03.pdf"),
          width = mm2in(W_TOTAL), height = mm2in(H_TOTAL), family = FONT_FAMILY)
render_vector_composite()
dev.off()

# PNG (raster preview, 600 DPI)
grDevices::png(file.path(OUT, "SuppFig_03.png"),
               width = W_TOTAL, height = H_TOTAL, units = "mm", res = ASSEMBLY_DPI,
               type = "cairo")
render_vector_composite()
dev.off()

# TIFF (600 DPI submission copy, LZW + Cairo for compact output)
grDevices::tiff(file.path(OUT, "SuppFig_03.tiff"),
                width = W_TOTAL, height = H_TOTAL, units = "mm", res = ASSEMBLY_DPI,
                compression = "lzw", type = "cairo")
render_vector_composite()
dev.off()

cat("  -> SuppFig_03.pdf (pure vector) / .png / .tiff\n")
cat("\n=== Supplementary Figure 3 rendering complete ===\n")
