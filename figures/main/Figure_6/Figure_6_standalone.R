#!/usr/bin/env Rscript
# =============================================================================
# Figure_6_standalone.R
# Figure 6: Multi-omics Integration and Consensus Molecular Subtyping
# Target: EBioMedicine (Lancet family)
# =============================================================================
# Panels:
#   (A) DIABLO supervised integration - top feature loadings lollipop
#   (B) MOFA2 variance decomposition across latent factors
#   (C) WGCNA cross-omics module overlap bubble plot
#   (D) Consensus clustering K evaluation metrics
#   (E) PCA projection coloured by molecular subtype
#   (F) Subtype multi-omics profile heatmap (TC + PR markers)
#   (G) Comprehensive subtype characterisation (pathway + immune)
# =============================================================================
# Layout (183x250mm, vector grid viewport assembly):
#   Row 1 (90mm):  A (95x90)  | B (88x38, top) + D (88x52, bottom)
#   Row 2 (58mm):  C (93x58)  | E (90x58)
#   Row 3 (102mm): F (100x102)| G (83x102)
# =============================================================================
# Usage: conda run -n multiomics Rscript Figure_6_standalone.R
# =============================================================================

cat("=== Figure 6: Multi-omics Integration and Consensus Molecular Subtyping ===\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggrepel)
  library(ggsci)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/main/Figure_6")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants — Gradient font hierarchy
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM <- 25.4
ASSEMBLY_DPI <- 600

# Gradient font sizes (strict hierarchy):
# Tag (a/b/c) 12pt > Panel title 10pt > Axis title 9pt > Body text 8pt > Min 7pt
FS_TAG <- 12
FS_TITLE <- 10
FS_AXIS_TITLE <- 9
FS_AXIS_TEXT <- 8
FS_LEGEND_T <- 8
FS_LEGEND_L <- 8
FS_GEOM <- 2.82

# Layout dimensions (mm) — Row 1 expanded so Panel A's 24 Y-axis labels have ample spacing
W_TOTAL <- 183; H_TOTAL <- 250
W_A <- 95; W_BD <- W_TOTAL - W_A      # 88
H1 <- 90; H_B <- 38; H_D <- H1 - H_B  # 52
W_C <- 93; W_E <- W_TOTAL - W_C       # 90
H2 <- 58
W_F <- 100; W_G <- W_TOTAL - W_F      # 83
H3 <- H_TOTAL - H1 - H2               # 102

# Colours
PAL_CAT <- pal_jco("default")(10)
COL_UP <- "#CD534CFF"; COL_DOWN <- "#0073C2FF"; COL_NS <- "#868686FF"
COL_TC <- "#0073C2FF"; COL_PR <- "#CD534CFF"; COL_MT <- "#EFC000FF"
COL_CS1 <- "#CD534CFF"; COL_CS2 <- "#0073C2FF"

col_div <- colorRamp2(c(-2, 0, 2), c(COL_DOWN, "#FFFFFF", COL_UP))
col_zscore <- colorRamp2(c(-2.5, 0, 2.5), c(COL_DOWN, "#FFFFFF", COL_UP))

# =============================================================================
# Theme
# =============================================================================
theme_nc <- theme_bw(base_size = FS_AXIS_TEXT, base_family = FONT_FAMILY) +
  theme(
    text = element_text(family = FONT_FAMILY, size = FS_AXIS_TEXT),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = FS_AXIS_TEXT, color = "black", family = FONT_FAMILY),
    axis.title = element_text(size = FS_AXIS_TITLE, face = "bold", family = FONT_FAMILY),
    plot.title = element_text(size = FS_TITLE, face = "bold", hjust = 0.5, family = FONT_FAMILY),
    plot.title.position = "panel",
    legend.text = element_text(size = FS_LEGEND_L, family = FONT_FAMILY),
    legend.title = element_text(size = FS_LEGEND_T, face = "bold", family = FONT_FAMILY),
    legend.key.size = unit(3, "mm"),
    legend.background = element_blank(),
    strip.text = element_text(size = FS_AXIS_TEXT, face = "bold", family = FONT_FAMILY),
    plot.margin = margin(0.5, 1.5, 0.5, 1.5, "mm")
  )
theme_set(theme_nc)

# ComplexHeatmap gpar helpers — all defaults respect 8pt minimum (gradient floor)
gp_rn <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID)
gp_cn <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
gp_lt <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
gp_ll <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID)
gp_rt <- function(sz = 9) gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
gp_an <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID)
ht_opt$message <- FALSE

std_lp <- function() list(title_gp = gp_lt(), labels_gp = gp_ll(),
                          legend_height = unit(18, "mm"), grid_width = unit(3, "mm"))

save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fp <- file.path(OUT, filename)
  cairo_pdf(fp, width = w_mm / MM, height = h_mm / MM, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat(sprintf("  -> %s (%.0f x %.0f mm)\n", filename, w_mm, h_mm))
}

scale_rows <- function(mat, lim = 2) {
  mat <- t(scale(t(as.matrix(mat)))); mat[is.na(mat)] <- 0
  mat[mat > lim] <- lim; mat[mat < -lim] <- -lim; mat
}

# ID mapping
ensg_to_symbol <- function(ids) {
  m <- tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    db <- AnnotationDbi::select(org.Hs.eg.db,
      keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
      columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
    db <- db[!is.na(db$SYMBOL) & !duplicated(db$ENSEMBL), ]
    setNames(db$SYMBOL, db$ENSEMBL)
  }, error = function(e) character(0))
  ids_clean <- gsub("\\.[0-9]+$", "", ids)
  out <- ifelse(ids_clean %in% names(m), m[ids_clean], ids)
  out
}

ensp_to_symbol <- function(ids) {
  fp <- file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv")
  if (!file.exists(fp)) return(ids)
  d <- read.csv(fp, stringsAsFactors = FALSE)
  id_col <- intersect(c("Protein", "protein_id"), colnames(d))[1]
  nm_col <- intersect(c("gene_name", "symbol", "Gene"), colnames(d))[1]
  if (is.na(id_col) || is.na(nm_col)) return(ids)
  m <- setNames(d[[nm_col]], gsub("\\.[0-9]+$", "", d[[id_col]]))
  ids_clean <- gsub("\\.[0-9]+$", "", ids)
  out <- ifelse(ids_clean %in% names(m) & !is.na(m[ids_clean]) & m[ids_clean] != "",
                m[ids_clean], ids)
  out
}

# =============================================================================
# Panel object placeholders (populated by panel blocks; used by vector assembly)
# =============================================================================
p6a <- p6b <- p6c <- p6d <- p6e <- NULL
grob_f <- grob_g <- NULL

# =============================================================================
# PANEL A: DIABLO supervised integration - feature loadings lollipop (95x90mm)
# =============================================================================
cat("\n--- Panel A: DIABLO feature loadings (95x90mm) ---\n")
tryCatch({
  diablo <- read.csv(file.path(RES, "phase5_integration/DIABLO_features_annotated.csv"),
                     stringsAsFactors = FALSE)
  perm <- read.csv(file.path(RES, "phase5_integration/DIABLO_permutation_test.csv"),
                   stringsAsFactors = FALSE)

  # Top 8 per omics layer on component 1 = 24 total
  d1 <- diablo[diablo$component == 1, ]
  d1$abs_loading <- abs(d1$value.var)
  top_per_view <- d1 %>% group_by(view) %>% slice_max(abs_loading, n = 8) %>% ungroup()
  top_per_view$view <- factor(top_per_view$view,
    levels = c("transcriptomics", "proteomics", "metabolomics"))
  top_per_view <- top_per_view %>% arrange(view, desc(abs_loading))
  top_per_view$feature_name[is.na(top_per_view$feature_name) | top_per_view$feature_name == ""] <-
    top_per_view$feature[is.na(top_per_view$feature_name) | top_per_view$feature_name == ""]
  top_per_view$feature_name <- factor(top_per_view$feature_name,
    levels = rev(top_per_view$feature_name))

  view_cols <- c(transcriptomics = COL_TC, proteomics = COL_PR, metabolomics = COL_MT)
  ber_text <- sprintf("BER = %.1f%% (null = %.1f%%, P < 0.001)",
                      perm$original_BER * 100, perm$mean_permuted_BER * 100)

  p6a <- ggplot(top_per_view, aes(x = value.var, y = feature_name, color = view)) +
    geom_segment(aes(x = 0, xend = value.var, yend = feature_name), linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_color_manual(values = view_cols, name = "Omics layer",
                       labels = c("Transcriptomics", "Proteomics", "Metabolomics")) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    annotate("text", x = max(top_per_view$value.var) * 0.55,
             y = 2.5, label = ber_text, size = FS_GEOM, family = FONT_FAMILY) +
    labs(title = "DIABLO component 1 loadings", x = "Loading", y = NULL) +
    theme(legend.position = c(0.12, 0.93),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
          axis.text.y = element_text(size = 8),
          plot.margin = margin(2, 6, 2, 2, "mm"))
  save_panel_pdf("Fig6a_DIABLO_loadings.pdf", W_A, H1, print(p6a))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL B: MOFA2 variance decomposition (88x38mm)
# =============================================================================
cat("\n--- Panel B: MOFA2 variance decomposition (88x38mm) ---\n")
tryCatch({
  mofa <- read.csv(file.path(RES, "phase5_integration/MOFA2_K_selection_scan.csv"),
                   stringsAsFactors = FALSE)
  var_df <- data.frame(
    Factor = paste0("F", mofa$K),
    Transcriptomics = c(mofa$cumR2_transcriptomics[1],
                        diff(mofa$cumR2_transcriptomics)),
    Proteomics = c(mofa$cumR2_proteomics[1], diff(mofa$cumR2_proteomics)),
    Metabolomics = c(mofa$cumR2_metabolomics[1], diff(mofa$cumR2_metabolomics))
  )
  var_df <- head(var_df, 5)
  var_long <- var_df %>%
    pivot_longer(-Factor, names_to = "View", values_to = "Variance") %>%
    mutate(Factor = factor(Factor, levels = paste0("F", 1:5)),
           View = factor(View, levels = c("Transcriptomics", "Proteomics", "Metabolomics")))

  view_fills <- c(Transcriptomics = COL_TC, Proteomics = COL_PR, Metabolomics = COL_MT)

  p6b <- ggplot(var_long, aes(x = Factor, y = Variance, fill = View)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = view_fills, name = "Omics") +
    labs(title = "MOFA2 variance per factor",
         x = "Latent factor", y = "Variance (%)") +
    theme(legend.position = "top",
          legend.key.size = unit(2.5, "mm"),
          plot.title = element_text(size = FS_TITLE),
          plot.margin = margin(1, 4, 1, 4, "mm"))
  save_panel_pdf("Fig6b_MOFA2_variance.pdf", W_BD, H_B, print(p6b))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL C: WGCNA cross-omics conserved modules barplot (93x58mm)
# =============================================================================
cat("\n--- Panel C: WGCNA cross-omics module overlap (93x58mm) ---\n")
tryCatch({
  wgcna <- read.csv(file.path(RES, "enhancement2_wgcna/cross_omics_module_overlap.csv"),
                    stringsAsFactors = FALSE)

  # WGCNA module colors (biological meaning — functional co-expression clusters)
  wgcna_cols <- c(blue = "#4682B4", brown = "#8B4513", yellow = "#FFD700",
                  green = "#228B22", turquoise = "#40E0D0", black = "#2D2D2D",
                  red = "#CD5C5C", magenta = "#8B008B", grey = "#808080")

  # Top enriched cross-omics module pairs (OR > 1, significant)
  sig_pairs <- wgcna[wgcna$OR > 1 & wgcna$padj < 0.05, ]
  sig_pairs <- sig_pairs[order(-sig_pairs$OR), ]
  sig_pairs <- head(sig_pairs, 10)

  # Create display labels
  sig_pairs$pair_label <- sprintf("TC-%s ~ PR-%s", sig_pairs$TC_module, sig_pairs$PR_module)
  sig_pairs$pair_label <- factor(sig_pairs$pair_label,
    levels = rev(sig_pairs$pair_label))
  # Significance stars
  sig_pairs$stars <- ifelse(sig_pairs$padj < 0.001, "***",
                    ifelse(sig_pairs$padj < 0.01, "**", "*"))
  # Color by TC module
  sig_pairs$fill_col <- wgcna_cols[sig_pairs$TC_module]

  p6c <- ggplot(sig_pairs, aes(x = OR, y = pair_label)) +
    geom_col(aes(fill = TC_module), width = 0.7) +
    geom_text(aes(label = sprintf("n=%d %s", overlap, stars)),
              hjust = -0.1, size = FS_GEOM, family = FONT_FAMILY) +
    scale_fill_manual(values = wgcna_cols, name = "TC module") +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "grey50") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title = "Cross-omics conserved modules",
         x = "Odds ratio (Fisher's exact test)", y = NULL) +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8),
          plot.margin = margin(1, 8, 1, 4, "mm"))
  save_panel_pdf("Fig6c_WGCNA_overlap.pdf", W_C, H2, print(p6c))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL D: Consensus clustering K evaluation (88x52mm)
# =============================================================================
cat("\n--- Panel D: Consensus clustering metrics (88x52mm) ---\n")
tryCatch({
  boot <- read.csv(file.path(RES, "subtype_robustness_v2/bootstrap_stability.csv"),
                   stringsAsFactors = FALSE)
  kmet <- read.csv(file.path(RES, "phase6_subtyping/K_selection_metrics_comprehensive.csv"),
                   stringsAsFactors = FALSE)

  plot_df <- data.frame(
    K = kmet$K,
    Silhouette = kmet$Silhouette_concat,
    PAC = kmet$PAC_concat
  )
  plot_df <- merge(plot_df, boot[, c("K", "Mean_Jaccard")], by = "K", all.x = TRUE)
  plot_long <- plot_df %>%
    pivot_longer(-K, names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value))
  plot_long$Metric <- factor(plot_long$Metric,
    levels = c("Mean_Jaccard", "Silhouette", "PAC"))

  p6d <- ggplot(plot_long, aes(x = factor(K), y = Value, color = Metric, group = Metric)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c(Mean_Jaccard = COL_UP, Silhouette = COL_TC, PAC = COL_MT),
                       labels = c("Bootstrap Jaccard", "Silhouette", "PAC")) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "grey40") +
    annotate("text", x = 1.3, y = max(plot_long$Value, na.rm = TRUE) * 0.95,
             label = "K=2", size = FS_GEOM, family = FONT_FAMILY, hjust = 0,
             fontface = "bold") +
    labs(title = "Consensus clustering K evaluation",
         x = "Number of clusters (K)", y = "Metric value", color = NULL) +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.9), color = "grey80", linewidth = 0.3),
          legend.key.size = unit(2.5, "mm"),
          plot.title = element_text(size = FS_TITLE),
          plot.margin = margin(1, 4, 1, 4, "mm"))
  save_panel_pdf("Fig6d_consensus_K_eval.pdf", W_BD, H_D, print(p6d))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL E: PCA projection by subtype (90x58mm)
# =============================================================================
cat("\n--- Panel E: Subtype PCA (90x58mm) ---\n")
tryCatch({
  sa <- read.csv(file.path(RES, "phase6_subtyping/subtype_K2.csv"), stringsAsFactors = FALSE)
  colnames(sa)[colnames(sa) == "subtype"] <- "cluster"
  tc_mat <- as.matrix(read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"),
                                row.names = 1, check.names = FALSE))
  shared <- intersect(sa$sample, colnames(tc_mat))
  pca <- prcomp(t(tc_mat[, shared]), scale. = TRUE)
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample = shared)
  pca_df <- merge(pca_df, sa[, c("sample", "cluster")], by = "sample")
  ve <- summary(pca)$importance[2, 1:2] * 100

  p6e <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
    stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.68, linewidth = 0.7) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_text_repel(aes(label = sample), size = FS_GEOM, family = FONT_FAMILY,
                    show.legend = FALSE, max.overlaps = 15, segment.size = 0.3) +
    scale_color_manual(values = c(CS1 = COL_CS1, CS2 = COL_CS2)) +
    scale_fill_manual(values = c(CS1 = COL_CS1, CS2 = COL_CS2)) +
    labs(title = "Subtype PCA (K=2)",
         x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
         color = "Subtype") +
    guides(fill = "none") +
    theme(legend.position = c(0.88, 0.15),
          plot.margin = margin(1, 4, 1, 4, "mm"))
  save_panel_pdf("Fig6e_subtype_PCA.pdf", W_E, H2, print(p6e))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL F: Subtype marker profile heatmap (100x102mm) — TC + PR markers
# — No dendrograms, clear CS1/CS2 grouping, direction annotation
# =============================================================================
cat("\n--- Panel F: Subtype profile heatmap (100x102mm) ---\n")
tryCatch({
  markers <- read.csv(file.path(RES, "phase7_characterization/subtype_top_markers.csv"),
                      stringsAsFactors = FALSE)
  sa <- read.csv(file.path(RES, "phase6_subtyping/subtype_K2.csv"), stringsAsFactors = FALSE)
  colnames(sa)[colnames(sa) == "subtype"] <- "cluster"
  sa$cluster <- gsub("^DS", "CS", sa$cluster)

  tc_mat <- as.matrix(read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"),
                                row.names = 1, check.names = FALSE))
  pr_mat <- as.matrix(read.csv(file.path(DATA, "proteomics_log2_norm.csv"),
                                row.names = 1, check.names = FALSE))
  rownames(pr_mat) <- gsub("\\.[0-9]+$", "", rownames(pr_mat))

  shared <- intersect(sa$sample, intersect(colnames(tc_mat), colnames(pr_mat)))
  sa_sub <- sa[sa$sample %in% shared, ]
  sa_sub <- sa_sub[order(sa_sub$cluster, sa_sub$sample), ]
  shared <- sa_sub$sample

  # TC markers — top 15 by |log2FC|
  tc_m <- markers[markers$omics == "TC", ]
  tc_m <- head(tc_m[order(-abs(tc_m$log2FC)), ], 20)
  tc_feat <- intersect(tc_m$feature, rownames(tc_mat))
  tc_heat <- NULL
  if (length(tc_feat) >= 5) {
    tc_heat <- scale_rows(tc_mat[tc_feat, shared])
    tc_labs <- ensg_to_symbol(rownames(tc_heat))
    keep <- !grepl("^ENS", tc_labs)
    tc_heat <- tc_heat[keep, , drop = FALSE]
    rownames(tc_heat) <- tc_labs[keep]
    tc_heat <- head(tc_heat, 15)
    # Direction annotation for TC
    tc_direction <- ifelse(
      rowMeans(tc_heat[, sa_sub$cluster == "CS1"]) > rowMeans(tc_heat[, sa_sub$cluster == "CS2"]),
      "CS1-high", "CS2-high")
  }

  # PR markers — top 15 by |log2FC|
  pr_m <- markers[markers$omics == "PR", ]
  pr_m <- head(pr_m[order(-abs(pr_m$log2FC)), ], 20)
  pr_feat <- intersect(pr_m$feature, rownames(pr_mat))
  pr_heat <- NULL
  if (length(pr_feat) < 5) {
    m_map <- ensp_to_symbol(rownames(pr_mat))
    names(m_map) <- rownames(pr_mat)
    rev_map <- setNames(names(m_map), m_map)
    pr_feat <- rev_map[intersect(pr_m$feature, m_map)]
    pr_feat <- pr_feat[!is.na(pr_feat)]
  }
  if (length(pr_feat) >= 5) {
    pr_heat <- scale_rows(pr_mat[pr_feat, shared])
    pr_labs <- ensp_to_symbol(rownames(pr_heat))
    keep <- !grepl("^ENS", pr_labs)
    pr_heat <- pr_heat[keep, , drop = FALSE]
    rownames(pr_heat) <- pr_labs[keep]
    pr_heat <- head(pr_heat, 15)
    # Direction annotation for PR
    pr_direction <- ifelse(
      rowMeans(pr_heat[, sa_sub$cluster == "CS1"]) > rowMeans(pr_heat[, sa_sub$cluster == "CS2"]),
      "CS1-high", "CS2-high")
  }

  col_split_raw <- sa_sub$cluster[match(shared, sa_sub$sample)]
  n_cs1 <- sum(col_split_raw == "CS1")
  n_cs2 <- sum(col_split_raw == "CS2")
  # Rename levels to include sample counts for prominent display
  col_split_labeled <- ifelse(col_split_raw == "CS1",
                              sprintf("CS1 (n=%d)", n_cs1),
                              sprintf("CS2 (n=%d)", n_cs2))
  col_split <- factor(col_split_labeled, levels = c(sprintf("CS1 (n=%d)", n_cs1),
                                                     sprintf("CS2 (n=%d)", n_cs2)))
  top_anno <- HeatmapAnnotation(
    Subtype = ifelse(col_split_raw == "CS1", "CS1", "CS2"),
    col = list(Subtype = c(CS1 = COL_CS1, CS2 = COL_CS2)),
    annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    simple_anno_size = unit(5, "mm"),
    annotation_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll()),
    show_annotation_name = TRUE,
    annotation_name_side = "left"
  )

  ht_list <- NULL
  if (!is.null(tc_heat) && nrow(tc_heat) >= 3) {
    # Direction bar annotation (left side) — prominent 5mm
    tc_dir_anno <- rowAnnotation(
      Direction = tc_direction,
      col = list(Direction = c("CS1-high" = COL_CS1, "CS2-high" = COL_CS2)),
      simple_anno_size = unit(5, "mm"),
      annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
      annotation_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll())
    )
    ht_tc <- Heatmap(tc_heat, name = "TC z-score", col = col_zscore,
      cluster_columns = FALSE, cluster_rows = FALSE,
      column_split = col_split,
      show_column_names = FALSE,
      row_names_gp = gp_rn(8), row_names_max_width = unit(28, "mm"),
      column_title_gp = gpar(fontsize = 10, fontfamily = FONT_GRID, fontface = "bold",
                             col = c(COL_CS1, COL_CS2)),
      top_annotation = top_anno,
      left_annotation = tc_dir_anno,
      row_title = "Transcriptomics", row_title_gp = gp_rt(9),
      heatmap_legend_param = std_lp(), rect_gp = gpar(col = "white", lwd = 0.3),
      height = unit(nrow(tc_heat) * 3.2, "mm"),
      column_gap = unit(3.5, "mm"))
    ht_list <- ht_tc
  }
  if (!is.null(pr_heat) && nrow(pr_heat) >= 3) {
    # Direction legend: show if PR is the only heatmap (TC absent), hide if TC already shows it
    pr_show_dir_legend <- is.null(ht_list)
    pr_dir_anno <- rowAnnotation(
      Direction = pr_direction,
      col = list(Direction = c("CS1-high" = COL_CS1, "CS2-high" = COL_CS2)),
      simple_anno_size = unit(5, "mm"),
      annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
      annotation_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll()),
      show_legend = pr_show_dir_legend
    )
    # If TC is absent, PR must show column_split titles; if TC exists, suppress to avoid duplication
    if (is.null(ht_list)) {
      ht_pr <- Heatmap(pr_heat, name = "PR z-score", col = col_zscore,
        cluster_columns = FALSE, cluster_rows = FALSE,
        column_split = col_split,
        show_column_names = FALSE,
        row_names_gp = gp_rn(8), row_names_max_width = unit(28, "mm"),
        column_title_gp = gpar(fontsize = 10, fontfamily = FONT_GRID, fontface = "bold",
                               col = c(COL_CS1, COL_CS2)),
        top_annotation = top_anno,
        left_annotation = pr_dir_anno,
        row_title = "Proteomics", row_title_gp = gp_rt(9),
        heatmap_legend_param = std_lp(), rect_gp = gpar(col = "white", lwd = 0.3),
        height = unit(nrow(pr_heat) * 3.2, "mm"),
        column_gap = unit(3.5, "mm"))
    } else {
      ht_pr <- Heatmap(pr_heat, name = "PR z-score", col = col_zscore,
        cluster_columns = FALSE, cluster_rows = FALSE,
        column_split = col_split,
        show_column_names = FALSE,
        row_names_gp = gp_rn(8), row_names_max_width = unit(28, "mm"),
        column_title = NULL,
        left_annotation = pr_dir_anno,
        row_title = "Proteomics", row_title_gp = gp_rt(9),
        heatmap_legend_param = std_lp(), rect_gp = gpar(col = "white", lwd = 0.3),
        height = unit(nrow(pr_heat) * 3.2, "mm"),
        column_gap = unit(3.5, "mm"))
    }
    if (is.null(ht_list)) ht_list <- ht_pr else ht_list <- ht_list %v% ht_pr
  }

  if (!is.null(ht_list)) {
    draw_f <- function() draw(ht_list, merge_legend = FALSE,
         heatmap_legend_side = "right", annotation_legend_side = "right",
         column_title = "Subtype marker profile",
         column_title_gp = gpar(fontsize = FS_TITLE, fontface = "bold", fontfamily = FONT_FAMILY),
         padding = unit(c(1, 2, 1, 2), "mm"))
    fp <- file.path(OUT, "Fig6f_subtype_profile.pdf")
    cairo_pdf(fp, width = W_F / MM, height = H3 / MM, family = FONT_FAMILY)
    draw_f()
    dev.off()
    cat(sprintf("  -> Fig6f_subtype_profile.pdf (%.0f x %.0f mm)\n", W_F, H3))
    grob_f <<- grid.grabExpr(draw_f())
  }
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# PANEL G: Subtype characterisation — Pathways + Immune cells (83x102mm)
# — No dendrograms, compact horizontal legends at bottom, maximise data area
# =============================================================================
cat("\n--- Panel G: Comprehensive characterisation (83x102mm) ---\n")
tryCatch({
  pw <- read.csv(file.path(RES, "phase7_characterization/subtype_pathway_scores.csv"),
                 stringsAsFactors = FALSE)
  imm <- read.csv(file.path(RES, "phase7_characterization/subtype_immune_scores.csv"),
                  stringsAsFactors = FALSE)

  # Top 10 pathways by |diff|
  pw$diff_val <- pw$diff_CS1_CS2
  pw_sig <- pw[pw$padj < 0.2, ]
  if (nrow(pw_sig) < 8) pw_sig <- head(pw[order(pw$pvalue), ], 10)
  pw_sig <- head(pw_sig[order(-abs(pw_sig$diff_val)), ], 10)

  pw_mat <- matrix(c(pw_sig$mean_CS1, pw_sig$mean_CS2), ncol = 2,
                   dimnames = list(pw_sig$pathway, c("CS1", "CS2")))
  pw_mat_s <- t(scale(t(pw_mat))); pw_mat_s[is.na(pw_mat_s)] <- 0
  pw_mat_s[pw_mat_s > 2] <- 2; pw_mat_s[pw_mat_s < -2] <- -2
  rownames(pw_mat_s) <- gsub("^HALLMARK_", "", rownames(pw_mat_s))
  rownames(pw_mat_s) <- gsub("_", " ", rownames(pw_mat_s))
  rownames(pw_mat_s) <- tools::toTitleCase(tolower(rownames(pw_mat_s)))

  # Significance annotation for pathways
  pw_stars <- ifelse(pw_sig$padj < 0.01, "***",
              ifelse(pw_sig$padj < 0.05, "**",
              ifelse(pw_sig$padj < 0.2, "*", "")))
  pw_dir <- ifelse(pw_sig$diff_val > 0, "CS1-high", "CS2-high")

  # Top 8 immune cell types (reduced from 10 to fit)
  imm$diff_val <- imm$diff
  imm_sig <- imm[imm$padj < 0.3, ]
  if (nrow(imm_sig) < 5) imm_sig <- head(imm[order(imm$pvalue), ], 8)
  imm_sig <- head(imm_sig[order(-abs(imm_sig$diff_val)), ], 8)

  imm_mat <- matrix(c(imm_sig$mean_CS1, imm_sig$mean_CS2), ncol = 2,
                    dimnames = list(imm_sig$cell_type, c("CS1", "CS2")))
  imm_mat_s <- t(scale(t(imm_mat))); imm_mat_s[is.na(imm_mat_s)] <- 0
  imm_mat_s[imm_mat_s > 2] <- 2; imm_mat_s[imm_mat_s < -2] <- -2
  rownames(imm_mat_s) <- gsub("_", " ", rownames(imm_mat_s))

  imm_stars <- ifelse(imm_sig$padj < 0.01, "***",
               ifelse(imm_sig$padj < 0.05, "**",
               ifelse(imm_sig$padj < 0.3, "*", "")))
  imm_dir <- ifelse(imm_sig$diff_val > 0, "CS1-high", "CS2-high")

  # Sort pathways by diff_val descending (CS1-high at top → CS2-high at bottom)
  pw_order <- order(-pw_sig$diff_val)
  pw_mat_s <- pw_mat_s[pw_order, , drop = FALSE]
  pw_dir <- pw_dir[pw_order]
  pw_stars <- pw_stars[pw_order]

  # Sort immune by diff_val descending
  imm_order <- order(-imm_sig$diff_val)
  imm_mat_s <- imm_mat_s[imm_order, , drop = FALSE]
  imm_dir <- imm_dir[imm_order]
  imm_stars <- imm_stars[imm_order]

  col_anno <- HeatmapAnnotation(
    Subtype = c("CS1", "CS2"),
    col = list(Subtype = c(CS1 = COL_CS1, CS2 = COL_CS2)),
    annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    simple_anno_size = unit(4, "mm"),
    annotation_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll(),
                                   nrow = 1, direction = "horizontal"),
    show_annotation_name = TRUE
  )

  # Pathway heatmap — Direction annotation bar (5mm) on right, sorted high→low
  pw_right_anno <- rowAnnotation(
    Direction = pw_dir,
    col = list(Direction = c("CS1-high" = COL_CS1, "CS2-high" = COL_CS2)),
    simple_anno_size = unit(5, "mm"),
    annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    annotation_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll(),
                                   nrow = 1, direction = "horizontal")
  )
  ht_pw <- Heatmap(pw_mat_s, name = "z-score", col = col_div,
    cluster_columns = FALSE, cluster_rows = FALSE,
    show_column_names = TRUE, show_heatmap_legend = TRUE,
    column_names_gp = gpar(fontsize = 9, fontfamily = FONT_GRID, fontface = "bold"),
    row_names_gp = gp_rn(8),
    row_names_max_width = unit(38, "mm"), top_annotation = col_anno,
    right_annotation = pw_right_anno,
    row_title = "Pathways", row_title_gp = gp_rt(9),
    heatmap_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll(),
                                direction = "horizontal", legend_width = unit(18, "mm")),
    rect_gp = gpar(col = "white", lwd = 0.5),
    width = unit(22, "mm"), height = unit(nrow(pw_mat_s) * 3.8, "mm"))

  # Immune heatmap — Direction annotation bar (5mm), no duplicate legend
  imm_right_anno <- rowAnnotation(
    Direction = imm_dir,
    col = list(Direction = c("CS1-high" = COL_CS1, "CS2-high" = COL_CS2)),
    simple_anno_size = unit(5, "mm"),
    annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    show_legend = FALSE
  )
  ht_imm <- Heatmap(imm_mat_s, name = "Immune z-score", col = col_div,
    cluster_columns = FALSE, cluster_rows = FALSE,
    show_column_names = TRUE, show_heatmap_legend = FALSE,
    column_names_gp = gpar(fontsize = 9, fontfamily = FONT_GRID, fontface = "bold"),
    row_names_gp = gp_rn(8),
    row_names_max_width = unit(38, "mm"),
    right_annotation = imm_right_anno,
    row_title = "Immune cells", row_title_gp = gp_rt(9),
    rect_gp = gpar(col = "white", lwd = 0.5),
    width = unit(22, "mm"), height = unit(nrow(imm_mat_s) * 3.8, "mm"))

  ht_g <- ht_pw %v% ht_imm

  draw_g <- function() draw(ht_g, merge_legend = FALSE,
       heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       column_title = "Subtype characterisation",
       column_title_gp = gpar(fontsize = FS_TITLE, fontface = "bold", fontfamily = FONT_FAMILY),
       padding = unit(c(1, 2, 1, 2), "mm"))
  fp <- file.path(OUT, "Fig6g_comprehensive_heatmap.pdf")
  cairo_pdf(fp, width = W_G / MM, height = H3 / MM, family = FONT_FAMILY)
  draw_g()
  dev.off()
  cat(sprintf("  -> Fig6g_comprehensive_heatmap.pdf (%.0f x %.0f mm)\n", W_G, H3))
  grob_g <<- grid.grabExpr(draw_g())
}, error = function(e) cat("  ERROR:", e$message, "\n"))

# =============================================================================
# Composite Assembly — pure VECTOR grid viewports (Adobe Illustrator-editable PDF)
# All panels (ggplot + ComplexHeatmap grobs) are placed via mm-based viewports.
# No rasterization for PDF; PNG/TIFF rasterized only at the cairo device level.
# =============================================================================
cat(sprintf("\n--- Assembling Figure_6 (vector grid, %.0fx%.0fmm, %d DPI) ---\n",
            W_TOTAL, H_TOTAL, ASSEMBLY_DPI))
tryCatch({
  # Panel layout (top-left origin in mm): list(x, y, w, h, content)
  layout_specs <- list(
    list(x = 0,    y = 0,        w = W_A,  h = H1,  obj = p6a),
    list(x = W_A,  y = 0,        w = W_BD, h = H_B, obj = p6b),
    list(x = W_A,  y = H_B,      w = W_BD, h = H_D, obj = p6d),
    list(x = 0,    y = H1,       w = W_C,  h = H2,  obj = p6c),
    list(x = W_C,  y = H1,       w = W_E,  h = H2,  obj = p6e),
    list(x = 0,    y = H1 + H2,  w = W_F,  h = H3,  obj = grob_f),
    list(x = W_F,  y = H1 + H2,  w = W_G,  h = H3,  obj = grob_g)
  )

  tag_labels <- c("A", "B", "C", "D", "E", "F", "G")
  tag_x_mm <- c(1, W_A + 1, 1, W_A + 1, W_C + 1, 1, W_F + 1)
  tag_y_mm <- c(1, 1, H1 + 1, H_B + 1, H1 + 1, H1 + H2 + 1, H1 + H2 + 1)

  render_final <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(width = unit(W_TOTAL, "mm"),
                                      height = unit(H_TOTAL, "mm"),
                                      xscale = c(0, W_TOTAL),
                                      yscale = c(0, H_TOTAL)))
    for (sp in layout_specs) {
      if (is.null(sp$obj)) next
      vp <- grid::viewport(
        x = unit(sp$x + sp$w / 2, "mm"),
        y = unit(H_TOTAL - sp$y - sp$h / 2, "mm"),
        width = unit(sp$w, "mm"),
        height = unit(sp$h, "mm")
      )
      if (inherits(sp$obj, "ggplot")) {
        print(sp$obj, vp = vp)
      } else {
        grid::pushViewport(vp)
        grid::grid.draw(sp$obj)
        grid::popViewport()
      }
    }
    for (i in seq_along(tag_labels)) {
      grid::grid.text(label = tag_labels[i],
        x = unit(tag_x_mm[i], "mm"),
        y = unit(H_TOTAL - tag_y_mm[i], "mm"),
        just = c("left", "top"),
        gp = gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY))
    }
    grid::popViewport()
  }

  # Pixel dimensions for raster outputs
  DPI <- ASSEMBLY_DPI
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)

  # Save PDF — fully vector, AI-editable
  cairo_pdf(file.path(OUT, "Figure_6.pdf"),
            width = W_TOTAL / MM, height = H_TOTAL / MM, family = FONT_FAMILY)
  render_final()
  dev.off()
  cat("  -> Figure_6.pdf (vector, AI-editable)\n")

  # Save PNG — vector rendered at 600 DPI by cairo
  grDevices::png(file.path(OUT, "Figure_6.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_final()
  dev.off()
  cat("  -> Figure_6.png\n")

  # Save TIFF — vector rendered at 600 DPI by cairo, LZW
  grDevices::tiff(file.path(OUT, "Figure_6.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_final()
  dev.off()
  cat("  -> Figure_6.tiff\n")

  cat(sprintf("  Figure_6 DONE (%.0fx%.0fmm, vector grid assembly, %d DPI)\n",
              W_TOTAL, H_TOTAL, DPI))
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Figure 6 rendering complete ===\n")
