#!/usr/bin/env Rscript
# =============================================================================
# panel_factories.R -- Reusable Factory Functions + Lazy Data Cache
# =============================================================================
# Extracted from nc_production_pipeline.R SECTION 1-3
# Dependencies: nc_theme.R (must be sourced first)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(RColorBrewer)
  library(patchwork)
  library(stringr)
})

# =============================================================================
# 1. theme_nc -- Nature Communications theme (theme_bw base, rectangular border)
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks         = element_line(linewidth = 0.3, color = "black"),
    axis.text          = element_text(size = 7, color = "black"),
    axis.title         = element_text(size = 8, face = "bold", color = "black"),
    plot.title         = element_text(size = 8.5, face = "bold", hjust = 0),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.key.size    = unit(2.5, "mm"),
    legend.background  = element_blank(),
    legend.box.background = element_blank(),
    legend.spacing.y   = unit(1, "mm"),
    strip.text         = element_text(size = 7.5, face = "bold"),
    strip.background   = element_rect(fill = "grey96", linewidth = 0.3, color = "grey80"),
    plot.margin        = margin(2, 2, 2, 2, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# 2. Factory Functions (12 reusable plot builders)
# =============================================================================

# --- 2.1 Volcano Plot Factory ---
factory_volcano <- function(df, title, x_col = "logFC", y_col = "P.Value",
                            fc_cut = 1.0, p_cut = 0.05, n_label = 15) {
  df <- as.data.frame(df)
  if ("significance" %in% colnames(df) && all(c("Up", "Down", "NS") %in% unique(df$significance))) {
    df$sig <- factor(df$significance, levels = c("Up", "Down", "NS"))
  } else {
    df$sig <- "NS"
    df$sig[df[[x_col]] >  fc_cut & df[[y_col]] < p_cut] <- "Up"
    df$sig[df[[x_col]] < -fc_cut & df[[y_col]] < p_cut] <- "Down"
    df$sig <- factor(df$sig, levels = c("Up", "Down", "NS"))
  }
  df$neg_log10_p <- -log10(df[[y_col]] + 1e-300)
  df_sig <- df[df$sig != "NS", ]
  label_genes <- if (nrow(df_sig) > 0) head(df_sig[order(df_sig[[y_col]]), ], n_label) else data.frame()
  gene_col <- find_col(df, c("gene_name", "gene", "Gene", "symbol", "name", "ID"))
  if (is.null(gene_col)) { df$rowname_id <- rownames(df); gene_col <- "rowname_id" }
  df[[gene_col]] <- ensg_to_symbol(df[[gene_col]])
  if (nrow(label_genes) > 0) {
    label_genes[[gene_col]] <- ensg_to_symbol(label_genes[[gene_col]])
    label_genes <- label_genes[!grepl("^ENSG", label_genes[[gene_col]]), ]
    # NC upgrade: adaptive label size based on gene name length
    label_genes$label_size <- ifelse(nchar(label_genes[[gene_col]]) > 8, 2, 2.5)
  }
  p <- ggplot(df, aes(x = .data[[x_col]], y = neg_log10_p, color = sig)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_manual(values = c(Up = COL_UP, Down = COL_DOWN, NS = COL_NS), drop = FALSE) +
    # NC upgrade: threshold lines with grey30 color and linewidth 0.4
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.4, color = "grey30") +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", linewidth = 0.4, color = "grey30") +
    labs(title = title, x = expression(log[2]~fold~change),
         y = expression(-log[10]~italic(P)), color = NULL) +
    theme(legend.position = c(0.85, 0.92))
  if (nrow(label_genes) > 0 && gene_col %in% colnames(label_genes)) {
    # NC upgrade: optimized ggrepel parameters for better label placement
    p <- p + geom_text_repel(
      data = label_genes, aes(label = .data[[gene_col]], size = label_size),
      fontface = "italic", family = FONT_FAMILY,
      force = 1.5, force_pull = 0.3,
      segment.size = 0.2, segment.color = "grey50",
      max.overlaps = 12, min.segment.length = 0.2,
      box.padding = 0.3, point.padding = 0.2) +
      scale_size_identity()
  }
  n_up   <- sum(df$sig == "Up", na.rm = TRUE)
  n_down <- sum(df$sig == "Down", na.rm = TRUE)
  p + annotate("text", x = Inf, y = Inf, label = paste0("Up: ", n_up),
               hjust = 1.1, vjust = 1.5, size = 3, color = COL_UP,
               fontface = "bold", family = FONT_FAMILY) +
    annotate("text", x = -Inf, y = Inf, label = paste0("Down: ", n_down),
             hjust = -0.1, vjust = 1.5, size = 3, color = COL_DOWN,
             fontface = "bold", family = FONT_FAMILY)
}

# --- 2.2 PCA Scatter Factory ---
factory_pca_scatter <- function(mat, title) {
  grp <- make_group_vector(colnames(mat))
  pca <- prcomp(t(mat), scale. = TRUE)
  ve <- summary(pca)$importance[2, 1:2] * 100
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = grp)
  ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 1.5, alpha = 0.8) +
    # NC upgrade: dashed ellipse with linewidth 0.5 for visual distinction
    stat_ellipse(type = "norm", level = 0.95, linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme(legend.position = c(0.82, 0.15))
}

# --- 2.3 Correlation Heatmap Factory ---
factory_cor_heatmap <- function(mat, title) {
  cor_mat <- cor(mat, use = "pairwise.complete.obs")
  grp <- make_group_vector(colnames(mat))
  ha <- HeatmapAnnotation(Group = grp,
    col = list(Group = c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)),
    annotation_name_gp = gp_anno_name(), show_legend = FALSE)
  # NC upgrade: explicit #FFFFFF for center color instead of "white"
  col_cor_qc <- colorRamp2(c(0.7, 0.85, 1), c(COL_DOWN, "#FFFFFF", COL_UP))
  Heatmap(cor_mat, name = "Pearson R", col = col_cor_qc,
    top_annotation = ha, show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    column_title = title,
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
    # NC upgrade: white grid lines with lwd 0.5
    border = TRUE, rect_gp = gpar(col = "#FFFFFF", lwd = 0.5),
    use_raster = FALSE,  # AI compatibility: vector cells for editing
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY, fontface = "bold"),
      labels_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
      grid_height = unit(3, "mm"), grid_width = unit(3, "mm")))
}

# --- 2.4 DE Expression Heatmap Factory ---
factory_expr_heatmap <- function(expr_mat, sig_csv, n_top = 30, title) {
  sig <- read.csv(sig_csv, stringsAsFactors = FALSE)
  fc_col <- find_col(sig, c("logFC", "log2FoldChange", "fc"))
  id_candidates_all <- c("gene_id", "Protein", "Compound_ID", "name", "gene_name",
                          "gene", "X", "metabolite", "metabolite_name", "ID")
  best_n <- 0; id_col <- NULL
  for (cand in id_candidates_all) {
    if (cand %in% colnames(sig)) {
      ids_try <- head(sig[order(abs(sig[[fc_col]]), decreasing = TRUE), cand], n_top * 2)
      n_match <- sum(ids_try %in% rownames(expr_mat))
      if (n_match > best_n) { best_n <- n_match; id_col <- cand }
    }
  }
  if (best_n < 5) return(NULL)
  sig_sorted <- sig[order(abs(sig[[fc_col]]), decreasing = TRUE), ]
  ids <- head(sig_sorted[[id_col]], n_top * 2)
  matched <- head(ids[ids %in% rownames(expr_mat)], n_top)
  mat <- scale_rows(expr_mat[matched, ], 2)
  name_col <- find_col(sig, c("gene_name", "metabolite_name", "name"))
  if (!is.null(name_col) && name_col != id_col) {
    name_map <- setNames(sig[[name_col]], sig[[id_col]])
    new_names <- name_map[rownames(mat)]
    new_names[is.na(new_names) | new_names == ""] <- rownames(mat)[is.na(new_names) | new_names == ""]
    rownames(mat) <- make.unique(new_names)
  }
  rownames(mat) <- ensg_to_symbol(rownames(mat))
  rownames(mat) <- ifelse(nchar(rownames(mat)) > 30,
                          paste0(substr(rownames(mat), 1, 27), "..."),
                          rownames(mat))
  ha <- make_group_annotation(colnames(mat))
  # NC upgrade: row_names fontsize >= 6.5, white grid lines, explicit legend params
  Heatmap(mat, name = "Z-score", col = col_div, top_annotation = ha,
    show_row_names = TRUE, row_names_gp = gp_row_names(6),
    show_column_names = FALSE, border = TRUE,
    rect_gp = gpar(col = "#FFFFFF", lwd = 0.5),
    use_raster = FALSE,  # AI compatibility: vector cells for editing
    show_row_dend = FALSE, show_column_dend = FALSE,
    column_title = title,
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
    width = unit(50, "mm"),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6.5, fontfamily = FONT_FAMILY),
      title_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY, fontface = "bold"),
      legend_height = unit(20, "mm"),
      grid_width = unit(3, "mm")))
}

# --- 2.5 GSEA Dotplot Factory ---
factory_gsea_dotplot <- function(csv_path, title, n = 15) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  name_col <- find_col(df, c("Description", "pathway_name", "pathway", "term"))
  if (is.null(name_col)) name_col <- find_col(df, c("ID"))
  p_col <- find_col(df, c("pval", "pvalue", "p.adjust", "padj", "qvalue", "P.Value"))
  size_col <- find_col(df, c("setSize", "Count", "count", "size"))
  if (is.null(name_col) || is.null(p_col)) return(NULL)
  df$pathway_name <- clean_pathway(df[[name_col]])
  # NC upgrade: truncate pathway names > 30 characters
  df$pathway_name <- ifelse(nchar(df$pathway_name) > 30,
                            paste0(substr(df$pathway_name, 1, 27), "..."),
                            df$pathway_name)
  df$neg_log10p <- -log10(df[[p_col]] + 1e-300)
  df <- head(df[order(df$neg_log10p, decreasing = TRUE), ], n)
  p <- ggplot(df, aes(x = neg_log10p, y = reorder(pathway_name, neg_log10p)))
  if (!is.null(size_col)) {
    p <- p + geom_point(aes(size = .data[[size_col]]), color = COL_TC, alpha = 0.8) +
      scale_size_continuous(range = c(2, 5), name = "Gene count")
  } else {
    p <- p + geom_point(size = 2.5, color = COL_TC, alpha = 0.8)
  }
  p + labs(title = title, x = expression(-log[10]~P), y = NULL) +
    theme(axis.text.y = element_text(size = 6, family = FONT_FAMILY))
}

# --- 2.6 GSEA Barplot Factory ---
factory_gsea_barplot <- function(csv_path, title, pattern = NULL, n = 15) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  name_col <- find_col(df, c("Description", "pathway_name", "pathway_label", "pathway", "term", "ID"))
  nes_col <- find_col(df, c("NES", "NES_TC", "enrichmentScore"))
  p_col <- find_col(df, c("pval", "pvalue", "padj", "padj_TC", "p.adjust", "P.Value"))
  if (is.null(name_col) || is.null(p_col)) return(NULL)
  if (!is.null(pattern)) df <- df[grepl(pattern, df[[name_col]], ignore.case = TRUE), ]
  if (nrow(df) == 0) return(NULL)
  df[[name_col]] <- clean_pathway(df[[name_col]])
  # NC upgrade: truncate pathway names > 30 characters
  df[[name_col]] <- ifelse(nchar(df[[name_col]]) > 30,
                           paste0(substr(df[[name_col]], 1, 27), "..."),
                           df[[name_col]])
  df$neg_log10p <- -log10(df[[p_col]] + 1e-300)
  df <- head(df[order(df$neg_log10p, decreasing = TRUE), ], n)
  df$direction <- if (!is.null(nes_col)) ifelse(df[[nes_col]] > 0, "Up", "Down") else "Up"
  ggplot(df, aes(x = neg_log10p, y = reorder(.data[[name_col]], neg_log10p), fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), guide = "none") +
    labs(title = title, x = expression(-log[10]~P), y = NULL) +
    theme(axis.text.y = element_text(size = 6, family = FONT_FAMILY))
}

# --- 2.7 WGCNA Module-Trait Heatmap Factory ---
factory_wgcna_trait <- function(cor_csv, pval_csv, title, outfile) {
  cor_mat <- as.matrix(read.csv(cor_csv, row.names = 1, check.names = FALSE))
  pv_mat <- as.matrix(read.csv(pval_csv, row.names = 1, check.names = FALSE))
  colnames(cor_mat) <- clean_varname(colnames(cor_mat))
  colnames(pv_mat) <- clean_varname(colnames(pv_mat))
  sig_txt <- matrix("", nrow(pv_mat), ncol(pv_mat))
  sig_txt[pv_mat < 0.001] <- "***"
  sig_txt[pv_mat >= 0.001 & pv_mat < 0.01] <- "**"
  sig_txt[pv_mat >= 0.01 & pv_mat < 0.05] <- "*"
  s <- sp(outfile, W_HALF, H_TALL)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  ht <- Heatmap(cor_mat, name = "Cor", col = col_cor,
    cell_fun = function(j, i, x, y, w, h, fill) grid.text(sig_txt[i, j], x, y, gp = gp_cell_text()),
    row_names_gp = gp_row_names(7), column_names_gp = gp_col_names(7, bold = FALSE),
    column_names_rot = 45, border = TRUE, rect_gp = gp_border(),
    use_raster = FALSE,  # AI compatibility: vector cells for editing
    show_row_dend = FALSE, show_column_dend = FALSE,
    column_title = title,
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
    heatmap_legend_param = std_legend_param())
  draw(ht, padding = unit(c(5, 15, 5, 5), "mm")); dev.off()
  cat("  ->", outfile, "\n")
}

# --- 2.8 Hub Gene Scatter Factory ---
factory_hub_scatter <- function(hub_csv, title, n_top = 15) {
  hubs <- read.csv(file.path(RES, hub_csv), stringsAsFactors = FALSE)
  hubs$gene <- ensg_to_symbol(hubs$gene)
  hubs$gene <- ensp_to_symbol(hubs$gene)
  hubs$abs_kME <- abs(hubs$kME); hubs$abs_GS <- abs(hubs$GS_group)
  top <- head(hubs[order(hubs$hub_score, decreasing = TRUE), ], n_top)
  top <- top[!grepl("^ENS[GP]", top$gene), ]
  # NC upgrade: add adaptive label size
  top$label_size <- ifelse(nchar(top$gene) > 8, 2, 2.5)
  ggplot(hubs, aes(x = abs_kME, y = abs_GS)) +
    geom_point(size = 0.8, alpha = 0.4, color = "grey60") +
    geom_point(data = top, size = 2, color = COL_TC) +
    geom_text_repel(data = top, aes(label = gene, size = label_size),
                    fontface = "italic", family = FONT_FAMILY,
                    force = 2, force_pull = 0.5,
                    max.overlaps = 15, segment.size = 0.2,
                    segment.color = "grey50") +
    scale_size_identity() +
    labs(title = title, x = "|Module membership|", y = "|Gene significance|")
}

# --- 2.9 Direction Bar Factory ---
factory_direction_bar <- function(df, name_col, value_col, title,
                                  colors = c(Up = COL_UP, Down = COL_DOWN),
                                  x_lab = NULL, name_size = 7) {
  dir_labels <- names(colors)
  df$direction <- ifelse(df[[value_col]] > 0, dir_labels[1], dir_labels[2])
  if (is.null(x_lab)) x_lab <- value_col
  ggplot(df, aes(x = .data[[value_col]], y = reorder(.data[[name_col]], .data[[value_col]]),
                 fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = colors, guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(title = title, x = x_lab, y = NULL) +
    theme(axis.text.y = element_text(size = name_size))
}

# --- 2.10 Horizontal Bar Factory ---
factory_barplot_horizontal <- function(df, name_col, value_col, title,
                                       fill_col = PAL_CAT[1], x_lab = NULL,
                                       name_size = 7) {
  if (is.null(x_lab)) x_lab <- value_col
  ggplot(df, aes(x = .data[[value_col]], y = reorder(.data[[name_col]], .data[[value_col]]))) +
    geom_col(fill = fill_col, width = 0.7) +
    labs(title = title, x = x_lab, y = NULL) +
    theme(axis.text.y = element_text(size = name_size))
}

# --- 2.11 Forest Plot Factory ---
factory_forest_plot <- function(df, title, n_genes = 6) {
  top_genes <- unique(df$gene)[1:min(n_genes, length(unique(df$gene)))]
  fs <- df[df$gene %in% top_genes, ]
  if (!"ci_lower" %in% colnames(fs)) {
    fs$ci_lower <- fs$logFC - 1.96 * fs$SE
    fs$ci_upper <- fs$logFC + 1.96 * fs$SE
  }
  ggplot(fs, aes(x = logFC, y = dataset, color = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2,
                  position = position_dodge(width = 0.6), linewidth = 0.5,
                  orientation = "y") +
    scale_color_manual(values = PAL_CAT[seq_along(top_genes)]) +
    labs(title = title, x = expression(log[2]~FC), y = NULL, color = "Gene") +
    theme(axis.text.y = element_text(size = 7.5),
          legend.text = element_text(size = 7, face = "italic"))
}

# --- 2.12 Dual NES Bar Factory ---
factory_dual_nes_bar <- function(csv_path, title) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  pw_col <- find_col(df, c("pathway", "Description", "term"))
  nes_cols <- grep("NES", colnames(df), value = TRUE)
  if (is.null(pw_col) || length(nes_cols) < 2) return(NULL)
  df$pw <- clean_pathway(df[[pw_col]])
  # NC upgrade: truncate pathway names > 30 characters
  df$pw <- ifelse(nchar(df$pw) > 30, paste0(substr(df$pw, 1, 27), "..."), df$pw)
  df_long <- data.frame(
    pathway = rep(df$pw, 2),
    NES = c(df[[nes_cols[1]]], df[[nes_cols[2]]]),
    Omics = rep(c("Transcriptomics", "Proteomics"), each = nrow(df)),
    stringsAsFactors = FALSE)
  df_long <- df_long[!is.na(df_long$NES), ]
  ggplot(df_long, aes(x = NES, y = reorder(pathway, NES), fill = Omics)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c(Transcriptomics = COL_TC, Proteomics = COL_PR)) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(title = title, x = "NES", y = NULL) +
    theme(axis.text.y = element_text(size = 6, family = FONT_FAMILY))
}

# =============================================================================
# 3. Lazy-Cached Data Loading
# =============================================================================
.data_cache <- new.env(parent = emptyenv())

cache_get <- function(key, loader) {
  if (!exists(key, envir = .data_cache)) {
    cat("  [cache] Loading", key, "...\n")
    assign(key, loader(), envir = .data_cache)
  }
  get(key, envir = .data_cache)
}

get_tc_mat <- function() cache_get("tc_mat", function()
  read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"), row.names = 1, check.names = FALSE))

get_pr_mat <- function() cache_get("pr_mat", function()
  read.csv(file.path(DATA, "proteomics_log2_norm.csv"), row.names = 1, check.names = FALSE))

get_mb_mat <- function() cache_get("mb_mat", function()
  read.csv(file.path(DATA, "metabolomics_log2_merged.csv"), row.names = 1, check.names = FALSE))

cat("[panel_factories.R] 12 factories + cache loaded.\n")
