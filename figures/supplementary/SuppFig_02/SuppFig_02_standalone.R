#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_02_standalone.R
# Supplementary Figure 2: Top DE Features + mRNA-Protein Discordance
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Nature Communications - vector PDF + PNG/TIFF (600 DPI) composite
# =============================================================================
# 8 panels (3-row layout, 183x268mm):
#   a = Volcano plot - Transcriptomics (DEGs, FC>log2(1.5), P<0.05)
#   b = Volcano plot - Proteomics (DEPs, |logFC|>log2(1.5), P<0.05)
#   c = Volcano plot - Metabolomics (DEMs, |logFC|>log2(1.5), P<0.05)
#   d = mRNA-protein Spearman correlation density (n=8,011 genes; 523 sig padj<0.05)
#   e = Per-gene mRNA(VST) vs Protein(log2) scatter for top 3 positive
#       (CYP1A2 / GSTA2 / EFHD1) and top 3 negative (CMC4 / LURAP1L / RPSA) rho
#   f = Top 20 DEGs horizontal logFC bar chart (rank: |logFC|, P<0.05)
#   g = Top 20 DEPs horizontal logFC bar chart
#   h = Top 20 DEMs horizontal logFC bar chart
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_02_standalone.R
# =============================================================================

cat("=== Supplementary Figure 2: Top DE Features + mRNA-Protein Discordance ===\n")
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
  library(patchwork)
  library(ggsci)
  library(ggpubr)
  library(magick)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_02")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4

W_SINGLE <- 89; W_DOUBLE <- 183; W_HALF <- 89
H_STD <- 85; H_TALL <- 100; H_MAX <- 240

ASSEMBLY_DPI <- 600   # pixel-exact rendering constant

mm2in <- function(mm) mm / MM_PER_INCH
FS_TITLE      <- 10
FS_AXIS_TITLE <- 9
FS_AXIS_TEXT  <- 8
FS_LEGEND_T   <- 8
FS_LEGEND_L   <- 8
FS_TAG        <- 12
FS_GEOM_TEXT  <- 2.82

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
COL_UP       <- "#CD534CFF"
COL_DOWN     <- "#0073C2FF"
COL_NS       <- "#868686FF"
COL_NA       <- "#F0F0F0"
COL_NORMAL   <- "#7AA6DCFF"
COL_ADJACENT <- "#CD534CFF"
COL_TC       <- "#0073C2FF"
COL_PR       <- "#CD534CFF"
COL_MT       <- "#EFC000FF"
PAL_CAT      <- pal_jco("default")(10)

col_div <- colorRamp2(c(-2, 0, 2), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 5: Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(2.5, "mm"),
    legend.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# SECTION 6: ComplexHeatmap gpar Factories
# =============================================================================
gp_anno_name    <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID)
gp_legend_title <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
gp_legend_labels <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID)
ht_opt$message <- FALSE

# =============================================================================
# SECTION 7: Helper Functions
# =============================================================================
sp <- function(filename, w_mm, h_mm) {
  list(path = file.path(OUT, filename),
       width = w_mm / MM_PER_INCH, height = h_mm / MM_PER_INCH)
}

save_pdf <- function(filename, w_mm, h_mm, expr) {
  s <- sp(filename, w_mm, h_mm)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR saving ", filename, ": ", e$message))
  dev.off()
  cat("  ->", basename(s$path), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
  invisible(s$path)
}

# --- Pixel-exact panel rendering (figure-assembly v2 pattern) ---
save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR saving ", filename, ": ", e$message))
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm pixel-exact)\n", w_mm, h_mm))
  invisible(fpath)
}

save_heatmap_panel_pdf <- function(filename, ht, w_mm, h_mm, title = NULL) {
  fpath <- file.path(OUT, filename)
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  if (!is.null(title)) {
    draw(ht, padding = unit(c(2, 2, 12, 2), "mm"))
    grid::grid.text(title, x = unit(0.5, "npc"), y = unit(1, "npc") - unit(4, "mm"),
                    just = c("center", "top"),
                    gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY))
  } else {
    draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  }
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
  invisible(fpath)
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else NULL
}

clamp_matrix <- function(mat, lim = 2) {
  mat[mat > lim] <- lim; mat[mat < -lim] <- -lim; mat
}

scale_rows <- function(mat, lim = 2) {
  mat <- t(scale(t(as.matrix(mat)))); mat[is.na(mat)] <- 0; clamp_matrix(mat, lim)
}

make_group_vector <- function(sample_names) {
  ifelse(grepl("Normal|^N\\d", sample_names), "Normal", "Adjacent")
}

make_group_annotation <- function(sample_names) {
  grp <- make_group_vector(sample_names)
  HeatmapAnnotation(
    Group = grp,
    col = list(Group = c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)),
    annotation_name_gp = gp_anno_name(),
    simple_anno_size = unit(3, "mm"),
    annotation_legend_param = list(title_gp = gp_legend_title(), labels_gp = gp_legend_labels())
  )
}

# --- ENSG -> gene symbol mapping ---
.ensg_map <- NULL
build_ensg_map <- function() {
  if (!is.null(.ensg_map)) return(.ensg_map)
  map <- data.frame(ensg = character(), symbol = character(), stringsAsFactors = FALSE)
  tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    db_map <- AnnotationDbi::select(org.Hs.eg.db,
      keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
      columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
    db_map <- db_map[!is.na(db_map$SYMBOL) & !duplicated(db_map$ENSEMBL), ]
    map <- rbind(map, data.frame(ensg = db_map$ENSEMBL, symbol = db_map$SYMBOL, stringsAsFactors = FALSE))
  }, error = function(e) cat("  ensg_map: org.Hs.eg.db unavailable\n"))
  for (f in c("phase1_diff/DEGs_Adjacent_vs_Normal.csv", "phase1_diff/DEPs_Adjacent_vs_Normal.csv")) {
    fp <- file.path(RES, f)
    if (file.exists(fp)) {
      d <- read.csv(fp, stringsAsFactors = FALSE)
      id_col <- find_col(d, c("gene_id", "Protein", "gene"))
      nm_col <- find_col(d, c("gene_name", "symbol", "Gene"))
      if (!is.null(id_col) && !is.null(nm_col)) {
        sub <- data.frame(ensg = d[[id_col]], symbol = d[[nm_col]], stringsAsFactors = FALSE)
        sub <- sub[!is.na(sub$symbol) & sub$symbol != "" & sub$symbol != sub$ensg, ]
        map <- rbind(map, sub)
      }
    }
  }
  map <- map[!duplicated(map$ensg), ]
  .ensg_map <<- map; map
}

ensg_to_symbol <- function(ids) {
  m <- build_ensg_map()
  ids_clean <- gsub("\\.[0-9]+$", "", ids); ids_clean <- gsub("_[0-9]+$", "", ids_clean)
  idx <- match(ids_clean, m$ensg)
  ifelse(is.na(idx), ids, m$symbol[idx])
}

# --- ENSP -> gene symbol mapping ---
.ensp_map <- NULL
build_ensp_map <- function() {
  if (!is.null(.ensp_map)) return(.ensp_map)
  map <- data.frame(ensp = character(), symbol = character(), stringsAsFactors = FALSE)
  for (f in c("phase1_diff/DEPs_Adjacent_vs_Normal.csv", "phase1_diff/DEPs_significant.csv")) {
    fp <- file.path(RES, f)
    if (file.exists(fp)) {
      d <- read.csv(fp, stringsAsFactors = FALSE)
      id_col <- find_col(d, c("Protein", "protein_id", "ENSP"))
      nm_col <- find_col(d, c("gene_name", "symbol", "Gene"))
      if (!is.null(id_col) && !is.null(nm_col)) {
        sub <- data.frame(ensp = d[[id_col]], symbol = d[[nm_col]], stringsAsFactors = FALSE)
        sub$ensp <- gsub("\\.[0-9]+$", "", sub$ensp)
        sub <- sub[!is.na(sub$symbol) & sub$symbol != "" & sub$symbol != "_", ]
        map <- rbind(map, sub)
      }
    }
  }
  tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    if ("ENSEMBLPROT" %in% keytypes(org.Hs.eg.db)) {
      db_map <- AnnotationDbi::select(org.Hs.eg.db,
        keys = keys(org.Hs.eg.db, keytype = "ENSEMBLPROT"),
        columns = c("ENSEMBLPROT", "SYMBOL"), keytype = "ENSEMBLPROT")
      db_map <- db_map[!is.na(db_map$SYMBOL) & !duplicated(db_map$ENSEMBLPROT), ]
      map <- rbind(map, data.frame(ensp = db_map$ENSEMBLPROT, symbol = db_map$SYMBOL, stringsAsFactors = FALSE))
    }
  }, error = function(e) NULL)
  map <- map[!duplicated(map$ensp), ]
  .ensp_map <<- map; map
}

ensp_to_symbol <- function(ids) {
  m <- build_ensp_map()
  ids_clean <- gsub("\\.[0-9]+$", "", ids); ids_clean <- gsub("_[0-9]+$", "", ids_clean)
  idx <- match(ids_clean, m$ensp)
  ifelse(is.na(idx), ids, m$symbol[idx])
}

# =============================================================================
# SECTION 8: Data Loading
# =============================================================================
cat("  Loading data...\n")

tc_vst <- read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"), row.names = 1, check.names = FALSE)
pr_log <- read.csv(file.path(DATA, "proteomics_log2_norm.csv"), row.names = 1, check.names = FALSE)
mb_log <- read.csv(file.path(DATA, "metabolomics_log2_merged.csv"), row.names = 1, check.names = FALSE)
gene_corr <- read.csv(file.path(RES, "phase8_biomarker/mRNA_protein_correlation.csv"), stringsAsFactors = FALSE)

cat("  Data: TC=", ncol(tc_vst), "samples, PR=", ncol(pr_log), ", MB=", ncol(mb_log), "\n")
cat("  Correlations:", nrow(gene_corr), "genes; Sig:", sum(gene_corr$padj < 0.05), "\n")

# Load DEG/DEP/DEM results for volcano plots
degs <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
deps <- read.csv(file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
dems <- read.csv(file.path(RES, "phase1_diff/DEMs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)

# =============================================================================
# SECTION 8b: Volcano Plot Rendering (Panels A-C) — Moved from Figure 1
# =============================================================================
cat("\n--- Rendering Panels A-C: Volcano Plots ---\n")

render_volcano <- function(df, label_col, fc_cut, p_cut, title, filename,
                           n_label = 5, max_label_len = 35, id_mapper = NULL,
                           w_mm = 61, h_mm = 80, use_italic = TRUE) {
  tryCatch({
    vdf <- data.frame(
      logFC = df$logFC, pval = df$P.Value,
      label = as.character(df[[label_col]]), stringsAsFactors = FALSE)
    vdf$neg_log10p <- -log10(vdf$pval + 1e-300)
    vdf$direction <- ifelse(abs(vdf$logFC) < fc_cut | vdf$pval >= p_cut, "NS",
                            ifelse(vdf$logFC > 0, "Up", "Down"))
    if (!is.null(id_mapper)) vdf$label <- id_mapper(vdf$label)
    vdf$label_clean <- ifelse(grepl("^ENS[GPT]", vdf$label), NA, vdf$label)
    vdf$label_clean <- ifelse(!is.na(vdf$label_clean) & nchar(vdf$label_clean) > max_label_len,
                              paste0(substr(vdf$label_clean, 1, max_label_len - 3), "..."),
                              vdf$label_clean)
    sig_up <- vdf[vdf$direction == "Up" & !is.na(vdf$label_clean), ]
    sig_down <- vdf[vdf$direction == "Down" & !is.na(vdf$label_clean), ]
    top_up <- head(sig_up[order(sig_up$pval), ], n_label)
    top_down <- head(sig_down[order(sig_down$pval), ], n_label)
    label_idx <- which(vdf$label_clean %in% c(top_up$label_clean, top_down$label_clean) &
                       vdf$direction != "NS")
    vdf$show_label <- ""; vdf$show_label[label_idx] <- vdf$label_clean[label_idx]
    vdf$is_labeled <- vdf$show_label != ""
    n_up <- sum(vdf$direction == "Up"); n_down <- sum(vdf$direction == "Down")
    label_face <- if (use_italic) "italic" else "plain"

    p <- ggplot(vdf, aes(x = logFC, y = neg_log10p)) +
      geom_point(data = vdf[!vdf$is_labeled, ], aes(color = direction), size = 0.4, alpha = 0.35) +
      geom_point(data = vdf[vdf$is_labeled & vdf$direction == "Up", ],
                 color = "#8B1A1AFF", size = 1.2, alpha = 1) +
      geom_point(data = vdf[vdf$is_labeled & vdf$direction == "Down", ],
                 color = "#003F72FF", size = 1.2, alpha = 1) +
      scale_color_manual(values = c(Up = COL_UP, Down = COL_DOWN, NS = COL_NS)) +
      geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", linewidth = 0.3, color = "grey40") +
      geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.3, color = "grey40") +
      ggrepel::geom_label_repel(
        data = vdf[vdf$show_label != "", ],
        aes(label = show_label), size = 2.2, family = FONT_FAMILY, fontface = label_face,
        max.overlaps = 15, force = 4, segment.size = 0.2, segment.color = "grey40",
        label.size = 0.12, label.padding = unit(0.1, "lines"),
        fill = alpha("white", 0.85), seed = 42) +
      labs(title = title, subtitle = sprintf("Up: %d  |  Down: %d", n_up, n_down),
           x = expression(log[2]~"fold change"), y = expression(-log[10]*(italic(P)))) +
      theme_nc + theme(legend.position = "none",
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5),
            plot.margin = margin(3, 3, 3, 3, "mm"))

    save_panel_pdf(filename, w_mm, h_mm, print(p))
    return(p)
  }, error = function(e) { cat("  ERROR", filename, ":", e$message, "\n"); return(NULL) })
}

suppressPackageStartupMessages(library(ggrepel))

p_va <- render_volcano(degs, "gene_name", fc_cut = log2(1.5), p_cut = 0.05,
                       title = "Transcriptomics", filename = "Supp02a_volcano_TC.pdf",
                       n_label = 5, id_mapper = ensg_to_symbol, w_mm = 61, h_mm = 78)

p_vb <- render_volcano(deps, "gene_name", fc_cut = log2(1.5), p_cut = 0.05,
                       title = "Proteomics", filename = "Supp02b_volcano_PR.pdf",
                       n_label = 5, id_mapper = ensp_to_symbol, w_mm = 61, h_mm = 78,
                       use_italic = FALSE)

p_vc <- render_volcano(dems, "metabolite_name", fc_cut = log2(1.5), p_cut = 0.05,
                       title = "Metabolomics", filename = "Supp02c_volcano_MT.pdf",
                       n_label = 5, max_label_len = 20, w_mm = 61, h_mm = 78,
                       use_italic = FALSE)


# =============================================================================
# SECTION 9: Heatmap Builder (Panels A-C)
# =============================================================================

make_top_heatmap <- function(expr_mat, sig_csv, n_top = 30, title, type,
                             body_w_mm = 38, row_max_mm = 42,
                             show_legend = TRUE) {
  sig <- read.csv(sig_csv, stringsAsFactors = FALSE)
  fc_col <- find_col(sig, c("logFC", "log2FoldChange", "fc"))
  id_candidates <- c("gene_id", "Protein", "Compound_ID", "name", "gene_name",
                     "gene", "X", "metabolite", "metabolite_name", "ID")
  best_n <- 0; id_col <- NULL
  for (cand in id_candidates) {
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

  if (type == "TC") {
    rownames(mat) <- ensg_to_symbol(rownames(mat))
    mat <- mat[!grepl("^ENSG", rownames(mat)), , drop = FALSE]
    row_gp <- gpar(fontsize = 8, fontface = "italic", fontfamily = FONT_GRID)
  } else if (type == "PR") {
    rownames(mat) <- ensp_to_symbol(rownames(mat))
    mat <- mat[!grepl("^ENSP", rownames(mat)), , drop = FALSE]
    row_gp <- gpar(fontsize = 8, fontface = "italic", fontfamily = FONT_GRID)
  } else {
    # MT: shorten long metabolite names per AGENTS.md rules
    rn <- rownames(mat)
    rn <- gsub("Phosphatidylcholine", "PC", rn, ignore.case = TRUE)
    rn <- gsub("Phosphatidylethanolamine", "PE", rn, ignore.case = TRUE)
    rn <- gsub("Lysophosphatidylcholine", "LPC", rn, ignore.case = TRUE)
    rn <- gsub("Lysophosphatidylethanolamine", "LPE", rn, ignore.case = TRUE)
    rn <- gsub("Sphingomyelin", "SM", rn, ignore.case = TRUE)
    rn <- gsub("Hydroxy", "OH-", rn, ignore.case = TRUE)
    rn <- gsub("alpha-", "\u03b1-", rn, ignore.case = TRUE)
    rn <- gsub("beta-",  "\u03b2-", rn, ignore.case = TRUE)
    rn <- gsub("gamma-", "\u03b3-", rn, ignore.case = TRUE)
    rn <- ifelse(nchar(rn) > 28, paste0(substr(rn, 1, 25), "..."), rn)
    rownames(mat) <- make.unique(rn)
    row_gp <- gpar(fontsize = 7, fontfamily = FONT_GRID)
  }

  if (nrow(mat) < 3) return(NULL)

  grp <- factor(make_group_vector(colnames(mat)), levels = c("Normal", "Adjacent"))
  legend_param <- list(title_gp = gp_legend_title(), labels_gp = gp_legend_labels(),
                       legend_height = unit(18, "mm"), grid_width = unit(2.5, "mm"))

  # Top annotation: Group color strip for clear visual group distinction
  ha_top <- HeatmapAnnotation(
    Group = grp,
    col = list(Group = c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)),
    annotation_name_gp = gp_anno_name(),
    simple_anno_size = unit(3, "mm"),
    show_legend = show_legend,
    annotation_legend_param = list(title_gp = gp_legend_title(), labels_gp = gp_legend_labels())
  )

  Heatmap(mat, name = "Z-score", col = col_div,
    show_row_names = TRUE, row_names_gp = row_gp,
    show_column_names = FALSE,
    top_annotation = ha_top,
    column_split = grp, cluster_column_slices = FALSE, cluster_columns = TRUE,
    column_title = c("Normal", "Adjacent"),
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_FAMILY),
    column_title_side = "top",
    column_gap = unit(2, "mm"),
    show_row_dend = FALSE, show_column_dend = FALSE,
    border = TRUE, rect_gp = gpar(col = NA),
    width = unit(body_w_mm, "mm"),
    show_heatmap_legend = show_legend,
    heatmap_legend_param = legend_param,
    row_names_max_width = unit(row_max_mm, "mm"))
}

# =============================================================================
# SECTION 10: Render Panels F-H (Top DE Features)
# =============================================================================
cat("\n--- Rendering panels F-H (Top DE Features) ---\n")

# Per-panel dimensions (3 panels side-by-side: 55+55+73=183mm):
HM_H <- 95
A_W  <- 55;  A_BODY <- 24;  A_RNW <- 22   # DEGs bar chart
B_W  <- 55;  B_BODY <- 24;  B_RNW <- 22   # DEPs bar chart
C_W  <- 73;  C_BODY <- 24;  C_RNW <- 30   # DEMs (heatmap with shared legend)

ht_A <- ht_B <- ht_C <- NULL

# --- Helper: Horizontal logFC bar chart for top DE features ---
make_logfc_bar <- function(sig_csv, type, n_top = 20, title_text) {
  sig <- read.csv(sig_csv, stringsAsFactors = FALSE)
  fc_col <- find_col(sig, c("logFC", "log2FoldChange", "fc"))
  sig <- sig[!is.na(sig[[fc_col]]) & sig$P.Value < 0.05, ]
  sig <- sig[order(abs(sig[[fc_col]]), decreasing = TRUE), ]

  if (type == "TC") {
    # For transcriptomics: use gene_id -> ensg_to_symbol for robust resolution
    id_col <- find_col(sig, c("gene_id", "Protein"))
    nm_col <- find_col(sig, c("gene_name", "symbol"))
    # First try gene_name, then resolve via ensg_to_symbol
    labels <- sig[[nm_col]]
    if (!is.null(id_col)) {
      mapped <- ensg_to_symbol(sig[[id_col]])
      # Replace ENS or empty labels with successfully mapped symbols
      for (i in seq_along(labels)) {
        if ((grepl("^ENS", labels[i]) || is.na(labels[i]) || labels[i] == "") &&
            !grepl("^ENS", mapped[i]) && !is.na(mapped[i]) && mapped[i] != "") {
          labels[i] <- mapped[i]
        }
      }
    }
    # Filter out any remaining ENSG IDs
    keep <- !grepl("^ENS", labels) & !is.na(labels) & labels != ""
    sig <- sig[keep, ]
    labels <- labels[keep]
    sig <- head(sig, n_top)
    labels <- head(labels, n_top)
    face_style <- "italic"
  } else if (type == "PR") {
    nm_col <- find_col(sig, c("gene_name", "symbol"))
    sig <- head(sig, n_top * 2)
    labels <- ensp_to_symbol(sig[[nm_col]])
    labels[grepl("^ENS", labels)] <- sig[[nm_col]][grepl("^ENS", labels)]
    keep <- !grepl("^ENS", labels) & !is.na(labels) & labels != ""
    sig <- sig[keep, ]
    labels <- labels[keep]
    sig <- head(sig, n_top)
    labels <- head(labels, n_top)
    face_style <- "italic"
  } else {
    # Metabolomics: use metabolite_name, abbreviate long names
    nm_col <- find_col(sig, c("metabolite_name", "name", "Compound_ID"))
    sig <- head(sig, n_top)
    labels <- sig[[nm_col]]
    # Abbreviate long metabolite names (per AGENTS.md rules)
    labels <- gsub("Phosphatidylcholine", "PC", labels, ignore.case = TRUE)
    labels <- gsub("Phosphatidylethanolamine", "PE", labels, ignore.case = TRUE)
    labels <- gsub("Lysophosphatidylcholine", "LPC", labels, ignore.case = TRUE)
    labels <- gsub("Lysophosphatidylethanolamine", "LPE", labels, ignore.case = TRUE)
    labels <- gsub("Sphingomyelin", "SM", labels, ignore.case = TRUE)
    labels <- gsub("Hydroxy", "OH-", labels, ignore.case = TRUE)
    labels <- gsub("alpha-", "\u03b1-", labels, ignore.case = TRUE)
    labels <- gsub("beta-",  "\u03b2-", labels, ignore.case = TRUE)
    labels <- gsub("gamma-", "\u03b3-", labels, ignore.case = TRUE)
    labels <- ifelse(nchar(labels) > 25, paste0(substr(labels, 1, 22), "..."), labels)
    face_style <- "plain"
  }

  labels[is.na(labels) | labels == ""] <- paste0("Feature_", which(is.na(labels) | labels == ""))

  df_bar <- data.frame(
    label = factor(make.unique(labels), levels = rev(make.unique(labels))),
    logFC = sig[[fc_col]],
    direction = ifelse(sig[[fc_col]] > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )

  p <- ggplot(df_bar, aes(x = logFC, y = label, fill = direction)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
    scale_fill_manual(values = c("Up" = COL_UP, "Down" = COL_DOWN)) +
    labs(title = title_text, x = expression(log[2]*FC), y = NULL) +
    theme_bw(base_size = 8, base_family = FONT_FAMILY) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.15, color = "grey90"),
      panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
      axis.text.y = element_text(size = 7, face = face_style, family = FONT_FAMILY),
      axis.text.x = element_text(size = 8, family = FONT_FAMILY),
      axis.title.x = element_text(size = 9, face = "bold", family = FONT_FAMILY),
      plot.title = element_text(size = 10, face = "bold", family = FONT_FAMILY, hjust = 0.5),
      plot.margin = margin(3, 4, 3, 3, "mm")
    )
  return(p)
}

# Panel F: Top 20 DEGs — horizontal bar chart
tryCatch({
  p_bar_f <- make_logfc_bar(
    file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"),
    type = "TC", n_top = 20, title_text = "Top 20 DEGs"
  )
  save_panel_pdf("Supp02f_top_DEGs.pdf", A_W, HM_H, print(p_bar_f))
}, error = function(e) cat("  ERROR Panel F:", e$message, "\n"))

# Panel G: Top 20 DEPs — horizontal bar chart
tryCatch({
  p_bar_g <- make_logfc_bar(
    file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"),
    type = "PR", n_top = 20, title_text = "Top 20 DEPs"
  )
  save_panel_pdf("Supp02g_top_DEPs.pdf", B_W, HM_H, print(p_bar_g))
}, error = function(e) cat("  ERROR Panel G:", e$message, "\n"))

# Panel H: Top 20 DEMs — horizontal bar chart (consistent with F-G)
tryCatch({
  p_bar_h <- make_logfc_bar(
    file.path(RES, "phase1_diff/DEMs_Adjacent_vs_Normal.csv"),
    type = "MT", n_top = 20, title_text = "Top 20 DEMs"
  )
  save_panel_pdf("Supp02h_top_DEMs.pdf", C_W, HM_H, print(p_bar_h))
}, error = function(e) cat("  ERROR Panel H:", e$message, "\n"))

# =============================================================================
# SECTION 11: Render Panels D-E (mRNA-Protein Discordance)
# =============================================================================
cat("\n--- Rendering panels D-E (mRNA-Protein Discordance) ---\n")

PD_W <- 88; PD_H <- 95    # density (Row 2 left: 88+95=183)
PE_W <- 95; PE_H <- 95    # scatter 3x2 grid (Row 2 right)

p_D <- NULL
p_E <- NULL

# Panel D: mRNA-protein correlation density
tryCatch({
  med_rho <- median(gene_corr$rho, na.rm = TRUE)
  n_sig   <- sum(gene_corr$padj < 0.05, na.rm = TRUE)
  p_D <<- ggplot(gene_corr, aes(x = rho)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = COL_TC, alpha = 0.7, color = "white", linewidth = 0.1) +
    geom_density(color = COL_UP, linewidth = 0.8) +
    geom_vline(xintercept = med_rho, linetype = "dashed", linewidth = 0.4, colour = "grey20") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey50") +
    annotate("text", x = med_rho + 0.03, y = Inf,
             label = sprintf("Median = %.3f", med_rho),
             vjust = 2, hjust = 0, size = FS_GEOM_TEXT, family = FONT_FAMILY) +
    labs(title = "mRNA-Protein Correlation Distribution",
         subtitle = sprintf("n = %d genes; %d significant (padj < 0.05)", nrow(gene_corr), n_sig),
         x = expression(Spearman~italic(rho)), y = "Density") +
    theme_nc + theme(plot.title = element_text(hjust = 0.5))
  save_panel_pdf("Supp02d_mRNA_prot_cor_density.pdf", PD_W, PD_H, print(p_D))
}, error = function(e) cat("  ERROR Panel D:", e$message, "\n"))

# Panel E: Top correlated gene scatter plots (3 pos + 3 neg = 6 in 3x2 grid)
tryCatch({
  shared_samples <- intersect(colnames(tc_vst), colnames(pr_log))
  degs <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
  deps <- read.csv(file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)

  ensg_map_local <- setNames(degs$gene_id, degs$gene_name)
  ensg_map_local <- ensg_map_local[!is.na(names(ensg_map_local)) & names(ensg_map_local) != ""]
  ensp_map_local <- setNames(deps$Protein, deps$gene_name)
  ensp_map_local <- ensp_map_local[!is.na(names(ensp_map_local)) & names(ensp_map_local) != ""]

  gene_corr$ensg <- ensg_map_local[gene_corr$gene]
  gene_corr$ensp <- ensp_map_local[gene_corr$gene]
  gene_corr_avail <- gene_corr[!is.na(gene_corr$ensg) & !is.na(gene_corr$ensp) &
                                gene_corr$ensg %in% rownames(tc_vst) &
                                gene_corr$ensp %in% rownames(pr_log), ]

  pos_sorted <- gene_corr_avail[order(-gene_corr_avail$rho), ]
  neg_sorted <- gene_corr_avail[order( gene_corr_avail$rho), ]
  top_pos <- head(pos_sorted[pos_sorted$rho > 0, ], 3)
  top_neg <- head(neg_sorted[neg_sorted$rho < 0, ], 3)
  top_ex  <- rbind(top_pos, top_neg)

  if (nrow(top_ex) >= 4) {
    grp <- ifelse(grepl("Normal|^N\\d", shared_samples), "Normal", "Adjacent")
    plist <- list()
    for (i in seq_len(nrow(top_ex))) {
      g_symbol <- top_ex$gene[i]
      g_ensg   <- top_ex$ensg[i]
      g_ensp   <- top_ex$ensp[i]
      rho_val  <- top_ex$rho[i]
      df_s <- data.frame(mRNA    = as.numeric(tc_vst[g_ensg, shared_samples]),
                         Protein = as.numeric(pr_log[g_ensp, shared_samples]),
                         Group   = grp)
      sub_title <- sprintf("%s  (\u03c1 = %.2f)", g_symbol, rho_val)
      plist[[i]] <- ggplot(df_s, aes(mRNA, Protein, color = Group)) +
        geom_smooth(method = "lm", se = TRUE, color = "grey30",
                    fill = "grey80", linewidth = 0.4, alpha = 0.4) +
        geom_point(size = 1.4, alpha = 0.85, stroke = 0) +
        scale_color_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT)) +
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 4) +
        labs(title = sub_title, x = "mRNA (VST)", y = "Protein (log2)") +
        theme_nc + theme(legend.position = "none",
                         plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
                         axis.text  = element_text(size = 6.5),
                         axis.title = element_text(size = 7.5),
                         plot.margin = margin(1, 0.5, 1, 0.5, "mm"))
    }
    p_E <<- ggpubr::ggarrange(plotlist = plist, ncol = 3, nrow = 2,
                              common.legend = TRUE, legend = "bottom")
    save_panel_pdf("Supp02e_mRNA_prot_scatter.pdf", PE_W, PE_H, print(p_E))
  } else {
    cat("  SKIP Panel E: insufficient matched genes\n")
  }
}, error = function(e) cat("  ERROR Panel E:", e$message, "\n"))

# =============================================================================
# SECTION 12: Full-Vector Composite Assembly (grid viewports — zero-gap)
# =============================================================================
cat("\n--- Assembling composite SuppFig_02 (grid viewports, 183x268mm, 3 rows) ---\n")

# Layout (W=183, H=268), 3 rows — all 8 panels, reading order a→h:
#   Row 1 (H=78):  a(61x78) + b(61x78) + c(61x78)     -- 3 Volcano plots
#   Row 2 (H=95):  d(88x95) + e(95x95)                 -- Density + Scatter 3x2
#   Row 3 (H=95):  f(55x95) + g(55x95) + h(73x95)     -- 3 Bar charts
W_TOTAL <- 183; H_TOTAL <- 268
FS_TAG  <- 12

tryCatch({
  # Ensure all panels exist
  panels_ok <- !is.null(p_va) && !is.null(p_vb) && !is.null(p_vc) &&
               !is.null(p_D) && !is.null(p_E) &&
               !is.null(p_bar_f) && !is.null(p_bar_g) && !is.null(p_bar_h)
  if (!panels_ok) stop("One or more panel ggplot objects are NULL")

  # Page margins (mm) — keep content within printable area
  M <- 3  # uniform 3mm margin on all sides
  AW <- W_TOTAL - 2 * M   # available width  = 177mm
  AH <- H_TOTAL - 2 * M   # available height = 262mm

  # Row heights (proportional to original 78:95:95, scaled to AH=262)
  rh1 <- round(AH * 78 / 268)   # 76mm
  rh2 <- round(AH * 95 / 268)   # 93mm
  rh3 <- AH - rh1 - rh2         # 93mm

  # Panel positioning: exact mm coordinates (x, y from top-left), with margins
  panels_layout <- list(
    list(x = M,            y = M,          w = round(AW/3),             h = rh1, p = p_va,    tag = "a"),
    list(x = M+round(AW/3),   y = M,      w = round(AW/3),             h = rh1, p = p_vb,    tag = "b"),
    list(x = M+2*round(AW/3), y = M,      w = AW-2*round(AW/3),       h = rh1, p = p_vc,    tag = "c"),
    list(x = M,            y = M+rh1,      w = round(AW*88/183),       h = rh2, p = p_D,     tag = "d"),
    list(x = M+round(AW*88/183), y = M+rh1, w = AW-round(AW*88/183), h = rh2, p = p_E,     tag = "e"),
    list(x = M,            y = M+rh1+rh2,  w = round(AW*55/183),       h = rh3, p = p_bar_f, tag = "f"),
    list(x = M+round(AW*55/183),   y = M+rh1+rh2, w = round(AW*55/183),       h = rh3, p = p_bar_g, tag = "g"),
    list(x = M+2*round(AW*55/183), y = M+rh1+rh2, w = AW-2*round(AW*55/183), h = rh3, p = p_bar_h, tag = "h")
  )

  # Render function: place each panel into exact mm viewport
  render_composite <- function() {
    grid::grid.newpage()
    for (pan in panels_layout) {
      vp <- grid::viewport(
        x      = unit(pan$x, "mm"),
        y      = unit(H_TOTAL - pan$y, "mm"),
        width  = unit(pan$w, "mm"),
        height = unit(pan$h, "mm"),
        just   = c("left", "top")
      )
      grid::pushViewport(vp)
      # Print ggplot with zero margin into this viewport
      p_stripped <- pan$p + theme(plot.margin = margin(0, 0, 0, 0))
      print(p_stripped, newpage = FALSE)
      # Draw panel tag (top-left corner)
      grid::grid.text(
        pan$tag,
        x    = unit(1.5, "mm"),
        y    = unit(pan$h - 1.5, "mm"),
        just = c("left", "top"),
        gp   = grid::gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY)
      )
      grid::popViewport()
    }
  }

  # --- PDF (full vector) ---
  cairo_pdf(file.path(OUT, "SuppFig_02.pdf"),
            width = W_TOTAL / 25.4, height = H_TOTAL / 25.4, family = FONT_FAMILY)
  render_composite()
  dev.off()
  cat("  -> SuppFig_02.pdf  (full-vector grid viewports)\n")

  # --- PNG (600 DPI) ---
  px_per_mm <- ASSEMBLY_DPI / 25.4
  px_W <- round(W_TOTAL * px_per_mm)
  px_H <- round(H_TOTAL * px_per_mm)
  grDevices::png(file.path(OUT, "SuppFig_02.png"),
                 width = px_W, height = px_H, res = ASSEMBLY_DPI, type = "cairo")
  render_composite()
  dev.off()
  cat(sprintf("  -> SuppFig_02.png  (%dx%dpx)\n", px_W, px_H))

  # --- TIFF (600 DPI, LZW) ---
  grDevices::tiff(file.path(OUT, "SuppFig_02.tiff"),
                  width = px_W, height = px_H, res = ASSEMBLY_DPI,
                  compression = "lzw", type = "cairo")
  render_composite()
  dev.off()
  cat("  -> SuppFig_02.tiff (LZW)\n")

  cat("  SuppFig_02 DONE — grid viewport zero-gap layout, all 8 panels\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 2 rendering complete ===\n")
