#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_01_standalone.R
# Complete, Self-Contained Code for Supplementary Figure 1
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: EBioMedicine (Lancet family) - vector PDF + PNG/TIFF (600 DPI) composite
# =============================================================================
# Supplementary Figure 1: QC - PCA & Sample Correlation
# 7 panels:
#   a-c) Per-omics paired PCA (Transcriptomics / Proteomics / Metabolomics)
#        with 90% confidence ellipse, paired-sample lines, PERMANOVA R^2/p,
#        and MAD-based PC outlier labelling
#   d)   Combined multi-omics PCA on shared samples (top-1000 TC + top-1000 PR
#        + top-300 MT, Z-scaled per layer, concatenated)
#   e-g) Sample-to-sample Pearson correlation heatmaps (TC / PR / MT)
#        split by Group, auto-floor color scale (5th-percentile of off-diagonal)
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript SuppFig_01_standalone.R
# =============================================================================

cat("=== Supplementary Figure 1: QC - PCA & Sample Correlation ===\n")
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
})

# --- Register Arial in R's PostScript/PDF font databases ---
# cairo_pdf natively embeds ArialMT TrueType from the system, but R's grid
# package uses the PostScript font database for text-width calculation. Without
# this registration, grid emits harmless but noisy "font family 'Arial' not
# found in PostScript font database" warnings.  Mapping Arial -> Helvetica
# metrics is safe because the two fonts are metrically identical.
pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
cat("  Arial registered in PDF/PostScript font databases\n")


# =============================================================================
# SECTION 2: Project Paths (EDIT THESE IF RUNNING ON A DIFFERENT MACHINE)
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_01")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants (Nature Communications Standard)
# =============================================================================
FONT_FAMILY <- "Arial"      # AI native font - better than Arial for AI editing
FONT_GRID   <- "Arial"      # unified for ComplexHeatmap
MM_PER_INCH <- 25.4

# NC page dimensions (mm)
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
FS_TAG        <- 12    # panel tag (A, B, C...)
FS_GEOM_TEXT  <- 2.82  # geom_text size (= 8pt in mm)

# =============================================================================
# SECTION 4: Color Palette (ggsci JCO / Journal of Clinical Oncology)
# =============================================================================
PAL_CAT <- pal_jco("default")(10)

# Semantic colors
COL_UP       <- "#CD534CFF"    # up-regulated (NPG red)
COL_DOWN     <- "#0073C2FF"    # down-regulated (NPG cyan)
COL_NS       <- "#868686FF"    # not significant (SCI standard grey)
COL_NA       <- "#F0F0F0"    # NA / background

# Omics layer colors
COL_TC       <- "#0073C2FF"    # transcriptomics (navy)
COL_PR       <- "#CD534CFF"    # proteomics (red)
COL_MT       <- "#EFC000FF"    # metabolomics (teal)

# Group colors
COL_NORMAL   <- "#7AA6DCFF"    # normal tissue (cyan)
COL_ADJACENT <- "#CD534CFF"    # adjacent/disease tissue (red)

# Subtype colors
COL_CS1      <- "#CD534CFF"    # consensus subtype 1 (red)
COL_CS2      <- "#0073C2FF"    # consensus subtype 2 (navy)

# =============================================================================
# SECTION 5: Color Scale Functions (circlize::colorRamp2)
# =============================================================================
col_div     <- colorRamp2(c(-2, 0, 2),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_nes     <- colorRamp2(c(-3, 0, 3),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_cor     <- colorRamp2(c(-1, 0, 1),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_zscore  <- colorRamp2(c(-2.5, 0, 2.5), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_immune  <- colorRamp2(c(-2, 0, 2),     c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 6: ggplot2 Theme (theme_bw base, Nature Communications style)
# =============================================================================
# Expert review requirement: ALL text elements must explicitly set family = "Arial"
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
    plot.title         = element_text(family = FONT_FAMILY, size = 10, face = "bold", hjust = 0),
    plot.title.position = "plot",
    legend.text        = element_text(family = FONT_FAMILY, size = 8),
    legend.title       = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    legend.key.size    = unit(3.5, "mm"),
    legend.background  = element_blank(),
    legend.box.background = element_blank(),
    legend.spacing.y   = unit(1, "mm"),
    strip.text         = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    strip.background   = element_rect(fill = "grey96", linewidth = 0.5, color = "grey80"),
    plot.margin        = margin(5, 5, 5, 5, "mm")
  )
theme_set(theme_nc)
# Alias for compatibility with nc_theme.R enhanced vis functions
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
  gpar(col = "grey90", lwd = 0.3)
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

# --- Save helper ---
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
  tryCatch(force(expr), error = function(e) message("  ERROR saving ", filename, ": ", e$message))
  dev.off()
  cat("  ->", basename(s$path), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
  invisible(s$path)
}

save_tiff <- function(filename, w_mm, h_mm, expr, out_dir = NULL) {
  fpath <- file.path(if (is.null(out_dir)) OUT else out_dir, sub("\\.pdf$", ".tiff", filename))
  invisible(fpath)
}

# --- CSV reader ---
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

# --- Matrix scaling ---
clamp_matrix <- function(mat, lim = 2) {
  mat[mat > lim] <- lim; mat[mat < -lim] <- -lim; mat
}

scale_rows <- function(mat, lim = 2) {
  mat <- t(scale(t(as.matrix(mat)))); mat[is.na(mat)] <- 0; clamp_matrix(mat, lim)
}

# --- Group helpers ---
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

# --- Clinical annotation bar for Figure 1 ---
get_nc_annotation <- function(subtypes, clinical_df, patient_ids = NULL) {
  if (!is.null(patient_ids)) {
    clin_sub <- clinical_df[match(patient_ids, clinical_df$patient_id), ]
  } else {
    clin_sub <- clinical_df
  }
  sex_vec <- ifelse(clin_sub$sex == "\u7537", "Male", "Female")
  bilirubin_vec <- ifelse(!is.na(clin_sub$total_bilirubin) & clin_sub$total_bilirubin > 20, "High", "Normal")
  bile_duct_vec <- ifelse(clin_sub$bile_duct_invasion == "1", "Invaded", "No")
  HeatmapAnnotation(
    Subtype = subtypes,
    Sex = sex_vec,
    Bilirubin = bilirubin_vec,
    `Bile Duct` = bile_duct_vec,
    col = list(
      Subtype = c("CS1" = COL_CS1, "CS2" = COL_CS2, "DS1" = COL_CS1, "DS2" = COL_CS2),
      Sex = c("Male" = COL_TC, "Female" = PAL_CAT[5]),
      Bilirubin = c("High" = COL_MT, "Normal" = COL_NA),
      `Bile Duct` = c("Invaded" = PAL_CAT[8], "No" = COL_NA)
    ),
    simple_anno_size = unit(4, "mm"),
    gap = unit(1, "mm"),
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_GRID),
    annotation_legend_param = list(
      title_gp  = gpar(fontsize = 8, fontfamily = FONT_GRID),
      labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID)
    )
  )
}

# --- Column finder ---
find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else NULL
}

# --- ID mapping: ENSG -> gene symbol ---
.ensg_map <- NULL
build_ensg_map <- function() {
  if (!is.null(.ensg_map)) return(.ensg_map)
  map <- data.frame(ensg = character(), symbol = character(), stringsAsFactors = FALSE)
  tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    db_map <- AnnotationDbi::select(org.Hs.eg.db,
      keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
      columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
    db_map <- db_map[!is.na(db_map$SYMBOL), ]
    db_map <- db_map[!duplicated(db_map$ENSEMBL), ]
    map <- rbind(map, data.frame(ensg = db_map$ENSEMBL, symbol = db_map$SYMBOL,
                                  stringsAsFactors = FALSE))
    cat("  ensg_map: loaded", nrow(db_map), "from org.Hs.eg.db\n")
  }, error = function(e) cat("  ensg_map: org.Hs.eg.db unavailable\n"))
  for (f in c("phase1_diff/DEGs_Adjacent_vs_Normal.csv",
              "phase1_diff/DEPs_Adjacent_vs_Normal.csv")) {
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
  .ensg_map <<- map
  cat("  ensg_map: total", nrow(map), "mappings\n")
  map
}

ensg_to_symbol <- function(ids) {
  m <- build_ensg_map()
  ids_clean <- gsub("\\.[0-9]+$", "", ids)
  ids_clean <- gsub("_[0-9]+$", "", ids_clean)
  idx <- match(ids_clean, m$ensg)
  out <- ifelse(is.na(idx), ids, m$symbol[idx])
  out
}

# --- ID mapping: ENSP -> gene symbol ---
.ensp_map <- NULL
build_ensp_map <- function() {
  if (!is.null(.ensp_map)) return(.ensp_map)
  map <- data.frame(ensp = character(), symbol = character(), stringsAsFactors = FALSE)
  for (f in c("phase1_diff/DEPs_Adjacent_vs_Normal.csv",
              "phase1_diff/DEPs_significant.csv")) {
    fp <- file.path(RES, f)
    if (file.exists(fp)) {
      d <- read.csv(fp, stringsAsFactors = FALSE)
      id_col <- find_col(d, c("Protein", "protein_id", "ENSP"))
      nm_col <- find_col(d, c("gene_name", "symbol", "Gene"))
      if (!is.null(id_col) && !is.null(nm_col)) {
        sub <- data.frame(ensp = d[[id_col]], symbol = d[[nm_col]], stringsAsFactors = FALSE)
        sub$ensp_clean <- gsub("\\.[0-9]+$", "", sub$ensp)
        sub <- sub[!is.na(sub$symbol) & sub$symbol != "" & sub$symbol != "_", ]
        map <- rbind(map, data.frame(ensp = sub$ensp_clean, symbol = sub$symbol, stringsAsFactors = FALSE))
      }
    }
  }
  tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    if ("ENSEMBLPROT" %in% keytypes(org.Hs.eg.db)) {
      db_map <- AnnotationDbi::select(org.Hs.eg.db,
        keys = keys(org.Hs.eg.db, keytype = "ENSEMBLPROT"),
        columns = c("ENSEMBLPROT", "SYMBOL"), keytype = "ENSEMBLPROT")
      db_map <- db_map[!is.na(db_map$SYMBOL), ]
      db_map <- db_map[!duplicated(db_map$ENSEMBLPROT), ]
      map <- rbind(map, data.frame(ensp = db_map$ENSEMBLPROT, symbol = db_map$SYMBOL,
                                    stringsAsFactors = FALSE))
      cat("  ensp_map: loaded", nrow(db_map), "from org.Hs.eg.db\n")
    }
  }, error = function(e) cat("  ensp_map: org.Hs.eg.db ENSEMBLPROT unavailable\n"))
  map <- map[!duplicated(map$ensp), ]
  .ensp_map <<- map
  cat("  ensp_map: total", nrow(map), "mappings\n")
  map
}

ensp_to_symbol <- function(ids) {
  m <- build_ensp_map()
  ids_clean <- gsub("\\.[0-9]+$", "", ids)
  ids_clean <- gsub("_[0-9]+$", "", ids_clean)
  idx <- match(ids_clean, m$ensp)
  out <- ifelse(is.na(idx), ids, m$symbol[idx])
  out
}

# --- Feature name formatting ---
format_feature_names <- function(ids, mapping_df, type = "MT", max_len = 30) {
  new_names <- mapping_df$name[match(ids, mapping_df$id)]
  new_names <- ifelse(is.na(new_names) | new_names == "", ids, new_names)
  if (type == "MT") {
    new_names <- ifelse(nchar(new_names) > max_len,
                        paste0(substr(new_names, 1, max_len - 3), "..."),
                        new_names)
  } else {
    unmapped_idx <- grepl("^ENS[GP]", new_names)
    new_names[unmapped_idx] <- NA
  }
  return(new_names)
}

# --- Clean cell type names ---
clean_celltype <- function(x) {
  x <- gsub("_", " ", x)
  x
}

# --- Top-variance feature selection (replaces noisy full-feature PCA) ---
select_top_var <- function(mat, n = 2000) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  v[is.na(v)] <- 0
  if (length(v) <= n) return(mat[v > 0, , drop = FALSE])
  keep <- head(order(v, decreasing = TRUE), n)
  mat[keep, , drop = FALSE]
}

# --- MAD outlier detection on PCA scores ---
detect_pca_outliers <- function(scores, threshold = 3) {
  m1 <- median(scores$PC1); s1 <- mad(scores$PC1, constant = 1.4826)
  m2 <- median(scores$PC2); s2 <- mad(scores$PC2, constant = 1.4826)
  if (s1 == 0) s1 <- sd(scores$PC1)
  if (s2 == 0) s2 <- sd(scores$PC2)
  z1 <- abs(scores$PC1 - m1) / s1
  z2 <- abs(scores$PC2 - m2) / s2
  scores$is_outlier <- (z1 > threshold) | (z2 > threshold)
  scores
}

# --- Self-contained PERMANOVA (Anderson 2001), Euclidean distance, two-group ---
permanova_two_group <- function(mat, grp, n_perm = 999, seed = 42) {
  set.seed(seed)
  X <- t(mat)
  D <- as.matrix(dist(X))
  N <- nrow(D)
  grp <- as.factor(grp)
  if (nlevels(grp) != 2) return(list(R2 = NA, p = NA, F = NA, n_perm = 0))
  ss_total <- sum(D[lower.tri(D)]^2) / N
  compute_ssw <- function(g) {
    sum(sapply(levels(g), function(lv) {
      idx <- which(g == lv)
      if (length(idx) < 2) return(0)
      sub <- D[idx, idx]
      sum(sub[lower.tri(sub)]^2) / length(idx)
    }))
  }
  ss_within_obs <- compute_ssw(grp)
  ss_between_obs <- ss_total - ss_within_obs
  F_obs <- (ss_between_obs / 1) / (ss_within_obs / (N - 2))
  R2 <- ss_between_obs / ss_total
  null_F <- replicate(n_perm, {
    g_perm <- sample(grp)
    ssw <- compute_ssw(g_perm)
    ssb <- ss_total - ssw
    (ssb / 1) / (ssw / (N - 2))
  })
  p <- (sum(null_F >= F_obs) + 1) / (n_perm + 1)
  list(R2 = R2, p = p, F = F_obs, n_perm = n_perm)
}

fmt_permanova <- function(res) {
  if (is.na(res$R2)) return("")
  p_str <- if (res$p < 0.001) "p<0.001"
           else if (res$p < 0.01) sprintf("p=%.3f", res$p)
           else sprintf("p=%.2f", res$p)
  sprintf("PERMANOVA: R\u00b2=%.2f, %s", res$R2, p_str)
}


# --- Factory Functions (from production pipeline) ---

# --- 2.2 PCA Scatter Factory (with paired lines, PERMANOVA, MAD outliers) ---
# Enhancements vs original:
#   1. Top-2000 variance feature selection (de-noise)
#   2. scale.=FALSE (data already VST/log2 normalised)
#   3. Group differentiation: filled circle (Normal) vs filled triangle (Adjacent)
#   4. Translucent ellipse FILL (alpha=0.15) + solid border
#   5. PERMANOVA R^2 + p-value annotated bottom-right
#   6. MAD outlier detection (|z|>3 on PC1 or PC2): label sample name in grey italic
#   7. Larger points (size=2.2) + visible paired lines
factory_pca_scatter <- function(mat, title, add_paired_lines = TRUE, top_n = 2000) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[vars > 0 & !is.na(vars), , drop = FALSE]
  mat <- select_top_var(mat, n = top_n)

  grp <- make_group_vector(colnames(mat))
  patient_id <- gsub("Normal|Adjacent", "", colnames(mat))

  pca <- prcomp(t(mat), scale. = FALSE, center = TRUE)
  ve  <- summary(pca)$importance[2, 1:2] * 100

  df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    Group = factor(grp, levels = c("Normal", "Adjacent")),
    patient_id = patient_id,
    sample = colnames(mat),
    stringsAsFactors = FALSE
  )
  df <- detect_pca_outliers(df, threshold = 3)

  perm_res <- tryCatch(permanova_two_group(mat, grp, n_perm = 999),
                       error = function(e) list(R2 = NA, p = NA, F = NA))
  perm_label <- fmt_permanova(perm_res)
  cat(sprintf("    %s -> %s\n", title, perm_label))

  fills <- c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)

  p <- ggplot(df, aes(x = PC1, y = PC2))
  if (add_paired_lines) {
    p <- p + geom_line(aes(group = patient_id),
                       colour = "grey45", linewidth = 0.35,
                       alpha = 0.55, linetype = "solid")
  }
  p <- p +
    stat_ellipse(aes(fill = Group, colour = Group),
                 geom = "polygon", level = 0.90,
                 alpha = 0.15, linewidth = 0.45) +
    geom_point(aes(fill = Group, shape = Group),
               size = 2.2, stroke = 0.35, colour = "black") +
    scale_fill_manual(values = fills) +
    scale_colour_manual(values = fills) +
    scale_shape_manual(values = c(Normal = 21, Adjacent = 24)) +
    labs(title = title,
         subtitle = perm_label,
         x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.box.spacing = unit(1, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          plot.subtitle = element_text(family = FONT_FAMILY, size = 7,
                                       colour = "grey25",
                                       margin = margin(0, 0, 1, 0))) +
    guides(colour = "none",
           fill = guide_legend(override.aes = list(shape = c(21, 24),
                                                   colour = "black",
                                                   alpha = 1)),
           shape = "none")

  if (any(df$is_outlier)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(df, is_outlier),
      aes(label = sample),
      family = FONT_FAMILY, size = 2.3, fontface = "italic",
      colour = "grey25",
      min.segment.length = 0, segment.size = 0.25, segment.colour = "grey55",
      box.padding = 0.3, point.padding = 0.2,
      max.overlaps = Inf
    )
  }
  p
}

# --- 2.2b Combined Multi-omics PCA Factory ---
factory_combined_pca <- function(tc_mat, pr_mat, mb_mat, title) {
  # Find shared samples (transcriptomics has only 12 pairs = 24 samples)
  shared_samples <- Reduce(intersect, list(colnames(tc_mat), colnames(pr_mat), colnames(mb_mat)))
  cat(sprintf("    Combined PCA: %d shared samples\n", length(shared_samples)))
  
  # Select top variable features from each layer
  select_top <- function(mat, n) {
    v <- apply(mat, 1, var, na.rm = TRUE)
    mat[head(order(v, decreasing = TRUE), min(n, nrow(mat))), ]
  }
  
  tc_top  <- select_top(tc_mat[, shared_samples], 1000)
  pr_top  <- select_top(pr_mat[, shared_samples], 1000)
  met_top <- select_top(mb_mat[, shared_samples], 300)
  
  # Z-score scale each layer independently
  tc_z  <- t(scale(t(tc_top)));  tc_z <- tc_z[complete.cases(tc_z), , drop = FALSE]
  pr_z  <- t(scale(t(pr_top)));  pr_z <- pr_z[complete.cases(pr_z), , drop = FALSE]
  met_z <- t(scale(t(met_top))); met_z <- met_z[complete.cases(met_z), , drop = FALSE]
  
  # Prefix feature names to avoid collisions
  rownames(tc_z)  <- paste0("TC_", rownames(tc_z))
  rownames(pr_z)  <- paste0("PR_", rownames(pr_z))
  rownames(met_z) <- paste0("MET_", rownames(met_z))
  
  # Combine matrices
  combined <- rbind(tc_z, pr_z, met_z)
  cat(sprintf("    Combined matrix: %d features\n", nrow(combined)))
  
  # Run PCA
  grp <- make_group_vector(colnames(combined))
  patient_id <- gsub("Normal|Adjacent", "", colnames(combined))
  
  pca <- prcomp(t(combined), scale. = FALSE, center = TRUE)
  ve <- summary(pca)$importance[2, 1:2] * 100
  
  df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    Group = factor(grp, levels = c("Normal", "Adjacent")),
    patient_id = patient_id,
    sample = colnames(combined),
    stringsAsFactors = FALSE
  )
  df <- detect_pca_outliers(df, threshold = 3)

  perm_res <- tryCatch(permanova_two_group(combined, grp, n_perm = 999),
                       error = function(e) list(R2 = NA, p = NA, F = NA))
  perm_label <- fmt_permanova(perm_res)
  cat(sprintf("    Combined PCA -> %s\n", perm_label))

  fills <- c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)

  p <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_line(aes(group = patient_id),
              colour = "grey45", linewidth = 0.35, alpha = 0.55) +
    stat_ellipse(aes(fill = Group, colour = Group),
                 geom = "polygon", level = 0.90,
                 alpha = 0.15, linewidth = 0.45) +
    geom_point(aes(fill = Group, shape = Group),
               size = 2.2, stroke = 0.35, colour = "black") +
    scale_fill_manual(values = fills) +
    scale_colour_manual(values = fills) +
    scale_shape_manual(values = c(Normal = 21, Adjacent = 24)) +
    labs(title = title,
         subtitle = perm_label,
         x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.box.spacing = unit(1, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          plot.subtitle = element_text(family = FONT_FAMILY, size = 7,
                                       colour = "grey25",
                                       margin = margin(0, 0, 1, 0))) +
    guides(colour = "none",
           fill = guide_legend(override.aes = list(shape = c(21, 24),
                                                   colour = "black",
                                                   alpha = 1)),
           shape = "none")

  if (any(df$is_outlier)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(df, is_outlier),
      aes(label = sample),
      family = FONT_FAMILY, size = 2.3, fontface = "italic",
      colour = "grey25",
      min.segment.length = 0, segment.size = 0.25, segment.colour = "grey55",
      box.padding = 0.3, point.padding = 0.2,
      max.overlaps = Inf
    )
  }
  p
}

# --- 2.3 Correlation Heatmap Factory ---
# Enhancements vs original:
#   1. Auto color scale: floor = 5th-percentile of off-diagonal, ceil = 1.0
#      (avoid the 0.7 hard-coded floor that turned MT panel uniformly deep blue)
#   2. Enable column dendrogram so reviewers can see whether samples cluster by group
#   3. Prominent 4mm group annotation bar (vs 3mm) + bold name + outlined block
#   4. White grid line lwd 0.2 (was 0.3) - cleaner for many-sample MT panel
#   5. Group split: enables visual block-diagonal pattern when QC is good
factory_cor_heatmap <- function(mat, title, split_by_group = TRUE) {
  cor_mat <- cor(mat, use = "pairwise.complete.obs")
  grp <- make_group_vector(colnames(mat))
  grp_factor <- factor(grp, levels = c("Normal", "Adjacent"))

  off_diag <- cor_mat[lower.tri(cor_mat)]
  lo <- max(0.5, min(0.85, round(quantile(off_diag, 0.05, na.rm = TRUE), 2)))
  mid <- (lo + 1) / 2
  col_cor_qc <- colorRamp2(c(lo, mid, 1), c(COL_DOWN, "#FFFFFF", COL_UP))

  ha_top <- HeatmapAnnotation(
    Group = grp_factor,
    col = list(Group = c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)),
    annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    show_legend = FALSE,
    simple_anno_size = unit(4, "mm"),
    border = TRUE
  )

  ha_left <- rowAnnotation(
    Group = grp_factor,
    col = list(Group = c(Normal = COL_NORMAL, Adjacent = COL_ADJACENT)),
    annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    show_legend = TRUE,
    simple_anno_size = unit(4, "mm"),
    border = TRUE,
    annotation_legend_param = list(
      title_gp  = gp_legend_title(),
      labels_gp = gp_legend_labels(),
      ncol = 1,
      grid_height = unit(3, "mm"),
      grid_width  = unit(3, "mm")
    )
  )

  ht <- Heatmap(
    cor_mat,
    name              = "Pearson R",
    col               = col_cor_qc,
    top_annotation    = ha_top,
    left_annotation   = ha_left,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    show_row_dend     = FALSE,
    show_column_dend  = FALSE,
    column_split      = if (split_by_group) grp_factor else NULL,
    row_split         = if (split_by_group) grp_factor else NULL,
    cluster_row_slices    = FALSE,
    cluster_column_slices = FALSE,
    column_title         = title,
    column_title_gp      = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
    row_title            = NULL,
    border               = TRUE,
    rect_gp              = gpar(col = "white", lwd = 0.2),
    use_raster           = FALSE,
    heatmap_legend_param = list(
      title_gp     = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
      labels_gp    = gpar(fontsize = 8, fontfamily = FONT_GRID),
      grid_height  = unit(3, "mm"),
      grid_width   = unit(3, "mm"),
      legend_height = unit(20, "mm"),
      at = c(lo, mid, 1),
      labels = sprintf("%.2f", c(lo, mid, 1))
    )
  )
  attr(ht, "auto_lo") <- lo
  ht
}

# =============================================================================
# SECTION 9: Data Loading
# =============================================================================
cat("  Loading data files...\n")

tc_mat <- read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"),
                   row.names = 1, check.names = FALSE)
pr_mat <- read.csv(file.path(DATA, "proteomics_log2_norm.csv"),
                   row.names = 1, check.names = FALSE)
mb_mat <- read.csv(file.path(DATA, "metabolomics_log2_merged.csv"),
                   row.names = 1, check.names = FALSE)

cat("  Data loaded: TC=", ncol(tc_mat), "samples,", nrow(tc_mat), "genes\n")
cat("               PR=", ncol(pr_mat), "samples,", nrow(pr_mat), "proteins\n")
cat("               MB=", ncol(mb_mat), "samples,", nrow(mb_mat), "metabolites\n")



# --- Panel object registry (for optimized assembly) ---
# Strategy A: ggplot objects kept directly (best quality - vector)
# Strategy B: ComplexHeatmap captured via grid.grabExpr (single rasterization)
# Strategy C: External/pre-rendered PDF read at 600 DPI (fallback)
obj_A <- obj_B <- obj_C <- obj_D <- obj_E <- obj_F <- obj_G <- NULL
# =============================================================================
# SECTION 10: Render Supplementary Figure 1 -- QC - PCA & Sample Correlation
# =============================================================================
cat("\n--- Rendering Supplementary Figure 1: QC - PCA & Sample Correlation ---\n")

cat("=== SUPP FIGURE 1 ===\n")
tc_vst <- tc_mat; pr_log <- pr_mat; mb_log <- mb_mat

# Supp1A-C: Per-omics PCA plots (with paired lines)
# Order: Single-platform exploration (A-C: PCA) → Multi-omics integration (D: combined PCA) → QC validation (E-G: correlation)
panel_abc_names <- c("A", "B", "C")
for (idx in seq_along(list(
  list(mat = tc_vst, t = "Transcriptomics", o = "Supp01A_PCA_TC.pdf"),
  list(mat = pr_log, t = "Proteomics",      o = "Supp01B_PCA_PR.pdf"),
  list(mat = mb_log, t = "Metabolomics",     o = "Supp01C_PCA_MT.pdf")
))) {
  info <- list(
    list(mat = tc_vst, t = "Transcriptomics", o = "Supp01A_PCA_TC.pdf"),
    list(mat = pr_log, t = "Proteomics",      o = "Supp01B_PCA_PR.pdf"),
    list(mat = mb_log, t = "Metabolomics",     o = "Supp01C_PCA_MT.pdf")
  )[[idx]]
  tryCatch({
    p <- factory_pca_scatter(info$mat, info$t, add_paired_lines = TRUE)
    save_pdf(info$o, 61, 80, print(p))
    assign(paste0("obj_", panel_abc_names[idx]), p, envir = .GlobalEnv)
  }, error = function(e) cat("  ERROR", info$o, ":", e$message, "\n"))
}

# Supp1D: Combined Multi-omics PCA (using shared samples only)
cat("  Generating combined multi-omics PCA...\n")
tryCatch({
  p_combined <- factory_combined_pca(tc_vst, pr_log, mb_log, "Combined Multi-omics PCA")
  save_pdf("Supp01D_PCA_combined.pdf", 61, 80, print(p_combined))
  obj_D <- p_combined  # Strategy A: keep ggplot object for direct assembly
}, error = function(e) cat("  ERROR Supp1D_PCA_combined:", e$message, "\n"))

# Supp1E-F: Correlation heatmaps (narrow panels)
panel_ef_names <- c("E", "F")
ef_list <- list(
  list(mat = tc_vst, t = "Sample correlation - TC", o = "Supp01E_cor_TC.pdf"),
  list(mat = pr_log, t = "Sample correlation - PR", o = "Supp01F_cor_PR.pdf")
)
for (idx in seq_along(ef_list)) {
  info <- ef_list[[idx]]
  tryCatch({
    s <- sp(info$o, 61, 80)
    cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
    ht <- factory_cor_heatmap(info$mat, info$t)
    draw(ht, padding = unit(c(5, 5, 5, 5), "mm")); dev.off()
    assign(paste0("ht_", panel_ef_names[idx]), ht, envir = .GlobalEnv)
    cat("  ->", info$o, "\n")
  }, error = function(e) cat("  ERROR", info$o, ":", e$message, "\n"))
}

# Supp1G: Correlation heatmap (full-width panel)
tryCatch({
  s <- sp("Supp01G_cor_MT.pdf", 183, 85)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  ht <- factory_cor_heatmap(mb_log, "Sample correlation - MT")
  draw(ht, padding = unit(c(5, 5, 5, 5), "mm")); dev.off()
  ht_G <- ht
  cat("  -> Supp01G_cor_MT.pdf\n")
}, error = function(e) cat("  ERROR Supp01G_cor_MT:", e$message, "\n"))


# =============================================================================
# =============================================================================
# SECTION 11: Composite Figure Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite SuppFig_01 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  # Row1(80mm): A(61)+B(61)+C(61), Row2(80mm): D(61)+E(61)+F(61), Row3(85mm): G(183)
  H1 <- 80; H2 <- 80; H3 <- H_TOTAL - H1 - H2
  W1 <- 61; W2 <- 61; W3 <- W_TOTAL - W1 - W2

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    # Panel A (Row1 col1)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H2+H3,"mm"),
      width=unit(W1,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(obj_A, newpage=FALSE); grid::popViewport()
    # Panel B (Row1 col2)
    grid::pushViewport(grid::viewport(x=unit(W1,"mm"), y=unit(H2+H3,"mm"),
      width=unit(W2,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(obj_B, newpage=FALSE); grid::popViewport()
    # Panel C (Row1 col3)
    grid::pushViewport(grid::viewport(x=unit(W1+W2,"mm"), y=unit(H2+H3,"mm"),
      width=unit(W3,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(obj_C, newpage=FALSE); grid::popViewport()
    # Panel D (Row2 col1)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H3,"mm"),
      width=unit(W1,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(obj_D, newpage=FALSE); grid::popViewport()
    # Panel E (Row2 col2) - heatmap
    grid::pushViewport(grid::viewport(x=unit(W1,"mm"), y=unit(H3,"mm"),
      width=unit(W2,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    draw(ht_E, padding=unit(c(5,5,5,5),"mm"), newpage=FALSE); grid::popViewport()
    # Panel F (Row2 col3) - heatmap
    grid::pushViewport(grid::viewport(x=unit(W1+W2,"mm"), y=unit(H3,"mm"),
      width=unit(W3,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    draw(ht_F, padding=unit(c(5,5,5,5),"mm"), newpage=FALSE); grid::popViewport()
    # Panel G (Row3 full-width) - heatmap
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_TOTAL,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    draw(ht_G, padding=unit(c(5,5,5,5),"mm"), newpage=FALSE); grid::popViewport()
    # Labels
    label_data <- data.frame(
      text=c("a","b","c","d","e","f","g"),
      x_mm=c(2,W1+2,W1+W2+2,2,W1+2,W1+W2+2,2),
      y_mm=c(H2+H3+H1-2,H2+H3+H1-2,H2+H3+H1-2,H3+H2-2,H3+H2-2,H3+H2-2,H3-2),
      stringsAsFactors=FALSE)
    for(i in seq_len(nrow(label_data))) {
      grid::grid.text(label=label_data$text[i],
        x=unit(label_data$x_mm[i],"mm"), y=unit(label_data$y_mm[i],"mm"),
        just=c("left","top"),
        gp=grid::gpar(fontsize=FS_TAG, fontface="bold", fontfamily=FONT_FAMILY))
    }
    grid::popViewport()
  }

  cairo_pdf(file.path(OUT,"SuppFig_01.pdf"),
            width=W_TOTAL/25.4, height=H_TOTAL/25.4, family=FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_01.pdf (vector)\n")

  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT,"SuppFig_01.png"),
                 width=px_W, height=px_H, res=DPI, type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_01.png\n")

  grDevices::tiff(file.path(OUT,"SuppFig_01.tiff"),
                  width=px_W, height=px_H, res=DPI, compression="lzw", type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_01.tiff\n")

  cat("  SuppFig_01 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))
# END
# =============================================================================
cat("\n=== Supplementary Figure 1 rendering complete ===\n")

