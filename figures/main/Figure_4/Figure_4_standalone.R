#!/usr/bin/env Rscript
# =============================================================================
# Figure_4_standalone.R
# Complete, Self-Contained Code for Figure 4
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: EBioMedicine (Lancet family) - TIFF 600dpi (optimized assembly) + PDF
# =============================================================================
# Figure 4: Immune Microenvironment
# 5 panels: ssGSEA heatmap + immune barplot + immune NES + checkpoint heatmap + immune-metab coupling
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript Figure_4_standalone.R
# =============================================================================

cat("=== Figure 4: Immune Microenvironment ===\n")
cat("  Loading libraries...\n")


# =============================================================================
# SECTION 1: Library Imports + Arial Font Registration
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggpubr)
  library(ggsci)
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
OUT  <- file.path(BASE, "04_figures/main/Figure_4")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants (EBioMedicine / Lancet family Standard)
# =============================================================================
FONT_FAMILY <- "Arial"      # AI native font - better than Arial for AI editing
FONT_GRID   <- "Arial"      # unified for ComplexHeatmap
MM_PER_INCH <- 25.4

# EBioMedicine page dimensions (mm)
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
FS_AXIS_TITLE <- 8     # axis.title (bold)
FS_AXIS_TEXT  <- 8     # axis.text, strip.text (minimum)
FS_LEGEND_T   <- 8     # legend.title (bold)
FS_LEGEND_L   <- 8     # legend.text
FS_ANNO       <- 8     # annotation name
FS_ROW_NAME   <- 8     # heatmap row/column names
FS_CELL       <- 8     # heatmap cell text
FS_TAG        <- 16    # panel tag (A, B, C...)
FS_GEOM_TEXT  <- 2.82  # geom_text size (= 8pt in mm)

# =============================================================================
# SECTION 4: Color Palette (ggsci JCO / Journal of Clinical Oncology)
# =============================================================================
PAL_CAT <- pal_jco("default")(10)

# Semantic colors
COL_UP       <- "#CD534CFF"    # up-regulated (NPG red)
COL_DOWN     <- "#0073C2FF"    # down-regulated (NPG cyan)
COL_NS       <- "#868686FF"    # not significant (darker grey for better contrast)
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
# SECTION 6: ggplot2 Theme (theme_bw base, EBioMedicine / Lancet family style)
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

# --- Save helper ---
sp <- function(filename, w_mm, h_mm, out_dir = NULL) {
  if (is.null(out_dir)) out_dir <- OUT
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  list(path = file.path(out_dir, filename),
       width = target_w_px / ASSEMBLY_DPI, height = target_h_px / ASSEMBLY_DPI)
}

save_pdf <- function(filename, w, h, expr) {
  fpath <- file.path(OUT, filename)
  target_w_px <- round(w * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h * ASSEMBLY_DPI / 25.4)
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  force(expr)
  dev.off()
  cat(sprintf("  Panel saved: %s (%d\u00d7%d px)\n", basename(fpath), target_w_px, target_h_px))
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
        # Clean ENSG IDs: strip version/chr suffixes to match ensg_to_symbol() cleaning
        sub$ensg <- gsub("\\.[0-9]+$", "", sub$ensg)
        sub$ensg <- gsub("_[0-9]+$", "", sub$ensg)
        # Remove rows where symbol is NA, empty, same as ensg, or still an ENSG ID
        sub <- sub[!is.na(sub$symbol) & sub$symbol != "" &
                   sub$symbol != sub$ensg & !grepl("^ENSG[0-9]", sub$symbol), ]
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

# --- Clean pathway / gene-set names for display ---
clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_|^REACTOME_|^KEGG_|^GOBP_|^GO_|^WP_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x
}

# --- Abbreviate metabolite names (lipid classes, greek letters, etc.) ---
abbreviate_metabolite <- function(names) {
  x <- names
  x <- gsub("glycerophosphoethanolamine", "-GPE", x, ignore.case = TRUE)
  x <- gsub("glycerophosphocholine", "-GPC", x, ignore.case = TRUE)
  x <- gsub("glycerophosphoserine", "-GPS", x, ignore.case = TRUE)
  x <- gsub("glycerophosphoinositol", "-GPI", x, ignore.case = TRUE)
  x <- gsub("lysophosphatidylcholine", "LPC", x, ignore.case = TRUE)
  x <- gsub("lysophosphatidylethanolamine", "LPE", x, ignore.case = TRUE)
  x <- gsub("phosphatidylcholine", "PC", x, ignore.case = TRUE)
  x <- gsub("phosphatidylethanolamine", "PE", x, ignore.case = TRUE)
  x <- gsub("sphingomyelin", "SM", x, ignore.case = TRUE)
  x <- gsub("([0-9]+[a-zA-Z]*)-[Hh]ydroxy-?", "\\1-OH-", x)
  x <- gsub("^[Hh]ydroxy-?", "OH-", x)
  x <- gsub("alpha", "\u03b1", x, ignore.case = TRUE)
  x <- gsub("beta", "\u03b2", x, ignore.case = TRUE)
  x <- gsub("gamma", "\u03b3", x, ignore.case = TRUE)
  x <- gsub("androst-4-ene-3,17-dione", "androstenedione", x, ignore.case = TRUE)
  x <- gsub("Nicotinamide adenine dinucleotide.*", "NAD", x, ignore.case = TRUE)
  x <- gsub("METHYL", "methyl", x)
  x <- gsub("PENTANEDIOL", "pentanediol", x)
  x <- gsub("HEXANEDIOL", "hexanediol", x)
  x <- gsub("--+", "-", x)
  x <- gsub("- ", " ", x)
  x <- gsub("-$", "", x)
  x <- gsub("^-", "", x)
  x
}


# --- Factory Functions (from production pipeline) ---

# --- 2.9 Direction Bar Factory (Up/Down or CS1/CS2 fill) ---
factory_direction_bar <- function(df, name_col, value_col, title,
                                  colors = c(Up = COL_UP, Down = COL_DOWN),
                                  x_lab = NULL, name_size = 8) {
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

# =============================================================================
# SECTION 9: Data Loading
# =============================================================================
cat("  Data loaded inline in render logic below.\n")



# --- Panel object registry (for optimized assembly) ---
# Strategy A: ggplot objects kept directly (best quality - vector)
# Strategy B: ComplexHeatmap captured via grid.grabExpr (single rasterization)
# Strategy C: External/pre-rendered PDF read at 600 DPI (fallback)
obj_A <- obj_B <- obj_C <- obj_D <- obj_E <- NULL
# =============================================================================
# SECTION 10: Render Figure 4 -- Immune Microenvironment
# =============================================================================
cat("\n--- Rendering Figure 4: Immune Microenvironment ---\n")

cat("\n=== FIGURE 4: Immune Microenvironment ===\n")

# Fig 4a: ssGSEA immune cell delta barplot — TRANSCRIPTOMICS (with Padj stars)
cat("  Fig4a (transcriptomic ssGSEA)\n")
tryCatch({
  ssgsea_im <- read.csv(file.path(RES, "phase4_immune/ssGSEA_transcriptomics.csv"),
                        stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  imd_tc <- read.csv(file.path(RES, "phase4_immune/immune_diff_transcriptomics.csv"),
                     stringsAsFactors = FALSE)
  # Beautify cell-type names (raw rownames -> publication labels)
  pretty_celltype <- function(x) {
    x <- gsub("_", " ", x)
    x <- gsub("Gamma delta T",            "\u03b3\u03b4 T cells", x)
    x <- gsub("Activated dendritic cells", "aDCs",   x)
    x <- gsub("Immature dendritic cells",  "iDCs",   x)
    x <- gsub("Plasmacytoid dendritic cells", "pDCs", x)
    x <- gsub("Natural killer cells",      "NK cells", x)
    x <- gsub("Central memory",            "CM",     x)
    x <- gsub("Effector memory",           "EM",     x)
    x <- gsub("regulatory T cells",        "Tregs",  x, ignore.case = TRUE)
    x <- gsub("Treg cells",                "Tregs",  x)
    x
  }
  rownames(ssgsea_im) <- pretty_celltype(rownames(ssgsea_im))
  imd_tc$cell_type    <- pretty_celltype(imd_tc$cell_type)

  mat_im <- scale_rows(ssgsea_im, 2.5)
  is_adj_a <- make_group_vector(colnames(mat_im)) == "Adjacent"
  delta_a  <- rowMeans(mat_im[,  is_adj_a, drop = FALSE]) -
              rowMeans(mat_im[, !is_adj_a, drop = FALSE])
  d4a <- data.frame(
    name  = names(delta_a),
    value = as.numeric(delta_a),
    stringsAsFactors = FALSE
  )
  d4a$padj <- imd_tc$padj[match(d4a$name, imd_tc$cell_type)]
  d4a$sig  <- ifelse(is.na(d4a$padj), "",
              ifelse(d4a$padj < 0.001, "***",
              ifelse(d4a$padj < 0.01,  "**",
              ifelse(d4a$padj < 0.05,  "*",
              ifelse(d4a$padj < 0.1,   "\u2020", "")))))
  d4a <- d4a[order(d4a$value, decreasing = TRUE), ]
  d4a$direction <- ifelse(d4a$value > 0, "Up", "Down")
  d4a$name <- factor(d4a$name, levels = d4a$name[order(d4a$value)])

  p4a <- ggplot(d4a, aes(x = value, y = name, fill = direction)) +
    geom_col(width = 0.72, alpha = 0.92) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_text(aes(label = sig,
                  x = ifelse(value > 0, value + 0.04, value - 0.04),
                  hjust = ifelse(value > 0, 0, 1)),
              size = FS_GEOM_TEXT, family = FONT_FAMILY) +
    scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    labs(title = "Immune Cells (Transcriptomics)",
         x = "Adjacent \u2212 Normal (\u0394Z-score)",
         y = NULL) +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5))
  obj_A <- p4a
  save_pdf("Fig4a_ssGSEA_heatmap.pdf", 80, 85, print(p4a))
  cat("    -> saved (TC; max Padj_min =", round(min(d4a$padj, na.rm = TRUE), 4), ")\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# Fig 4b: PROTEOMIC ssGSEA immune cell delta barplot (mast-cell enrichment, Padj=0.032)
# Replaces previous transcriptomic redundant panel (was redundant with Fig4a)
cat("  Fig4b (proteomic ssGSEA)\n")
tryCatch({
  imd_pr <- read.csv(file.path(RES, "phase4_immune/immune_diff_proteomics.csv"),
                     stringsAsFactors = FALSE)
  pretty_celltype <- function(x) {
    x <- gsub("_", " ", x)
    x <- gsub("Mast cells",        "Mast cells",  x)
    x <- gsub("Plasma cells",      "Plasma cells", x)
    x <- gsub("Treg cells",        "Tregs",        x)
    x <- gsub("Th17 cells",        "Th17 cells",   x)
    x <- gsub("Th2 cells",         "Th2 cells",    x)
    x <- gsub("CD4 T cells",       "CD4 T cells",  x)
    x <- gsub("CD8 T cells",       "CD8 T cells",  x)
    x <- gsub("NK cells",          "NK cells",     x)
    x <- gsub("Macrophages M1",    "Macrophages M1", x)
    x <- gsub("Macrophages M2",    "Macrophages M2", x)
    x <- gsub("Dendritic cells",   "DCs",          x)
    x <- gsub("Endothelial",       "Endothelial",  x)
    x <- gsub("Fibroblasts",       "Fibroblasts",  x)
    x <- gsub("Eosinophils",       "Eosinophils",  x)
    x <- gsub("Neutrophils",       "Neutrophils",  x)
    x <- gsub("B cells",           "B cells",      x)
    x
  }
  imd_pr$cell_type <- pretty_celltype(imd_pr$cell_type)
  imd_pr <- imd_pr[order(imd_pr$diff, decreasing = TRUE), ]
  imd_pr$direction <- ifelse(imd_pr$diff > 0, "Up", "Down")
  imd_pr$sig <- ifelse(imd_pr$padj < 0.001, "***",
                ifelse(imd_pr$padj < 0.01,  "**",
                ifelse(imd_pr$padj < 0.05,  "*",
                ifelse(imd_pr$padj < 0.1,   "\u2020", ""))))
  imd_pr$cell_type <- factor(imd_pr$cell_type, levels = imd_pr$cell_type[order(imd_pr$diff)])

  p4b <- ggplot(imd_pr, aes(x = diff, y = cell_type, fill = direction)) +
    geom_col(width = 0.72, alpha = 0.92) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_text(aes(label = sig,
                  x = ifelse(diff > 0, diff + 0.005, diff - 0.005),
                  hjust = ifelse(diff > 0, 0, 1)),
              size = FS_GEOM_TEXT, family = FONT_FAMILY) +
    scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    labs(title = "Immune Cells (Proteomics)",
         x = "Score difference (Adjacent \u2212 Normal)",
         y = NULL) +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5))
  p4c <- p4b   # alias for legacy composite reference
  save_pdf("Fig4b_immune_barplot.pdf", 80, 85, print(p4b))
  cat("    -> saved (PR; mast cells Padj =",
      round(imd_pr$padj[imd_pr$cell_type == "Mast cells"], 4), ")\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# Fig 3c: Immune hallmark NES
cat("  Fig4c\n")
tryCatch({
  ihn <- read.csv(file.path(RES, "phase4_immune/immune_hallmark_NES.csv"), stringsAsFactors = FALSE)
  # Explicit pathway name mapping to preserve abbreviations
  pathway_name_map <- c(
    "ALLOGRAFT REJECTION" = "Allograft Rejection",
    "COAGULATION" = "Coagulation",
    "COMPLEMENT" = "Complement",
    "IL2 STAT5 SIGNALING" = "IL2-STAT5 Signaling",
    "IL6 JAK STAT3 SIGNALING" = "IL6-JAK-STAT3 Signaling",
    "INFLAMMATORY RESPONSE" = "Inflammatory Response",
    "INTERFERON ALPHA RESPONSE" = "Interferon Alpha Response",
    "INTERFERON GAMMA RESPONSE" = "Interferon Gamma Response",
    "KRAS SIGNALING UP" = "KRAS Signaling Up",
    "TNFA SIGNALING VIA NFKB" = "TNF\u03b1 Signaling via NF-\u03baB"
  )
  ihn$pw_name <- ifelse(ihn$pathway_label %in% names(pathway_name_map),
                        pathway_name_map[ihn$pathway_label],
                        clean_pathway(ihn$pathway_label))
  p4d <- factory_direction_bar(ihn, "pw_name", "NES_TC", "Immune-Related Hallmark Pathways",
                               x_lab = "NES (Transcriptomics)") +
    theme(plot.title.position = "plot")
  save_pdf("Fig4c_immune_NES.pdf", 83, 85, print(p4d))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# Fig 3d: Checkpoint markers delta barplot
cat("  Fig4d\n")
tryCatch({
  tc_mat <- as.matrix(read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"),
                               row.names = 1, check.names = FALSE))
  rownames(tc_mat) <- ensg_to_symbol(rownames(tc_mat))
  if (any(duplicated(rownames(tc_mat)))) {
    rv <- apply(tc_mat, 1, var, na.rm = TRUE)
    tc_mat <- tc_mat[order(-rv), ]
    tc_mat <- tc_mat[!duplicated(rownames(tc_mat)), ]
  }
  checkpoint_genes <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "CD274",
                        "PDCD1LG2", "BTLA", "VSIR", "CD28", "CD80", "CD86",
                        "ICOS", "CD40", "CD40LG", "TNFRSF4", "TNFRSF9",
                        "ADORA2A", "TOX", "ENTPD1", "LAYN", "CXCL13",
                        "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1")
  found <- intersect(checkpoint_genes, rownames(tc_mat))
  cat("    Found", length(found), "of", length(checkpoint_genes), "checkpoint genes\n")
  mat_ck <- tc_mat[found, , drop = FALSE]
  mat_ck <- scale_rows(mat_ck, 2.5)
  is_adj_d <- make_group_vector(colnames(mat_ck)) == "Adjacent"
  delta_d   <- rowMeans(mat_ck[,  is_adj_d, drop=FALSE]) -
               rowMeans(mat_ck[, !is_adj_d, drop=FALSE])
  ord_d     <- order(delta_d, decreasing = TRUE)
  delta_d   <- delta_d[ord_d]
  nm_d      <- rownames(mat_ck)[ord_d]
  p4d_ck <- factory_direction_bar(
    data.frame(name = nm_d, value = delta_d, stringsAsFactors = FALSE),
    "name", "value",
    title  = "Checkpoint Markers (Transcriptomics)",
    x_lab  = "Adjacent \u2212 Normal (\u0394Z-score)",
    name_size = 7.5
  ) + theme(axis.text.y = element_text(face = "italic"),
            legend.position = "none")
  obj_D <- p4d_ck
  save_pdf("Fig4d_checkpoint_heatmap.pdf", 80, 110, print(p4d_ck))
  cat("    -> saved\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# --- Panel E: Immune-Metabolic Coupling Heatmap ---
cat("  Fig4e\n")
tryCatch({
  cor_e <- as.matrix(read.csv(file.path(RES, "enhancement10_mechanism/immune_metabolic_correlation.csv"),
                              row.names = 1, check.names = FALSE))
  pval_e <- as.matrix(read.csv(file.path(RES, "enhancement10_mechanism/immune_metabolic_pvalues.csv"),
                               row.names = 1, check.names = FALSE))
  # Clean names with proper capitalization (preserve abbreviations)
  row_name_map <- c(
    "ICP_INHIBITORY" = "Inhibitory Checkpoints",
    "ICP_STIMULATORY" = "Stimulatory Checkpoints",
    "T_CELL_EXHAUSTION" = "T Cell Exhaustion",
    "MYELOID_CHECKPOINT" = "Myeloid Checkpoint"
  )
  col_name_map <- c(
    "GLYCOLYSIS" = "Glycolysis",
    "TCA_CYCLE" = "TCA Cycle",
    "FAT_ACID_OXIDATION" = "Fatty Acid Oxidation",
    "FAT_ACID_SYNTHESIS" = "Fatty Acid Synthesis",
    "BILE_ACID_SYNTH" = "Bile Acid Synthesis",
    "UREA_CYCLE" = "Urea Cycle",
    "OXPHOS" = "OXPHOS",
    "PENTOSE_PHOSPHATE" = "Pentose Phosphate",
    "GLUTAMINE_METAB" = "Glutamine Metabolism"
  )
  rownames(cor_e) <- ifelse(rownames(cor_e) %in% names(row_name_map),
                            row_name_map[rownames(cor_e)], rownames(cor_e))
  colnames(cor_e) <- ifelse(colnames(cor_e) %in% names(col_name_map),
                            col_name_map[colnames(cor_e)], colnames(cor_e))
  rownames(pval_e) <- rownames(cor_e)
  colnames(pval_e) <- colnames(cor_e)
  
  sig_mat_e <- ifelse(pval_e < 0.001, "***", ifelse(pval_e < 0.01, "**", ifelse(pval_e < 0.05, "*", "")))
  
  col_fn_3e <- colorRamp2(c(-1, 0, 1), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
  ht_3e <- Heatmap(cor_e, name = "Correlation", col = col_fn_3e,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = gp_row_names(8),
    column_names_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(sig_mat_e[i, j], x, y, gp = gpar(fontsize = 12, fontfamily = FONT_GRID))
    },
    heatmap_legend_param = std_legend_param(),
    row_names_max_width = unit(50, "mm"))

  s <- sp("Fig4e_immune_metab_coupling.pdf", 183, 75)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  draw(ht_3e, padding = unit(c(3, 3, 5, 15), "mm"),
    column_title = "Immune-Metabolic Coupling",
    column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_GRID))
  dev.off()
  cat(sprintf("    -> saved (%.0fx%.0fmm)\n", s$width * 25.4, s$height * 25.4))
}, error = function(e) cat("    ERROR:", e$message, "\n"))


# =============================================================================

# =============================================================================
# SECTION 11: Composite Figure Assembly (Pure Vector Grid Viewport — AI-editable)
# =============================================================================
cat("\n--- Assembling composite Figure_4 (VECTOR, 183x245mm) ---\n")

tryCatch({
  DPI <- 600
  W_TOTAL <- 183; H_TOTAL <- 265

  # Row 1 ( 85mm): A(103mm) + B(80mm) = 183mm
  # Row 2 (105mm): C(83mm) + D(100mm) = 183mm  — taller for 28 checkpoint genes
  # Row 3 ( 75mm): E(183mm) = full width
  W_A <- 103; W_B <- W_TOTAL - W_A
  W_C <- 83;  W_D <- W_TOTAL - W_C
  H1 <- 85; H2 <- 105; H3 <- H_TOTAL - H1 - H2

  # --- Render function: places all panels + labels into current device ---
  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL),
      clip = "off"
    ))

    # Row 1 (top): A + B — y = H2 + H3 from bottom
    y_row1 <- H2 + H3

    # Panel a (ggplot - ssGSEA delta barplot)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row1, "mm"),
      width = unit(W_A, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(p4a, newpage = FALSE)
    grid::popViewport()

    # Panel b (ggplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_A, "mm"), y = unit(y_row1, "mm"),
      width = unit(W_B, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(p4c, newpage = FALSE)
    grid::popViewport()

    # Row 2 (middle): C + D — y = H3 from bottom
    y_row2 <- H3

    # Panel c (ggplot)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_C, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(p4d, newpage = FALSE)
    grid::popViewport()

    # Panel d (ggplot - checkpoint delta barplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_C, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_D, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(p4d_ck, newpage = FALSE)
    grid::popViewport()

    # Row 3 (bottom): E full-width — y = 0
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W_TOTAL, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom"),
      clip = "off"
    ))
    draw(ht_3e, padding = unit(c(3, 3, 8, 15), "mm"),
      column_title = "Immune-Metabolic Coupling",
      column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_GRID),
      newpage = FALSE)
    grid::popViewport()

    # --- Panel labels (bold, FS_TAG pt) ---
    label_data <- data.frame(
      text  = c("A", "B", "C", "D", "E"),
      x_mm  = c(1, W_A + 1, 1, W_C + 1, 1),
      y_mm  = c(H_TOTAL - 1, H_TOTAL - 1, y_row1 - 1, y_row1 - 1, H3 - 1),
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
  cairo_pdf(file.path(OUT, "Figure_4.pdf"),
            width = W_TOTAL / 25.4, height = H_TOTAL / 25.4, family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> Figure_4.pdf (VECTOR, AI-editable)\n")

  # --- Save PNG (600 DPI for review) ---
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "Figure_4.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_4.png\n")

  # --- Save TIFF (600 DPI for journal submission) ---
  grDevices::tiff(file.path(OUT, "Figure_4.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_4.tiff\n")

  cat("  Figure_4 DONE (183x265mm, VECTOR PDF + 600DPI PNG/TIFF)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

# =============================================================================
cat("\n=== Figure 4 rendering complete ===\n")

