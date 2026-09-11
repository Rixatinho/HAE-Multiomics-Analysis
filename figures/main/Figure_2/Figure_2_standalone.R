#!/usr/bin/env Rscript
# =============================================================================
# Figure_2_standalone.R
# Complete, Self-Contained Code for Figure 2
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: EBioMedicine (Lancet family) - TIFF 600dpi (optimized assembly) + PDF
# =============================================================================
# Figure 2: Metabolic Reprogramming
# 6 panels: A=metabolic bubble, B=enzyme heatmap, C=DEMs heatmap,
#           D=MSEA dotplot, E=prot-metab network, F=GSEA barplot
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript Figure_2_standalone.R
# =============================================================================

cat("=== Figure 2: Metabolic Reprogramming ===\n")
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
  library(ggsci)
  library(ggrepel)
  library(igraph)
  library(ggnewscale)
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
OUT  <- file.path(BASE, "04_figures/main/Figure_2")
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
FS_AXIS_TITLE <- 9     # axis.title (bold)
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
    plot.title         = element_text(family = FONT_FAMILY, size = 10, face = "bold", hjust = 0),
    plot.title.position = "plot",
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

# --- Clean pathway / gene-set names for display ---
clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_|^REACTOME_|^KEGG_|^GOBP_|^GO_|^WP_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x
}

# --- Abbreviate metabolite names for compact display ---
abbreviate_metabolite <- function(names) {
  x <- names
  # --- Polyunsaturated fatty-acid common-name shorthand (resolved BEFORE wrap) ---
  # These long IUPAC names exceed 30 chars and would otherwise fall back to
  # Compound_ID. Replace with biochemist-recognized common names.
  fa_map <- c(
    # Source name (lowercased contains-match) -> Display name
    "cis-7,10,13,16-docosic acidtraenoic acid" = "Adrenic acid (22:4n-6)",
    "cis-7,10,13,16-docosatetraenoic acid"     = "Adrenic acid (22:4n-6)",
    "cis-7,10,13,16,19-docosapentaenoic acid"  = "DPA (22:5n-3)",
    "cis-4,7,10,13,16-docosapentaenoic acid"   = "DPA (22:5n-6)",
    "cis-11,14,17-eicosatrienoic acid"         = "ETE (20:3n-3)",
    "cis-8,11,14-eicosatrienoic acid"          = "DGLA (20:3n-6)",
    "cis-5,8,11,14-eicosatetraenoic acid"      = "Arachidonic acid (20:4n-6)",
    "cis-5,8,11,14,17-eicosapentaenoic acid"   = "EPA (20:5n-3)",
    "cis-4,7,10,13,16,19-docosahexaenoic acid" = "DHA (22:6n-3)",
    "alpha-linolenic acid"                     = "\u03b1-Linolenic acid (18:3n-3)",
    "gamma-linolenic acid"                     = "\u03b3-Linolenic acid (18:3n-6)"
  )
  for (k in names(fa_map)) {
    hit <- grepl(k, tolower(x), fixed = TRUE)
    if (any(hit)) x[hit] <- fa_map[[k]]
  }
  # Long acyl-carnitine: (13Z,16Z,19Z)-Docosa-13,16,19-trienoylcarnitine -> 22:3-Carnitine
  x <- gsub("\\(?[0-9]{1,2}[EZ],?\\s*[0-9]{1,2}[EZ],?\\s*[0-9]{1,2}[EZ]\\)?-?\\s*Docosa-?[0-9,]*-trienoylcarnitine",
            "22:3-Carnitine", x, ignore.case = TRUE)
  x <- gsub("\\(?[0-9]{1,2}[EZ],?\\s*[0-9]{1,2}[EZ],?\\s*[0-9]{1,2}[EZ],?\\s*[0-9]{1,2}[EZ]\\)?-?\\s*Docosa-?[0-9,]*-tetraenoylcarnitine",
            "22:4-Carnitine", x, ignore.case = TRUE)
  # Lipid class abbreviations (longest patterns first)
  x <- gsub("glycerophosphoethanolamine", "-GPE", x, ignore.case = TRUE)
  x <- gsub("glycerophosphocholine", "-GPC", x, ignore.case = TRUE)
  x <- gsub("glycerophosphoserine", "-GPS", x, ignore.case = TRUE)
  x <- gsub("glycerophosphoinositol", "-GPI", x, ignore.case = TRUE)
  x <- gsub("lysophosphatidylcholine", "LPC", x, ignore.case = TRUE)
  x <- gsub("lysophosphatidylethanolamine", "LPE", x, ignore.case = TRUE)
  x <- gsub("phosphatidylcholine", "PC", x, ignore.case = TRUE)
  x <- gsub("phosphatidylethanolamine", "PE", x, ignore.case = TRUE)
  x <- gsub("sphingomyelin", "SM", x, ignore.case = TRUE)
  # Positional hydroxy -> OH (before greek letter substitution)
  x <- gsub("([0-9]+[a-zA-Z]*)-[Hh]ydroxy-?", "\\1-OH-", x)
  x <- gsub("^[Hh]ydroxy-?", "OH-", x)
  # Greek letters
  x <- gsub("alpha", "\u03b1", x, ignore.case = TRUE)
  x <- gsub("beta", "\u03b2", x, ignore.case = TRUE)
  x <- gsub("gamma", "\u03b3", x, ignore.case = TRUE)
  # Steroid core abbreviations
  x <- gsub("androst-4-ene-3,17-dione", "androstenedione", x, ignore.case = TRUE)
  # Common biochemical abbreviations
  x <- gsub("Nicotinamide adenine dinucleotide.*", "NAD", x, ignore.case = TRUE)
  # Normalize ALL-CAPS segments to lowercase
  x <- gsub("METHYL", "methyl", x)
  x <- gsub("PENTANEDIOL", "pentanediol", x)
  x <- gsub("HEXANEDIOL", "hexanediol", x)
  # Clean up artifacts (double hyphens, trailing hyphens before space/end)
  x <- gsub("--+", "-", x)
  x <- gsub("- ", " ", x)
  x <- gsub("-$", "", x)
  x <- gsub("^-", "", x)
  x
}

# --- Format metabolite display names (>25 chars: wrap; >40 chars: use ID) ---
format_metabolite_display <- function(abbr_names, compound_ids = NULL,
                                      width = 25, max_len = 40) {
  x <- abbr_names
  id_replacements <- character(0)  # track ID replacements for figure legend
  for (i in seq_along(x)) {
    nc <- nchar(x[i])
    if (nc > max_len && !is.null(compound_ids)) {
      # Very long name: replace with compound ID
      id_replacements <- c(id_replacements,
                           paste0(compound_ids[i], " = ", x[i]))
      x[i] <- compound_ids[i]
    } else if (nc > width) {
      # Moderately long: wrap at natural break point near 'width'
      left <- substr(x[i], 1, min(width + 2, nc))
      brk_pos <- gregexpr("[- ,/(]", left)[[1]]
      brk_pos <- brk_pos[brk_pos >= 10]  # don't break too early
      if (length(brk_pos) > 0) {
        brk <- max(brk_pos)
        x[i] <- paste0(substr(x[i], 1, brk), "\n", substring(x[i], brk + 1))
      } else {
        x[i] <- paste0(substr(x[i], 1, width), "\n", substring(x[i], width + 1))
      }
    }
  }
  if (length(id_replacements) > 0) {
    cat("    [Legend note] Metabolite ID replacements:\n")
    for (r in id_replacements) cat("      ", r, "\n")
  }
  x
}
# --- Enhanced Visualization Functions (from nc_theme.R) ---

# --- make_nc_group_anno() -- Standardized group annotation bar ---
make_nc_group_anno <- function(meta, group_col = "Group") {
  groups <- meta[[group_col]]
  HeatmapAnnotation(
    Group = groups,
    col = list(Group = c("Normal" = COL_DOWN, "Adjacent" = COL_UP)),
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm"),
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
  )
}

# =============================================================================
# SECTION 9: Data Loading
# =============================================================================
cat("  Data loaded inline in render logic below.\n")

# Load metabolomics matrix for Fig 3d
mb_mat <- tryCatch({
  read.csv(file.path(DATA, "metabolomics_log2_merged.csv"), row.names = 1, check.names = FALSE)
}, error = function(e) {
  cat("  Warning: Could not load metabolomics data\n")
  NULL
})



# --- Panel object registry (for optimized assembly) ---
# Strategy A: ggplot objects kept directly (best quality - vector)
# Strategy B: ComplexHeatmap captured via grid.grabExpr (single rasterization)
# Strategy C: External/pre-rendered PDF read at 600 DPI (fallback)
obj_A <- obj_B <- obj_C <- obj_D <- obj_E <- obj_F <- NULL

# Helper: standalone horizontal barplot for group-difference panels
make_delta_bar <- function(names_vec, values_vec, title,
                           x_lab = "Adjacent \u2212 Normal (\u0394Z-score)",
                           name_size = 7.5, italic = FALSE) {
  df <- data.frame(
    name  = factor(names_vec, levels = rev(names_vec)),
    value = values_vec,
    dir   = ifelse(values_vec > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )
  face_val <- if (italic) "italic" else "plain"
  ggplot(df, aes(x = value, y = name, fill = dir)) +
    geom_col(width = 0.74, alpha = 0.92) +
    scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
    labs(title = title, x = x_lab, y = NULL) +
    theme_classic(base_size = 8, base_family = FONT_GRID) +
    theme(
      plot.title         = element_text(face = "bold", size = FS_TITLE,
                                        hjust = 0.5, family = FONT_FAMILY,
                                        margin = margin(0, 0, 3, 0)),
      axis.text.y        = element_text(size = name_size, face = face_val,
                                        family = FONT_FAMILY, color = "black"),
      axis.text.x        = element_text(size = 7, family = FONT_FAMILY, color = "black"),
      axis.title.x       = element_text(size = 8, family = FONT_FAMILY, color = "black"),
      axis.ticks.y       = element_blank(),
      axis.line.y        = element_blank(),
      axis.line.x        = element_line(linewidth = 0.4, color = "black"),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      plot.margin        = margin(4, 5, 3, 4, "mm")
    )
}

# =============================================================================
# SECTION 10: Render Figure 2 -- Metabolic Reprogramming
# =============================================================================
cat("\n--- Rendering Figure 2: Metabolic Reprogramming ---\n")

cat("\n=== FIGURE 2: Metabolic Reprogramming ===\n")

# Fig 2a: Metabolic pathway bubble
cat("  Fig2a\n")
tryCatch({
  mpw <- read.csv(file.path(RES, "phase3_metabolic/metabolic_pathway_summary.csv"), stringsAsFactors = FALSE)
  mpw$pw_name <- clean_pathway(mpw$pathway)
  mpw$n_val <- mpw$n_genes_defined
  p3a <- ggplot(mpw, aes(x = TC_mean_logFC, y = reorder(pw_name, TC_mean_logFC))) +
    geom_segment(aes(xend = 0, yend = reorder(pw_name, TC_mean_logFC)),
                 color = "grey75", linewidth = 0.45) +
    geom_point(aes(size = n_val, color = PR_mean_logFC)) +
    scale_color_gradient2(low = COL_DOWN, mid = "grey92", high = COL_UP,
                          midpoint = 0, name = "PR logFC") +
    scale_size_continuous(range = c(2, 7), name = "Genes") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey55") +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
    labs(title = "Metabolic pathway dysregulation",
         x = "TC mean log\u2082 FC", y = NULL) +
    theme_nc +
    theme(axis.text.y         = element_text(size = 8, family = FONT_FAMILY,
                                             color = "black"),
          plot.title          = element_text(hjust = 0.5, margin = margin(0, 0, 3, 0)),
          plot.title.position = "plot",
          panel.grid.major.x  = element_line(color = "grey94", linewidth = 0.3),
          plot.margin         = margin(4, 5, 3, 4, "mm"))
  save_pdf("Fig2a_metabolic_bubble.pdf", 113, 80, print(p3a))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# Fig 2b: Enzyme mean logFC barplot (top10 UP + top10 DOWN)
cat("  Fig2b\n")
tryCatch({
  enzyme <- read.csv(file.path(RES, "phase3_metabolic/metabolic_enzyme_expression.csv"), stringsAsFactors = FALSE)
  obj_A <- p3a  # Strategy A: keep ggplot object for direct assembly
  enz_tc <- aggregate(logFC ~ gene_name, data = enzyme[enzyme$omics == "Transcriptomics", ], FUN = mean)
  enz_pr <- aggregate(logFC ~ gene_name, data = enzyme[enzyme$omics == "Proteomics", ], FUN = mean)
  colnames(enz_tc)[2] <- "TC"; colnames(enz_pr)[2] <- "PR"
  enz_wide <- merge(enz_tc, enz_pr, by = "gene_name", all = TRUE)
  rownames(enz_wide) <- make.unique(enz_wide$gene_name)
  enz_mat  <- as.matrix(enz_wide[, c("TC", "PR")])
  enz_mat  <- enz_mat[complete.cases(enz_mat), ]
  enz_mean <- rowMeans(enz_mat)
  # Top10 UP + Top10 DOWN by |mean logFC|
  up_idx   <- order(enz_mean[enz_mean > 0], decreasing = TRUE)
  dn_idx   <- order(enz_mean[enz_mean < 0], decreasing = FALSE)
  top_up   <- head(names(sort(enz_mean[enz_mean > 0], decreasing = TRUE)), 10)
  top_dn   <- head(names(sort(enz_mean[enz_mean < 0], decreasing = FALSE)),  10)
  sel_nm   <- c(top_up, top_dn)
  sel_val  <- enz_mean[sel_nm]
  # Sort combined set: highest first → displayed at chart top
  ord      <- order(sel_val, decreasing = TRUE)
  sel_nm   <- sel_nm[ord];  sel_val <- sel_val[ord]
  p_b <- make_delta_bar(sel_nm, sel_val,
    title  = "Enzyme logFC",
    x_lab  = "Mean logFC (TC + PR)",
    name_size = 7.5, italic = TRUE)
  obj_B <- p_b
  save_pdf("Fig2b_enzyme_heatmap.pdf", 62, 78, print(p_b))
  cat("    -> saved\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# Fig 2c: Top DEMs heatmap
cat("  Fig2c\n")
tryCatch({
  mt_expr <- mb_mat
  mt_sig <- read.csv(file.path(RES, "phase1_diff/DEMs_significant.csv"), stringsAsFactors = FALSE)
  id_col <- find_col(mt_sig, c("Compound_ID", "name", "gene_name", "X", "metabolite", "ID", "metabolite_name"))
  fc_col <- find_col(mt_sig, c("logFC", "log2FoldChange"))
  top_ids <- head(mt_sig[order(abs(mt_sig[[fc_col]]), decreasing = TRUE), ], 20)
  matched <- top_ids[[id_col]][top_ids[[id_col]] %in% rownames(mt_expr)]
  mat <- scale_rows(mt_expr[matched, ], 3)
  nm_col <- find_col(mt_sig, c("metabolite_name", "name"))
  if (!is.null(nm_col)) {
    nm_map <- setNames(mt_sig[[nm_col]], mt_sig[[id_col]])
    new_rn <- nm_map[rownames(mat)]
    new_rn[is.na(new_rn) | new_rn == ""] <- rownames(mat)[is.na(new_rn) | new_rn == ""]
    new_rn <- abbreviate_metabolite(new_rn)
    new_rn <- format_metabolite_display(new_rn, compound_ids = rownames(mat), width = 35)
    rownames(mat) <- make.unique(new_rn)
  }
  # Adj − Normal delta-z, sorted descending
  is_adj_c <- make_group_vector(colnames(mat)) == "Adjacent"
  delta_c   <- rowMeans(mat[,  is_adj_c, drop=FALSE]) -
               rowMeans(mat[, !is_adj_c, drop=FALSE])
  ord_c     <- order(delta_c, decreasing = TRUE)
  delta_c   <- delta_c[ord_c]
  nm_c      <- rownames(mat)[ord_c]
  p_c <- make_delta_bar(nm_c, delta_c,
    title = "Top DEMs \u0394Z-score",
    x_lab = "Adjacent \u2212 Normal (\u0394Z-score)",
    name_size = 7)
  obj_C <- p_c
  save_pdf("Fig2c_top_DEMs_heatmap.pdf", 88, 82, print(p_c))
  cat("    -> saved\n")
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# =============================================================================
# SECTION 10d-f: Panels D, E, F (previously pre-rendered, now rendered in-script)
# =============================================================================

# --- Panel D: MSEA Enrichment Dotplot ---
cat("\n  Fig2d: MSEA enrichment dotplot\n")
tryCatch({
  msea_df <- read.csv(file.path(RES, "auxiliary_metabolomics_enrichment/msea_results.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  msea_top <- msea_df %>%
    filter(!is.na(NES)) %>%
    arrange(pval) %>%
    head(20) %>%
    mutate(
      pw_clean = gsub("_", " ", pathway_name),
      pw_clean = tools::toTitleCase(tolower(pw_clean)),
      neg_log10p = -log10(pval + 1e-300)
    )

  p2d <- ggplot(msea_top, aes(x = NES, y = reorder(pw_clean, NES))) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey55") +
    geom_point(aes(size = size, color = neg_log10p),
               alpha = 0.92, stroke = 0.25) +
    scale_color_gradient(low = "#FCD7D2", high = COL_UP,
                         name = "\u2013log\u2081\u2080 P") +
    scale_size_continuous(range = c(2, 6.5), name = "Set Size") +
    labs(title = "MSEA Enrichment (KEGG)", x = "NES", y = NULL) +
    theme_nc +
    theme(axis.text.y         = element_text(size = 8, family = FONT_FAMILY,
                                             color = "black"),
          plot.title          = element_text(hjust = 0.5, margin = margin(0, 0, 3, 0)),
          plot.title.position = "plot",
          panel.grid.major.x  = element_line(color = "grey94", linewidth = 0.3),
          plot.margin         = margin(4, 5, 3, 4, "mm"))

  save_pdf("Fig2d_msea_dotplot.pdf", 90, 88, print(p2d))
  cat("    -> Fig2d_msea_dotplot.pdf\n")
}, error = function(e) cat("    ERROR Panel D:", e$message, "\n"))

# --- Panel E: Protein-Metabolite Network ---
cat("\n  Fig2e: Protein-metabolite network\n")
tryCatch({
  net_edges <- read.csv(file.path(RES, "enhancement4_prot_metab_network/prot_metab_edges.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  hub_nodes <- read.csv(file.path(RES, "enhancement4_prot_metab_network/network_hub_nodes.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)

  # Filter to top edges by absolute correlation
  net_top <- net_edges %>%
    mutate(abs_rho = abs(rho)) %>%
    arrange(desc(abs_rho)) %>%
    head(150)

  # Build edge data for igraph
  el <- data.frame(from = net_top$protein, to = net_top$metabolite_name_full,
                   stringsAsFactors = FALSE)
  g <- graph_from_data_frame(el, directed = FALSE)

  # Node type
  all_proteins <- unique(net_top$protein)
  V(g)$type <- ifelse(V(g)$name %in% all_proteins, "Protein", "Metabolite")
  V(g)$color <- ifelse(V(g)$type == "Protein", COL_PR, COL_MT)

  # Node size by degree
  deg <- igraph::degree(g)
  V(g)$size <- sqrt(deg) * 2.5 + 2

  # Labels: top hub nodes (increased density for richer panel)
  top_hubs <- head(hub_nodes$node[hub_nodes$node %in% V(g)$name], 14)
  V(g)$label <- ifelse(V(g)$name %in% top_hubs, V(g)$name, "")
  # Abbreviate metabolite labels and apply display formatting
  is_metab <- V(g)$name %in% top_hubs & V(g)$type == "Metabolite"
  metab_labels <- V(g)$label[is_metab]
  metab_labels <- abbreviate_metabolite(metab_labels)
  # Map metabolite names to compound IDs for >30 char fallback
  metab_id_map <- setNames(net_top$metabolite_id, net_top$metabolite_name_full)
  metab_ids <- metab_id_map[V(g)$name[is_metab]]
  V(g)$label[is_metab] <- format_metabolite_display(metab_labels, compound_ids = metab_ids)

  # Edge color by direction
  E(g)$color <- ifelse(net_top$direction == "positive",
                        adjustcolor(COL_UP, alpha.f = 0.4),
                        adjustcolor(COL_DOWN, alpha.f = 0.4))
  E(g)$width <- abs(net_top$rho) * 1.5

  # --- Convert igraph to ggplot2 + ggrepel for clear labels with leader lines ---
  set.seed(42)
  lo <- layout_with_fr(g)
  
  # Build node data.frame
  node_df <- data.frame(
    x = lo[, 1], y = lo[, 2],
    name = V(g)$name,
    type = V(g)$type,
    size = V(g)$size,
    label = V(g)$label,
    stringsAsFactors = FALSE
  )
  node_df$show_label <- nchar(node_df$label) > 0
  
  # Build edge data.frame
  el_idx <- igraph::as_edgelist(g, names = FALSE)
  edge_df <- data.frame(
    x = lo[el_idx[, 1], 1], y = lo[el_idx[, 1], 2],
    xend = lo[el_idx[, 2], 1], yend = lo[el_idx[, 2], 2],
    direction = net_top$direction[1:nrow(el_idx)],
    rho = abs(net_top$rho[1:nrow(el_idx)]),
    stringsAsFactors = FALSE
  )
  
  # ggplot network — legend placed BELOW the plot panel to avoid overlap with nodes
  p2e <- ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend, color = direction),
                 alpha = 0.3, linewidth = edge_df$rho * 0.8, show.legend = TRUE) +
    scale_color_manual(name = NULL,
                       values = c("positive" = COL_UP, "negative" = COL_DOWN),
                       labels = c("positive" = "Positive corr.",
                                  "negative" = "Negative corr.")) +
    ggnewscale::new_scale_color() +
    geom_point(data = node_df,
               aes(x = x, y = y, fill = type, size = size),
               shape = 21, color = "grey40", stroke = 0.3, show.legend = TRUE) +
    scale_fill_manual(name = NULL,
                      values = c("Protein" = COL_PR, "Metabolite" = COL_MT)) +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    ggrepel::geom_label_repel(
      data = node_df[node_df$show_label, ],
      aes(x = x, y = y, label = label),
      size = 2.5, fontface = "bold", family = FONT_FAMILY,
      box.padding = 0.5, point.padding = 0.3,
      segment.color = "grey30", segment.size = 0.4,
      min.segment.length = 0.2, max.overlaps = 20,
      fill = alpha("white", 0.85), label.size = 0.2,
      force = 8, force_pull = 0.5, seed = 42
    ) +
    labs(title = "Protein-Metabolite Network") +
    coord_cartesian(clip = "off") +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title           = element_text(size = 9, face = "bold", hjust = 0.5,
                                          family = FONT_FAMILY,
                                          margin = margin(0, 0, 2, 0)),
      legend.position      = "bottom",
      legend.box           = "vertical",
      legend.box.just      = "center",
      legend.text          = element_text(size = 7, family = FONT_FAMILY),
      legend.title         = element_blank(),
      legend.key.size      = unit(2.8, "mm"),
      legend.spacing.x     = unit(2, "mm"),
      legend.spacing.y     = unit(0.5, "mm"),
      legend.box.spacing   = unit(0.5, "mm"),
      legend.background    = element_blank(),
      legend.margin        = margin(0, 1, 0, 1),
      plot.margin          = margin(2, 2, 1, 2, "mm")
    ) +
    guides(fill  = guide_legend(order = 1, nrow = 1,
                                override.aes = list(size = 3)),
           color = guide_legend(order = 2, nrow = 1,
                                override.aes = list(linewidth = 1, alpha = 0.8)))
  
  save_pdf("Fig2e_prot_metab_network.pdf", 80, 81, print(p2e))
  # Store ggplot object for vector assembly (no grob_E needed)
  grob_E <- NULL  # flag: use p2e directly
  cat("    -> Fig2e_prot_metab_network.pdf\n")
}, error = function(e) cat("    ERROR Panel E:", e$message, "\n"))

# --- Panel F: Leading-edge metabolites of top 2 enriched pathways ---
# Replaces the previous redundant MSEA barplot. Now shows pathway-level
# leading-edge metabolite contributions, complementing Panel D.
cat("\n  Fig2f: Leading-edge metabolites barplot\n")
tryCatch({
  msea_df2 <- read.csv(file.path(RES, "auxiliary_metabolomics_enrichment/msea_results.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  dem_anno <- read.csv(file.path(RES, "phase1_diff/DEMs_significant.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  # Full metabolite annotation (covers all 3300+ compounds, not just DEMs)
  full_anno <- tryCatch(
    read.csv(file.path(DATA, "metabolomics_annotation.csv"),
             stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )

  # Choose top-2 pathways by p-value (these are the strongest signals from Panel D)
  msea_top2 <- msea_df2 %>%
    filter(!is.na(NES)) %>%
    arrange(pval) %>%
    head(2)

  # Build name lookup: prefer full annotation, fall back to DEMs annotation
  name_map <- setNames(dem_anno$metabolite_name, dem_anno$Compound_ID)
  if (!is.null(full_anno) && all(c("Compound_ID", "Name") %in% colnames(full_anno))) {
    full_map <- setNames(full_anno$Name, full_anno$Compound_ID)
    # Full annotation provides primary names; DEMs map fills any gaps
    combined_map <- full_map
    missing_in_full <- setdiff(names(name_map), names(combined_map))
    combined_map[missing_in_full] <- name_map[missing_in_full]
    name_map <- combined_map
  }

  # Sample-level Adjacent - Normal log2 fold change
  is_adj_f <- make_group_vector(colnames(mb_mat)) == "Adjacent"

  rows <- list()
  for (i in seq_len(nrow(msea_top2))) {
    le_ids <- strsplit(msea_top2$leadingEdge_str[i], ";")[[1]]
    le_ids <- le_ids[le_ids %in% rownames(mb_mat)]
    if (length(le_ids) == 0) next
    sub_mat <- mb_mat[le_ids, , drop = FALSE]
    delta <- rowMeans(sub_mat[,  is_adj_f, drop = FALSE]) -
             rowMeans(sub_mat[, !is_adj_f, drop = FALSE])
    # Keep only top 8 by |delta| per pathway for visual clarity
    if (length(delta) > 8) {
      keep_idx <- order(abs(delta), decreasing = TRUE)[1:8]
      delta <- delta[keep_idx]
      le_ids <- le_ids[keep_idx]
    }
    nms <- name_map[le_ids]
    nms[is.na(nms) | nms == ""] <- le_ids[is.na(nms) | nms == ""]
    nms <- abbreviate_metabolite(nms)
    nms <- format_metabolite_display(nms, compound_ids = le_ids,
                                     width = 24, max_len = 42)
    pw_label <- tools::toTitleCase(tolower(msea_top2$pathway_name[i]))
    # Display pathway name in full — Panel F strip is wide enough (103mm)
    rows[[i]] <- data.frame(
      pathway = pw_label,
      metabolite = nms,
      delta = unname(delta),
      compound_id = le_ids,
      stringsAsFactors = FALSE
    )
  }
  le_df <- do.call(rbind, rows)
  # Order metabolites within each pathway by delta (descending)
  le_df <- le_df %>%
    group_by(pathway) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    ungroup() %>%
    as.data.frame()
  le_df$dir <- ifelse(le_df$delta > 0, "Up", "Down")
  # Build a unique factor ordering across the two facets
  le_df$metab_id <- paste0(le_df$pathway, "__", le_df$metabolite)
  le_df$metab_id <- factor(le_df$metab_id, levels = rev(le_df$metab_id))
  # Preserve pathway order from msea_top2 (top1 first)
  le_df$pathway <- factor(le_df$pathway, levels = unique(le_df$pathway))

  p2f <- ggplot(le_df, aes(x = delta, y = metab_id, fill = dir)) +
    geom_col(width = 0.72, alpha = 0.92) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
    facet_wrap(~ pathway, ncol = 1, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
    scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), guide = "none") +
    labs(title = "Leading-edge metabolites",
         x = "Adjacent \u2212 Normal (log\u2082 abundance)", y = NULL) +
    theme_classic(base_size = 8, base_family = FONT_GRID) +
    theme(
      plot.title          = element_text(face = "bold", size = FS_TITLE,
                                         hjust = 0.5, family = FONT_FAMILY),
      axis.text.y         = element_text(size = 7, family = FONT_FAMILY),
      axis.text.x         = element_text(size = 7, family = FONT_FAMILY),
      axis.title.x        = element_text(size = 8, family = FONT_FAMILY),
      axis.ticks.y        = element_blank(),
      panel.grid.major.x  = element_line(colour = "grey90", linewidth = 0.3),
      strip.background    = element_rect(fill = "grey94", color = "grey80",
                                         linewidth = 0.4),
      strip.text          = element_text(size = 7.5, face = "bold",
                                         family = FONT_FAMILY,
                                         margin = margin(2, 2, 2, 2)),
      plot.margin         = margin(4, 4, 2, 4, "mm")
    )

  save_pdf("Fig2f_GSEA_metab.pdf", 103, 77, print(p2f))
  cat("    -> Fig2f_GSEA_metab.pdf (leading-edge metabolites)\n")
}, error = function(e) cat("    ERROR Panel F:", e$message, "\n"))

# =============================================================================
# SECTION 11: Composite Figure Assembly (Pure Vector Grid Viewport — AI-editable)
# Layout: Row1(82mm): A(91)+C(92); Row2(82mm): B(70)+D(113); Row3(81mm): E(80)+F(103)
# =============================================================================
cat("\n--- Assembling composite Figure_2 (VECTOR, 183x245mm) ---\n")

tryCatch({
  DPI <- 600
  W_TOTAL <- 183; H_TOTAL <- 245
  # New layout: Row1=A+C, Row2=B+D, Row3=E+F
  W_A <- 91;  W_C <- W_TOTAL - W_A   # 92mm
  W_B <- 70;  W_D <- W_TOTAL - W_B   # 113mm
  W_E <- 80;  W_F <- W_TOTAL - W_E   # 103mm
  H_R1 <- 82; H_R2 <- 82; H_R3 <- H_TOTAL - H_R1 - H_R2  # 81mm

  # --- Render function: places all panels + labels into current device ---
  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))

    # Row 1 (top): A + C — y = H_R2 + H_R3 from bottom
    y_row1 <- H_R2 + H_R3

    # Panel A (ggplot - metabolic bubble)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row1, "mm"),
      width = unit(W_A, "mm"), height = unit(H_R1, "mm"),
      just = c("left", "bottom")
    ))
    print(obj_A, newpage = FALSE)
    grid::popViewport()

    # Panel C (ggplot - DEMs delta barplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_A, "mm"), y = unit(y_row1, "mm"),
      width = unit(W_C, "mm"), height = unit(H_R1, "mm"),
      just = c("left", "bottom")
    ))
    print(p_c, newpage = FALSE)
    grid::popViewport()

    # Row 2 (middle): B + D — y = H_R3 from bottom
    y_row2 <- H_R3

    # Panel B (ggplot - enzyme logFC barplot)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_B, "mm"), height = unit(H_R2, "mm"),
      just = c("left", "bottom")
    ))
    print(p_b, newpage = FALSE)
    grid::popViewport()

    # Panel D (ggplot - MSEA dotplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_B, "mm"), y = unit(y_row2, "mm"),
      width = unit(W_D, "mm"), height = unit(H_R2, "mm"),
      just = c("left", "bottom")
    ))
    print(p2d, newpage = FALSE)
    grid::popViewport()

    # Row 3 (bottom): E + F — y = 0
    # Panel E (igraph network — grid grob)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W_E, "mm"), height = unit(H_R3, "mm"),
      just = c("left", "bottom")
    ))
    print(p2e, newpage = FALSE)
    grid::popViewport()

    # Panel F (ggplot - MSEA barplot)
    grid::pushViewport(grid::viewport(
      x = unit(W_E, "mm"), y = unit(0, "mm"),
      width = unit(W_F, "mm"), height = unit(H_R3, "mm"),
      just = c("left", "bottom")
    ))
    print(p2f, newpage = FALSE)
    grid::popViewport()

    # --- Panel labels (bold, FS_TAG pt) ---
    label_data <- data.frame(
      text = c("A", "B", "C", "D", "E", "F"),
      x_mm = c(1, 1, W_A + 1, W_B + 1, 1, W_E + 1),
      y_mm = c(H_TOTAL - 1, y_row1 - 1, H_TOTAL - 1, y_row1 - 1, H_R3 - 1, H_R3 - 1),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(label_data))) {
      grid::grid.text(
        label = label_data$text[i],
        x = unit(label_data$x_mm[i], "mm"),
        y = unit(label_data$y_mm[i], "mm"),
        just = c("left", "top"),
        gp = gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY)
      )
    }

    grid::popViewport()  # pop full-page viewport
  }

  # --- Save vector PDF (AI-editable) ---
  cairo_pdf(file.path(OUT, "Figure_2.pdf"), width = W_TOTAL / 25.4,
            height = H_TOTAL / 25.4, family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> Figure_2.pdf (VECTOR, AI-editable)\n")

  # --- Save PNG (600 DPI for review) ---
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "Figure_2.png"), width = px_W, height = px_H,
                 res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_2.png\n")

  # --- Save TIFF (600 DPI for journal submission) ---
  grDevices::tiff(file.path(OUT, "Figure_2.tiff"), width = px_W, height = px_H,
                  res = DPI, compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_2.tiff\n")

  cat("  Figure_2 DONE (183x245mm, VECTOR PDF + 600DPI PNG/TIFF)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))
# END
# =============================================================================
cat("\n=== Figure 2 rendering complete ===\n")
