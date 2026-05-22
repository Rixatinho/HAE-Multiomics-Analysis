#!/usr/bin/env Rscript
# =============================================================================
# nc_theme.R -- Single Source of Truth for Nature Communications Styling
# =============================================================================
# Design reference: ggsci NPG palette + theme_classic() (L-shaped axes)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(survival)
  library(survminer)
})

# =============================================================================
# 1. FONT CONSTANTS
# =============================================================================
FONT_FAMILY <- "Helvetica"      # AI native font - better than Arial for AI editing
FONT_GRID   <- "Helvetica"      # unified for ComplexHeatmap
FONT        <- "Helvetica"      # single constant for assembly labels
MM_PER_INCH <- 25.4

# Force ComplexHeatmap global font family to Arial
ht_opt$TITLE_PADDING    <- unit(c(4, 4), "pt")
ht_opt$message          <- FALSE

# =============================================================================
# 2. NATURE COMMUNICATIONS DIMENSION CONSTANTS (mm)
# =============================================================================
W_SINGLE  <- 89
W_DOUBLE  <- 183
W_HALF    <- 88
H_STD     <- 85
H_TALL    <- 100
H_MAX     <- 240

mm2in <- function(mm) mm / MM_PER_INCH
in2mm <- function(inch) inch * MM_PER_INCH

# =============================================================================
# 3. COLOR PALETTE (ggsci NPG / Nature Publishing Group)
# =============================================================================

# NPG categorical palette (Nature standard)
PAL_CAT <- c(
  "#E64B35",  # red (NPG)
  "#4DBBD5",  # cyan (NPG)
  "#00A087",  # teal (NPG)
  "#3C5488",  # navy (NPG)
  "#F39B7F",  # salmon (NPG)
  "#8491B4",  # lavender (NPG)
  "#91D1C2",  # mint (NPG)
  "#B09C85"   # tan (NPG)
)

# Semantic colors -- NPG-derived
COL_UP       <- "#E64B35"    # up-regulated (NPG red)
COL_DOWN     <- "#4DBBD5"    # down-regulated (NPG cyan)
COL_NS       <- "#999999"    # not significant (WCAG 2.0 AA compliant)
COL_NA       <- "#F0F0F0"    # NA / background

# Omics layer colors (NPG-derived, high contrast)
COL_TC       <- "#3C5488"    # transcriptomics (navy)
COL_PR       <- "#E64B35"    # proteomics (red)
COL_MT       <- "#00A087"    # metabolomics (teal)

# Group colors
COL_NORMAL   <- "#4DBBD5"    # normal tissue (cyan)
COL_ADJACENT <- "#E64B35"    # adjacent/disease tissue (red)

# Subtype colors
COL_CS1      <- "#E64B35"    # consensus subtype 1 (red)
COL_CS2      <- "#3C5488"    # consensus subtype 2 (navy)

# =============================================================================
# 4. COLOR SCALE FUNCTIONS (circlize::colorRamp2)
# =============================================================================

col_div     <- colorRamp2(c(-2, 0, 2),     c("#3C5488", "#FFFFFF", "#E64B35"))
col_nes     <- colorRamp2(c(-3, 0, 3),     c("#3C5488", "#FFFFFF", "#E64B35"))
col_cor     <- colorRamp2(c(-1, 0, 1),     c("#3C5488", "#FFFFFF", "#E64B35"))
col_zscore  <- colorRamp2(c(-2.5, 0, 2.5), c("#3C5488", "#FFFFFF", "#E64B35"))
col_immune  <- colorRamp2(c(-2, 0, 2),     c("#00A087", "#FFFFFF", "#8491B4"))

# =============================================================================
# 5. UNIFIED ggplot2 THEME (theme_classic base = L-shaped axes)
# =============================================================================

theme_pub <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    axis.line          = element_line(linewidth = 0.5, color = "black"),
    axis.ticks         = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length  = unit(1.5, "mm"),
    axis.text          = element_text(size = 8, family = FONT_FAMILY, color = "black"),
    axis.title         = element_text(size = 9, family = FONT_FAMILY, face = "bold", color = "black"),
    plot.title         = element_text(size = 10, family = FONT_FAMILY, face = "bold", hjust = 0),
    legend.text        = element_text(size = 7, family = FONT_FAMILY),
    legend.title       = element_text(size = 8, family = FONT_FAMILY, face = "bold"),
    legend.key.size    = unit(2.5, "mm"),
    legend.background  = element_blank(),
    legend.box.background = element_blank(),
    legend.spacing.y   = unit(1, "mm"),
    strip.text         = element_text(size = 9, family = FONT_FAMILY, face = "bold"),
    strip.background   = element_rect(fill = "grey96", linewidth = 0.3, color = "grey80"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    plot.margin        = margin(5, 5, 5, 5, "mm")
  )
theme_set(theme_pub)

# =============================================================================
# 6. ComplexHeatmap gpar FACTORY FUNCTIONS
# =============================================================================

gp_row_names <- function(size = 7, italic = FALSE) {

  gpar(fontsize = size, fontfamily = FONT_GRID,
       fontface = if (italic) "italic" else "plain")
}

gp_col_names <- function(size = 7, bold = TRUE) {
  gpar(fontsize = size, fontfamily = FONT_GRID,
       fontface = if (bold) "bold" else "plain")
}

gp_legend_title <- function(size = 8) {
  gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
}

gp_legend_labels <- function(size = 7) {
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

# Heat map cell border (white grid lines for NC standard)
gp_rect <- function() {
  gpar(col = "#FFFFFF", lwd = 0.5)  # Explicit hex code for AI compatibility
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

# Standard heatmap parameters for Nature Communications
std_heatmap_param <- function() {
  list(
    rect_gp              = gp_rect(),
    row_names_gp         = gp_row_names(),
    column_names_gp      = gp_col_names(),
    column_title_gp      = gpar(fontsize = 9, fontfamily = FONT_GRID, fontface = "bold"),
    row_title_gp         = gp_row_title(),
    show_row_dend        = FALSE,
    show_column_dend     = FALSE,
    use_raster           = FALSE,    # AI compatibility: vector cells for editing
    heatmap_legend_param = std_legend_param()
  )
}

# Panel label gpar (for assembly scripts)
gp_panel_label <- function(size = 10) {
  gpar(fontsize = size, fontfamily = FONT_FAMILY, fontface = "bold")
}

# =============================================================================
# 6.5 AI-COMPATIBLE PDF OUTPUT HELPER
# =============================================================================

# AI-compatible PDF output helper (for Adobe Illustrator editing)
ai_pdf <- function(filename, width_mm, height_mm) {
  w_in <- width_mm / 25.4
  h_in <- height_mm / 25.4
  cairo_pdf(filename, width = w_in, height = h_in, 
            family = "Helvetica")
}

# =============================================================================
# 7. HELPER FUNCTIONS
# =============================================================================

# Map panel filename to per-figure subdirectory (e.g. "Fig3a_*.pdf" -> "Figure_3")
filename_to_figdir <- function(filename) {
  if (grepl("^Fig(\\d+)", filename))
    return(paste0("Figure_", sub("^Fig(\\d+).*", "\\1", filename)))
  if (grepl("^Supp(\\d+)", filename))
    return(sprintf("SuppFig_%02d", as.integer(sub("^Supp(\\d+).*", "\\1", filename))))
  return("")
}

sp <- function(filename, w_mm, h_mm, out_dir = NULL) {
  if (is.null(out_dir)) {
    subdir <- filename_to_figdir(filename)
    out_dir <- file.path(OUT, subdir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  list(path = file.path(out_dir, filename),
       width = w_mm / MM_PER_INCH, height = h_mm / MM_PER_INCH)
}

save_cairo <- function(filename, w_mm, h_mm, expr, out_dir = NULL) {
  s <- sp(filename, w_mm, h_mm, out_dir)
  cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR saving ", filename, ": ", e$message))
  dev.off()
  cat("  ->", basename(s$path), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
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

# --- Standardized Figure 1 clinical annotation bar ---
# @param subtypes: vector of subtype assignments (CS1/CS2)
# @param clinical_df: data.frame with patient clinical info
# @param patient_ids: vector of patient IDs to subset clinical_df
get_nc_annotation <- function(subtypes, clinical_df, patient_ids = NULL) {
  # Subset clinical data if patient_ids provided
  if (!is.null(patient_ids)) {
    clin_sub <- clinical_df[match(patient_ids, clinical_df$patient_id), ]
  } else {
    clin_sub <- clinical_df
  }
  
  # Derive annotation vectors
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
      Sex = c("Male" = "#3C5488", "Female" = "#F39B7F"),
      Bilirubin = c("High" = "#FFD700", "Normal" = COL_NA),
      `Bile Duct` = c("Invaded" = "#7E6148", "No" = COL_NA)
    ),
    simple_anno_size = unit(4, "mm"),
    gap = unit(1, "mm"),
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
  )
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else NULL
}

# --- Clean pathway / gene-set names for display ---
clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_|^REACTOME_|^KEGG_|^GOBP_|^GO_|^WP_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x
}

# --- Clean clinical variable names ---
clean_varname <- function(x) {
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(x)
  x <- gsub(" Pct$", " (%)", x)
  x <- gsub("^Wbc$", "WBC", x); x <- gsub("^Alt$", "ALT", x)
  x <- gsub("^Ast$", "AST", x); x <- gsub("^Ggt$", "GGT", x)
  x
}

# --- Map ENSG IDs to gene symbols using org.Hs.eg.db ---
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
  ids_clean <- gsub("\\.[0-9]+$", "", ids)       # strip version
  ids_clean <- gsub("_[0-9]+$", "", ids_clean)    # strip suffix like _3
  idx <- match(ids_clean, m$ensg)
  out <- ifelse(is.na(idx), ids, m$symbol[idx])
  out
}

# --- Map ENSP protein IDs to gene symbols ---
.ensp_map <- NULL
build_ensp_map <- function() {
  if (!is.null(.ensp_map)) return(.ensp_map)
  map <- data.frame(ensp = character(), symbol = character(), stringsAsFactors = FALSE)

  # Source 1: DEPs file (has Protein ENSP and gene_name)
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

  # Source 2: org.Hs.eg.db ENSEMBLPROT
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
  ids_clean <- gsub("\\.[0-9]+$", "", ids)       # strip version
  ids_clean <- gsub("_[0-9]+$", "", ids_clean)    # strip suffix like _3
  idx <- match(ids_clean, m$ensp)
  out <- ifelse(is.na(idx), ids, m$symbol[idx])
  out
}

# --- Format feature names with ID mapping and truncation ---
# @param ids: vector of feature IDs to map
# @param mapping_df: data.frame with 'id' and 'name' columns
# @param type: "MT" for metabolites (truncate long names), "TC"/"PR" for genes/proteins (mark ENS* as NA)
# @param max_len: maximum length before truncation (for MT type)
format_feature_names <- function(ids, mapping_df, type = "MT", max_len = 30) {
  # Perform ID -> name mapping
  new_names <- mapping_df$name[match(ids, mapping_df$id)]
  
  # Fallback: unmapped IDs keep original value
  new_names <- ifelse(is.na(new_names) | new_names == "", ids, new_names)
  
  if (type == "MT") {
    # Metabolites: truncate long names
    new_names <- ifelse(nchar(new_names) > max_len,
                        paste0(substr(new_names, 1, max_len - 3), "..."),
                        new_names)
  } else {
    # Transcriptomics/Proteomics: mark unmapped ENS* IDs as NA
    # Caller can filter out NA entries as needed
    unmapped_idx <- grepl("^ENS[GP]", new_names)
    new_names[unmapped_idx] <- NA
  }
  
  return(new_names)
}

# --- Clean cell type / signature names ---
clean_celltype <- function(x) {
  x <- gsub("_", " ", x)
  x
}

# --- Clean metabolite / pathway names with truncation ---
clean_label_name <- function(x, max_chars = 25) {
  x <- gsub("\\s*\\(.*?\\)$", "", x)  # 去除末尾括号内容
  x <- gsub("^(HALLMARK_|KEGG_|GO_|REACTOME_)", "", x)  # 去除数据库前缀
  x <- gsub("_", " ", x)  # 下划线转空格
  x <- tools::toTitleCase(tolower(x))  # 标题大小写
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

# =============================================================================
# 7.5 ENHANCED VISUALIZATION FUNCTIONS
# =============================================================================

# --- plot_nc_volcano() -- Enhanced Volcano Plot ---
plot_nc_volcano <- function(df, title, fc_col = "logFC", p_col = "P.Value",
                            symbol_col = "symbol", fc_cutoff = 1.0, 
                            p_cutoff = 0.05, top_n = 15) {
  # 确保列存在
  if (!symbol_col %in% colnames(df)) {
    # 尝试常见列名
    for (cand in c("gene_name", "metabolite_name", "name", "Name")) {
      if (cand %in% colnames(df)) { symbol_col <- cand; break }
    }
  }
  
  df$neg_log10p <- -log10(df[[p_col]])
  df$fc <- df[[fc_col]]
  
  # 分类标记
  df$direction <- "NS"
  df$direction[df$fc > fc_cutoff & df[[p_col]] < p_cutoff] <- "Up"
  df$direction[df$fc < -fc_cutoff & df[[p_col]] < p_cutoff] <- "Down"
  
  # 标注 top N 显著基因
  sig <- df[df$direction != "NS", ]
  sig <- sig[order(sig[[p_col]]), ]
  top_genes <- head(sig[[symbol_col]], top_n)
  df$label <- ifelse(df[[symbol_col]] %in% top_genes, df[[symbol_col]], "")
  
  # 过滤掉未映射的 ENS ID 标签
  df$label <- ifelse(grepl("^ENS[GP]", df$label), "", df$label)
  
  ggplot(df, aes(x = fc, y = neg_log10p)) +
    geom_point(aes(color = fc), size = 0.4, alpha = 0.4) +
    scale_color_gradient2(low = COL_DOWN, mid = "grey90", high = COL_UP, midpoint = 0,
                          guide = "none") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", 
               linewidth = 0.3, color = "grey40") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", 
               linewidth = 0.3, color = "grey40") +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5, segment.size = 0.2, 
                              fontface = "italic", family = FONT_FAMILY,
                              max.overlaps = 15, force = 1.5, box.padding = 0.3) +
    labs(title = title, x = expression(log[2]~"(Fold Change)"), 
         y = expression(-log[10]~"(P-value)")) +
    theme_pub
}

# --- plot_nc_enrichment() -- Enhanced Enrichment Bubble Plot ---
plot_nc_enrichment <- function(en_df, title, n = 15, desc_col = NULL, 
                                nes_col = NULL, pval_col = NULL, 
                                size_col = NULL) {
  # 自动探测列名（与现有 make_gsea_dotplot 类似的逻辑）
  if (is.null(desc_col)) {
    for (cand in c("Description", "pathway", "Term", "name")) {
      if (cand %in% colnames(en_df)) { desc_col <- cand; break }
    }
  }
  if (is.null(pval_col)) {
    for (cand in c("p.adjust", "padj", "pvalue", "P.Value", "FDR")) {
      if (cand %in% colnames(en_df)) { pval_col <- cand; break }
    }
  }
  if (is.null(size_col)) {
    for (cand in c("Count", "setSize", "size", "n_genes")) {
      if (cand %in% colnames(en_df)) { size_col <- cand; break }
    }
  }
  
  # 过滤（符合元数据 Fig3c 过滤逻辑）
  en_df <- en_df[!grepl("^map[0-9]|^hsa[0-9]", en_df[[desc_col]], ignore.case = TRUE), ]
  en_df <- en_df[nchar(en_df[[desc_col]]) > 5, ]
  en_df[[desc_col]] <- clean_pathway(en_df[[desc_col]])
  
  # 排序取 top n
  en_df <- en_df[order(en_df[[pval_col]]), ]
  en_df <- head(en_df, n)
  
  en_df$neg_log10p <- -log10(en_df[[pval_col]])
  
  p <- ggplot(en_df, aes(x = neg_log10p, y = reorder(en_df[[desc_col]], neg_log10p)))
  
  if (!is.null(size_col) && size_col %in% colnames(en_df)) {
    p <- p + geom_point(aes(size = .data[[size_col]], color = neg_log10p)) +
      scale_size_continuous(range = c(2, 5))
  } else {
    p <- p + geom_point(aes(color = neg_log10p), size = 2.5)
  }
  
  p + scale_color_gradient(low = "#F0C1B8", high = COL_UP) +
    labs(title = title, x = expression(-log[10]~"(P-adj)"), y = NULL) +
    theme_pub +
    theme(axis.text.y = element_text(size = 6.5))
}

# --- clean_pathway_nc() -- Enhanced pathway name cleaning ---
clean_pathway_nc <- function(df, desc_col = "Description") {
  # Auto-detect description column
  if (!desc_col %in% colnames(df)) {
    for (cand in c("Description", "pathway", "Term", "name", "pw_name")) {
      if (cand %in% colnames(df)) { desc_col <- cand; break }
    }
  }
  
  df %>%
    # Exclude rows with map/hsa numeric IDs
    dplyr::filter(!grepl("^map[0-9]|^hsa[0-9]", .data[[desc_col]], ignore.case = TRUE)) %>%
    # Exclude too short names
    dplyr::filter(nchar(as.character(.data[[desc_col]])) > 5) %>%
    # Format names: remove prefix -> underscore to space -> sentence case
    dplyr::mutate(!!desc_col := stringr::str_to_sentence(gsub("_", " ", .data[[desc_col]])))
}

# --- make_nc_gsea_barplot() -- Production-grade GSEA bar plot ---
make_nc_gsea_barplot <- function(gsea_res, title, n = 15, 
                                  desc_col = NULL, nes_col = NULL, pval_col = NULL) {
  # Auto-detect column names
  if (is.null(desc_col)) {
    for (cand in c("Description", "pathway", "Term", "name", "pw_name")) {
      if (cand %in% colnames(gsea_res)) { desc_col <- cand; break }
    }
  }
  if (is.null(nes_col)) {
    for (cand in c("NES", "enrichmentScore", "nes")) {
      if (cand %in% colnames(gsea_res)) { nes_col <- cand; break }
    }
  }
  if (is.null(pval_col)) {
    for (cand in c("pvalue", "p.adjust", "padj", "P.Value", "FDR")) {
      if (cand %in% colnames(gsea_res)) { pval_col <- cand; break }
    }
  }
  
  # Clean pathway names
  gsea_res[[desc_col]] <- clean_pathway(gsea_res[[desc_col]])
  
  # Filter invalid pathways
  gsea_res <- gsea_res[!grepl("^map[0-9]|^hsa[0-9]", gsea_res[[desc_col]], ignore.case = TRUE), ]
  gsea_res <- gsea_res[nchar(gsea_res[[desc_col]]) > 5, ]
  
  # Sort and take top n
  gsea_res <- gsea_res[order(-abs(gsea_res[[nes_col]])), ]
  df_clean <- head(gsea_res, n)
  
  # Auto-wrap pathway names longer than 40 characters
  df_clean[[desc_col]] <- stringr::str_wrap(df_clean[[desc_col]], width = 40)
  
  ggplot(df_clean, aes(x = .data[[nes_col]], 
                        y = reorder(.data[[desc_col]], .data[[nes_col]]))) +
    geom_col(aes(fill = .data[[nes_col]] > 0), width = 0.8) +
    scale_fill_manual(values = c("TRUE" = COL_UP, "FALSE" = COL_DOWN), guide = "none") +
    geom_text(aes(label = sprintf("p=%.3f", .data[[pval_col]]),
                  hjust = ifelse(.data[[nes_col]] > 0, -0.1, 1.1)), 
              size = 2.8, family = FONT_FAMILY) +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_pub +
    theme(axis.text.y = element_text(size = 8))
}

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

# --- plot_nc_immune_barplot() -- Enhanced immune cell boxplot with stats ---
plot_nc_immune_barplot <- function(immune_df, cell_col = NULL, value_col = NULL, 
                                    group_col = "Group", title = "Immune Infiltration") {
  # Auto-detect column names
  if (is.null(cell_col)) {
    for (cand in c("cell_type", "Cell_Type", "celltype", "name")) {
      if (cand %in% colnames(immune_df)) { cell_col <- cand; break }
    }
  }
  if (is.null(value_col)) {
    for (cand in c("diff", "score", "Proportion", "value", "enrichment")) {
      if (cand %in% colnames(immune_df)) { value_col <- cand; break }
    }
  }
  
  # Clean cell type names
  immune_df[[cell_col]] <- clean_celltype(immune_df[[cell_col]])
  
  ggplot(immune_df, aes(x = .data[[cell_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.3) +
    ggpubr::stat_compare_means(aes(group = .data[[group_col]]), 
                                label = "p.signif", size = 2.8, 
                                label.y = max(immune_df[[value_col]], na.rm = TRUE) * 1.1) +
    scale_fill_manual(values = c("Normal" = COL_DOWN, "Adjacent" = COL_UP)) +
    labs(x = NULL, y = "Estimated Proportion", title = title) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "top")
}

# --- plot_nc_checkpoint_heatmap() -- Immune checkpoint heatmap with clinical annotation ---
plot_nc_checkpoint_heatmap <- function(mat, meta, title = "Immune Checkpoints") {
  # Standardized Group annotation bar
  top_anno <- HeatmapAnnotation(
    Group = meta$Group,
    col = list(Group = c("Normal" = COL_DOWN, "Adjacent" = COL_UP)),
    simple_anno_size = unit(3, "mm"),
    show_annotation_name = FALSE,
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
  )
  
  Heatmap(mat, 
          name = "Z-score",
          top_annotation = top_anno,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          show_row_dend = FALSE,
          show_column_dend = FALSE,
          show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = FONT_GRID),
          col = col_zscore,
          rect_gp = gp_rect(),
          column_title = title,
          column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
          heatmap_legend_param = std_legend_param())
}

# =============================================================================
# 7.6 FIGURE 5 / SUPPFIG 16 OPTIMIZED FUNCTIONS
# =============================================================================

# --- plot_nc_keygene_scatter() -- Key Gene Bubble Scatter Plot ---
plot_nc_keygene_scatter <- function(df, symbol_col = NULL, fc_col = NULL, 
                                     pval_col = NULL, title = "Key Gene Differentiation",
                                     top_n = 20, point_color = "#3C5488") {
  # Auto-detect column names
  if (is.null(symbol_col)) {
    for (cand in c("symbol", "gene_name", "name", "Name")) {
      if (cand %in% colnames(df)) { symbol_col <- cand; break }
    }
  }
  if (is.null(fc_col)) {
    for (cand in c("logFC", "log2FoldChange", "log2FC", "fc")) {
      if (cand %in% colnames(df)) { fc_col <- cand; break }
    }
  }
  if (is.null(pval_col)) {
    for (cand in c("P.Value", "pvalue", "padj", "p.adjust", "fdr", "p_val")) {
      if (cand %in% colnames(df)) { pval_col <- cand; break }
    }
  }
  
  df$neg_log10p <- -log10(df[[pval_col]] + 1e-300)
  
  # Take top N by p-value
  df <- df[order(df[[pval_col]]), ]
  top_df <- head(df, top_n)
  
  # Filter out unmapped ENS IDs from labels
  top_df$label <- top_df[[symbol_col]]
  top_df$label <- ifelse(grepl("^ENS[GP]", top_df$label), "", top_df$label)
  
  ggplot(df, aes(x = .data[[fc_col]], y = neg_log10p)) +
    geom_point(aes(size = neg_log10p), color = point_color, alpha = 0.6) +
    ggrepel::geom_text_repel(data = top_df, aes(label = label),
                              size = 2.8, fontface = "italic", 
                              segment.size = 0.2, max.overlaps = 20,
                              family = FONT_FAMILY) +
    scale_size_continuous(range = c(0.5, 3), guide = "none") +
    labs(title = title, 
         x = expression(log[2]~"(Fold Change)"), 
         y = expression(-log[10]~"(P-value)")) +
    theme_pub
}

# --- plot_nc_cor_heatmap() -- Cross-module Correlation Heatmap ---
plot_nc_cor_heatmap <- function(cor_mat, title = "Cross-module Correlation") {
  # Clean names
  colnames(cor_mat) <- clean_celltype(colnames(cor_mat))
  rownames(cor_mat) <- clean_celltype(rownames(cor_mat))
  
  # Clamp to [-1, 1]
  cor_mat[cor_mat > 1] <- 1
  cor_mat[cor_mat < -1] <- -1
  
  Heatmap(cor_mat, 
          name = "Cor",
          col = col_cor,
          rect_gp = gp_rect(),
          column_title = title,
          column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
          column_names_rot = 45,
          show_row_dend = FALSE,
          show_column_dend = FALSE,
          row_names_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
          column_names_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
          heatmap_legend_param = std_legend_param())
}

# --- plot_nc_direction_bar() -- Generic Bidirectional Bar Plot (Genes/Pathways) ---
plot_nc_direction_bar <- function(df, name_col = NULL, value_col = NULL, 
                                   title = "", n = 15, 
                                   is_gene = TRUE) {
  # Auto-detect column names
  if (is.null(name_col)) {
    for (cand in c("name", "gene_name", "symbol", "Description", "pathway")) {
      if (cand %in% colnames(df)) { name_col <- cand; break }
    }
  }
  if (is.null(value_col)) {
    for (cand in c("logFC", "diff", "value", "log2FoldChange", "NES", "fc")) {
      if (cand %in% colnames(df)) { value_col <- cand; break }
    }
  }
  
  # Take top n by absolute value
  df <- df[order(-abs(df[[value_col]])), ]
  df <- head(df, n)
  
  # Clean names based on type
  if (is_gene) {
    df[[name_col]] <- ensg_to_symbol(df[[name_col]])
    df <- df[!grepl("^ENS[GP]", df[[name_col]]), ]
  } else {
    df[[name_col]] <- clean_celltype(df[[name_col]])
    # Wrap long names (>40 chars)
    df[[name_col]] <- stringr::str_wrap(df[[name_col]], width = 40)
  }
  
  df$direction <- ifelse(df[[value_col]] > 0, "Up", "Down")
  
  face_style <- if (is_gene) "italic" else "plain"
  font_size <- 8
  
  ggplot(df, aes(x = .data[[value_col]], 
                  y = reorder(.data[[name_col]], .data[[value_col]]), 
                  fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Up" = COL_UP, "Down" = COL_DOWN), guide = "none") +
    labs(title = title, x = expression(log[2]~"(Fold Change)"), y = NULL) +
    theme_pub +
    theme(axis.text.y = element_text(face = face_style, size = font_size))
}

# =============================================================================
# 7.7 FIGURE 8 / SUPPFIG 17 EXTERNAL VALIDATION FUNCTIONS
# =============================================================================

# --- plot_nc_concordance_heatmap() -- Directional Concordance Heatmap ---
plot_nc_concordance_heatmap <- function(mat, title = "Directional Concordance") {
  col_fun <- colorRamp2(c(-1, 0, 1), c(COL_DOWN, "#FFFFFF", COL_UP))
  
  Heatmap(mat, 
          name = "log2FC",
          col = col_fun,
          rect_gp = gp_rect(),
          column_title = title,
          column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
          row_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = FONT_GRID),
          column_names_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
          column_names_rot = 45,
          show_row_dend = FALSE,
          show_column_dend = FALSE,
          cluster_columns = FALSE,
          heatmap_legend_param = std_legend_param())
}

# --- plot_nc_pathway_heatmap() -- External Pathway NES Heatmap ---
plot_nc_pathway_heatmap <- function(nes_mat, title = "External Pathway Validation") {
  # Clean pathway names
  rownames(nes_mat) <- clean_pathway(rownames(nes_mat))
  
  Heatmap(nes_mat,
          name = "NES",
          col = col_nes,
          rect_gp = gp_rect(),
          column_title = title,
          column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_GRID),
          column_names_rot = 45,
          row_names_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
          column_names_gp = gpar(fontsize = 7, fontfamily = FONT_GRID),
          show_row_dend = FALSE,
          show_column_dend = FALSE,
          heatmap_legend_param = std_legend_param())
}

# --- plot_nc_sample_bar() -- Effective Sample Size Bar Plot ---
plot_nc_sample_bar <- function(sample_df, n_col = NULL, name_col = NULL,
                                title = "Effective Sample Size") {
  # Auto-detect column names
  if (is.null(n_col)) {
    for (cand in c("effective_n", "n_total", "n", "sample_size", "Count")) {
      if (cand %in% colnames(sample_df)) { n_col <- cand; break }
    }
  }
  if (is.null(name_col)) {
    for (cand in c("dataset", "Dataset", "study", "name")) {
      if (cand %in% colnames(sample_df)) { name_col <- cand; break }
    }
  }
  
  ggplot(sample_df, aes(x = .data[[n_col]], y = reorder(.data[[name_col]], .data[[n_col]]))) +
    geom_col(fill = COL_UP, width = 0.7) +
    geom_text(aes(label = .data[[n_col]]), hjust = -0.2, size = 2.8, family = FONT_FAMILY) +
    labs(x = "Effective Sample Size (n)", y = NULL, title = title) +
    theme_pub +
    theme(axis.text.y = element_text(size = 8))
}

# =============================================================================
# 7.8 FIGURE 6 / SURVIVAL & CLINICAL FUNCTIONS
# =============================================================================

# --- plot_nc_survival() -- Production-grade KM Survival Curve ---
plot_nc_survival <- function(data, time_col = "OS_time", status_col = "OS_status",
                              group_col = "Risk_Group", title = "",
                              group_colors = c("High" = "#E64B35", "Low" = "#4DBBD5"),
                              group_labels = NULL) {
  # Build survival formula
  surv_formula <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col))
  fit <- survival::survfit(surv_formula, data = data)
  
  # Auto-generate labels if not provided
  if (is.null(group_labels)) {
    group_labels <- names(group_colors)
  }
  
  survminer::ggsurvplot(
    fit,
    data = data,
    palette = unname(group_colors),
    pval = TRUE,
    pval.size = 3,
    pval.coord = c(0, 0.1),
    conf.int = TRUE,
    conf.int.alpha = 0.1,
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 2.8,
    legend.labs = group_labels,
    xlab = "Time (Months)",
    ggtheme = theme_pub,
    title = title
  )
}

# --- plot_nc_clinical_box() -- Clinical Feature Boxplot with Stats ---
plot_nc_clinical_box <- function(df, feature_col, group_col = "Subtype", 
                                  title = NULL,
                                  group_colors = c("CS1" = COL_CS1, "CS2" = COL_CS2)) {
  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[feature_col]], 
                       fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    scale_fill_manual(values = group_colors) +
    ggpubr::stat_compare_means(comparisons = list(names(group_colors)),
                                label = "p.format", size = 2.8) +
    labs(y = feature_col, x = NULL, title = title) +
    theme_pub +
    theme(legend.position = "none")
  
  return(p)
}

# =============================================================================
# 8. PROJECT PATHS
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "analysis/results")
DATA <- file.path(BASE, "analysis/data/processed")
OUT  <- file.path(BASE, "图片")

# =============================================================================
# Self-test
# =============================================================================
if (sys.nframe() == 0) {
  cat("nc_theme.R self-test:\n")
  cat("  FONT_FAMILY:", FONT_FAMILY, "\n")
  cat("  COL_UP:", COL_UP, " COL_DOWN:", COL_DOWN, "\n")
  cat("  PAL_CAT (NPG):", paste(PAL_CAT, collapse=", "), "\n")
  cat("  All checks passed.\n")
}
