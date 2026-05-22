#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_07_standalone.R
# Supplementary Figure 7: HAE-Distinctive Pathway Signatures
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# =============================================================================
# Panels:
#   A = HAE-distinctive pathway lollipop (8 pathways with opposite direction)
#   B = Cross-disease Spearman correlation to HAE (rho + significance)
#   C = 8 distinctive pathways x 5 diseases delta heatmap with HAE highlighted
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_07_standalone.R
# =============================================================================

cat("=== Supplementary Figure 7: HAE-Distinctive Pathway Signatures ===\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
CROSS_DIR <- file.path(RES, "cross_disease_positioning")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_07")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants & theme
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4
ASSEMBLY_DPI <- 600
mm2in <- function(mm) mm / MM_PER_INCH
FS_TAG <- 12

COL_HAE      <- "#CD534CFF"
COL_HCC      <- "#0073C2FF"
COL_CCA      <- "#EFC000FF"
COL_FIBROSIS <- "#20854EFF"
COL_NAFLD    <- "#7876B1FF"
DISEASE_COLS <- c("HAE" = COL_HAE, "HCC" = COL_HCC, "CCA" = COL_CCA,
                  "Fibrosis" = COL_FIBROSIS, "NAFLD" = COL_NAFLD)

theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.text  = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                              margin = margin(b = 1)),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5,
                                 margin = margin(b = 2)),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(2.5, "mm"),
    plot.margin = margin(4, 4, 4, 4, "mm")
  )
theme_set(theme_nc)
ht_opt$message <- FALSE

save_pdf <- function(fn, w, h, expr) {
  cairo_pdf(file.path(OUT, fn), width = mm2in(w), height = mm2in(h),
            family = FONT_FAMILY)
  force(expr); dev.off()
  cat(sprintf("  -> %s (%dx%dmm)\n", fn, w, h))
}

clean_pathway <- function(x) {
  # Curated display map preserving biological acronyms (HUGO/IUPHAR conventions)
  display_map <- c(
    "TNFA SIGNALING VIA NFKB"            = "TNF-\u03b1 signaling via NF-\u03baB",
    "HYPOXIA"                            = "Hypoxia",
    "CHOLESTEROL HOMEOSTASIS"            = "Cholesterol homeostasis",
    "MITOTIC SPINDLE"                    = "Mitotic spindle",
    "WNT BETA CATENIN SIGNALING"         = "Wnt/\u03b2-catenin signaling",
    "TGF BETA SIGNALING"                 = "TGF-\u03b2 signaling",
    "IL6 JAK STAT3 SIGNALING"            = "IL-6/JAK/STAT3 signaling",
    "DNA REPAIR"                         = "DNA repair",
    "G2M CHECKPOINT"                     = "G2/M checkpoint",
    "APOPTOSIS"                          = "Apoptosis",
    "NOTCH SIGNALING"                    = "Notch signaling",
    "ADIPOGENESIS"                       = "Adipogenesis",
    "ESTROGEN RESPONSE EARLY"            = "Estrogen response (early)",
    "ESTROGEN RESPONSE LATE"             = "Estrogen response (late)",
    "ANDROGEN RESPONSE"                  = "Androgen response",
    "MYOGENESIS"                         = "Myogenesis",
    "PROTEIN SECRETION"                  = "Protein secretion",
    "INTERFERON ALPHA RESPONSE"          = "Interferon-\u03b1 response",
    "INTERFERON GAMMA RESPONSE"          = "Interferon-\u03b3 response",
    "APICAL JUNCTION"                    = "Apical junction",
    "APICAL SURFACE"                     = "Apical surface",
    "HEDGEHOG SIGNALING"                 = "Hedgehog signaling",
    "COMPLEMENT"                         = "Complement",
    "UNFOLDED PROTEIN RESPONSE"          = "Unfolded protein response",
    "PI3K AKT MTOR SIGNALING"            = "PI3K/AKT/mTOR signaling",
    "MTORC1 SIGNALING"                   = "mTORC1 signaling",
    "E2F TARGETS"                        = "E2F targets",
    "MYC TARGETS V1"                     = "MYC targets V1",
    "MYC TARGETS V2"                     = "MYC targets V2",
    "EPITHELIAL MESENCHYMAL TRANSITION"  = "Epithelial-mesenchymal transition",
    "INFLAMMATORY RESPONSE"              = "Inflammatory response",
    "XENOBIOTIC METABOLISM"              = "Xenobiotic metabolism",
    "FATTY ACID METABOLISM"              = "Fatty acid metabolism",
    "OXIDATIVE PHOSPHORYLATION"          = "Oxidative phosphorylation",
    "GLYCOLYSIS"                         = "Glycolysis",
    "REACTIVE OXYGEN SPECIES PATHWAY"    = "Reactive oxygen species pathway",
    "P53 PATHWAY"                        = "p53 pathway",
    "UV RESPONSE UP"                     = "UV response (up)",
    "UV RESPONSE DN"                     = "UV response (down)",
    "ANGIOGENESIS"                       = "Angiogenesis",
    "HEME METABOLISM"                    = "Heme metabolism",
    "COAGULATION"                        = "Coagulation",
    "BILE ACID METABOLISM"               = "Bile acid metabolism",
    "PEROXISOME"                         = "Peroxisome",
    "ALLOGRAFT REJECTION"                = "Allograft rejection",
    "SPERMATOGENESIS"                    = "Spermatogenesis",
    "PANCREAS BETA CELLS"                = "Pancreas \u03b2 cells",
    "KRAS SIGNALING UP"                  = "KRAS signaling (up)",
    "KRAS SIGNALING DN"                  = "KRAS signaling (down)"
  )
  # Normalize input: strip HALLMARK_ prefix, replace _ with space, upper-case
  key <- toupper(gsub("_", " ", sub("^HALLMARK_", "", x, ignore.case = TRUE)))
  out <- ifelse(key %in% names(display_map), display_map[key], NA_character_)
  unmapped <- is.na(out)
  if (any(unmapped)) {
    fb <- tolower(key[unmapped])
    fb <- sub("^(\\w)", "\\U\\1", fb, perl = TRUE)
    out[unmapped] <- fb
  }
  unname(out)
}

sig_star <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", "ns")))
}

# =============================================================================
# Data
# =============================================================================
cross_json <- fromJSON(file.path(CROSS_DIR, "cross_disease_results.json"),
                       simplifyDataFrame = TRUE)
unique_paths <- cross_json$unique_details
unique_df <- data.frame(
  pathway = unique_paths$pathway,
  delta   = unique_paths$delta,
  stringsAsFactors = FALSE
) %>%
  mutate(pathway_clean = clean_pathway(pathway),
         direction = ifelse(delta > 0, "Up in HAE", "Down in HAE")) %>%
  arrange(delta)
unique_df$pathway_clean <- factor(unique_df$pathway_clean,
                                  levels = unique_df$pathway_clean)

# correlations to HAE
cor_list <- cross_json$correlations
cor_df <- data.frame(
  Disease = names(cor_list),
  rho     = vapply(cor_list, function(x) x$rho, numeric(1)),
  P       = vapply(cor_list, function(x) x$P,   numeric(1)),
  stringsAsFactors = FALSE
)
cor_df$Sig <- sig_star(cor_df$P)
cor_df <- cor_df %>% arrange(desc(rho))
cor_df$Disease <- factor(cor_df$Disease, levels = cor_df$Disease)

# distances
dist_list <- cross_json$distances
dist_df <- data.frame(
  Disease = names(dist_list),
  Distance = unlist(dist_list),
  stringsAsFactors = FALSE
) %>% arrange(Distance)
dist_df$Disease <- factor(dist_df$Disease, levels = dist_df$Disease)

# pathway comparison matrix (28 hallmarks x 5 diseases)
pathmat <- read.csv(file.path(CROSS_DIR, "pathway_comparison_matrix.csv"),
                    row.names = 1, stringsAsFactors = FALSE)

# subset distinctive 8 pathways
distinctive_keys <- toupper(unique_df$pathway)
distinctive_keys <- intersect(distinctive_keys, rownames(pathmat))
heat_mat <- as.matrix(pathmat[distinctive_keys, c("HAE","HCC","CCA","Fibrosis","NAFLD")])
rownames(heat_mat) <- clean_pathway(rownames(heat_mat))
# order by HAE delta
heat_mat <- heat_mat[order(heat_mat[,"HAE"]), , drop = FALSE]

# =============================================================================
# Panel A: HAE-distinctive lollipop
# =============================================================================
cat("  Panel A: HAE-distinctive pathways...\n")

pA <- ggplot(unique_df, aes(x = delta, y = pathway_clean, color = direction)) +
  geom_segment(aes(x = 0, xend = delta, y = pathway_clean, yend = pathway_clean),
               linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_color_manual(values = c("Up in HAE" = COL_HAE, "Down in HAE" = COL_HCC),
                     name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  labs(x = "HAE GSVA delta",
       y = NULL,
       title = "HAE-distinctive pathways",
       subtitle = "Opposite to all reference diseases (n=8)") +
  theme(legend.position = c(0.98, 0.04),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                         color = NA),
        legend.key.size = unit(3, "mm"),
        axis.text.y = element_text(size = 8))

save_pdf("Supp07a_distinctive_paths.pdf", 91, 110, print(pA))

# =============================================================================
# Panel B: Cross-disease Spearman correlations
# =============================================================================
cat("  Panel B: Cross-disease correlations...\n")

pB <- ggplot(cor_df, aes(x = Disease, y = rho, fill = Disease)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.35) +
  geom_text(aes(label = sprintf("\u03c1 = %.2f", rho)),
            vjust = -2.0, size = 2.85, family = FONT_FAMILY,
            fontface = "bold") +
  geom_text(aes(label = Sig), vjust = -0.5, size = 3.0,
            family = FONT_FAMILY, color = "grey25") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = DISEASE_COLS) +
  scale_y_continuous(limits = c(0, 0.95), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Spearman \u03c1 to HAE pathway profile",
       title = "Pathway profile concordance",
       subtitle = "HAE vs four reference liver diseases (Spearman)") +
  theme(legend.position = "none")

save_pdf("Supp07b_correlations.pdf", 92, 110, print(pB))

# =============================================================================
# Panel C: 8 distinctive pathways x 5 diseases delta heatmap
# =============================================================================
cat("  Panel C: Distinctive pathways heatmap...\n")

col_delta <- colorRamp2(c(-0.5, 0, 0.5), c("#0073C2FF","#FFFFFF","#CD534CFF"))

# Column annotation highlighting HAE
group_vec <- ifelse(colnames(heat_mat) == "HAE", "HAE", "Reference")
ha_col <- HeatmapAnnotation(
  which = "column",
  Group = group_vec,
  col = list(Group = c("HAE" = COL_HAE, "Reference" = "grey75")),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
  show_legend = TRUE,
  annotation_legend_param = list(
    Group = list(
      title = "Group",
      title_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
      grid_height = unit(3.2, "mm"),
      grid_width  = unit(3.2, "mm")
    )
  ),
  simple_anno_size = unit(2.5, "mm")
)

ht_C <- Heatmap(
  heat_mat,
  name = "Delta",
  col  = col_delta,
  cluster_rows = FALSE, cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  top_annotation = ha_col,
  row_names_gp    = gpar(fontsize = 8, fontfamily = FONT_GRID),
  column_names_gp = gpar(fontsize = 8.5, fontfamily = FONT_GRID, fontface = "bold"),
  column_names_rot = 0, column_names_centered = TRUE,
  heatmap_legend_param = list(
    title = "GSVA delta",
    title_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
    legend_height = unit(2.4, "cm"),
    title_position = "leftcenter-rot"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", heat_mat[i, j]), x, y,
              gp = gpar(fontsize = 7.95, fontfamily = FONT_GRID,
                        col = ifelse(abs(heat_mat[i,j]) > 0.35, "white", "black")))
  },
  column_title = "Distinctive pathway direction across diseases",
  column_title_gp = gpar(fontsize = 10, fontfamily = FONT_GRID, fontface = "bold"),
  row_names_max_width = max_text_width(rownames(heat_mat),
                                       gp = gpar(fontsize = 8))
)

save_pdf("Supp07c_distinctive_heatmap.pdf", 183, 125, draw(ht_C))

# =============================================================================
# Composite Assembly
# =============================================================================
cat("\n--- Assembling composite SuppFig_07 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  H1 <- 110; H2 <- H_TOTAL - H1   # 110 + 135 = 245
  W_L <- 91; W_R <- W_TOTAL - W_L

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    # Panel A (Row 1 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H2,"mm"),
      width=unit(W_L,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(pA, newpage=FALSE); grid::popViewport()
    # Panel B (Row 1 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H2,"mm"),
      width=unit(W_R,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(pB, newpage=FALSE); grid::popViewport()
    # Panel C (Row 2 full width)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_TOTAL,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    draw(ht_C, newpage=FALSE); grid::popViewport()
    # Labels
    label_data <- data.frame(
      text=c("a","b","c"),
      x_mm=c(2, W_L+2, 2),
      y_mm=c(H_TOTAL-2, H_TOTAL-2, H2-2), stringsAsFactors=FALSE)
    for(i in seq_len(nrow(label_data))) {
      grid::grid.text(label=label_data$text[i],
        x=unit(label_data$x_mm[i],"mm"), y=unit(label_data$y_mm[i],"mm"),
        just=c("left","top"),
        gp=grid::gpar(fontsize=FS_TAG, fontface="bold", fontfamily=FONT_FAMILY))
    }
    grid::popViewport()
  }

  cairo_pdf(file.path(OUT,"SuppFig_07.pdf"),
            width=W_TOTAL/25.4, height=H_TOTAL/25.4, family=FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_07.pdf (vector)\n")

  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT,"SuppFig_07.png"),
                 width=px_W, height=px_H, res=DPI, type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_07.png\n")

  grDevices::tiff(file.path(OUT,"SuppFig_07.tiff"),
                  width=px_W, height=px_H, res=DPI, compression="lzw", type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_07.tiff\n")

  cat("  SuppFig_07 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 7 rendering complete ===\n")
