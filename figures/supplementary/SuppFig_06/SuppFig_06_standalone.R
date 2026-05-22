#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_06_standalone.R — Deep Optimization v2
# Supplementary Figure 6: Cross-Disease Hallmark Benchmarking & Molecular Positioning
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Nature Communications | 183x245mm | 600DPI | Zero-Distortion Assembly
# =============================================================================
# Panels:
#   a = Cross-disease hallmark radar chart (10 hallmarks x 4 diseases)
#   b = Cross-disease similarity index barplot (cosine similarity)
#   c = Hallmark classification with bootstrap CI
#   d = Proliferation gap lollipop (E2F, G2M, MYC)
#   e = Cross-disease Spearman correlation heatmap
#   f = MDS positioning plot
#   g = Euclidean distance barplot
#   h = Pathway comparison heatmap (top 18)
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_06_standalone.R
# =============================================================================

cat("=== Supplementary Figure 6: Cross-Disease Hallmark Benchmarking ===\n")
cat("  [v2 Deep Optimization: Content-Adaptive Layout + Zero-Distortion Assembly]\n\n")

# =============================================================================
# SECTION 1: Library Imports
# =============================================================================
cat("  Loading libraries...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(magick)
  library(ggsci)
  library(ggrepel)
  library(jsonlite)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
cat("  Arial registered\n")

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
MIM_DIR   <- file.path(RES, "cancer_hallmark_mimicry")
CROSS_DIR <- file.path(RES, "cross_disease_positioning")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_06")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Assembly Constants (Content-Adaptive Layout)
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
FS_TAG <- 12
MM_PER_INCH <- 25.4
ASSEMBLY_DPI <- 600

# Total canvas
W_TOTAL <- 183; H_TOTAL <- 245

# Content-adaptive row heights (sum = 245)
H1 <- 68; H2 <- 56; H3 <- 58; H4 <- 63
stopifnot(H1 + H2 + H3 + H4 == H_TOTAL)

# Content-adaptive column widths per row (each row sums to 183)
W1_L <- 102; W1_R <- W_TOTAL - W1_L  # Row 1: radar(102) + bar(81)
W2_L <- 108; W2_R <- W_TOTAL - W2_L  # Row 2: classification(108) + lollipop(75)
W3_L <-  89; W3_R <- W_TOTAL - W3_L  # Row 3: cor heatmap(89) + MDS(94)
W4_L <-  72; W4_R <- W_TOTAL - W4_L  # Row 4: distance(72) + pathway heatmap(111)

cat(sprintf("  Layout: Row heights = %d+%d+%d+%d = %dmm\n", H1, H2, H3, H4, H_TOTAL))
cat(sprintf("  Row widths: %d+%d | %d+%d | %d+%d | %d+%d\n",
            W1_L, W1_R, W2_L, W2_R, W3_L, W3_R, W4_L, W4_R))

# =============================================================================
# SECTION 4: Refined Color Palette (NC Premium)
# =============================================================================
# Elevated disease colors - deeper saturation, harmonious
COL_HAE      <- "#C44E52"  # warm terracotta red
COL_HCC      <- "#4C72B0"  # slate blue
COL_CCA      <- "#DD8452"  # amber orange
COL_FIBROSIS <- "#55A868"  # sage green
COL_NAFLD    <- "#8172B3"  # muted violet

DISEASE_COLS <- c("HAE" = COL_HAE, "HCC" = COL_HCC, "CCA" = COL_CCA,
                  "Fibrosis" = COL_FIBROSIS, "NAFLD" = COL_NAFLD)
CAT_COLS <- c("Paradoxical" = "#C44E52", "Partially mimicked" = "#DD8452",
              "Not mimicked" = "#AAAAAA")

col_nes <- colorRamp2(c(-0.8, 0, 0.8), c("#4C72B0", "#FAFAFA", "#C44E52"))

# =============================================================================
# SECTION 5: Premium Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text               = element_text(family = FONT_FAMILY, size = 8),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(linewidth = 0.4, color = "black", fill = NA),
    axis.line          = element_blank(),
    axis.ticks         = element_line(linewidth = 0.3, color = "black"),
    axis.ticks.length  = unit(1.2, "mm"),
    axis.text          = element_text(size = 8, color = "black"),
    axis.title         = element_text(size = 8, face = "bold"),
    plot.title         = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.title.position = "panel",
    plot.subtitle      = element_text(size = 8, hjust = 0.5, color = "grey35",
                                       margin = margin(b = 2)),
    legend.text        = element_text(size = 8),
    legend.title       = element_text(size = 8, face = "bold"),
    legend.key.size    = unit(3, "mm"),
    legend.background  = element_blank(),
    legend.margin      = margin(0, 0, 0, 0),
    plot.margin        = margin(5, 5, 4, 5, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# SECTION 6: Helper Functions
# =============================================================================
ht_opt$message <- FALSE

# Pixel-exact PDF rendering (skill v2 standard)
save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat(sprintf("  -> %s (%d x %d mm)\n", basename(fpath), w_mm, h_mm))
}

save_heatmap_panel_pdf <- function(filename, ht_obj, w_mm, h_mm, title = NULL) {
  fpath <- file.path(OUT, filename)
  target_w_px <- round(w_mm * ASSEMBLY_DPI / 25.4)
  target_h_px <- round(h_mm * ASSEMBLY_DPI / 25.4)
  w_in <- target_w_px / ASSEMBLY_DPI
  h_in <- target_h_px / ASSEMBLY_DPI
  cairo_pdf(fpath, width = w_in, height = h_in, family = FONT_FAMILY)
  if (!is.null(title)) {
    draw(ht_obj, padding = unit(c(2, 2, 8, 2), "mm"))
    grid::grid.text(title, x = unit(0.5, "npc"), y = unit(1, "npc") - unit(3, "mm"),
                    just = c("center", "top"),
                    gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY))
  } else {
    draw(ht_obj, padding = unit(c(2, 2, 2, 2), "mm"))
  }
  dev.off()
  cat(sprintf("  -> %s (%d x %d mm)\n", basename(fpath), w_mm, h_mm))
}

clean_pathway <- function(x) {
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  # Scientific capitalization corrections
  x <- gsub("\\bIl6\\b", "IL-6", x)
  x <- gsub("\\bIl2\\b", "IL-2", x)
  x <- gsub("\\bJak\\b", "JAK", x)
  x <- gsub("\\bStat3\\b", "STAT3", x)
  x <- gsub("\\bStat5\\b", "STAT5", x)
  x <- gsub("\\bE2f\\b", "E2F", x)
  x <- gsub("\\bG2m\\b", "G2M", x)
  x <- gsub("\\bMyc\\b", "MYC", x)
  x <- gsub("\\bDna\\b", "DNA", x)
  x <- gsub("\\bRna\\b", "RNA", x)
  x <- gsub("\\bTnfa\\b", "TNF-a", x)
  x <- gsub("\\bNfkb\\b", "NF-kB", x)
  x <- gsub("\\bNf-kb\\b", "NF-kB", x)
  x <- gsub("\\bWnt\\b", "WNT", x)
  x <- gsub("\\bTgf\\b", "TGF", x)
  x <- gsub("\\bMtor\\b", "mTOR", x)
  x <- gsub("\\bMtorc1\\b", "mTORC1", x)
  x <- gsub("\\bP53\\b", "p53", x)
  x <- gsub("\\bPi3k\\b", "PI3K", x)
  x <- gsub("\\bAkt\\b", "AKT", x)
  x <- gsub("\\bKras\\b", "KRAS", x)
  x <- gsub("\\bUv\\b", "UV", x)
  x <- gsub("\\bRos\\b", "ROS", x)
  x
}

# =============================================================================
# SECTION 7: Data Loading
# =============================================================================
cat("  Loading data...\n")

radar <- read.csv(file.path(MIM_DIR, "hallmark_radar_data.csv"), stringsAsFactors = FALSE)
mimicry <- read.csv(file.path(MIM_DIR, "cancer_mimicry_index.csv"), stringsAsFactors = FALSE)
hclass <- read.csv(file.path(MIM_DIR, "hallmark_classification.csv"), stringsAsFactors = FALSE)
bootci <- read.csv(file.path(MIM_DIR, "hallmark_bootstrap_ci.csv"), stringsAsFactors = FALSE)
cormat <- read.csv(file.path(CROSS_DIR, "correlation_matrix.csv"), row.names = 1, stringsAsFactors = FALSE)
pathmat <- read.csv(file.path(CROSS_DIR, "pathway_comparison_matrix.csv"), row.names = 1, stringsAsFactors = FALSE)
cross_json <- fromJSON(file.path(CROSS_DIR, "cross_disease_results.json"))

cat("  Data loaded.\n\n")

# =============================================================================
# PANEL A: Radar Chart (Manual Cartesian Polar - Vector Safe)
# Target: 102 x 68 mm
# =============================================================================
cat("  Panel A: Cross-disease hallmark radar (102x68mm)...\n")

abbrev_hallmark <- function(s) {
  s <- gsub("Sustaining Proliferative Signaling", "Proliferative\nSignaling", s)
  s <- gsub("Evading Growth Suppressors",         "Growth\nSuppressors", s)
  s <- gsub("Resisting Cell Death",                "Cell Death", s)
  s <- gsub("Enabling Replicative Immortality",    "Replicative\nImmortality", s)
  s <- gsub("Inducing Angiogenesis",               "Angiogenesis", s)
  s <- gsub("Activating Invasion & Metastasis",    "Invasion &\nMetastasis", s)
  s <- gsub("Genome Instability & Mutation",       "Genome\nInstability", s)
  s <- gsub("Tumor-promoting Inflammation",        "Tumor\nInflammation", s)
  s <- gsub("Deregulating Cellular Energetics",    "Cellular\nEnergetics", s)
  s <- gsub("Avoiding Immune Destruction",         "Immune\nEvasion", s)
  s
}

# Drop hallmarks with no representative pathway (all-zero across diseases)
# "Enabling Replicative Immortality" has no mapped MSigDB pathway -> all 0 -> meaningless axis
radar_vals  <- radar[, -1]
keep_axes   <- rowSums(abs(radar_vals)) > 1e-9
radar       <- radar[keep_axes, ]

hallmarks   <- radar$hallmark
n_axes      <- length(hallmarks)
axis_angles <- pi/2 - 2*pi*(seq_len(n_axes) - 1)/n_axes
axis_labels <- abbrev_hallmark(hallmarks)

delta_min <- min(unlist(radar[, -1]), na.rm = TRUE)
delta_max <- max(unlist(radar[, -1]), na.rm = TRUE)
pad       <- 0.08 * (delta_max - delta_min)
r_min     <- delta_min - pad
r_max     <- delta_max + pad
to_r      <- function(v) (v - r_min) / (r_max - r_min)

# Concentric grid (3 levels for cleaner look)
grid_levels <- c(0.25, 0.5, 0.75, 1.0)
grid_circles <- do.call(rbind, lapply(grid_levels, function(rr) {
  ang <- seq(0, 2*pi, length.out = 200)
  data.frame(r = rr, x = rr*cos(ang), y = rr*sin(ang), grp = sprintf("g%.2f", rr))
}))

spokes <- data.frame(x0 = 0, y0 = 0, x1 = cos(axis_angles), y1 = sin(axis_angles))

zero_r <- to_r(0)
zero_circle <- data.frame(
  x = zero_r*cos(seq(0, 2*pi, length.out = 200)),
  y = zero_r*sin(seq(0, 2*pi, length.out = 200))
)

axis_label_df <- data.frame(
  x = 1.22*cos(axis_angles), y = 1.22*sin(axis_angles),
  label = axis_labels
)

build_poly <- function(disease_name) {
  vals <- to_r(radar[[disease_name]])
  data.frame(
    Disease = disease_name,
    x = c(vals, vals[1])*cos(c(axis_angles, axis_angles[1])),
    y = c(vals, vals[1])*sin(c(axis_angles, axis_angles[1]))
  )
}
poly_df <- do.call(rbind, lapply(c("HAE", "HCC", "CCA", "Fibrosis"), build_poly))
poly_df$Disease <- factor(poly_df$Disease, levels = c("HAE", "HCC", "CCA", "Fibrosis"))

pA <- ggplot() +
  # Outer boundary
  annotate("path", x = cos(seq(0, 2*pi, length.out = 200)),
           y = sin(seq(0, 2*pi, length.out = 200)),
           color = "grey70", linewidth = 0.4) +
  # Grid circles
  geom_path(data = grid_circles, aes(x = x, y = y, group = grp),
            color = "grey85", linewidth = 0.2) +
  # Zero-reference
  geom_path(data = zero_circle, aes(x = x, y = y),
            color = "grey50", linewidth = 0.35, linetype = "21") +
  # Spokes
  geom_segment(data = spokes, aes(x = x0, y = y0, xend = x1, yend = y1),
               color = "grey82", linewidth = 0.2) +
  # Disease polygons (HAE emphasized)
  geom_polygon(data = poly_df %>% filter(Disease != "HAE"),
               aes(x = x, y = y, color = Disease, fill = Disease, group = Disease),
               alpha = 0.05, linewidth = 0.5, linetype = "dashed") +
  geom_polygon(data = poly_df %>% filter(Disease == "HAE"),
               aes(x = x, y = y, color = Disease, fill = Disease, group = Disease),
               alpha = 0.15, linewidth = 0.8) +
  geom_point(data = poly_df %>% filter(Disease == "HAE"),
             aes(x = x, y = y, color = Disease), size = 1.5, shape = 16) +
  geom_point(data = poly_df %>% filter(Disease != "HAE"),
             aes(x = x, y = y, color = Disease), size = 0.9, shape = 16) +
  # Axis labels
  geom_text(data = axis_label_df, aes(x = x, y = y, label = label),
            family = FONT_FAMILY, size = 2.85, lineheight = 0.82, color = "grey20") +
  scale_color_manual(values = DISEASE_COLS[c("HAE", "HCC", "CCA", "Fibrosis")]) +
  scale_fill_manual(values = DISEASE_COLS[c("HAE", "HCC", "CCA", "Fibrosis")]) +
  coord_fixed(xlim = c(-1.55, 1.55), ylim = c(-1.45, 1.45), clip = "off") +
  labs(title = "Cross-disease Hallmark profiles",
       subtitle = "GSVA delta (Adjacent \u2013 Normal)") +
  theme_void(base_family = FONT_FAMILY) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                              margin = margin(b = 1)),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey35",
                                 margin = margin(b = 1)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8, family = FONT_FAMILY),
    legend.key.size = unit(3.5, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.margin = margin(t = -4),
    legend.spacing.x = unit(1, "mm"),
    plot.margin = margin(3, 3, 2, 3, "mm")
  ) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

save_panel_pdf("Supp06a_hallmark_radar.pdf", W1_L, H1, print(pA))

# =============================================================================
# PANEL B: Similarity Index Barplot
# Target: 81 x 68 mm
# =============================================================================
cat("  Panel B: Similarity index (81x68mm)...\n")

mimicry <- mimicry %>% arrange(desc(cosine_similarity))
mimicry$disease <- factor(mimicry$disease, levels = mimicry$disease)

pB <- ggplot(mimicry, aes(x = disease, y = cosine_similarity, fill = disease)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", cosine_similarity)),
            vjust = ifelse(mimicry$cosine_similarity >= 0, -0.6, 1.3),
            size = 2.85, family = FONT_FAMILY, color = "grey20") +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = DISEASE_COLS) +
  scale_y_continuous(limits = c(-0.22, 0.65), expand = c(0, 0)) +
  labs(x = NULL, y = "Cosine similarity\nto HAE",
       title = "Cross-disease similarity index",
       subtitle = "HAE vs. reference liver diseases") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, face = "bold"))

save_panel_pdf("Supp06b_similarity_index.pdf", W1_R, H1, print(pB))

# =============================================================================
# PANEL C: Hallmark Classification
# Target: 108 x 56 mm
# =============================================================================
cat("  Panel C: Hallmark classification (108x56mm)...\n")

hclass <- hclass %>% arrange(HAE_delta)
# Short single-line labels to prevent overlap
hclass$hallmark_short <- c(
  "Sustaining Proliferative Signaling" = "Proliferative Signaling",
  "Evading Growth Suppressors" = "Growth Suppressors",
  "Resisting Cell Death" = "Cell Death Resistance",
  "Enabling Replicative Immortality" = "Replicative Immortality",
  "Inducing Angiogenesis" = "Angiogenesis",
  "Activating Invasion & Metastasis" = "Invasion & Metastasis",
  "Genome Instability & Mutation" = "Genome Instability",
  "Tumor-promoting Inflammation" = "Tumor Inflammation",
  "Deregulating Cellular Energetics" = "Cellular Energetics",
  "Avoiding Immune Destruction" = "Immune Evasion"
)[hclass$hallmark]
hclass$hallmark_short <- factor(hclass$hallmark_short, levels = hclass$hallmark_short)

pC <- ggplot(hclass, aes(x = HAE_delta, y = hallmark_short, fill = category)) +
  geom_col(width = 0.55) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey30") +
  scale_fill_manual(values = CAT_COLS, name = "Classification") +
  labs(x = "HAE GSVA delta", y = NULL,
       title = "Hallmark classification",
       subtitle = "Paradoxical / Partially mimicked / Not mimicked") +
  theme(legend.position = c(0.75, 0.12),
        legend.direction = "vertical",
        legend.key.size = unit(2.5, "mm"),
        legend.background = element_rect(fill = alpha("white", 0.85), linewidth = 0),
        axis.text.y = element_text(size = 8))

save_panel_pdf("Supp06c_hallmark_class.pdf", W2_L, H2, print(pC))

# =============================================================================
# PANEL D: Proliferation Gap Lollipop
# Target: 75 x 56 mm
# =============================================================================
cat("  Panel D: Proliferation gap (75x56mm)...\n")

prolif_paths <- c("E2F TARGETS", "G2M CHECKPOINT", "MYC TARGETS V1")
prolif_data <- pathmat[prolif_paths, c("HAE", "HCC")] %>%
  as.data.frame() %>%
  mutate(pathway = rownames(.)) %>%
  pivot_longer(cols = c(HAE, HCC), names_to = "Disease", values_to = "Delta") %>%
  mutate(pathway = clean_pathway(pathway))

prolif_data$pathway <- factor(prolif_data$pathway,
                              levels = rev(clean_pathway(prolif_paths)))

pD <- ggplot(prolif_data, aes(x = Delta, y = pathway, color = Disease)) +
  geom_segment(data = prolif_data %>% pivot_wider(names_from = Disease, values_from = Delta),
               aes(x = HAE, xend = HCC, y = pathway, yend = pathway),
               color = "grey70", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(size = 3.5, shape = 16) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("HAE" = COL_HAE, "HCC" = COL_HCC)) +
  labs(x = "Pathway delta (GSVA)", y = NULL,
       title = "Proliferation gap",
       subtitle = "HAE suppressed vs. HCC activated") +
  theme(legend.position = c(0.80, 0.20),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.85), linewidth = 0),
        legend.key.size = unit(3, "mm"))

save_panel_pdf("Supp06d_proliferation_gap.pdf", W2_R, H2, print(pD))

# =============================================================================
# PANEL E: Cross-Disease Correlation Heatmap
# Target: 89 x 58 mm
# =============================================================================
cat("  Panel E: Cross-disease correlation heatmap (89x58mm)...\n")

cor_mat <- as.matrix(cormat)
col_cor <- colorRamp2(c(-0.1, 0, 0.4, 0.7, 1.0),
                      c("#4C72B0", "#E8EDF5", "#FAFAFA", "#F5D5D5", "#C44E52"))

ht_E <- Heatmap(cor_mat,
  name = "\u03c1",
  col = col_cor,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
  column_names_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 1.5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
    legend_width = unit(20, "mm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- cor_mat[i, j]
    txt_col <- ifelse(v >= 0.99, "white", "black")
    grid.text(sprintf("%.2f", v), x, y,
              gp = gpar(fontsize = 8, fontfamily = FONT_GRID, col = txt_col))
  }
)

save_heatmap_panel_pdf("Supp06e_cross_disease_cor.pdf", ht_E, W3_L, H3,
                       title = "Cross-disease correlation")

# =============================================================================
# PANEL F: MDS Positioning
# Target: 94 x 58 mm
# =============================================================================
cat("  Panel F: MDS positioning (94x58mm)...\n")

dist_mat <- as.dist(1 - cor_mat)
mds <- cmdscale(dist_mat, k = 2)
mds_df <- data.frame(Disease = rownames(mds), Dim1 = mds[, 1], Dim2 = mds[, 2])

pF <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = Disease)) +
  geom_point(size = 4.5, shape = 16) +
  geom_text_repel(aes(label = Disease), size = 2.85, family = FONT_FAMILY,
                  show.legend = FALSE, box.padding = 0.5, point.padding = 0.3,
                  segment.color = "grey60", segment.size = 0.3,
                  force = 2, max.overlaps = 20) +
  scale_color_manual(values = DISEASE_COLS) +
  labs(x = "MDS Dimension 1", y = "MDS Dimension 2",
       title = "Molecular positioning (MDS)",
       subtitle = "Based on pathway correlation distance") +
  theme(legend.position = "none",
        panel.grid.major = element_line(linewidth = 0.15, color = "grey92"))

save_panel_pdf("Supp06f_MDS.pdf", W3_R, H3, print(pF))

# =============================================================================
# PANEL G: Euclidean Distance
# Target: 72 x 63 mm
# =============================================================================
cat("  Panel G: Euclidean distance (72x63mm)...\n")

dist_df <- mimicry %>%
  select(disease, euclidean_distance) %>%
  arrange(euclidean_distance)
dist_df$disease <- factor(dist_df$disease, levels = dist_df$disease)

pG <- ggplot(dist_df, aes(x = disease, y = euclidean_distance, fill = disease)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", euclidean_distance)),
            vjust = -0.6, size = 2.85, family = FONT_FAMILY, color = "grey20") +
  scale_fill_manual(values = DISEASE_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Euclidean distance\nfrom HAE",
       title = "Distance to HAE",
       subtitle = "Pathway profile distance") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, face = "bold"))

save_panel_pdf("Supp06g_euclidean_dist.pdf", W4_L, H4, print(pG))

# =============================================================================
# PANEL H: Pathway Comparison Heatmap (Top 18)
# Target: 111 x 63 mm
# =============================================================================
cat("  Panel H: Pathway comparison heatmap (111x63mm)...\n")

# Optimal row ordering: maximize cross-disease visual contrast
# Strategy: group by HAE direction vs others, sort by divergence magnitude
hae_abs <- abs(pathmat$HAE)
top18_idx <- order(hae_abs, decreasing = TRUE)[1:18]
top18_mat <- as.matrix(pathmat[top18_idx, ])
rownames(top18_mat) <- clean_pathway(rownames(top18_mat))

# Compute divergence score: HAE vs mean of others (directional difference)
hae_vals <- top18_mat[, "HAE"]
other_mean <- rowMeans(top18_mat[, colnames(top18_mat) != "HAE"])
divergence <- hae_vals - other_mean

# Group: negative HAE (paradoxical, blue in HAE vs red in others) first,
# then positive HAE; within each group sort by divergence magnitude
grp_neg <- which(hae_vals < 0)
grp_pos <- which(hae_vals >= 0)
# Negative group: most divergent first (most negative divergence = strongest contrast)
ord_neg <- grp_neg[order(divergence[grp_neg])]
# Positive group: most divergent first (most positive divergence = HAE uniquely activated)
ord_pos <- grp_pos[order(divergence[grp_pos], decreasing = TRUE)]
optimal_order <- c(ord_neg, ord_pos)
top18_mat <- top18_mat[optimal_order, ]

ht_H <- Heatmap(top18_mat,
  name = "Delta",
  col = col_nes,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
  column_names_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
  column_names_rot = 0,
  column_names_centered = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.8),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
    legend_height = unit(25, "mm")
  ),
  row_names_max_width = max_text_width(rownames(top18_mat), gp = gpar(fontsize = 8))
)

save_heatmap_panel_pdf("Supp06h_pathway_heatmap.pdf", ht_H, W4_R, H4,
                       title = "Top 18 HAE-dysregulated pathways")

# =============================================================================
# COMPOSITE ASSEMBLY: Pure Vector (cairo_pdf + grid viewports)
# =============================================================================
cat("\n--- Assembling composite SuppFig_06 (Pure Vector, 183x245mm) ---\n")

tryCatch({
  DPI <- ASSEMBLY_DPI
  px_per_mm <- DPI / 25.4
  px_W <- round(W_TOTAL * px_per_mm)
  px_H <- round(H_TOTAL * px_per_mm)

  # --- Pure vector PDF via grid viewports ---
  render_vector_composite <- function() {
    # Master viewport: mm coordinate system, origin bottom-left
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))

    # Row 1: Panel A (left) + Panel B (right)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(H2 + H3 + H4, "mm"),
      width = unit(W1_L, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")))
    print(pA, newpage = FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x = unit(W1_L, "mm"), y = unit(H2 + H3 + H4, "mm"),
      width = unit(W1_R, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")))
    print(pB, newpage = FALSE)
    grid::popViewport()

    # Row 2: Panel C (left) + Panel D (right)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(H3 + H4, "mm"),
      width = unit(W2_L, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")))
    print(pC, newpage = FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x = unit(W2_L, "mm"), y = unit(H3 + H4, "mm"),
      width = unit(W2_R, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")))
    print(pD, newpage = FALSE)
    grid::popViewport()

    # Row 3: Panel E (left, heatmap) + Panel F (right)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(H4, "mm"),
      width = unit(W3_L, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")))
    draw(ht_E, newpage = FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x = unit(W3_L, "mm"), y = unit(H4, "mm"),
      width = unit(W3_R, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")))
    print(pF, newpage = FALSE)
    grid::popViewport()

    # Row 4: Panel G (left) + Panel H (right, heatmap)
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W4_L, "mm"), height = unit(H4, "mm"),
      just = c("left", "bottom")))
    print(pG, newpage = FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      x = unit(W4_L, "mm"), y = unit(0, "mm"),
      width = unit(W4_R, "mm"), height = unit(H4, "mm"),
      just = c("left", "bottom")))
    draw(ht_H, newpage = FALSE)
    grid::popViewport()

    # Panel labels: 12pt bold, 2mm inset from each panel top-left
    labels <- c("a", "b", "c", "d", "e", "f", "g", "h")
    lx_mm <- c(2, W1_L + 2, 2, W2_L + 2, 2, W3_L + 2, 2, W4_L + 2)
    ly_mm <- c(H2+H3+H4+H1-2, H2+H3+H4+H1-2,
               H3+H4+H2-2, H3+H4+H2-2,
               H4+H3-2, H4+H3-2,
               H4-2, H4-2)
    for (i in seq_along(labels)) {
      grid::grid.text(labels[i],
        x = unit(lx_mm[i], "mm"), y = unit(ly_mm[i], "mm"),
        just = c("left", "top"),
        gp = gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY))
    }

    grid::popViewport()
  }

  # PDF (pure vector, AI-editable)
  cairo_pdf(file.path(OUT, "SuppFig_06.pdf"),
            width = W_TOTAL/25.4, height = H_TOTAL/25.4, family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> SuppFig_06.pdf (PURE VECTOR, AI-editable)\n")

  # PNG (raster for review)
  grDevices::png(file.path(OUT, "SuppFig_06.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> SuppFig_06.png\n")

  # TIFF (raster for submission)
  grDevices::tiff(file.path(OUT, "SuppFig_06.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> SuppFig_06.tiff\n")

  # Verify pure vector
  cat("\n  Verifying vector quality...\n")
  pdf_size <- file.info(file.path(OUT, "SuppFig_06.pdf"))$size / 1024
  cat(sprintf("  PDF size: %.1f KB (pure vector expected <300KB)\n", pdf_size))

  cat("\n  SuppFig_06 COMPLETE (pure vector, content-adaptive)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 6 Deep Optimization Complete ===\n")
