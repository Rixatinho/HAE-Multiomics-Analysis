#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_10_standalone.R
# Supplementary Figure 10: ULM TF Activity Inference + L-R Interaction Network
# HAE Multi-omics Study | Journal of Hepatology
# =============================================================================
# Panels (A-J):
#   A = ULM TF activity waterfall (all 24 TFs, FDR threshold)
#   B = TF regulatory circuit summary (dumbbell: z-score vs mean target logFC)
#   C = TF activity vs protein logFC concordance (scatter)
#   D = TF evidence integration heatmap (top 10 TFs × metrics, by category)
#   E = HNF4A/HNF1A target gene expression (DoRothEA targets × DEG logFC)
#   F = Mean TF activity by functional category (grouped dot plot)
#   G = Pathway-level L-R communication (bubble plot)
#   H = Top L-R pairs by effect size (lollipop)
#   I = L-R expression heatmap by pathway (per sample)
#   J = L-R direction by biological function (diverging stacked bar)
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_10_standalone.R
# =============================================================================

cat("=== Supplementary Figure 10: ULM TF Activity + L-R Network ===\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggsci)
  library(ggrepel)
  library(dorothea)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "02_analysis/results")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_10")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants
# =============================================================================
FONT <- "Arial"
MM   <- 25.4
FS_TAG   <- 12
COL_UP   <- "#CD534CFF"
COL_DOWN <- "#0073C2FF"
COL_NS   <- "#868686FF"
PAL_CAT  <- pal_jco("default")(10)

# Functional category colour palette
COL_LIVER   <- "#2166AC"
COL_IMMUNE  <- "#B2182B"
COL_FIBRO   <- "#E08214"
COL_OTHER   <- "#999999"

# =============================================================================
# Theme
# =============================================================================
theme_pub <- function(base_sz = 8) {
  theme_bw(base_size = base_sz, base_family = FONT) %+replace%
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(linewidth = 0.45, color = "black", fill = NA),
      axis.ticks       = element_line(linewidth = 0.3),
      axis.text        = element_text(size = base_sz, color = "black"),
      axis.title       = element_text(size = base_sz + 1, face = "bold"),
      plot.title       = element_text(size = base_sz + 2, face = "bold", hjust = 0),
      legend.text      = element_text(size = base_sz - 1),
      legend.title     = element_text(size = base_sz - 1, face = "bold"),
      legend.key.size  = unit(2.5, "mm"),
      legend.background = element_blank(),
      plot.margin      = margin(4, 4, 4, 4, "mm")
    )
}
theme_set(theme_pub())
ht_opt$message <- FALSE

# =============================================================================
# Helpers
# =============================================================================
save_panel <- function(p, fname, w_mm, h_mm) {
  fpath <- file.path(OUT, fname)
  cairo_pdf(fpath, width = w_mm / MM, height = h_mm / MM, family = FONT)
  if (inherits(p, "ggplot")) print(p) else force(p)
  dev.off()
  cat("  ->", fname, "\n")
}

gp_rn <- function(sz = 7) gpar(fontsize = sz, fontfamily = FONT)
gp_lt <- function(sz = 7) gpar(fontsize = sz, fontfamily = FONT, fontface = "bold")

# =============================================================================
# Load data
# =============================================================================
cat("Loading data...\n")
tf_scores  <- read.csv(file.path(RES, "tf_activity_inference/tf_activity_scores.csv"),
                       stringsAsFactors = FALSE)
tf_circuits <- read.csv(file.path(RES, "tf_activity_inference/tf_regulatory_circuits.csv"),
                        stringsAsFactors = FALSE)
pw_comm     <- read.csv(file.path(RES, "ligand_receptor_network/pathway_communication_summary.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
lr_data     <- read.csv(file.path(RES, "ligand_receptor_network/lr_interaction_results.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
lr_scores   <- read.csv(file.path(RES, "enhancement11_cell_communication/lr_pair_scores.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
deg_data    <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"),
                        stringsAsFactors = FALSE)

# DoRothEA regulons (confidence A/B/C)
data("dorothea_hs")
dorothea_abc <- dorothea_hs[dorothea_hs$confidence %in% c("A", "B", "C"), ]

# Functional category assignment
assign_category <- function(tf_names) {
  liver_tfs   <- c("HNF1A", "HNF4A", "NFE2L2")
  immune_tfs  <- c("STAT1", "IRF7", "NFKB1", "RELA", "IRF1", "STAT6")
  fibrosis_tfs <- c("MYOCD", "SMAD2", "SMAD3", "SRF", "HIF1A")
  ifelse(tf_names %in% liver_tfs, "Liver Identity",
  ifelse(tf_names %in% immune_tfs, "Immune/Inflammatory",
  ifelse(tf_names %in% fibrosis_tfs, "Fibrosis/Remodeling", "Other")))
}

cat("Data loaded successfully.\n")

# =============================================================================
# Panel A: ULM TF activity waterfall
# =============================================================================
cat("  Panel A: TF activity waterfall\n")
tf_a <- tf_scores %>%
  arrange(z_score) %>%
  mutate(
    TF  = factor(TF, levels = TF),
    sig = ifelse(ulm_fdr < 0.05, "FDR < 0.05", "NS"),
    cat = assign_category(TF)
  )

p_a <- ggplot(tf_a, aes(x = TF, y = z_score, fill = sig)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
             linewidth = 0.25, color = "grey55") +
  scale_fill_manual(values = c("FDR < 0.05" = PAL_CAT[1], "NS" = COL_NS),
                    name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "ULM z-score", title = "TF Activity (DoRothEA/ULM)") +
  theme_pub(7) +
  theme(axis.text.y = element_text(face = "italic", size = 5.5),
        legend.position = c(0.22, 0.88),
        legend.key.size = unit(2, "mm"))
save_panel(p_a, "Supp10a_TF_waterfall.pdf", 84, 80)


# =============================================================================
# Panel B: TF regulatory circuit dumbbell (z-score vs mean target logFC)
# =============================================================================
cat("  Panel B: TF regulatory circuit dumbbell\n")
tf_b <- tf_scores %>%
  filter(ulm_fdr < 0.10) %>%
  arrange(z_score) %>%
  mutate(
    TF  = factor(TF, levels = TF),
    cat = assign_category(TF)
  )

p_b <- ggplot(tf_b) +
  # Segment connecting z_score and mean_logFC (scaled for visual)
  geom_segment(aes(x = TF, xend = TF, y = z_score, yend = mean_logFC * 10),
               linewidth = 0.5, color = "grey65") +
  geom_point(aes(x = TF, y = z_score, color = "ULM z-score"),
             size = 2.2, shape = 16) +
  geom_point(aes(x = TF, y = mean_logFC * 10, color = "Target logFC (\u00d710)"),
             size = 2.2, shape = 17) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_color_manual(values = c("ULM z-score" = COL_DOWN,
                                "Target logFC (\u00d710)" = COL_UP),
                     name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Score", title = "TF Regulatory Circuit Activity") +
  theme_pub(7) +
  theme(axis.text.y = element_text(face = "italic", size = 6.5),
        legend.position = c(0.22, 0.88),
        legend.key.size = unit(2, "mm"))
save_panel(p_b, "Supp10b_TF_dumbbell.pdf", 84, 66)


# =============================================================================
# Panel C: TF activity vs protein logFC scatter
# =============================================================================
cat("  Panel C: TF activity vs protein logFC scatter\n")
tf_c <- tf_scores %>%
  filter(!is.na(mean_logFC) & mean_logFC != 0) %>%
  mutate(cat = assign_category(TF))

cor_test <- cor.test(tf_c$z_score, tf_c$mean_logFC, method = "pearson")
cor_lab  <- sprintf("r = %.2f, P = %.2e", cor_test$estimate, cor_test$p.value)

p_c <- ggplot(tf_c, aes(x = z_score, y = mean_logFC)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey55", fill = "grey90",
              linewidth = 0.5, alpha = 0.3) +
  geom_point(aes(color = cat), size = 2.5, alpha = 0.85) +
  geom_text_repel(aes(label = TF), size = 2.4, family = FONT,
                  fontface = "italic", max.overlaps = 15,
                  segment.size = 0.2, segment.color = "grey70") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  annotate("text", x = -Inf, y = Inf, label = cor_lab,
           hjust = -0.05, vjust = 1.5, size = 2.3, family = FONT, color = "grey30") +
  scale_color_manual(values = c("Liver Identity" = COL_LIVER,
                                "Immune/Inflammatory" = COL_IMMUNE,
                                "Fibrosis/Remodeling" = COL_FIBRO,
                                "Other" = COL_OTHER),
                     name = "Category") +
  labs(x = "ULM Activity z-score", y = "Mean Target logFC",
       title = "TF Activity vs Target Expression") +
  theme_pub(7) +
  theme(legend.position = c(0.82, 0.2),
        legend.key.size = unit(2, "mm"))
save_panel(p_c, "Supp10c_TF_vs_protein.pdf", 84, 66)


# =============================================================================
# Panel D: TF evidence integration heatmap (top 10 TFs × metrics)
# =============================================================================
cat("  Panel D: TF evidence heatmap\n")
tf_d <- tf_scores %>%
  arrange(ulm_fdr) %>%
  slice_head(n = 12) %>%
  mutate(cat = assign_category(TF))

# Build matrix: z_score, mean_logFC, n_targets_ratio
mat_d <- data.frame(
  row.names     = tf_d$TF,
  `ULM z`       = tf_d$z_score,
  `Target logFC` = tf_d$mean_logFC,
  `log10(FDR)`  = -log10(tf_d$ulm_fdr + 1e-110),
  check.names   = FALSE
)
# Scale each column to [-1, 1] for balanced colour
mat_d_scaled <- apply(mat_d, 2, function(x) {
  mx <- max(abs(x), na.rm = TRUE)
  if (mx == 0) return(x)
  x / mx
})
rownames(mat_d_scaled) <- rownames(mat_d)

cat_colors <- c("Liver Identity" = COL_LIVER, "Immune/Inflammatory" = COL_IMMUNE,
                "Fibrosis/Remodeling" = COL_FIBRO, "Other" = COL_OTHER)
ra <- rowAnnotation(
  Category = tf_d$cat,
  col = list(Category = cat_colors),
  annotation_name_gp = gpar(fontsize = 6, fontfamily = FONT),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 6, fontfamily = FONT, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5, fontfamily = FONT),
    grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm")
  ),
  width = unit(3, "mm")
)

col_fun_d <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
ht_d <- Heatmap(mat_d_scaled,
  name = "Scaled\nvalue",
  col = col_fun_d,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp  = gpar(fontsize = 7, fontfamily = FONT, fontface = "italic"),
  column_names_gp = gpar(fontsize = 5.5, fontfamily = FONT, fontface = "bold"),
  column_names_rot = 45,
  left_annotation = ra,
  rect_gp = gpar(col = "white", lwd = 0.8),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.1f", mat_d[i, j]),
              x, y, gp = gpar(fontsize = 5.5, fontfamily = FONT,
                               col = ifelse(abs(mat_d_scaled[i, j]) > 0.65, "white", "black")))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontfamily = FONT, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5, fontfamily = FONT),
    legend_height = unit(25, "mm"), grid_width = unit(3, "mm")
  ),
  width = unit(32, "mm"), height = unit(50, "mm")
)
save_panel({
  draw(ht_d, padding = unit(c(4, 4, 4, 4), "mm"))
}, "Supp10d_TF_evidence_heatmap.pdf", 84, 70)


# =============================================================================
# Panel E: HNF4A/HNF1A target gene expression (DoRothEA targets × DEG logFC)
# =============================================================================
cat("  Panel E: HNF4A/HNF1A target gene expression\n")

# Get DoRothEA A/B/C targets for HNF1A and HNF4A
hnf1a_targets <- unique(dorothea_abc$target[dorothea_abc$tf == "HNF1A"])
hnf4a_targets <- unique(dorothea_abc$target[dorothea_abc$tf == "HNF4A"])

# Match with DEG data
hnf1a_deg <- deg_data[deg_data$gene_name %in% hnf1a_targets, c("gene_name", "logFC", "P.Value")]
hnf4a_deg <- deg_data[deg_data$gene_name %in% hnf4a_targets, c("gene_name", "logFC", "P.Value")]

hnf1a_deg$TF <- "HNF1A"
hnf4a_deg$TF <- "HNF4A"
# Remove shared targets from HNF4A to avoid duplicates
shared_genes <- intersect(hnf1a_deg$gene_name, hnf4a_deg$gene_name)
if (length(shared_genes) > 0) hnf4a_deg <- hnf4a_deg[!hnf4a_deg$gene_name %in% shared_genes, ]
e_data <- rbind(hnf1a_deg, hnf4a_deg) %>%
  arrange(logFC) %>%
  mutate(sig = ifelse(P.Value < 0.05, "P < 0.05", "NS"))

# Keep top/bottom 15 for readability
if (nrow(e_data) > 30) {
  e_data <- bind_rows(
    e_data %>% slice_head(n = 15),
    e_data %>% slice_tail(n = 15)
  ) %>% distinct()
}
e_data$gene_name <- factor(e_data$gene_name, levels = unique(e_data$gene_name))

mean_lfc <- mean(e_data$logFC, na.rm = TRUE)
p_e <- ggplot(e_data, aes(x = logFC, y = gene_name, fill = TF)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = mean_lfc, linetype = "dashed", linewidth = 0.3, color = "red3") +
  annotate("text", x = mean_lfc, y = Inf, label = sprintf("mean = %.2f", mean_lfc),
           hjust = -0.1, vjust = 2, size = 2.2, family = FONT, color = "red3") +
  scale_fill_manual(values = c("HNF1A" = "#4393C3", "HNF4A" = "#D6604D"), name = "Regulon") +
  labs(x = "log\u2082FC (Adjacent / Normal)", y = NULL,
       title = "HNF4A/HNF1A Target Gene Expression") +
  theme_pub(7) +
  theme(axis.text.y = element_text(face = "italic", size = 5),
        legend.position = c(0.85, 0.15),
        legend.key.size = unit(2, "mm"))
save_panel(p_e, "Supp10e_HNF_target_genes.pdf", 84, 70)


# =============================================================================
# Panel F: TF activity by functional category (faceted horizontal lollipop)
# =============================================================================
cat("  Panel F: TF by functional category\n")
tf_f <- tf_scores %>%
  mutate(cat = assign_category(TF)) %>%
  filter(cat != "Other") %>%
  arrange(cat, z_score) %>%
  mutate(
    TF  = factor(TF, levels = TF),
    dir = ifelse(z_score > 0, "Activated", "Repressed"),
    cat = factor(cat, levels = c("Liver Identity", "Immune/Inflammatory", "Fibrosis/Remodeling"))
  )

p_f <- ggplot(tf_f, aes(x = TF, y = z_score)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
  geom_segment(aes(xend = TF, y = 0, yend = z_score, color = cat),
               linewidth = 0.5) +
  geom_point(aes(color = cat), size = 2) +
  scale_color_manual(values = c("Liver Identity" = COL_LIVER,
                                "Immune/Inflammatory" = COL_IMMUNE,
                                "Fibrosis/Remodeling" = COL_FIBRO)) +
  facet_wrap(~ cat, scales = "free_y", ncol = 1, strip.position = "right") +
  coord_flip() +
  labs(x = NULL, y = "ULM z-score",
       title = "TF activity by functional category") +
  theme_pub(7) +
  theme(legend.position = "none",
        strip.text.y.right = element_text(size = 6, angle = 0, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA),
        axis.text.y = element_text(face = "italic", size = 6),
        panel.spacing = unit(2, "mm"))
save_panel(p_f, "Supp10f_TF_category_lollipop.pdf", 84, 70)


# =============================================================================
# Panel G: Pathway-level L-R communication (bubble plot)
# =============================================================================
cat("  Panel G: Pathway L-R bubble plot\n")
pw_g <- pw_comm
pw_g$pathway <- rownames(pw_g)
pw_g <- pw_g %>%
  mutate(
    dir = ifelse(mean_LR_score > 0, "Enhanced", "Reduced"),
    abs_score = abs(mean_LR_score),
    n_total = enhanced + reduced + discordant
  ) %>%
  arrange(mean_LR_score) %>%
  mutate(pathway = factor(pathway, levels = pathway))

p_g <- ggplot(pw_g, aes(x = mean_LR_score, y = pathway)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_segment(aes(xend = 0, yend = pathway), linewidth = 0.2, color = "grey75") +
  geom_point(aes(size = n_with_data, color = dir), alpha = 0.8) +
  scale_size_continuous(range = c(1.5, 5), name = "L-R pairs", breaks = c(2, 4, 6, 9)) +
  scale_color_manual(values = c("Enhanced" = COL_UP, "Reduced" = COL_DOWN), name = "Direction") +
  labs(x = "Mean L-R Communication Score", y = NULL,
       title = "Pathway-Level L-R Signalling") +
  theme_pub(7) +
  theme(axis.text.y = element_text(size = 6),
        legend.position = c(0.18, 0.82),
        legend.key.size = unit(2, "mm"),
        legend.spacing.y = unit(1, "mm"))
save_panel(p_g, "Supp10g_pathway_LR_bubble.pdf", 84, 70)


# =============================================================================
# Panel H: Top L-R pairs by effect size (lollipop)
# =============================================================================
cat("  Panel H: Top LR pairs lollipop\n")
lr_h <- lr_data %>%
  mutate(
    pair_name = paste0(Ligand, " \u2192 ", Receptor),
    abs_score = abs(LR_score_RNA)
  ) %>%
  arrange(desc(abs_score)) %>%
  slice_head(n = 20) %>%
  mutate(pair_name = factor(pair_name, levels = rev(pair_name)))

p_h <- ggplot(lr_h, aes(x = pair_name, y = LR_score_RNA)) +
  geom_segment(aes(xend = pair_name, y = 0, yend = LR_score_RNA),
               linewidth = 0.35, color = "grey70") +
  geom_point(aes(color = Direction), size = 2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_color_manual(values = c("Both Up (Enhanced)" = COL_UP,
                                "Both Down (Reduced)" = COL_DOWN,
                                "Discordant" = COL_NS),
                     name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "L-R Score (RNA)",
       title = "Top 20 L-R Pairs by Effect Size") +
  theme_pub(7) +
  theme(axis.text.y = element_text(size = 5.5),
        legend.position = c(0.82, 0.12),
        legend.key.size = unit(2, "mm"))
save_panel(p_h, "Supp10h_top_LR_pairs.pdf", 84, 70)


# =============================================================================
# Panel I: L-R expression heatmap by pathway (per sample)
# =============================================================================
cat("  Panel I: L-R heatmap by pathway\n")

# lr_scores: rows = LR pairs, columns = samples (Normal + Adjacent)
mat_i <- as.matrix(lr_scores)
# Remove rows with all zeros
mat_i <- mat_i[rowSums(mat_i != 0) > 0, ]
# Scale by row for visualisation
mat_i_z <- t(scale(t(mat_i)))
mat_i_z[is.nan(mat_i_z)] <- 0

# Split by Normal vs Adjacent for column annotation
sample_group <- ifelse(grepl("^Normal", colnames(mat_i_z)), "Normal", "Adjacent")
col_ann <- HeatmapAnnotation(
  Group = sample_group,
  col = list(Group = c("Normal" = COL_DOWN, "Adjacent" = COL_UP)),
  annotation_name_gp = gpar(fontsize = 6, fontfamily = FONT),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 6, fontfamily = FONT, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5, fontfamily = FONT),
    grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm")
  ),
  simple_anno_size = unit(3, "mm")
)

# Match LR pair names to pathways
pair_names <- gsub("_", " \u2192 ", rownames(mat_i_z))
# Try to match with lr_data for pathway annotation
lr_pathway_map <- setNames(lr_data$Pathway, paste0(lr_data$Ligand, "_", lr_data$Receptor))
pair_pathways <- lr_pathway_map[rownames(mat_i_z)]
pair_pathways[is.na(pair_pathways)] <- "Other"

# Take top variable pairs if too many
if (nrow(mat_i_z) > 30) {
  rv <- apply(mat_i_z, 1, var, na.rm = TRUE)
  keep <- names(sort(rv, decreasing = TRUE))[1:30]
  mat_i_z <- mat_i_z[keep, ]
  pair_pathways <- pair_pathways[keep]
}

col_fun_i <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
ht_i <- Heatmap(mat_i_z,
  name = "z-score",
  col = col_fun_i,
  top_annotation = col_ann,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  column_split = factor(sample_group, levels = c("Normal", "Adjacent")),
  column_title = NULL,
  column_gap = unit(1.5, "mm"),
  row_names_gp  = gpar(fontsize = 5, fontfamily = FONT),
  column_names_gp = gpar(fontsize = 5, fontfamily = FONT),
  column_names_rot = 60,
  show_column_names = FALSE,
  row_names_max_width = unit(35, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.3),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontfamily = FONT, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5, fontfamily = FONT),
    legend_height = unit(22, "mm"), grid_width = unit(3, "mm")
  ),
  width = unit(50, "mm")
)
save_panel({
  draw(ht_i, padding = unit(c(4, 4, 4, 4), "mm"))
}, "Supp10i_LR_heatmap.pdf", 94, 75)


# =============================================================================
# Panel J: L-R direction by biological function (diverging stacked bar)
# =============================================================================
cat("  Panel J: L-R direction by function\n")
pw_j <- pw_comm
pw_j$pathway <- rownames(pw_j)

# Create diverging data: enhanced goes right, reduced goes left
long_j <- data.frame(
  Pathway = rep(pw_j$pathway, 2),
  Category = rep(c("Enhanced", "Reduced"), each = nrow(pw_j)),
  Count = c(pw_j$enhanced, -pw_j$reduced),
  stringsAsFactors = FALSE
)
long_j$Pathway <- factor(long_j$Pathway,
                         levels = pw_j$pathway[order(pw_j$enhanced - pw_j$reduced)])

p_j <- ggplot(long_j, aes(x = Pathway, y = Count, fill = Category)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = c("Enhanced" = COL_UP, "Reduced" = COL_DOWN), name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Number of L-R Pairs (Enhanced \u2192 | \u2190 Reduced)",
       title = "L-R Direction by Pathway") +
  theme_pub(7) +
  theme(axis.text.y = element_text(size = 5.5),
        legend.position = c(0.85, 0.12),
        legend.key.size = unit(2, "mm"))
save_panel(p_j, "Supp10j_LR_direction_diverging.pdf", 84, 70)


# =============================================================================
# Clean up stale files
# =============================================================================
cat("  Cleaning stale files...\n")
old_files <- list.files(OUT,
  pattern = "^Supp09[a-z]_|^Supp10[a-z]_CellChat_LR_diff\\.pdf$|^Supp10f_TF_category_dotplot\\.pdf$",
  full.names = TRUE)
if (length(old_files) > 0) {
  file.remove(old_files)
  cat("  Removed", length(old_files), "stale file(s)\n")
}


# =============================================================================
# Composite Assembly — 5 rows × 2 cols (183 × 280 mm)
# =============================================================================
cat("\n--- Assembling composite SuppFig_10 (vector, 183x280mm, 600DPI) ---\n")

W_TOTAL <- 183; H_TOTAL <- 280; DPI <- 600
H_ROW <- 56   # 5 rows × 56 = 280
W_L   <- 91; W_R <- W_TOTAL - W_L

render_composite <- function() {
  grid::pushViewport(grid::viewport(
    width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
    xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
  ))

  # Row 1: A (left) + B (right)
  y1 <- 4 * H_ROW
  grid::pushViewport(grid::viewport(x = unit(0, "mm"), y = unit(y1, "mm"),
    width = unit(W_L, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_a, newpage = FALSE); grid::popViewport()

  grid::pushViewport(grid::viewport(x = unit(W_L, "mm"), y = unit(y1, "mm"),
    width = unit(W_R, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_b, newpage = FALSE); grid::popViewport()

  # Row 2: C (left) + D (right)
  y2 <- 3 * H_ROW
  grid::pushViewport(grid::viewport(x = unit(0, "mm"), y = unit(y2, "mm"),
    width = unit(W_L, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_c, newpage = FALSE); grid::popViewport()

  grid::pushViewport(grid::viewport(x = unit(W_L, "mm"), y = unit(y2, "mm"),
    width = unit(W_R, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  draw(ht_d, newpage = FALSE, padding = unit(c(3, 3, 3, 3), "mm")); grid::popViewport()

  # Row 3: E (left) + F (right)
  y3 <- 2 * H_ROW
  grid::pushViewport(grid::viewport(x = unit(0, "mm"), y = unit(y3, "mm"),
    width = unit(W_L, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_e, newpage = FALSE); grid::popViewport()

  grid::pushViewport(grid::viewport(x = unit(W_L, "mm"), y = unit(y3, "mm"),
    width = unit(W_R, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_f, newpage = FALSE); grid::popViewport()

  # Row 4: G (left) + H (right)
  y4 <- 1 * H_ROW
  grid::pushViewport(grid::viewport(x = unit(0, "mm"), y = unit(y4, "mm"),
    width = unit(W_L, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_g, newpage = FALSE); grid::popViewport()

  grid::pushViewport(grid::viewport(x = unit(W_L, "mm"), y = unit(y4, "mm"),
    width = unit(W_R, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_h, newpage = FALSE); grid::popViewport()

  # Row 5: I (left, wider) + J (right)
  y5 <- 0
  grid::pushViewport(grid::viewport(x = unit(0, "mm"), y = unit(y5, "mm"),
    width = unit(W_L + 5, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  draw(ht_i, newpage = FALSE, padding = unit(c(3, 2, 3, 2), "mm")); grid::popViewport()

  grid::pushViewport(grid::viewport(x = unit(W_L + 5, "mm"), y = unit(y5, "mm"),
    width = unit(W_R - 5, "mm"), height = unit(H_ROW, "mm"), just = c("left", "bottom")))
  print(p_j, newpage = FALSE); grid::popViewport()

  # Panel labels (a-j)
  labs <- data.frame(
    lab = letters[1:10],
    x = c(2, W_L+2, 2, W_L+2, 2, W_L+2, 2, W_L+2, 2, W_L+7),
    y = c(y1+H_ROW-2, y1+H_ROW-2, y2+H_ROW-2, y2+H_ROW-2,
          y3+H_ROW-2, y3+H_ROW-2, y4+H_ROW-2, y4+H_ROW-2,
          y5+H_ROW-2, y5+H_ROW-2),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(labs))) {
    grid::grid.text(label = labs$lab[i],
      x = unit(labs$x[i], "mm"), y = unit(labs$y[i], "mm"),
      just = c("left", "top"),
      gp = gpar(fontsize = FS_TAG, fontface = "bold", fontfamily = FONT))
  }
  grid::popViewport()
}

# Output with /tmp intermediary for iCloud safety
tryCatch({
  # PDF
  tmp_pdf <- file.path(tempdir(), "SuppFig_10_tmp.pdf")
  cairo_pdf(tmp_pdf, width = W_TOTAL / MM, height = H_TOTAL / MM, family = FONT)
  render_composite(); dev.off()
  file.copy(tmp_pdf, file.path(OUT, "SuppFig_10.pdf"), overwrite = TRUE)
  cat("  -> SuppFig_10.pdf\n")

  # PNG
  px_W <- round(W_TOTAL * DPI / MM); px_H <- round(H_TOTAL * DPI / MM)
  tmp_png <- file.path(tempdir(), "SuppFig_10_tmp.png")
  grDevices::png(tmp_png, width = px_W, height = px_H, res = DPI, type = "cairo")
  render_composite(); dev.off()
  file.copy(tmp_png, file.path(OUT, "SuppFig_10.png"), overwrite = TRUE)
  cat("  -> SuppFig_10.png\n")

  # TIFF
  tmp_tiff <- file.path(tempdir(), "SuppFig_10_tmp.tiff")
  grDevices::tiff(tmp_tiff, width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_composite(); dev.off()
  file.copy(tmp_tiff, file.path(OUT, "SuppFig_10.tiff"), overwrite = TRUE)
  cat("  -> SuppFig_10.tiff\n")

  cat("  SuppFig_10 DONE (10 panels, vector)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 10 rendering complete ===\n")
