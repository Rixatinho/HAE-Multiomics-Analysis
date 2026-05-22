#!/usr/bin/env Rscript
# =============================================================================
# fix_consensus_panels.R
# Regenerate consensus heatmap + silhouette panels using cairo_pdf + Arial
# Replaces old Helvetica-based panels (Fig5a, Fig5b, Supp3c, Supp3d)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ConsensusClusterPlus)
  library(cluster)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# --- Constants ---
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
DATA <- file.path(BASE, "analysis/data/processed")
RES  <- file.path(BASE, "analysis/results/phase6_subtyping")

FONT_FAMILY <- "Arial"
MM_PER_INCH <- 25.4
W_HALF <- 88
H_STD  <- 85

COL_CS1 <- "#E64B35"
COL_CS2 <- "#3C5488"

mm2in <- function(mm) mm / MM_PER_INCH

cat("=== Regenerating consensus panels (cairo_pdf + Arial) ===\n")

# --- 1. Load & prepare data (same as phase6_subtyping.R) ---
cat("  Loading expression data...\n")
tc_mat <- read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"),
                   row.names = 1, check.names = FALSE)
pr_mat <- read.csv(file.path(DATA, "proteomics_log2_norm.csv"),
                   row.names = 1, check.names = FALSE)
met_mat <- read.csv(file.path(DATA, "metabolomics_log2_merged.csv"),
                    row.names = 1, check.names = FALSE)

adj_samples <- grep("Adjacent|^A\\d", colnames(tc_mat), value = TRUE)
shared <- Reduce(intersect, list(
  adj_samples,
  grep("Adjacent|^A\\d", colnames(pr_mat), value = TRUE),
  grep("Adjacent|^A\\d", colnames(met_mat), value = TRUE)
))

# Scale and concatenate
tc_scaled  <- t(scale(t(tc_mat[, shared])))
pr_scaled  <- t(scale(t(pr_mat[, shared])))
met_scaled <- t(scale(t(met_mat[, shared])))

# Top variable features
tc_var  <- head(order(apply(tc_scaled, 1, var, na.rm = TRUE), decreasing = TRUE), 2000)
pr_var  <- head(order(apply(pr_scaled, 1, var, na.rm = TRUE), decreasing = TRUE), 2000)
met_var <- seq_len(nrow(met_scaled))

rownames(tc_scaled)  <- paste0("TC_", rownames(tc_scaled))
rownames(pr_scaled)  <- paste0("PR_", rownames(pr_scaled))
rownames(met_scaled) <- paste0("MET_", rownames(met_scaled))

concat_mat <- rbind(tc_scaled[tc_var, ], pr_scaled[pr_var, ], met_scaled[met_var, ])
concat_mat[is.na(concat_mat)] <- 0
cat(sprintf("  Concat matrix: %d features x %d samples\n", nrow(concat_mat), ncol(concat_mat)))

# --- 2. Run ConsensusClusterPlus (no native plots) ---
cat("  Running ConsensusClusterPlus (K=2-6, 1000 reps)...\n")
cc_dir <- tempdir()
cc <- ConsensusClusterPlus(
  d = concat_mat, maxK = 6, reps = 1000,
  pItem = 0.8, pFeature = 1,
  clusterAlg = "hc", distance = "pearson",
  innerLinkage = "ward.D2", finalLinkage = "ward.D2",
  seed = 42, title = cc_dir, plot = NULL
)
cat("  ConsensusClusterPlus complete.\n")

# --- 3. Consensus Heatmap (K=2) with ComplexHeatmap + cairo_pdf ---
cat("  Generating consensus heatmap (K=2)...\n")
cons_mat <- cc[[2]]$consensusMatrix
cl <- cc[[2]]$consensusClass

# Use display names (remove prefix)
display_names <- gsub("^Adjacent_", "A", colnames(concat_mat))
display_names <- gsub("^A(\\d+)$", "A\\1", display_names)
rownames(cons_mat) <- display_names
colnames(cons_mat) <- display_names

subtype_vec <- paste0("CS", cl)
names(subtype_vec) <- display_names

col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "#B8CCE4", "#3C5488"))
ha <- HeatmapAnnotation(
  Subtype = subtype_vec,
  col = list(Subtype = c("CS1" = COL_CS1, "CS2" = COL_CS2)),
  simple_anno_size = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_FAMILY),
    labels_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY)
  )
)

ht <- Heatmap(
  cons_mat, name = "Consensus",
  col = col_fun,
  top_annotation = ha,
  column_title = "Consensus Matrix (K=2, 1000 iterations)",
  column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = FONT_FAMILY),
  row_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
  column_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_FAMILY),
    labels_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY)
  )
)

# Save to both Figure_5 and SuppFig_03 (updated figure numbers)
for (out_info in list(
  list(dir = "Figure_5",  name = "Fig5a_consensus_K2.pdf"),
  list(dir = "SuppFig_03", name = "Supp3c_consensus_heatmap.pdf")
)) {
  out_dir <- file.path(BASE, "图片", out_info$dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, out_info$name)
  cairo_pdf(out_path, width = mm2in(W_HALF), height = mm2in(H_STD), family = FONT_FAMILY)
  draw(ht)
  dev.off()
  cat(sprintf("  -> %s/%s\n", out_info$dir, out_info$name))
}

# --- 4. Silhouette Plot (K=2) with ggplot2 + cairo_pdf ---
cat("  Generating silhouette plot (K=2)...\n")
d <- as.dist(1 - cons_mat)
sil <- silhouette(cl, d)
sil_df <- data.frame(
  sample = display_names[sil[, 1]],
  cluster = factor(paste0("CS", sil[, "cluster"])),
  sil_width = sil[, "sil_width"]
)

# Sort by cluster then by silhouette width (descending)
sil_df <- sil_df[order(sil_df$cluster, -sil_df$sil_width), ]
sil_df$order <- seq_len(nrow(sil_df))
mean_sil <- mean(sil_df$sil_width)

theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7, face = "bold"),
    legend.key.size = unit(2.5, "mm"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

p_sil <- ggplot(sil_df, aes(x = order, y = sil_width, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_hline(yintercept = mean_sil, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  scale_fill_manual(values = c("CS1" = COL_CS1, "CS2" = COL_CS2)) +
  coord_flip() +
  labs(
    title = sprintf("Silhouette Plot (K=2, mean=%.3f)", mean_sil),
    x = "Samples", y = "Silhouette Width", fill = "Subtype"
  ) +
  scale_x_continuous(breaks = sil_df$order, labels = sil_df$sample) +
  theme_nc +
  theme(axis.text.y = element_text(size = 6))

for (out_info in list(
  list(dir = "Figure_5",  name = "Fig5b_silhouette_K2.pdf"),
  list(dir = "SuppFig_03", name = "Supp3d_silhouette_plot.pdf")
)) {
  out_dir <- file.path(BASE, "图片", out_info$dir)
  out_path <- file.path(out_dir, out_info$name)
  cairo_pdf(out_path, width = mm2in(W_HALF), height = mm2in(H_STD), family = FONT_FAMILY)
  print(p_sil)
  dev.off()
  cat(sprintf("  -> %s/%s\n", out_info$dir, out_info$name))
}

cat("\n=== All 4 consensus panels regenerated with cairo_pdf + Arial ===\n")
