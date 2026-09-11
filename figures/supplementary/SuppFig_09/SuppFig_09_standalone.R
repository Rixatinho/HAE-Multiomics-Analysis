# ============================================================
# SuppFig_09_standalone.R
# Gene regulatory network, WGCNA, and extended deconvolution
# analysis — 12 panels (a–l)
# Target journal : EBioMedicine (Lancet family)  183 × 245 mm
# Run: conda run -n multiomics Rscript SuppFig_09_standalone.R
# ============================================================
FS_TAG <- 12   # font size for panel tags

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
  library(ggrepel)
  library(cowplot)
  library(grid)
  library(gridExtra)
  library(igraph)
  library(scales)
})

# ── 路径 ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
if (length(file_arg)) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", file_arg)))
} else {
  SCRIPT_DIR <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学/04_figures/supplementary/SuppFig_09"
}
OUT  <- SCRIPT_DIR
DATA <- file.path(SCRIPT_DIR, "../../../02_analysis/results")

# ── 配色常量 ───────────────────────────────────────────────────
COL_UP   <- "#CD534CFF"
COL_DOWN <- "#0073C2FF"
COL_NS   <- "#868686FF"
PAL_CAT  <- pal_jco("default")(10)

# WGCNA module colors
PR_MOD_COLS <- c("0"="#BEBEBE","1"="#00BCD4","2"="#2196F3",
                 "3"="#795548","4"="#FFEB3B","5"="#4CAF50")
TC_MOD_COLS <- c("black"="#212121","red"="#E53935","turquoise"="#00BCD4")
PR_HUB_COLS <- c("blue"="#2196F3","brown"="#795548","green"="#4CAF50",
                 "turquoise"="#00BCD4","yellow"="#FFC107")

# ── 基础主题 ───────────────────────────────────────────────────
# macOS: "Arial" = TrueType font embedded by cairo_pdf
FONT <- "Arial"

theme_nc <- function(base_size = 8) {
  theme_classic(base_size = base_size) %+replace% theme(
    text            = element_text(family = FONT),
    axis.text       = element_text(size = base_size, color = "black"),
    axis.title      = element_text(size = base_size + 1),
    legend.text     = element_text(size = base_size - 1),
    legend.title    = element_text(size = base_size, face = "bold"),
    legend.key.size = unit(3, "mm"),
    plot.title      = element_text(size = base_size + 1, face = "bold", hjust = 0),
    strip.text      = element_text(size = base_size, face = "bold"),
    panel.border    = element_blank(),
    axis.line       = element_line(color = "black", linewidth = 0.4),
    plot.background = element_blank(),
    panel.background = element_blank()
  )
}

mm2in <- function(x) x / 25.4

# ── Helper: save ggplot panel ─────────────────────────────────
save_panel <- function(p, filename, w_mm, h_mm) {
  ggsave(file.path(OUT, filename), plot = p,
         width = mm2in(w_mm), height = mm2in(h_mm),
         device = cairo_pdf)
  invisible(p)
}

# ── Helper: network builder (ggplot2 hub-spoke) ───────────────
make_network_gg <- function(top_tfs, all_links, n_per_tf = 4,
                             node_col_tf = COL_UP, node_col_tg = COL_DOWN,
                             title = "") {
  sub_e <- all_links %>%
    filter(TF %in% top_tfs) %>%
    group_by(TF) %>%
    arrange(desc(Weight)) %>%
    slice_head(n = n_per_tf) %>%
    ungroup()

  targets <- unique(sub_e$Target)
  n_tf  <- length(top_tfs)
  n_tg  <- length(targets)

  # TF nodes on inner ring
  ang_tf <- seq(0, 2 * pi * (1 - 1/n_tf), length.out = n_tf)
  tf_nodes <- data.frame(
    name  = top_tfs,
    type  = "TF",
    x     = 0.35 * cos(ang_tf),
    y     = 0.35 * sin(ang_tf),
    stringsAsFactors = FALSE
  )
  # Target nodes on outer ring
  ang_tg <- seq(0, 2 * pi * (1 - 1/n_tg), length.out = n_tg)
  tg_nodes <- data.frame(
    name  = targets,
    type  = "Target",
    x     = 0.85 * cos(ang_tg),
    y     = 0.85 * sin(ang_tg),
    stringsAsFactors = FALSE
  )
  nodes <- bind_rows(tf_nodes, tg_nodes)

  edge_df <- sub_e %>%
    left_join(nodes[, c("name","x","y")], by = c("TF" = "name")) %>%
    rename(x0 = x, y0 = y) %>%
    left_join(nodes[, c("name","x","y")], by = c("Target" = "name")) %>%
    rename(x1 = x, y1 = y)

  ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x0, y = y0, xend = x1, yend = y1, alpha = Weight),
                 color = "grey60", linewidth = 0.3,
                 arrow = arrow(length = unit(0.12, "cm"), type = "open")) +
    geom_point(data = nodes,
               aes(x, y, color = type, size = type), show.legend = FALSE) +
    geom_text(data = tf_nodes,
              aes(x, y, label = name), fontface = "bold",
              size = 2.2, family = FONT, color = "black") +
    geom_text_repel(data = tg_nodes,
                    aes(x, y, label = name),
                    size = 1.8, family = FONT, color = "grey20",
                    max.overlaps = 20, segment.size = 0.2,
                    box.padding = 0.1) +
    scale_color_manual(values = c("TF" = node_col_tf, "Target" = node_col_tg)) +
    scale_size_manual(values  = c("TF" = 4.5,         "Target" = 1.8)) +
    scale_alpha_continuous(range = c(0.3, 0.9)) +
    ggtitle(title) +
    theme_void(base_size = 8) %+replace%
    theme(text = element_text(family = FONT),
          plot.title = element_text(size = 8, face = "bold",
                                    hjust = 0.5, margin = margin(b = 2)),
          plot.background = element_blank(),
          legend.position = "none")
}

# =============================================================
# DATA LOADING
# =============================================================
message("Loading data ...")

pr_cor   <- read.csv(file.path(DATA, "enhancement2_wgcna/PR_module_trait_cor.csv"),    row.names = 1, check.names = FALSE)
pr_pval  <- read.csv(file.path(DATA, "enhancement2_wgcna/PR_module_trait_pval_BH.csv"), row.names = 1, check.names = FALSE)
tc_hub   <- read.csv(file.path(DATA, "enhancement2_wgcna/TC_hub_genes.csv"),   check.names = FALSE)
pr_hub   <- read.csv(file.path(DATA, "enhancement2_wgcna/PR_hub_genes.csv"),   check.names = FALSE)
ml_evid  <- read.csv(file.path(DATA, "enhancement12_integrative/multi_layer_evidence_matrix.csv"), check.names = FALSE)
grn_links  <- read.csv(file.path(DATA, "enhancement15_grn/grn_significant_links.csv"),          check.names = FALSE)
grn_tf     <- read.csv(file.path(DATA, "enhancement15_grn/grn_tf_differential_activity.csv"),   check.names = FALSE)
grn_cmp    <- read.csv(file.path(DATA, "enhancement15_grn/grn_condition_comparison.csv"),        check.names = FALSE)
grn_tf_sub <- read.csv(file.path(DATA, "enhancement15_grn/grn_tf_subtype_activity.csv"),        check.names = FALSE)
cellchat   <- read.csv(file.path(DATA, "enhancement16_deconvolution/cellchat_LR_comparison.csv"), check.names = FALSE)
deconv_val <- read.csv(file.path(DATA, "enhancement16_deconvolution/bayesprism_validation_vs_ssGSEA.csv"), check.names = FALSE)
bp_props   <- read.csv(file.path(DATA, "enhancement16_deconvolution/bayesprism_cell_proportions.csv"), row.names = 1, check.names = FALSE)
gsva_hm    <- read.csv(file.path(DATA, "celltype_pathway_correlation/gsva_hallmark_scores_per_sample.csv"), row.names = 1, check.names = FALSE)

message("  Data loaded. Starting panel generation ...")

# =============================================================
# PANEL a — PR WGCNA Module-Trait Correlation Heatmap
# =============================================================
message("Panel a: PR WGCNA module-trait heatmap")

# Rename rows: module number → color name
pr_mod_map <- c("0"="Grey","1"="Turquoise","2"="Blue","3"="Brown","4"="Yellow","5"="Green")
mat_a  <- as.matrix(pr_cor)
pmat_a <- as.matrix(pr_pval)
rownames(mat_a)  <- pr_mod_map[rownames(mat_a)]
rownames(pmat_a) <- pr_mod_map[rownames(pmat_a)]

# Significance label function
sig_label <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

col_fun_a <- colorRamp2(c(-0.4, 0, 0.4),
                         c("#1565C0", "white", "#C62828"))

# Row annotation: module color bar
row_ann_col <- list(Module = setNames(
  c("#BEBEBE","#00BCD4","#2196F3","#795548","#FFEB3B","#4CAF50"),
  c("Grey","Turquoise","Blue","Brown","Yellow","Green")))
row_ann <- rowAnnotation(
  Module = rownames(mat_a),
  col    = row_ann_col,
  annotation_legend_param = list(Module = list(title_gp = gpar(fontsize = 7, fontface = "bold"),
                                               labels_gp = gpar(fontsize = 6))),
  show_annotation_name = FALSE,
  width = unit(3, "mm"))

ht_a <- Heatmap(
  mat_a,
  name          = "Pearson r",
  col           = col_fun_a,
  cell_fun      = function(j, i, x, y, w, h, fill) {
    lbl <- sig_label(pmat_a[i, j])
    if (nchar(lbl) > 0)
      grid.text(lbl, x, y, gp = gpar(fontsize = 7))
  },
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  show_row_dend     = FALSE,
  show_column_dend  = FALSE,
  row_names_gp      = gpar(fontsize = 7, fontfamily = FONT),
  column_names_gp   = gpar(fontsize = 6.5, fontfamily = FONT),
  column_names_rot  = 45,
  heatmap_legend_param = list(
    title_gp         = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT),
    labels_gp        = gpar(fontsize = 6, fontfamily = FONT),
    legend_direction = "horizontal",
    legend_width     = unit(25, "mm"),
    legend_height    = unit(3,  "mm")),
  border            = TRUE,
  rect_gp           = gpar(col = "white", lwd = 0.5),
  width             = unit(52, "mm"),
  height            = unit(30, "mm")
)

# Save individual panel PDF
cairo_pdf(file.path(OUT, "Supp09a_PR_module_trait.pdf"),
          width = mm2in(61), height = mm2in(65))
draw(ht_a, padding = unit(c(4,4,4,4), "mm"), heatmap_legend_side = "bottom")
dev.off()

# =============================================================
# PANEL b — Multi-Layer Evidence Matrix
# =============================================================
message("Panel b: Multi-layer evidence matrix")

evid_b <- ml_evid %>%
  mutate(
    pathway_short = gsub("_", " ", pathway),
    # Shorter category labels to fit in 61mm panel
    category      = factor(category,
                           levels = c("Cell Death","Immune Checkpoint","Metabolic Enzyme"),
                           labels = c("Cell Death","Immune Chkpt.","Metabolic Enzyme")),
    sig_label     = ifelse(ssGSEA_P < 0.05, "*", "")
  ) %>%
  arrange(category, desc(abs(ssGSEA_ES)))

evid_b$pathway_short <- factor(evid_b$pathway_short,
                                levels = rev(evid_b$pathway_short))

# Long format for faceted tiles
evid_long <- evid_b %>%
  select(pathway_short, category, direction, ssGSEA_ES,
         cross_immune_cor, cross_metab_cor, lr_enriched) %>%
  pivot_longer(
    cols     = c(ssGSEA_ES, cross_immune_cor, cross_metab_cor, lr_enriched),
    names_to = "evidence",
    values_to = "value"
  ) %>%
  mutate(
    evidence = factor(evidence,
                      levels  = c("ssGSEA_ES","cross_immune_cor",
                                  "cross_metab_cor","lr_enriched"),
                      labels  = c("ssGSEA ES","Immune cor.",
                                  "Metab. cor.","LR enr."))
  )

# Category color strip — updated labels to match short versions
cat_cols <- c("Cell\nDeath"      = PAL_CAT[1],
              "Immune\nChkpt."   = PAL_CAT[3],
              "Metabolic\nEnzyme"= PAL_CAT[5])

p_b <- ggplot(evid_long, aes(x = evidence, y = pathway_short)) +
  geom_tile(aes(fill = value), color = "white", linewidth = 0.3) +
  geom_text(data = evid_long %>% filter(evidence == "ssGSEA ES"),
            aes(label = ifelse(value > 0, "↑", "↓"), color = value > 0),
            size = 2, fontface = "bold") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP,
                       midpoint = 0, name = "Score",
                       breaks = c(-0.5, 0, 0.5),
                       guide = guide_colorbar(barheight = unit(12, "mm"),
                                              barwidth  = unit(2, "mm"))) +
  scale_color_manual(values = c("TRUE" = COL_UP, "FALSE" = COL_DOWN),
                     guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_nc(8) +
  theme(
    axis.text.x          = element_text(size = 7, angle = 45,
                                        hjust = 1, vjust = 1),
    axis.text.y          = element_text(size = 6),
    strip.text.y.left    = element_text(size = 6.5, angle = 90, hjust = 0.5,
                                        face = "bold"),
    strip.placement      = "outside",
    panel.border         = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
    axis.line            = element_blank(),
    plot.margin          = margin(2, 2, 2, 2)
  )

save_panel(p_b, "Supp09b_multi_layer_evidence.pdf", 61, 65)

# =============================================================
# PANEL c — TC WGCNA Hub Gene Scatter
# =============================================================
message("Panel c: TC hub gene scatter")

top_tc <- tc_hub %>%
  group_by(module) %>%
  slice_max(hub_score, n = 3) %>%
  ungroup()

p_c <- ggplot(tc_hub, aes(x = abs_kME, y = abs_GS, color = module)) +
  geom_point(alpha = 0.75, size = 1.5) +
  geom_text_repel(data = top_tc,
                  aes(label = gene),
                  size = 2.2, fontface = "italic",
                  family = FONT, max.overlaps = 15,
                  segment.size = 0.25, box.padding = 0.3) +
  scale_color_manual(values = TC_MOD_COLS, name = "Module") +
  labs(x = "Module membership (|kME|)",
       y = "Gene significance (|GS|)",
       title = "TC hub genes") +
  theme_nc(8) +
  guides(color = guide_legend(override.aes = list(size = 2)))

save_panel(p_c, "Supp09c_TC_hub_scatter.pdf", 61, 65)

# =============================================================
# PANEL d — PR WGCNA Hub Gene Scatter
# =============================================================
message("Panel d: PR hub gene scatter")

top_pr <- pr_hub %>%
  group_by(module) %>%
  slice_max(hub_score, n = 3) %>%
  ungroup()

p_d <- ggplot(pr_hub, aes(x = abs_kME, y = abs_GS, color = module)) +
  geom_point(alpha = 0.75, size = 1.5) +
  geom_text_repel(data = top_pr,
                  aes(label = gene),
                  size = 2.2, fontface = "italic",
                  family = FONT, max.overlaps = 15,
                  segment.size = 0.25, box.padding = 0.3) +
  scale_color_manual(values = PR_HUB_COLS, name = "Module") +
  labs(x = "Module membership (|kME|)",
       y = "Gene significance (|GS|)",
       title = "PR hub genes") +
  theme_nc(8) +
  guides(color = guide_legend(override.aes = list(size = 2)))

save_panel(p_d, "Supp09d_PR_hub_scatter.pdf", 61, 60)

# =============================================================
# PANEL e — Peri-lesional Condition-Specific Network
# =============================================================
message("Panel e: Peri-lesional regulatory network")

# Top TFs enriched in Adjacent (Diff > 0, highest Degree_Adjacent)
tf_adj <- grn_tf %>%
  filter(Diff > 0) %>%
  arrange(desc(Degree_Adjacent)) %>%
  slice_head(n = 6) %>%
  pull(TF)

p_e <- make_network_gg(
  top_tfs    = tf_adj,
  all_links  = grn_links,
  n_per_tf   = 4,
  node_col_tf = COL_UP,
  title      = "Peri-lesional network"
)

save_panel(p_e, "Supp09e_peri_lesional_network.pdf", 61, 60)

# =============================================================
# PANEL f — Normal Condition-Specific Network
# =============================================================
message("Panel f: Normal regulatory network")

# Top TFs enriched in Normal (Diff < 0, highest Degree_Normal)
tf_nor <- grn_tf %>%
  filter(Diff < 0) %>%
  arrange(desc(Degree_Normal)) %>%
  slice_head(n = 6) %>%
  pull(TF)

p_f <- make_network_gg(
  top_tfs    = tf_nor,
  all_links  = grn_links,
  n_per_tf   = 4,
  node_col_tf = COL_DOWN,
  title      = "Normal tissue network"
)

save_panel(p_f, "Supp09f_normal_network.pdf", 61, 60)

# =============================================================
# PANEL g — Regulatory Edge Overlap Venn Diagram
# =============================================================
message("Panel g: Regulatory edge overlap Venn")

venn_vals <- setNames(as.numeric(grn_cmp$Value), grn_cmp$Metric)
adj_only  <- venn_vals["Adjacent_Only"]
nor_only  <- venn_vals["Normal_Only"]
shared    <- venn_vals["Shared"]

# Draw two overlapping circles using geom_path
make_circle_df <- function(x0, y0, r, n = 200, id) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(theta), y = y0 + r * sin(theta), id = id)
}

r      <- 1.2
offset <- 0.9
circ <- bind_rows(
  make_circle_df(-offset/2, 0, r, id = "A"),
  make_circle_df( offset/2, 0, r, id = "N")
)

label_df <- data.frame(
  x     = c(-offset/2 - r/2,   0,            offset/2 + r/2),
  y     = c(0,                  0,             0),
  label = c(format(adj_only, big.mark=","),
            format(shared, big.mark=","),
            format(nor_only, big.mark=",")),
  sub   = c("Peri-lesional\nonly","Shared","Normal\nonly")
)

p_g <- ggplot() +
  geom_polygon(data = circ %>% filter(id == "A"),
               aes(x, y), fill = COL_UP,   alpha = 0.25) +
  geom_polygon(data = circ %>% filter(id == "N"),
               aes(x, y), fill = COL_DOWN, alpha = 0.25) +
  geom_path(data = circ %>% filter(id == "A"),
            aes(x, y), color = COL_UP,   linewidth = 0.8) +
  geom_path(data = circ %>% filter(id == "N"),
            aes(x, y), color = COL_DOWN, linewidth = 0.8) +
  geom_text(data = label_df,
            aes(x, y + 0.2, label = label),
            size = 3.5, fontface = "bold", family = FONT) +
  geom_text(data = label_df,
            aes(x, y - 0.25, label = sub),
            size = 2.2, family = FONT, color = "grey30",
            lineheight = 0.9) +
  annotate("text", x = -offset/2, y = r + 0.35,
           label = "Peri-lesional",
           size = 2.5, color = COL_UP, fontface = "bold") +
  annotate("text", x =  offset/2, y = r + 0.35,
           label = "Normal",
           size = 2.5, color = COL_DOWN, fontface = "bold") +
  ggtitle("Regulatory edge overlap") +
  coord_fixed() +
  theme_void(base_size = 8) %+replace%
  theme(text = element_text(family = FONT),
        plot.title = element_text(size = 8, face = "bold",
                                  hjust = 0.5, margin = margin(b = 2)),
        plot.background = element_blank())

save_panel(p_g, "Supp09g_edge_venn.pdf", 61, 60)

# =============================================================
# PANEL h — Differential TF Activity Heatmap
# =============================================================
message("Panel h: Differential TF activity heatmap")

# Top 16 TFs by |Diff|
top_tf_h <- grn_tf %>%
  arrange(desc(abs(Diff))) %>%
  slice_head(n = 16)

# Pivot to long format: two "conditions"
tf_long <- top_tf_h %>%
  select(TF, Degree_Adjacent, Degree_Normal, Diff) %>%
  pivot_longer(cols = c(Degree_Adjacent, Degree_Normal),
               names_to = "Condition",
               values_to = "Degree") %>%
  mutate(
    Condition = factor(Condition,
                       levels = c("Degree_Adjacent","Degree_Normal"),
                       labels = c("Peri-lesional","Normal")),
    TF = factor(TF, levels = rev(top_tf_h$TF))
  )

p_h <- ggplot(tf_long, aes(x = Condition, y = TF, fill = Degree)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Degree), size = 2.2, family = FONT) +
  scale_fill_gradient(low = "white", high = COL_UP,
                      name = "Degree",
                      guide = guide_colorbar(barheight = unit(12,"mm"),
                                             barwidth  = unit(2,"mm"))) +
  labs(x = NULL, y = NULL,
       title = "TF regulatory degree") +
  theme_nc(8) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 7.5),
    axis.text.y  = element_text(size = 6.5, face = "italic"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
    axis.line    = element_blank(),
    plot.margin  = margin(2, 4, 2, 2)
  )

save_panel(p_h, "Supp09h_TF_differential_activity.pdf", 61, 60)

# =============================================================
# PANEL i — GRN TF Differential Regulatory Activity (Adjacent vs Normal)
# =============================================================
message("Panel i: GRN TF differential regulatory activity")

top_tf_diff <- grn_tf %>%
  arrange(desc(abs(Diff))) %>%
  slice_head(n = 20) %>%
  mutate(
    Direction = ifelse(Diff > 0, "Up in peri-lesional", "Down in peri-lesional"),
    TF        = factor(TF, levels = rev(TF))
  )

COL_UP_I   <- "#E07B54"   # warm orange
COL_DOWN_I <- "#4A90A4"   # teal-blue

p_i <- ggplot(top_tf_diff, aes(x = Diff, y = TF, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.65, show.legend = TRUE) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
  geom_text(aes(label = TF,
                x     = ifelse(Diff > 0, -1.5, 1.5),
                hjust = ifelse(Diff > 0, 1, 0)),
            size = 2.1, family = FONT, color = "grey15") +
  scale_fill_manual(
    values = c("Up in peri-lesional"   = COL_UP_I,
               "Down in peri-lesional" = COL_DOWN_I),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  labs(x = "Regulatory degree change (Adjacent \u2212 Normal)", y = NULL) +
  theme_nc(8) +
  theme(
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    legend.position = c(0.28, 0.08),
    legend.key.size = unit(3, "mm"),
    legend.text     = element_text(size = 6, family = FONT)
  )

save_panel(p_i, "Supp09i_GRN_TF_differential_activity.pdf", 61, 60)

# =============================================================
# PANEL j — Cell Communication: Normal vs Peri-lesional slope chart
# =============================================================
message("Panel j: Cell communication strength comparison")

cc_path <- cellchat %>%
  group_by(pathway) %>%
  summarise(Normal   = mean(Normal,   na.rm = TRUE),
            Adjacent = mean(Adjacent, na.rm = TRUE),
            .groups  = "drop") %>%
  filter(!is.na(Normal) & !is.na(Adjacent)) %>%
  arrange(desc(abs(Normal - Adjacent))) %>%
  mutate(
    Direction = ifelse(Adjacent > Normal, "Up in peri-lesional", "Down in peri-lesional"),
    pathway   = factor(pathway, levels = pathway)
  )

cc_long <- cc_path %>%
  pivot_longer(cols = c(Normal, Adjacent),
               names_to  = "Condition",
               values_to = "Strength") %>%
  mutate(Condition = factor(Condition,
                            levels = c("Normal", "Adjacent"),
                            labels = c("Normal", "Peri-lesional")))

n_cc    <- nrow(cc_path)
cc_cols <- setNames(colorRampPalette(PAL_CAT)(n_cc), levels(cc_path$pathway))

p_j <- ggplot(cc_long,
              aes(x = Condition, y = Strength, group = pathway)) +
  geom_line(aes(color = pathway), linewidth = 0.8, alpha = 0.85) +
  geom_point(aes(color = pathway, shape = Condition), size = 2.5) +
  geom_text_repel(
    data = cc_long %>% filter(Condition == "Normal"),
    aes(label = pathway, color = pathway),
    hjust = 1.1, size = 2.2, direction = "y",
    family = FONT, segment.size = 0.2,
    max.overlaps = 20, show.legend = FALSE
  ) +
  scale_color_manual(values = cc_cols, guide = "none") +
  scale_shape_manual(values = c("Normal" = 16, "Peri-lesional" = 17),
                     name = "Condition") +
  scale_x_discrete(expand = expansion(add = c(1.5, 0.3))) +
  labs(x = NULL, y = "Mean communication score",
       title = "Cell communication by pathway") +
  theme_nc(8) +
  theme(legend.position    = "bottom",
        legend.key.size    = unit(3, "mm"),
        axis.text.x        = element_text(size = 8, face = "bold"))

save_panel(p_j, "Supp09j_cell_communication_chord.pdf", 61, 60)
message("  Panel j saved.")

# =============================================================
# PANEL k — Deconvolution Validation Heatmap
# =============================================================
message("Panel k: Deconvolution validation")

dv <- deconv_val %>%
  mutate(
    sig      = ifelse(pvalue < 0.05, "p < 0.05", "ns"),
    lab_pair = paste0(bayesprism_type, " / ", ssgsea_type),
    r_label  = sprintf("%.2f", correlation)
  ) %>%
  arrange(desc(correlation))

dv$lab_pair <- factor(dv$lab_pair, levels = rev(dv$lab_pair))

p_k <- ggplot(dv, aes(x = correlation, y = lab_pair,
                       fill = correlation, color = sig)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  geom_point(aes(x = correlation + 0.02, shape = sig),
             size = 2, show.legend = TRUE) +
  geom_text(aes(x = 0.01, label = r_label),
            hjust = 0, size = 2.2, color = "white", family = FONT) +
  scale_fill_gradient2(low = COL_DOWN, mid = "#F5F5F5", high = COL_UP,
                       midpoint = 0.5) +
  scale_color_manual(values = c("p < 0.05" = "black", "ns" = COL_NS),
                     guide = "none") +
  scale_shape_manual(values = c("p < 0.05" = 8, "ns" = 4),
                     name = "Significance") +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(x = "Pearson correlation", y = NULL,
       title = "Deconvolution validation") +
  theme_nc(8) +
  theme(
    axis.text.y   = element_text(size = 6),
    legend.position = c(0.82, 0.18),
    legend.key.size = unit(2.5, "mm"),
    legend.background = element_blank()
  )

save_panel(p_k, "Supp09k_deconvolution_validation.pdf", 61, 60)

# =============================================================
# PANEL l — Stellate Cell Proportion vs Myogenesis Score
# =============================================================
message("Panel l: Stellate cell vs Myogenesis pathway scatter")

# Extract stellate cell proportion and Myogenesis GSVA score
stell_vec  <- setNames(bp_props[, "Stellate_cell"], rownames(bp_props))
myog_row   <- "MYOGENESIS"
myog_vec   <- setNames(as.numeric(gsva_hm[myog_row, ]), colnames(gsva_hm))
common_l   <- intersect(names(stell_vec), names(myog_vec))
df_l <- data.frame(
  sample   = common_l,
  Stellate = stell_vec[common_l],
  Myog     = myog_vec[common_l],
  Group    = ifelse(grepl("^Normal", common_l), "Normal", "Peri-lesional")
)

cr_l <- cor.test(df_l$Myog, df_l$Stellate, method = "pearson")
r_lab_l <- sprintf("r = %.2f, p = %.1e", cr_l$estimate, cr_l$p.value)

p_l <- ggplot(df_l, aes(x = Myog, y = Stellate, color = Group)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1),
              color = "grey40", fill = "grey85", linewidth = 0.6, alpha = 0.3) +
  geom_point(size = 2.5, alpha = 0.85) +
  annotate("text", x = Inf, y = Inf, label = r_lab_l,
           hjust = 1.05, vjust = 1.4, size = 2.5, family = FONT, color = "grey20") +
  scale_color_manual(values = c(Normal = COL_DOWN, "Peri-lesional" = COL_UP),
                     name = "Tissue") +
  labs(x = "Myogenesis score (GSVA)",
       y = "Stellate cell proportion",
       title = "Stellate cell vs Myogenesis") +
  theme_nc(8) +
  theme(legend.position = c(0.82, 0.15))

save_panel(p_l, "Supp09l_stellate_myogenesis_scatter.pdf", 61, 60)

# =============================================================
# COMPOSITE ASSEMBLY — 183 × 245 mm  (pure-vector cairo_pdf)
# Layout: 4 rows × 3 cols
#   Row heights: 65 / 60 / 60 / 60 mm  (total = 245 mm)
#   Col widths:  61 / 61 / 61 mm        (total = 183 mm)
# All panels rendered natively via grid viewports → zero raster
# Verification: pdfimages -list SuppFig_09.pdf → empty (no images)
# =============================================================
message("Assembling composite SuppFig_09 (pure-vector grid viewport) ...")

FONT_FAMILY <- FONT
W_mm  <- 183; H_mm  <- 245
ROW_H <- c(65, 60, 60, 60)   # mm, must sum to H_mm
COL_W <- c(61, 61, 61)       # mm, must sum to W_mm

# Row 4 uses asymmetric widths: j=47mm, k=75mm, l=61mm (total=183mm)
ROW4_X <- c(0, 47, 122)
ROW4_W <- c(47, 75, 61)

# Column left edges and row bottom edges (mm from bottom of outer viewport)
x_left <- c(0, COL_W[1], COL_W[1]+COL_W[2])
y_bot  <- c(
  H_mm - ROW_H[1],                                    # row 1 bottom = 180 mm
  H_mm - ROW_H[1] - ROW_H[2],                         # row 2 bottom = 120 mm
  H_mm - ROW_H[1] - ROW_H[2] - ROW_H[3],              # row 3 bottom =  60 mm
  0                                                    # row 4 bottom =   0 mm
)

# Panel label positions: rows 1-3 uniform cols; row 4 asymmetric
lbl_x <- c(rep(x_left, times = 3) + 2, ROW4_X + 2)
lbl <- data.frame(
  lab  = letters[1:12],
  x_mm = lbl_x,
  y_mm = rep(y_bot + ROW_H - 2, each = 3),
  stringsAsFactors = FALSE
)

# ── Master render function (called for every output format) ──
render_composite <- function() {
  grid.newpage()

  # Outer viewport — provides mm coordinate space for absolute positioning
  grid::pushViewport(grid::viewport(
    width  = unit(W_mm, "mm"),
    height = unit(H_mm, "mm"),
    xscale = c(0, W_mm),
    yscale = c(0, H_mm)
  ))

  # Helper: push a panel viewport with explicit mm position and size
  pvp <- function(r, c, x_mm = x_left[c], w_mm = COL_W[c])
    grid::pushViewport(grid::viewport(
      x      = unit(x_mm,      "mm"),
      y      = unit(y_bot[r],  "mm"),
      width  = unit(w_mm,      "mm"),
      height = unit(ROW_H[r],  "mm"),
      just   = c("left", "bottom")
    ))

  # ── Row 1 ────────────────────────────────────────────────────
  pvp(1,1)
  draw(ht_a, newpage = FALSE, padding = unit(c(5, 3, 3, 3), "mm"), heatmap_legend_side = "bottom")
  grid::popViewport()

  pvp(1,2); print(p_b, newpage = FALSE); grid::popViewport()
  pvp(1,3); print(p_c, newpage = FALSE); grid::popViewport()

  # ── Row 2 ────────────────────────────────────────────────────
  pvp(2,1); print(p_d, newpage = FALSE); grid::popViewport()
  pvp(2,2); print(p_e, newpage = FALSE); grid::popViewport()
  pvp(2,3); print(p_f, newpage = FALSE); grid::popViewport()

  # ── Row 3 ────────────────────────────────────────────────────
  pvp(3,1); print(p_g, newpage = FALSE); grid::popViewport()
  pvp(3,2); print(p_h, newpage = FALSE); grid::popViewport()
  pvp(3,3); print(p_i, newpage = FALSE); grid::popViewport()

  # ── Row 4 — asymmetric widths (j=47, k=75, l=61 mm) ─────────
  pvp(4,1, x_mm=ROW4_X[1], w_mm=ROW4_W[1]); print(p_j, newpage=FALSE); grid::popViewport()
  pvp(4,2, x_mm=ROW4_X[2], w_mm=ROW4_W[2]); print(p_k, newpage=FALSE); grid::popViewport()
  pvp(4,3, x_mm=ROW4_X[3], w_mm=ROW4_W[3]); print(p_l, newpage=FALSE); grid::popViewport()

  # ── Panel labels (drawn in outer viewport mm coordinates) ────
  for (i in seq_len(nrow(lbl))) {
    grid::grid.text(
      lbl$lab[i],
      x    = unit(lbl$x_mm[i], "mm"),
      y    = unit(lbl$y_mm[i], "mm"),
      just = c("left", "top"),
      gp   = grid::gpar(fontsize   = 12,
                        fontface   = "bold",
                        fontfamily = FONT_FAMILY,
                        col        = "black")
    )
  }

  grid::popViewport()  # outer viewport
}

# ── PDF — pure vector (write to /tmp, then copy to avoid iCloud lock) ──
{
  tmp_pdf <- file.path(tempdir(), "SuppFig_09_tmp.pdf")
  cairo_pdf(tmp_pdf, width = mm2in(W_mm), height = mm2in(H_mm), family = FONT_FAMILY)
  render_composite()
  dev.off()
  file.copy(tmp_pdf, file.path(OUT, "SuppFig_09.pdf"), overwrite = TRUE)
  file.remove(tmp_pdf)
  message("  SuppFig_09.pdf saved (pure vector).")
}

# ── PNG — 600 DPI preview (write to /tmp first, then copy to avoid iCloud write error) ──
{
  px_W    <- round(W_mm * 600 / 25.4)
  px_H    <- round(H_mm * 600 / 25.4)
  tmp_png <- file.path(tempdir(), "SuppFig_09_tmp.png")
  grDevices::png(tmp_png, width = px_W, height = px_H, res = 600, type = "cairo")
  render_composite()
  dev.off()
  file.copy(tmp_png, file.path(OUT, "SuppFig_09.png"), overwrite = TRUE)
  file.remove(tmp_png)
  message("  SuppFig_09.png saved.")
}

# ── TIFF — 600 DPI, LZW (write to /tmp, then copy) ───────────
{
  tmp_tiff <- file.path(tempdir(), "SuppFig_09_tmp.tiff")
  grDevices::tiff(tmp_tiff, width = W_mm, height = H_mm,
                  units = "mm", res = 600, compression = "lzw", type = "cairo")
  render_composite()
  dev.off()
  file.copy(tmp_tiff, file.path(OUT, "SuppFig_09.tiff"), overwrite = TRUE)
  file.remove(tmp_tiff)
  message("  SuppFig_09.tiff saved.")
}





# =============================================================
# CLEAN UP old mis-named panel files (only SuppFig_08 files that
# may have been copied here from the previous wrong version)
# =============================================================
old_files <- list.files(OUT, pattern = "^Supp08[a-l]_.*\\.pdf$|^Supp09i_LR_category_barplot\\.pdf$|^Supp09l_hepatocyte_IFNg_scatter\\.pdf$",
                        full.names = TRUE)
if (length(old_files)) {
  file.remove(old_files)
  message(sprintf("  Removed %d old mis-named panel files.", length(old_files)))
}

message("=== SuppFig_09 complete — FS_TAG = ", FS_TAG, " panels ===")
