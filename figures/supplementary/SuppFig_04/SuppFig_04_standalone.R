#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_04_standalone.R
# Supplementary Figure 4: Immune Pathway NES Discordance + Hepatic Zonation
#                          Collapse + CYP Enzyme Suppression
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: EBioMedicine (Lancet family) - vector PDF + PNG/TIFF (600 DPI) composite
# =============================================================================
# 6 panels (3-row layout, 183x245mm canvas):
#   a = Immune Hallmark NES grouped bar (TC vs PR; 8 immune pathways) — 91x75
#   b = Zone scores boxplot (Periportal | Pericentral) x (Normal | Adjacent)
#       with 1-sided paired Wilcoxon P + Cohen's d_z annotation — 92x75
#   c = Zone marker mRNA logFC lollipop (top 12 / zone, faceted) — 91x95
#   d = Paired periportal score lines (Normal -> Adjacent) coloured by
#       collapse (7/12) vs non-collapse (5/12); 1-sided P + d_z annotation — 92x95
#   e = CYP enzyme multi-omics grouped bar (mRNA + protein log2FC, 8 detected
#       CYPs incl. CYP3A4 mRNA = -1.486, protein = -0.671) — 91x75
#   f = Peri-lesional periportal score vs total bilirubin scatter
#       with linear fit; Spearman rho = -0.706, P = 0.010, n = 12 — 92x75
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_04_standalone.R
# =============================================================================

cat("=== Supplementary Figure 4: Zonation Collapse & CYP Suppression ===\n")
cat("  Loading libraries...\n")

# =============================================================================
# SECTION 1: Library Imports
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(grid)
  library(ggpubr)
  library(ggrepel)
  library(ggsci)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_04")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants
# =============================================================================
FONT_FAMILY <- "Arial"
MM_PER_INCH <- 25.4
W_SINGLE <- 89; W_DOUBLE <- 183
ASSEMBLY_DPI <- 600   # pixel-exact rendering constant

mm2in <- function(mm) mm / MM_PER_INCH

FS_TITLE <- 10; FS_AXIS_TITLE <- 9; FS_AXIS_TEXT <- 8
FS_LEGEND_T <- 8; FS_LEGEND_L <- 8; FS_TAG <- 12

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
COL_UP       <- "#CD534CFF"
COL_DOWN     <- "#0073C2FF"
COL_NS       <- "#868686FF"
COL_NORMAL   <- "#7AA6DCFF"
COL_ADJACENT <- "#CD534CFF"
COL_TC       <- "#0073C2FF"
COL_PR       <- "#CD534CFF"
PAL_CAT      <- pal_jco("default")(10)

# =============================================================================
# SECTION 5: Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.45, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.ticks.length = unit(0.8, "mm"),
    axis.text = element_text(size = 7.5, color = "black", family = FONT_FAMILY),
    axis.title = element_text(size = 8.5, face = "bold", family = FONT_FAMILY),
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0,
                              family = FONT_FAMILY,
                              margin = margin(b = 1.5, unit = "mm")),
    plot.subtitle = element_text(size = 7.5, color = "grey40", family = FONT_FAMILY),
    legend.text = element_text(size = 7, family = FONT_FAMILY),
    legend.title = element_text(size = 7.5, face = "bold", family = FONT_FAMILY),
    legend.key.size = unit(2.5, "mm"),
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0, "mm"),
    legend.box.spacing = unit(0.8, "mm"),
    strip.text = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.margin = margin(2, 2.5, 2, 2.5, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# SECTION 6: Helper Functions
# =============================================================================
save_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  cairo_pdf(fpath, width = w_mm / MM_PER_INCH, height = h_mm / MM_PER_INCH, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# =============================================================================
# PANEL A: Immune Hallmark NES (TC vs PR discordance)
# =============================================================================
cat("\n  Panel A: Immune Hallmark NES discordance\n")
tryCatch({
  gsea_tc <- read.csv(file.path(RES, "phase2_enrichment/GSEA_transcriptomics_Hallmark.csv"),
                      stringsAsFactors = FALSE)
  gsea_pr <- read.csv(file.path(RES, "phase2_enrichment/GSEA_proteomics_Hallmark.csv"),
                      stringsAsFactors = FALSE)

  immune_pathways <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_COMPLEMENT",
                       "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",
                       "HALLMARK_COAGULATION", "HALLMARK_ALLOGRAFT_REJECTION",
                       "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

  tc_imm <- gsea_tc %>% filter(pathway %in% immune_pathways) %>%
    select(pathway, NES) %>% mutate(layer = "Transcriptomic")
  pr_imm <- gsea_pr %>% filter(pathway %in% immune_pathways) %>%
    select(pathway, NES) %>% mutate(layer = "Proteomic")

  df_A <- bind_rows(tc_imm, pr_imm) %>%
    mutate(pw_clean = clean_pathway(pathway),
           pw_clean = factor(pw_clean, levels = rev(unique(clean_pathway(immune_pathways)))))

  p_A <- ggplot(df_A, aes(x = NES, y = pw_clean, fill = layer)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65,
             color = "black", linewidth = 0.15, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.45, color = "black") +
    scale_fill_manual(values = c("Transcriptomic" = COL_TC, "Proteomic" = COL_PR),
                      name = NULL) +
    labs(title = "Immune pathways: TC vs PR NES",
         x = "Normalised enrichment score (NES)", y = NULL) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.height = unit(2, "mm"),
          legend.key.width  = unit(3.5, "mm"),
          axis.text.y = element_text(size = 7.5, lineheight = 0.9))

  save_pdf("Supp04a_immune_NES.pdf", 91, 75, print(p_A))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# =============================================================================
# PANEL B: Zone Scores Boxplot
# =============================================================================
cat("\n  Panel B: Zone scores boxplot\n")
tryCatch({
  zon <- read.csv(file.path(RES, "zonation_collapse/zonation_scores_per_sample.csv"),
                  stringsAsFactors = FALSE)

  df_B <- zon %>%
    select(patient, periportal_normal, periportal_adjacent,
           pericentral_normal, pericentral_adjacent) %>%
    pivot_longer(-patient, names_to = "key", values_to = "score") %>%
    mutate(
      zone = ifelse(grepl("periportal", key), "Periportal", "Pericentral"),
      group = ifelse(grepl("normal", key), "Normal", "Adjacent"),
      group = factor(group, levels = c("Normal", "Adjacent"))
    )

  # Cohen's dz + one-sided P (biological prior: zonation collapse)
  pp_d <- zon$periportal_normal - zon$periportal_adjacent
  pc_d <- zon$pericentral_normal - zon$pericentral_adjacent
  dz_pp <- mean(pp_d) / sd(pp_d); dz_pc <- mean(pc_d) / sd(pc_d)
  p1_pp <- wilcox.test(zon$periportal_normal, zon$periportal_adjacent,
                       paired = TRUE, alternative = "greater")$p.value
  p1_pc <- wilcox.test(zon$pericentral_normal, zon$pericentral_adjacent,
                       paired = TRUE, alternative = "greater")$p.value
  ann_B <- data.frame(
    zone  = c("Periportal", "Pericentral"),
    label = c(sprintf("P[1-sided]=='%.2f'~~italic(d[z])=='%.2f'", p1_pp, dz_pp),
              sprintf("P[1-sided]=='%.2f'~~italic(d[z])=='%.2f'", p1_pc, dz_pc)))

  p_B <- ggplot(df_B, aes(x = group, y = score, fill = group)) +
    geom_boxplot(width = 0.55, outlier.size = 0.8, outlier.alpha = 0.6,
                 alpha = 0.75, linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.55, color = "grey25") +
    geom_text(data = ann_B, aes(x = 1.5, y = Inf, label = label),
              parse = TRUE, inherit.aes = FALSE, size = 2.3, vjust = 1.6,
              family = FONT_FAMILY) +
    facet_wrap(~zone) +
    scale_fill_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT),
                      guide = "none") +
    labs(title = "Hepatic zonation scores",
         x = NULL, y = "Zone score (mean log\u2082 FC)") +
    coord_cartesian(clip = "off")

  save_pdf("Supp04b_zone_scores.pdf", 92, 75, print(p_B))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# =============================================================================
# PANEL C: Zone Marker Gene Expression (Forest Plot)
# =============================================================================
cat("\n  Panel C: Zone marker logFC\n")
tryCatch({
  markers <- read.csv(file.path(RES, "zonation_collapse/zone_marker_gene_expression.csv"),
                      stringsAsFactors = FALSE)

  markers_plot <- markers %>%
    filter(detected_mRNA == "True") %>%
    mutate(zone = factor(zone, levels = c("Periportal", "Pericentral"))) %>%
    group_by(zone) %>%
    arrange(desc(abs(mRNA_logFC)), .by_group = TRUE) %>%
    slice_head(n = 12) %>%
    ungroup() %>%
    arrange(zone, mRNA_logFC) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  p_C <- ggplot(markers_plot, aes(x = mRNA_logFC, y = gene, color = zone)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_segment(aes(x = 0, xend = mRNA_logFC, y = gene, yend = gene),
                 linewidth = 0.45, alpha = 0.55) +
    geom_point(size = 1.7) +
    facet_grid(zone ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(values = c("Periportal" = COL_TC, "Pericentral" = "#EFC000FF"),
                       guide = "none") +
    labs(title = "Zone marker mRNA expression (top 12 / zone)",
         x = "log\u2082 fold change (Adjacent vs Normal)", y = NULL) +
    theme(axis.text.y = element_text(size = 7, face = "italic", lineheight = 0.85),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, face = "bold", size = 8),
          strip.background = element_rect(fill = "grey95", color = NA),
          panel.spacing.y = unit(1.5, "mm"))

  save_pdf("Supp04c_zone_markers.pdf", 91, 95, print(p_C))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# =============================================================================
# PANEL D: Paired Periportal Score Changes
# =============================================================================
cat("\n  Panel D: Paired zone score changes\n")
tryCatch({
  zon <- read.csv(file.path(RES, "zonation_collapse/zonation_scores_per_sample.csv"),
                  stringsAsFactors = FALSE)

  df_D <- zon %>%
    select(patient, periportal_normal, periportal_adjacent) %>%
    mutate(direction = ifelse(periportal_adjacent < periportal_normal,
                              "Collapse", "Non-collapse")) %>%
    pivot_longer(c(periportal_normal, periportal_adjacent),
                 names_to = "condition", values_to = "score") %>%
    mutate(condition = ifelse(grepl("normal", condition), "Normal", "Adjacent"),
           condition = factor(condition, levels = c("Normal", "Adjacent")))

  # One-sided paired Wilcoxon (collapse hypothesis: Normal > Adjacent) + Cohen's dz
  d_diff <- zon$periportal_normal - zon$periportal_adjacent
  dz_D <- mean(d_diff) / sd(d_diff)
  p1_D <- wilcox.test(zon$periportal_normal, zon$periportal_adjacent,
                      paired = TRUE, alternative = "greater")$p.value
  n_collapse <- sum(d_diff > 0); n_total <- length(d_diff)
  ann_D <- sprintf("P[1-sided]=='%.2f'~~italic(d[z])=='%.2f'~~(%d/%d~collapse)",
                   p1_D, dz_D, n_collapse, n_total)

  p_D <- ggplot(df_D, aes(x = condition, y = score)) +
    geom_line(aes(group = patient, color = direction),
              linewidth = 0.5, alpha = 0.85) +
    geom_point(aes(fill = condition), shape = 21, size = 2.4,
               color = "grey20", stroke = 0.3) +
    scale_color_manual(values = c("Collapse" = COL_DOWN, "Non-collapse" = "grey60"),
                       name = "Direction") +
    scale_fill_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT),
                      guide = "none") +
    annotate("text", x = 1.5, y = Inf, label = ann_D, parse = TRUE,
             vjust = 1.6, size = 2.4, family = FONT_FAMILY) +
    labs(title = "Periportal score (paired, Normal \u2192 Adjacent)",
         x = NULL, y = "Periportal zone score") +
    coord_cartesian(clip = "off") +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.height = unit(2, "mm"),
          legend.key.width  = unit(3.5, "mm"))

  save_pdf("Supp04d_zone_changes.pdf", 92, 95, print(p_D))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# =============================================================================
# PANEL E: CYP Enzyme Multi-omics Comparison
# =============================================================================
cat("\n  Panel E: CYP enzyme multi-omics\n")
tryCatch({
  cyp <- read.csv(file.path(RES, "zonation_collapse/cyp_enzyme_multiomics.csv"),
                  stringsAsFactors = FALSE)

  cyp_plot <- cyp %>%
    filter(detected_mRNA == "True" | detected_protein == "True") %>%
    select(gene, mRNA_logFC, protein_logFC) %>%
    pivot_longer(-gene, names_to = "layer", values_to = "logFC") %>%
    mutate(layer = ifelse(grepl("mRNA", layer), "mRNA", "Protein"),
           gene = factor(gene, levels = rev(unique(cyp$gene[cyp$detected_mRNA == "True"]))))

  p_E <- ggplot(cyp_plot, aes(x = logFC, y = gene, fill = layer)) +
    geom_col(position = position_dodge(width = 0.65), width = 0.55,
             color = "black", linewidth = 0.15, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.45, color = "black") +
    scale_fill_manual(values = c("mRNA" = COL_TC, "Protein" = COL_PR), name = NULL) +
    labs(title = "CYP enzyme suppression",
         x = "log\u2082 fold change (Adjacent vs Normal)", y = NULL) +
    theme(axis.text.y = element_text(size = 7.5, face = "italic"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.height = unit(2, "mm"),
          legend.key.width  = unit(3.5, "mm"))

  save_pdf("Supp04e_CYP_multiomics.pdf", 91, 75, print(p_E))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# =============================================================================
# PANEL F: Periportal Zone Score vs Total Bilirubin
# =============================================================================
cat("\n  Panel F: Bilirubin correlation\n")
tryCatch({
  clin <- read.csv(file.path(RES, "zonation_collapse/zonation_clinical_merged.csv"),
                   stringsAsFactors = FALSE)

  # Use peri-lesional periportal score (matches manuscript ¶23: rho = -0.706, P = 0.010)
  clin$pp_adj <- clin$periportal_adjacent

  # Spearman correlation
  ct <- cor.test(clin$pp_adj, clin$total_bilirubin, method = "spearman", exact = FALSE)
  lbl <- sprintf("rho = %.3f\nP = %.3f", ct$estimate, ct$p.value)

  p_F <- ggplot(clin, aes(x = pp_adj, y = total_bilirubin)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey25",
                fill = "grey80", linewidth = 0.55, alpha = 0.35) +
    geom_point(size = 2.2, color = COL_ADJACENT, alpha = 0.85) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
             label = lbl, size = 2.6, family = FONT_FAMILY,
             label.size = 0.2, label.r = unit(0.4, "mm"),
             fill = "white", color = "black") +
    labs(title = "Zonation collapse vs hepatic dysfunction",
         x = "Periportal score (Peri-lesional)",
         y = "Total bilirubin (\u00b5mol/L)")

  save_pdf("Supp04f_bilirubin_corr.pdf", 92, 75, print(p_F))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

# =============================================================================
# =============================================================================
# Composite Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite SuppFig_04 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  H1 <- 75; H2 <- 95; H3 <- 75   # rows: A/B | C/D (taller for zone markers) | E/F
  W_L <- 91; W_R <- W_TOTAL - W_L

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    # Panel A (Row1 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H2+H3,"mm"),
      width=unit(W_L,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_A, newpage=FALSE); grid::popViewport()
    # Panel B (Row1 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H2+H3,"mm"),
      width=unit(W_R,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_B, newpage=FALSE); grid::popViewport()
    # Panel C (Row2 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H3,"mm"),
      width=unit(W_L,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(p_C, newpage=FALSE); grid::popViewport()
    # Panel D (Row2 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H3,"mm"),
      width=unit(W_R,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(p_D, newpage=FALSE); grid::popViewport()
    # Panel E (Row3 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_L,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(p_E, newpage=FALSE); grid::popViewport()
    # Panel F (Row3 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(0,"mm"),
      width=unit(W_R,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(p_F, newpage=FALSE); grid::popViewport()
    # Labels
    label_data <- data.frame(
      text=c("a","b","c","d","e","f"),
      x_mm=c(2,W_L+2,2,W_L+2,2,W_L+2),
      y_mm=c(H2+H3+H1-2,H2+H3+H1-2,H3+H2-2,H3+H2-2,H3-2,H3-2),
      stringsAsFactors=FALSE)
    for(i in seq_len(nrow(label_data))) {
      grid::grid.text(label=label_data$text[i],
        x=unit(label_data$x_mm[i],"mm"), y=unit(label_data$y_mm[i],"mm"),
        just=c("left","top"),
        gp=grid::gpar(fontsize=FS_TAG, fontface="bold", fontfamily=FONT_FAMILY))
    }
    grid::popViewport()
  }

  cairo_pdf(file.path(OUT,"SuppFig_04.pdf"),
            width=W_TOTAL/25.4, height=H_TOTAL/25.4, family=FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_04.pdf (vector)\n")

  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT,"SuppFig_04.png"),
                 width=px_W, height=px_H, res=DPI, type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_04.png\n")

  grDevices::tiff(file.path(OUT,"SuppFig_04.tiff"),
                  width=px_W, height=px_H, res=DPI, compression="lzw", type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_04.tiff\n")

  cat("  SuppFig_04 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 4 rendering complete ===\n")

