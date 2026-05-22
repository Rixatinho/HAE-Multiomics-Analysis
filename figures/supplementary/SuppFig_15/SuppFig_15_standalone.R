#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_15_standalone.R
# Supplementary Figure 15: MR Sensitivity, Drug Repurposing & Pan-liver Comparison
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Nature Communications
# =============================================================================
# 8 Panels:
#   A = MR multi-method evidence matrix (inference strength + n_methods + pleiotropy)
#   B = MR heterogeneity (Cochran's Q) for the 4 sensitivity-tested pairs
#   C = MR-Egger pleiotropy intercept for the 4 sensitivity-tested pairs
#   D = Drug candidates ranked by integrated evidence score (lollipop)
#   E = Network proximity z-score distribution
#   F = Drug target evidence heatmap (raw counts, consistent encoding)
#   G = Cross-cohort DEG overlap (echinococcosis comparator cohorts)
#   H = HAE robust DEGs (top by |log2FC|, gene-symbol-only)
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_15_standalone.R
# =============================================================================

cat("=== Supplementary Figure 15: MR Sensitivity, Drug Repurposing & Pan-liver ===\n")
cat("  Loading libraries...\n")

# =============================================================================
# SECTION 1: Library Imports
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(ggsci)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_15")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FS_TAG <- 12
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4

W_SINGLE <- 89; W_DOUBLE <- 183; W_HALF <- 89
H_STD <- 85; H_TALL <- 100

ASSEMBLY_DPI <- 600   # pixel-exact rendering constant

mm2in <- function(mm) mm / MM_PER_INCH

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
PAL_CAT  <- pal_jco("default")(10)
COL_UP   <- "#CD534CFF"
COL_DOWN <- "#0073C2FF"
COL_NS   <- "#868686FF"
COL_MT   <- "#EFC000FF"

# =============================================================================
# SECTION 5: Theme (tightened per S14 lessons to prevent title truncation)
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = 8, color = "black", family = FONT_FAMILY),
    axis.title = element_text(size = 9, face = "bold", family = FONT_FAMILY),
    plot.title = element_text(size = 10, face = "bold", hjust = 0, family = FONT_FAMILY),
    plot.subtitle = element_text(size = 8, color = "grey25", family = FONT_FAMILY),
    legend.text = element_text(size = 8, family = FONT_FAMILY),
    legend.title = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    legend.key.size = unit(3, "mm"),
    legend.background = element_blank(),
    plot.margin = margin(1.5, 4, 1.5, 4, "mm")
  )
theme_set(theme_nc)

# =============================================================================
# SECTION 6: Helper Functions
# =============================================================================
save_pdf <- function(filename, w_mm, h_mm, expr) {
  fpath <- file.path(OUT, filename)
  cairo_pdf(fpath, width = w_mm/MM_PER_INCH, height = h_mm/MM_PER_INCH, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat("  ->", basename(fpath), sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

gp_rn <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID)
gp_cn <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
gp_lt <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID, fontface = "bold")
gp_ll <- function(size = 8) gpar(fontsize = size, fontfamily = FONT_GRID)

# Initialise plot containers (prevents stale-p reuse if a tryCatch fails)
p_A <- p_B <- p_C <- p_D <- p_E <- p_G <- p_H <- ggplot() + theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "Panel render error", size = 3)
ht_F <- NULL

# =============================================================================
# SECTION 7: Data Loading
# =============================================================================
cat("  Loading MR, drug repurposing, and pan-liver data...\n")

mr_summary  <- read.csv(file.path(RES, "enhancement18_mr/MR_sensitivity_summary_table.csv"), stringsAsFactors = FALSE)
mr_hetero   <- read.csv(file.path(RES, "enhancement18_mr/mr_heterogeneity.csv"), stringsAsFactors = FALSE)
mr_pleio    <- read.csv(file.path(RES, "enhancement18_mr/mr_pleiotropy.csv"), stringsAsFactors = FALSE)
drug_rank   <- read.csv(file.path(RES, "optimization_drug_repurposing/drug_candidates_ranked_optimized.csv"), stringsAsFactors = FALSE)
net_prox    <- read.csv(file.path(RES, "optimization_drug_repurposing/network_proximity_fullPPI.csv"), stringsAsFactors = FALSE)
pan_overlap <- read.csv(file.path(RES, "enhancement20_pan_liver/cross_dataset_overlap_stats.csv"), stringsAsFactors = FALSE)
hae_degs    <- read.csv(file.path(RES, "enhancement20_pan_liver/hae_robust_degs.csv"), stringsAsFactors = FALSE)
hae_sig     <- read.csv(file.path(RES, "enhancement20_pan_liver/hae_signature_summary.csv"), stringsAsFactors = FALSE)
ext_inv     <- read.csv(file.path(RES, "enhancement20_pan_liver/external_dataset_inventory.csv"), stringsAsFactors = FALSE)

cat("  Data loaded.\n")

# =============================================================================
# Panel A: MR Multi-method Evidence Matrix (REDESIGNED)
# Was: degenerate single-column tile heatmap with x="Evidence"
# Now: 3-column evidence matrix (strength category, # methods significant,
#      # methods direction-consistent), with dagger marker for pleiotropy
# =============================================================================
cat("  Panel A: MR multi-method evidence matrix\n")
p_A <- tryCatch({
  mr_summary$pair <- paste0(mr_summary$exposure, " -> ", mr_summary$outcome)
  # Order pairs by strength then n_methods_significant
  strength_order <- c("Strong" = 4, "Moderate" = 3, "Weak" = 2, "Inconclusive" = 1)
  mr_summary$strength_num <- strength_order[mr_summary$inference_strength]
  mr_summary <- mr_summary[order(-mr_summary$strength_num,
                                 -mr_summary$n_methods_significant), ]
  mr_summary$pair_f <- factor(mr_summary$pair, levels = rev(mr_summary$pair))

  # Short labels for in-tile text to avoid horizontal clipping inside narrow tiles
  strength_short <- c("Strong" = "Strong", "Moderate" = "Mod.",
                      "Weak" = "Weak", "Inconclusive" = "Inc.")

  # Long-format evidence matrix
  ev_df <- data.frame(
    pair = rep(mr_summary$pair_f, 3),
    metric = factor(rep(c("Strength", "n.sig", "n.conc"),
                        each = nrow(mr_summary)),
                    levels = c("Strength", "n.sig", "n.conc")),
    value = c(strength_short[mr_summary$inference_strength],
              as.character(mr_summary$n_methods_significant),
              as.character(mr_summary$n_methods_consistent_direction)),
    fill_class = c(mr_summary$inference_strength,
                   ifelse(mr_summary$n_methods_significant >= 3, "Strong",
                          ifelse(mr_summary$n_methods_significant >= 2, "Moderate",
                                 ifelse(mr_summary$n_methods_significant >= 1, "Weak", "Inconclusive"))),
                   ifelse(mr_summary$n_methods_consistent_direction >= 3, "Strong",
                          ifelse(mr_summary$n_methods_consistent_direction >= 2, "Moderate",
                                 "Weak"))),
    pleio = c(mr_summary$pleiotropy_detected, rep(FALSE, 2*nrow(mr_summary))),
    stringsAsFactors = FALSE)

  ggplot(ev_df, aes(x = metric, y = pair, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = value), size = 2.6, family = FONT_FAMILY,
              colour = "black") +
    geom_text(data = subset(ev_df, metric == "Strength" & pleio),
              aes(label = "\u2020"),
              size = 2.6, vjust = -1.6, hjust = 1.0,
              family = FONT_FAMILY, colour = "black") +
    scale_fill_manual(values = c("Strong" = COL_UP, "Moderate" = COL_MT,
                                  "Weak" = COL_DOWN, "Inconclusive" = COL_NS),
                      name = "Evidence",
                      breaks = c("Strong", "Moderate", "Weak", "Inconclusive")) +
    scale_x_discrete(position = "top",
                     labels = c("Strength" = "Strength",
                                "n.sig"    = "n.sig",
                                "n.conc"   = "n.conc")) +
    coord_cartesian(clip = "off") +
    labs(title = "MR causal evidence matrix",
         subtitle = "12 exposure-outcome pairs; \u2020 = pleiotropy detected",
         x = NULL, y = NULL) +
    theme(axis.text.x.top = element_text(angle = 35, hjust = 0, vjust = 0,
                                         size = 8, lineheight = 0.9,
                                         margin = margin(b = 1)),
          legend.position = "right",
          legend.key.height = unit(3, "mm"),
          plot.margin = margin(2, 4, 1.5, 4, "mm"))
}, error = function(e) { cat("    ERROR Panel A:", e$message, "\n"); p_A })
save_pdf("Supp15a_MR_strength.pdf", W_HALF, H_TALL, print(p_A))

# =============================================================================
# Panel B: MR Heterogeneity (Cochran's Q) — IVW only, 4 sensitivity-tested pairs
# =============================================================================
cat("  Panel B: MR heterogeneity\n")
p_B <- tryCatch({
  hetero_ivw <- mr_hetero[mr_hetero$method == "Inverse variance weighted", ]
  hetero_ivw$pair <- paste0(hetero_ivw$exposure, " -> ", hetero_ivw$outcome)
  hetero_ivw$sig <- ifelse(hetero_ivw$Q_pval < 0.05, "Significant", "NS (pass)")
  hetero_ivw <- hetero_ivw[order(-hetero_ivw$Q), ]
  hetero_ivw$pair_f <- factor(hetero_ivw$pair, levels = rev(hetero_ivw$pair))
  # Per-pair chi-square 95% threshold (Q_df differs per pair)
  hetero_ivw$threshold <- qchisq(0.95, df = hetero_ivw$Q_df)

  ggplot(hetero_ivw, aes(x = Q, y = pair_f)) +
    geom_segment(aes(x = 0, xend = Q, yend = pair_f, color = sig), linewidth = 0.4) +
    geom_point(aes(color = sig), size = 3) +
    geom_point(aes(x = threshold), shape = 4, color = "grey40",
               size = 2.5, stroke = 0.5) +
    geom_text(aes(label = sprintf("P=%.3f", Q_pval)),
              vjust = -1.2, hjust = 0.5,
              size = 2.83, family = FONT_FAMILY, color = "black") +
    scale_color_manual(values = c("Significant" = COL_UP, "NS (pass)" = COL_DOWN), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    scale_y_discrete(expand = expansion(mult = c(0.30, 0.40))) +
    coord_cartesian(clip = "off") +
    labs(title = "Cochran's Q heterogeneity (IVW)",
         subtitle = "All four pairs pass; \u2715 = \u03c7\u00b2 95% threshold",
         x = "Q statistic", y = NULL) +
    theme(legend.position = "bottom")
}, error = function(e) { cat("    ERROR Panel B:", e$message, "\n"); p_B })
save_pdf("Supp15b_heterogeneity.pdf", W_HALF, H_TALL, print(p_B))

# =============================================================================
# Panel C: MR-Egger Pleiotropy Intercept — 4 sensitivity-tested pairs
# =============================================================================
cat("  Panel C: MR-Egger intercept\n")
p_C <- tryCatch({
  mr_pleio$pair <- paste0(mr_pleio$exposure, " -> ", mr_pleio$outcome)
  mr_pleio$sig <- ifelse(mr_pleio$pval < 0.05, "Pleiotropy", "No pleiotropy")
  mr_pleio <- mr_pleio[order(mr_pleio$egger_intercept), ]
  mr_pleio$pair_f <- factor(mr_pleio$pair, levels = mr_pleio$pair)

  ggplot(mr_pleio, aes(x = egger_intercept, y = pair_f, color = sig)) +
    geom_errorbarh(aes(xmin = egger_intercept - 1.96 * se,
                       xmax = egger_intercept + 1.96 * se),
                   height = 0.25, linewidth = 0.4) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_text(aes(label = sprintf("P=%.2f", pval)),
              vjust = -1.2, hjust = 0.5,
              size = 2.83, family = FONT_FAMILY,
              color = "black", show.legend = FALSE) +
    scale_color_manual(values = c("Pleiotropy" = COL_UP, "No pleiotropy" = COL_DOWN), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    scale_y_discrete(expand = expansion(mult = c(0.30, 0.40))) +
    coord_cartesian(clip = "off") +
    labs(title = "MR-Egger pleiotropy intercept",
         subtitle = "No directional pleiotropy (P > 0.05)",
         x = "Egger intercept (95% CI)", y = NULL) +
    theme(legend.position = "bottom")
}, error = function(e) { cat("    ERROR Panel C:", e$message, "\n"); p_C })
save_pdf("Supp15c_pleiotropy.pdf", W_HALF, H_TALL, print(p_C))

# =============================================================================
# Panel D: Drug Candidates — ranked by integrated evidence score (lollipop)
# Was: bar plot with x = combined_rank (lower=better, counter-intuitive)
# Now: lollipop with x = evidence_score (higher=better), proximity-z annotated
# =============================================================================
cat("  Panel D: Drug candidates ranked by evidence score\n")
p_D <- tryCatch({
  top_drugs <- head(drug_rank[order(drug_rank$combined_rank), ], 12)
  top_drugs <- top_drugs[order(top_drugs$evidence_score, top_drugs$z_score,
                               decreasing = c(FALSE, TRUE)), ]
  top_drugs$drug <- factor(top_drugs$drug, levels = top_drugs$drug)
  top_drugs$sig_prox <- ifelse(top_drugs$p_value < 0.05,
                               "Proximity P < 0.05", "Proximity NS")

  ggplot(top_drugs, aes(x = evidence_score, y = drug)) +
    geom_segment(aes(x = 0, xend = evidence_score, yend = drug,
                     color = sig_prox), linewidth = 0.4) +
    geom_point(aes(color = sig_prox, size = abs(z_score))) +
    geom_text(aes(label = sprintf("z=%.2f", z_score)),
              hjust = -0.25, size = 2.83, family = FONT_FAMILY, color = "black") +
    scale_color_manual(values = c("Proximity P < 0.05" = COL_UP,
                                  "Proximity NS" = COL_NS), name = NULL) +
    scale_size_continuous(range = c(1.5, 3.8), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.40))) +
    scale_y_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    coord_cartesian(clip = "off") +
    labs(title = "Top 12 drug candidates",
         subtitle = "Ordered by evidence score (proximity + multi-omics)",
         x = "Evidence score", y = NULL) +
    theme(legend.position = "bottom")
}, error = function(e) { cat("    ERROR Panel D:", e$message, "\n"); p_D })
save_pdf("Supp15d_drug_ranking.pdf", W_HALF, H_TALL, print(p_D))

# =============================================================================
# Panel E: Network Proximity Z-score Distribution (annotation positioning fixed)
# =============================================================================
cat("  Panel E: Proximity z-score distribution\n")
p_E <- tryCatch({
  n_drugs <- nrow(net_prox)
  n_sig   <- sum(net_prox$p_value < 0.05, na.rm = TRUE)
  z_min   <- min(net_prox$z_score, na.rm = TRUE)
  z_max   <- max(net_prox$z_score, na.rm = TRUE)

  ggplot(net_prox, aes(x = z_score)) +
    geom_histogram(bins = 22, fill = "grey80", color = "grey55",
                   linewidth = 0.25, alpha = 0.95) +
    geom_vline(xintercept = -2, linetype = "dashed",
               color = COL_UP, linewidth = 0.5) +
    geom_rug(data = net_prox[net_prox$p_value < 0.05, ],
             color = COL_UP, linewidth = 0.4) +
    annotate("text", x = -1.85, y = Inf, vjust = 1.6, hjust = 0,
             label = "z = \u22122",
             size = 2.83, family = FONT_FAMILY, color = COL_UP, fontface = "bold") +
    annotate("text", x = z_max, y = Inf, vjust = 1.6, hjust = 1,
             label = sprintf("n = %d drugs\n%d significant", n_drugs, n_sig),
             size = 2.83, family = FONT_FAMILY) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(clip = "off") +
    labs(title = "Drug-disease proximity distribution",
         subtitle = "Network proximity z-score; rug marks significant candidates",
         x = "Network proximity z-score", y = "Drug count")
}, error = function(e) { cat("    ERROR Panel E:", e$message, "\n"); p_E })
save_pdf("Supp15e_zscore_dist.pdf", W_HALF, H_STD, print(p_E))

# =============================================================================
# Panel F: Drug Target Evidence Heatmap (consistent encoding)
# Was: column-normalised colour vs raw integer labels (mismatch)
# Now: 4 informative columns, value-encoded colour intensity, raw counts as labels
# =============================================================================
cat("  Panel F: Drug target evidence heatmap\n")
ht_F <- tryCatch({
  top15 <- head(drug_rank[order(drug_rank$combined_rank), ], 15)
  # Trim to 4 high-information columns; drop columns with all zero
  candidate_cols <- c("n_targets", "targets_in_network", "targets_DEG",
                      "targets_in_module")
  candidate_cols <- intersect(candidate_cols, colnames(top15))
  candidate_cols <- candidate_cols[sapply(candidate_cols, function(cc)
    sum(top15[[cc]], na.rm = TRUE) > 0)]
  if (length(candidate_cols) == 0) stop("No non-zero evidence columns")

  mat_raw <- as.matrix(top15[, candidate_cols, drop = FALSE])
  rownames(mat_raw) <- top15$drug
  colnames(mat_raw) <- gsub("targets_", "T:", candidate_cols)
  colnames(mat_raw) <- gsub("n_T:", "N targets", colnames(mat_raw))
  colnames(mat_raw) <- gsub("T:in_network", "In PPI", colnames(mat_raw))
  colnames(mat_raw) <- gsub("T:DEG",        "DEG", colnames(mat_raw))
  colnames(mat_raw) <- gsub("T:in_module",  "In module", colnames(mat_raw))

  # Single shared colour scale across whole matrix → label and colour agree
  mx <- max(mat_raw, na.rm = TRUE)
  col_evid <- colorRamp2(c(0, mx/2, mx), c("#FFFFFF", COL_MT, COL_UP))

  s_path <- file.path(OUT, "Supp15f_evidence_heatmap.pdf")
  cairo_pdf(s_path, width = mm2in(W_HALF), height = mm2in(H_TALL), family = FONT_FAMILY)
  ht <- Heatmap(mat_raw, name = "Targets",
    col = col_evid, border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    row_names_gp = gp_rn(7.5), column_names_gp = gp_cn(7.5),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%d", as.integer(mat_raw[i, j])), x, y,
                gp = gpar(fontsize = 7, fontfamily = FONT_GRID))
    },
    column_title = "Drug-target evidence (raw counts)",
    column_title_gp = gpar(fontsize = 8.5, fontface = "bold", fontfamily = FONT_FAMILY),
    heatmap_legend_param = list(title_gp = gp_lt(), labels_gp = gp_ll(),
                                 grid_width = unit(3, "mm")),
    show_row_dend = FALSE, show_column_dend = FALSE)
  draw(ht, padding = unit(c(3, 12, 3, 3), "mm"))
  dev.off()
  cat("  -> Supp15f_evidence_heatmap.pdf\n")
  ht
}, error = function(e) { cat("    ERROR Panel F:", e$message, "\n"); ht_F })

# =============================================================================
# Panel G: Cross-cohort DEG overlap (echinococcosis comparator cohorts)
# Was: misleading "Cross-disease" framing, hardcoded fallback
# Now: clarified subtitle naming the actual comparator cohorts (echinococcosis only)
# =============================================================================
cat("  Panel G: Cross-cohort overlap statistics\n")
p_G <- tryCatch({
  pan_sub <- pan_overlap[pan_overlap$dataset != "HAE_Human", ]
  if (nrow(pan_sub) == 0) stop("No comparator cohorts in cross_dataset_overlap_stats.csv")

  # Map dataset IDs to descriptive labels using external_dataset_inventory.csv
  # (best-effort; fallback to original ID)
  desc_map <- setNames(ext_inv$disease, ext_inv$gse_id)
  pan_sub$gse <- gsub("_.*", "", pan_sub$dataset)
  pan_sub$disease <- desc_map[pan_sub$gse]
  pan_sub$disease[is.na(pan_sub$disease)] <- "Comparator cohort"
  pan_sub$pretty <- paste0(pan_sub$gse, "\n(", pan_sub$disease, ")")
  pan_sub <- pan_sub[order(pan_sub$jaccard, decreasing = FALSE), ]
  pan_sub$pretty <- factor(pan_sub$pretty, levels = pan_sub$pretty)

  mean_j <- mean(pan_sub$jaccard, na.rm = TRUE)

  ggplot(pan_sub, aes(x = jaccard, y = pretty)) +
    geom_segment(aes(x = 0, xend = jaccard, yend = pretty),
                 color = COL_DOWN, linewidth = 0.4) +
    geom_point(color = COL_DOWN, size = 3) +
    geom_text(aes(label = sprintf("J=%.4f (n=%d)", jaccard, overlap_with_HAE)),
              hjust = -0.15, size = 2.6, family = FONT_FAMILY) +
    geom_vline(xintercept = mean_j, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.85))) +
    coord_cartesian(clip = "off") +
    labs(title = "Cross-cohort DEG overlap",
         subtitle = sprintf("HAE vs comparator cohorts; mean J = %.4f", mean_j),
         x = "Jaccard similarity", y = NULL) +
    theme(axis.text.y = element_text(size = 8, lineheight = 0.85),
          plot.margin = margin(2, 6, 1.5, 4, "mm"))
}, error = function(e) { cat("    ERROR Panel G:", e$message, "\n"); p_G })
save_pdf("Supp15g_overlap.pdf", W_HALF, H_STD, print(p_G))

# =============================================================================
# Panel H: HAE Robust DEGs — gene-symbol-only (filter ENSG-only IDs)
# =============================================================================
cat("  Panel H: HAE robust DEGs (gene-symbol-only)\n")
p_H <- tryCatch({
  hae_degs$abs_lfc <- abs(hae_degs$log2FC)
  # Use gene_name; filter rows where the symbol is just an ENSG ID
  nm_col <- if ("gene_name" %in% colnames(hae_degs)) "gene_name" else "gene"
  hae_degs$display <- hae_degs[[nm_col]]
  hae_named <- hae_degs[!grepl("^ENSG\\d+$", hae_degs$display), ]
  if (nrow(hae_named) == 0) hae_named <- hae_degs   # safety fallback
  top_n <- min(20, nrow(hae_named))
  top20 <- head(hae_named[order(-hae_named$abs_lfc), ], top_n)
  top20$display <- make.unique(top20$display)
  top20 <- top20[order(top20$log2FC), ]
  top20$display <- factor(top20$display, levels = top20$display)

  ggplot(top20, aes(x = log2FC, y = display, fill = direction)) +
    geom_col(width = 0.6, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c("up" = COL_UP, "down" = COL_DOWN), name = "Direction") +
    labs(title = "HAE robust DEGs (named, top by |log\u2082 FC|)",
         subtitle = sprintf("%d gene-symbol-resolved entries; ENSG-only IDs filtered", top_n),
         x = "log\u2082 FC", y = NULL) +
    theme(axis.text.y = element_text(size = 7, face = "italic"),
          legend.position = "top",
          legend.margin = margin(0, 0, 1, 0, "mm"))
}, error = function(e) { cat("    ERROR Panel H:", e$message, "\n"); p_H })
save_pdf("Supp15h_DEG_barplot.pdf", W_HALF, H_TALL, print(p_H))

# =============================================================================
# Composite Assembly (Pure Vector Grid Viewport, 183x245mm, 600DPI)
# =============================================================================
cat("\n--- Assembling composite SuppFig_15 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  H1 <- 61; H2 <- 61; H3 <- 61; H4 <- H_TOTAL - H1 - H2 - H3
  W_L <- 91; W_R <- W_TOTAL - W_L

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    # Panel A (Row1 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_A, newpage=FALSE); grid::popViewport()
    # Panel B (Row1 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(p_B, newpage=FALSE); grid::popViewport()
    # Panel C (Row2 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(p_C, newpage=FALSE); grid::popViewport()
    # Panel D (Row2 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(p_D, newpage=FALSE); grid::popViewport()
    # Panel E (Row3 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(p_E, newpage=FALSE); grid::popViewport()
    # Panel F (Row3 right)
    if (!is.null(ht_F)) {
      grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H4,"mm"),
        width=unit(W_R,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
      draw(ht_F, padding=unit(c(3,12,3,3),"mm"), newpage=FALSE)
      grid::popViewport()
    }
    # Panel G (Row4 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_L,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(p_G, newpage=FALSE); grid::popViewport()
    # Panel H (Row4 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(0,"mm"),
      width=unit(W_R,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(p_H, newpage=FALSE); grid::popViewport()
    # Labels
    label_data <- data.frame(
      text=c("a","b","c","d","e","f","g","h"),
      x_mm=c(2,W_L+2,2,W_L+2,2,W_L+2,2,W_L+2),
      y_mm=c(H2+H3+H4+H1-2,H2+H3+H4+H1-2,H3+H4+H2-2,H3+H4+H2-2,
             H4+H3-2,H4+H3-2,H4-2,H4-2), stringsAsFactors=FALSE)
    for(i in seq_len(nrow(label_data))) {
      grid::grid.text(label=label_data$text[i],
        x=unit(label_data$x_mm[i],"mm"), y=unit(label_data$y_mm[i],"mm"),
        just=c("left","top"),
        gp=grid::gpar(fontsize=FS_TAG, fontface="bold", fontfamily=FONT_FAMILY))
    }
    grid::popViewport()
  }

  cairo_pdf(file.path(OUT,"SuppFig_15.pdf"),
            width=W_TOTAL/25.4, height=H_TOTAL/25.4, family=FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_15.pdf (vector)\n")

  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT,"SuppFig_15.png"),
                 width=px_W, height=px_H, res=DPI, type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_15.png\n")

  grDevices::tiff(file.path(OUT,"SuppFig_15.tiff"),
                  width=px_W, height=px_H, res=DPI, compression="lzw", type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_15.tiff\n")

  cat("  SuppFig_15 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 15 rendering complete ===\n")
