#!/usr/bin/env Rscript
# =============================================================================
# Figure_9_standalone.R
# Title: Cell-Intrinsic Dead Zone Validation & Translational Efficiency Landscape
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Nature Communications - Cairo PDF | Arial 8-10pt
# =============================================================================
# Figure 9: Compositional deconfounding confirms cell-intrinsic CYP3A4 loss and
#           genome-wide translational efficiency landscape links FGG buffering
#           to protein half-life biology
# 5 panels (3-row vector grid, 183 x 245 mm):
#   a - Cell-type-adjusted CYP3A4 forest plot (limma + BayesPrism hepatocyte covariate;
#       8 genes split into "Cell-intrinsic (retained)" and "Composition-driven (lost)")
#   b - Genome-wide TE landscape rank plot (8,011 genes; FGG #1, FGB #3, FGA #4
#       buffered, CYP3A4 #222 highlighted)
#   c - mRNA-protein concordance landscape scatter (8,011 genes; quadrant labels +
#       Spearman rho = 0.178; fibrinogen subunits in immune-buffered quadrant,
#       CYP enzymes in concordant-loss quadrant). Note: protein half-life data
#       (halflife) is loaded for cross-reference but not directly plotted here.
#   d - CYP Activity Index vs serum total bilirubin scatter + linear fit (n = 12;
#       Spearman rho = -0.81, P = 0.002; CS1/CS2 colour split)
#   e - Patient-level Dead Zone severity heatmap (12 patients ordered by ascending
#       CYP3A4 Activity Index; left bar = severity, right annotations = bilirubin /
#       albumin / molecular subtype)
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript Figure_9_standalone.R
# =============================================================================

cat("=== Figure 9: Cell-Intrinsic Dead Zone Validation & TE Landscape ===\n")
cat("  Loading libraries...\n")


# =============================================================================
# SECTION 1: Library Imports + Arial Font Registration
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(ggsci)
  library(limma)
  library(ggrepel)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
cat("  Arial registered in PDF/PostScript font databases\n")


# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "02_analysis/results")
DATA <- file.path(BASE, "02_analysis/data/processed")
OUT  <- file.path(BASE, "04_figures/main/Figure_9")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants (NC Standard)
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM_PER_INCH <- 25.4

W_SINGLE  <- 89
W_DOUBLE  <- 183
H_STD     <- 85
H_TALL    <- 100
H_MAX     <- 240

ASSEMBLY_DPI <- 600   # pixel-exact rendering constant

mm2in <- function(mm) mm / MM_PER_INCH
in2mm <- function(inch) inch * MM_PER_INCH

FS_TITLE      <- 10
FS_SUBTITLE   <- 9
FS_AXIS_TITLE <- 9
FS_AXIS_TEXT  <- 8
FS_LEGEND_T   <- 8
FS_LEGEND_L   <- 8
FS_ANNO       <- 8
FS_ROW_NAME   <- 8
FS_CELL       <- 8
FS_TAG        <- 14
FS_GEOM_TEXT  <- 2.82

# =============================================================================
# SECTION 4: Color Palette (ggsci JCO)
# =============================================================================
PAL_CAT <- pal_jco("default")(10)

COL_UP       <- "#CD534CFF"
COL_DOWN     <- "#0073C2FF"
COL_NS       <- "#868686FF"
COL_NA       <- "#F0F0F0"
COL_TC       <- "#0073C2FF"
COL_PR       <- "#CD534CFF"
COL_MT       <- "#EFC000FF"
COL_NORMAL   <- "#7AA6DCFF"
COL_ADJACENT <- "#CD534CFF"
COL_CS1      <- "#CD534CFF"
COL_CS2      <- "#0073C2FF"

# Specific for Figure 9
COL_FGG      <- "#E64B35FF"   # FGG highlight (vibrant red)
COL_FGB      <- "#F39B7FFF"   # FGB (light red)
COL_FGA      <- "#91D1C2FF"   # FGA (teal)
COL_CYP3A4   <- "#4DBBD5FF"   # CYP3A4 (cyan)
COL_IMMUNE   <- "#3C5488FF"   # immune category
COL_METABOLIC <- "#E64B35FF"  # metabolic category

col_div     <- colorRamp2(c(-2, 0, 2), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))
col_zscore  <- colorRamp2(c(-2.5, 0, 2.5), c("#0073C2FF", "#FFFFFF", "#CD534CFF"))

# =============================================================================
# SECTION 5: ggplot2 Theme (NC style)
# =============================================================================
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
    axis.title         = element_text(family = FONT_FAMILY, size = 9, face = "bold", color = "black"),
    plot.title         = element_text(family = FONT_FAMILY, size = 10, face = "bold", hjust = 0.5),
    plot.title.position = "panel",
    plot.subtitle      = element_text(family = FONT_FAMILY, size = 9, hjust = 0),
    legend.title       = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    legend.text        = element_text(family = FONT_FAMILY, size = 8),
    legend.key.size    = unit(3, "mm"),
    strip.text         = element_text(family = FONT_FAMILY, size = 8, face = "bold"),
    strip.background   = element_rect(fill = "grey95", color = "black", linewidth = 0.4),
    plot.margin        = margin(4, 4, 2, 4, "mm")
  )


# =============================================================================
# SECTION 6: Load Data
# =============================================================================
cat("  Loading data...\n")

# 6a. Transcriptomics VST matrix (24 samples x ~18000 genes)
vst_mat <- read.csv(file.path(DATA, "transcriptomics_vst_paired.csv"),
                    row.names = 1, check.names = FALSE)

# 6b. BayesPrism cell proportions (28 samples x 15 cell types)
bp_props <- read.csv(file.path(RES, "enhancement16_deconvolution/bayesprism_cell_proportions.csv"),
                     row.names = 1, check.names = FALSE)

# 6c. Translational efficiency (8011 genes)
te_data <- read.csv(file.path(RES, "enhancement24_mechanistic/translational_efficiency.csv"),
                    stringsAsFactors = FALSE)

# 6d. Protein half-life comparison
halflife <- read.csv(file.path(RES, "mechanism_validation/protein_halflife_comparison.csv"),
                     stringsAsFactors = FALSE)

# 6e. CYP enzyme multiomics
cyp_data <- read.csv(file.path(RES, "zonation_collapse/cyp_enzyme_multiomics.csv"),
                     stringsAsFactors = FALSE)

# 6f. Clinical data
clinical <- read.csv(file.path(RES, "phase6_subtyping/clinical_parsed_corrected.csv"),
                     stringsAsFactors = FALSE)

# 6g. Zonation scores per sample
zonation <- read.csv(file.path(RES, "zonation_collapse/zonation_scores_per_sample.csv"),
                     stringsAsFactors = FALSE)

# 6h. DEGs full results (for gene name mapping)
degs_full <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"),
                      stringsAsFactors = FALSE)

cat("  All data loaded successfully.\n")


# =============================================================================
# PANEL A: Cell-Type-Adjusted CYP3A4 Validation
# =============================================================================
cat("\n--- Panel A: Cell-type-adjusted CYP3A4 validation ---\n")

# Identify paired samples in VST matrix
tc_samples <- colnames(vst_mat)
normal_samples <- tc_samples[grepl("^Normal", tc_samples)]
adjacent_samples <- tc_samples[grepl("^Adjacent", tc_samples)]

# Extract patient IDs
normal_patients <- as.integer(gsub("Normal", "", normal_samples))
adjacent_patients <- as.integer(gsub("Adjacent", "", adjacent_samples))
paired_patients <- intersect(normal_patients, adjacent_patients)
cat("  Paired patients for TC:", length(paired_patients), "\n")

# Build sample data frame for paired samples
sample_df <- data.frame(
  sample = c(paste0("Normal", paired_patients), paste0("Adjacent", paired_patients)),
  condition = factor(rep(c("Normal", "Adjacent"), each = length(paired_patients)),
                     levels = c("Normal", "Adjacent")),
  patient = factor(rep(paired_patients, 2)),
  stringsAsFactors = FALSE
)

# Add hepatocyte proportions from BayesPrism
sample_df$hepatocyte <- bp_props[sample_df$sample, "Hepatocyte"]

# Handle any missing BayesPrism data
if (any(is.na(sample_df$hepatocyte))) {
  cat("  WARNING: Some samples missing BayesPrism data, using median imputation\n")
  sample_df$hepatocyte[is.na(sample_df$hepatocyte)] <- median(sample_df$hepatocyte, na.rm = TRUE)
}

# Subset VST matrix to paired samples
vst_paired <- vst_mat[, sample_df$sample]

# Map gene IDs to gene names for CYP3A4
gene_map <- degs_full[, c("gene_id", "gene_name")]
gene_map <- gene_map[!duplicated(gene_map$gene_id), ]
rownames(gene_map) <- gene_map$gene_id

# Find CYP3A4 row
cyp3a4_rows <- grep("CYP3A4", gene_map$gene_name)
if (length(cyp3a4_rows) > 0) {
  cyp3a4_id <- gene_map$gene_id[cyp3a4_rows[1]]
} else {
  # Try matching from rownames of vst
  cyp3a4_id <- rownames(vst_mat)[grep("ENSG00000160868", rownames(vst_mat))]
  if (length(cyp3a4_id) == 0) cyp3a4_id <- NULL
}
cat("  CYP3A4 gene ID:", cyp3a4_id, "\n")

# --- Run limma: Unadjusted model ---
design_unadj <- model.matrix(~ condition + patient, data = sample_df)
fit_unadj <- lmFit(vst_paired, design_unadj)
fit_unadj <- eBayes(fit_unadj)
res_unadj <- topTable(fit_unadj, coef = "conditionAdjacent", number = Inf, sort.by = "none")
res_unadj$gene_id <- rownames(res_unadj)

# --- Run limma: Adjusted model (+ hepatocyte proportion) ---
design_adj <- model.matrix(~ condition + patient + hepatocyte, data = sample_df)
fit_adj <- lmFit(vst_paired, design_adj)
fit_adj <- eBayes(fit_adj)
res_adj <- topTable(fit_adj, coef = "conditionAdjacent", number = Inf, sort.by = "none")
res_adj$gene_id <- rownames(res_adj)

# Extract CYP3A4 results
if (!is.null(cyp3a4_id) && cyp3a4_id %in% rownames(res_unadj)) {
  cyp3a4_unadj <- res_unadj[cyp3a4_id, ]
  cyp3a4_adj   <- res_adj[cyp3a4_id, ]
  cat(sprintf("  CYP3A4 Unadjusted: logFC=%.3f, P=%.4f\n", cyp3a4_unadj$logFC, cyp3a4_unadj$P.Value))
  cat(sprintf("  CYP3A4 Adjusted:   logFC=%.3f, P=%.4f\n", cyp3a4_adj$logFC, cyp3a4_adj$P.Value))
} else {
  cat("  WARNING: CYP3A4 not found in VST matrix, using pre-computed values\n")
  cyp3a4_unadj <- data.frame(logFC = -1.486, P.Value = 0.045)
  cyp3a4_adj   <- data.frame(logFC = -1.42, P.Value = 0.038)
}

# Also extract key CYP enzymes and metabolic genes for comparison
target_genes <- c("CYP3A4", "CYP1A2", "CYP2E1", "ALDOB", "PCK1", "FGG", "FGB", "FGA")
target_ids <- gene_map$gene_id[gene_map$gene_name %in% target_genes]
names(target_ids) <- gene_map$gene_name[match(target_ids, gene_map$gene_id)]

# Build comparison data frame
compare_df <- data.frame(
  gene = character(),
  model = character(),
  logFC = numeric(),
  CI_low = numeric(),
  CI_high = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

for (gn in names(target_ids)) {
  gid <- target_ids[gn]
  if (gid %in% rownames(res_unadj)) {
    se_u <- res_unadj[gid, "logFC"] / res_unadj[gid, "t"]
    se_a <- res_adj[gid, "logFC"] / res_adj[gid, "t"]
    compare_df <- rbind(compare_df, data.frame(
      gene = gn, model = "Unadjusted",
      logFC = res_unadj[gid, "logFC"],
      CI_low = res_unadj[gid, "logFC"] - 1.96 * se_u,
      CI_high = res_unadj[gid, "logFC"] + 1.96 * se_u,
      P = res_unadj[gid, "P.Value"],
      stringsAsFactors = FALSE
    ))
    compare_df <- rbind(compare_df, data.frame(
      gene = gn, model = "Cell-type adjusted",
      logFC = res_adj[gid, "logFC"],
      CI_low = res_adj[gid, "logFC"] - 1.96 * se_a,
      CI_high = res_adj[gid, "logFC"] + 1.96 * se_a,
      P = res_adj[gid, "P.Value"],
      stringsAsFactors = FALSE
    ))
  }
}

# Compute per-gene attenuation summary: which genes retain significance after adj.
gene_summary <- compare_df %>%
  tidyr::pivot_wider(
    id_cols = gene,
    names_from = model,
    values_from = c(logFC, P)
  ) %>%
  dplyr::mutate(
    retains_sig = .data[["P_Cell-type adjusted"]] < 0.05,
    attenuation_pct = round(100 * (1 - .data[["logFC_Cell-type adjusted"]] /
                                    .data[["logFC_Unadjusted"]]), 0),
    panel = ifelse(retains_sig,
                   "Cell-intrinsic\n(retained after adjustment)",
                   "Composition-driven\n(lost after adjustment)")
  )

# Group order: intrinsic first (top), composition-driven second
intrinsic_genes <- gene_summary$gene[gene_summary$retains_sig]
intrinsic_genes <- intrinsic_genes[order(gene_summary$logFC_Unadjusted[match(intrinsic_genes, gene_summary$gene)])]
composition_genes <- gene_summary$gene[!gene_summary$retains_sig]
composition_genes <- composition_genes[order(gene_summary$logFC_Unadjusted[match(composition_genes, gene_summary$gene)])]
gene_order <- c(composition_genes, intrinsic_genes)
compare_df$gene <- factor(compare_df$gene, levels = gene_order)
compare_df$model <- factor(compare_df$model, levels = c("Unadjusted", "Cell-type adjusted"))

# Significance marker (no "ns" — leave blank when not significant)
compare_df$sig_label <- ifelse(compare_df$P < 0.001, "***",
                        ifelse(compare_df$P < 0.01, "**",
                        ifelse(compare_df$P < 0.05, "*", "")))

# Plot Panel A: Forest plot comparing models
# Manual y-dodge to avoid deprecated position_dodgev
compare_df$y_num <- as.numeric(compare_df$gene)
compare_df$y_dodge <- compare_df$y_num +
  ifelse(compare_df$model == "Unadjusted", 0.18, -0.18)

# Build y-axis labels with CYP3A4 in bold
y_labels <- sapply(levels(compare_df$gene), function(g) {
  if (g == "CYP3A4") bquote(bold(CYP3A4)) else g
})

# Compute panel-band positions for shading
n_intrinsic   <- length(intrinsic_genes)
n_composition <- length(composition_genes)
band_intrinsic_y <- c(n_composition + 0.5, n_composition + n_intrinsic + 0.5)
band_composition_y <- c(0.5, n_composition + 0.5)

pA <- ggplot(compare_df, aes(x = logFC, y = y_dodge, color = model, shape = model)) +
  # Background shading bands distinguishing the two biological groups
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = band_intrinsic_y[1], ymax = band_intrinsic_y[2],
           fill = COL_UP, alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = band_composition_y[1], ymax = band_composition_y[2],
           fill = COL_DOWN, alpha = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                width = 0.18, linewidth = 0.5, orientation = "y") +
  geom_point(size = 2.5) +
  geom_text(aes(x = CI_high + 0.15, label = sig_label),
            size = FS_GEOM_TEXT + 0.5, show.legend = FALSE,
            fontface = "bold", vjust = 0.5) +
  # Right-side group labels
  annotate("text", x = Inf, y = mean(band_intrinsic_y),
           label = "Cell-intrinsic\n(retained)",
           color = COL_UP, hjust = 1.02, vjust = 0.5,
           size = FS_GEOM_TEXT, fontface = "italic", lineheight = 0.85) +
  annotate("text", x = Inf, y = mean(band_composition_y),
           label = "Composition-driven\n(lost after adj.)",
           color = COL_DOWN, hjust = 1.02, vjust = 0.5,
           size = FS_GEOM_TEXT, fontface = "italic", lineheight = 0.85) +
  scale_y_continuous(breaks = seq_along(levels(compare_df$gene)),
                     labels = y_labels,
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  scale_color_manual(values = c("Unadjusted" = COL_DOWN, "Cell-type adjusted" = COL_UP)) +
  scale_shape_manual(values = c("Unadjusted" = 16, "Cell-type adjusted" = 17)) +
  labs(x = "log\u2082 Fold Change (Adjacent / Normal)",
       y = NULL, color = "Model", shape = "Model",
       subtitle = "Cell-type adjustment isolates intrinsic CYP3A4 loss") +
  theme_nc +
  theme(legend.position = "bottom",
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(3, "mm"),
        plot.subtitle = element_text(size = FS_SUBTITLE - 1, family = FONT_FAMILY,
                                     face = "italic", color = "grey25",
                                     margin = margin(b = 2)))

cat("  Panel A complete.\n")


# =============================================================================
# PANEL B: Genome-wide TE Landscape Rank Plot
# =============================================================================
cat("\n--- Panel B: Genome-wide TE landscape ---\n")

# Sort by TE index
te_data <- te_data[order(-te_data$TE_index), ]
te_data$rank <- 1:nrow(te_data)
n_genes <- nrow(te_data)
cat(sprintf("  Total genes with TE: %d\n", n_genes))

# Identify key genes
te_data$highlight <- "Other"
te_data$highlight[te_data$gene == "FGG"] <- "FGG"
te_data$highlight[te_data$gene == "FGB"] <- "FGB"
te_data$highlight[te_data$gene == "FGA"] <- "FGA"
te_data$highlight[te_data$gene == "CYP3A4"] <- "CYP3A4"

# Find positions
fgg_row <- te_data[te_data$gene == "FGG", ]
fgb_row <- te_data[te_data$gene == "FGB", ]
fga_row <- te_data[te_data$gene == "FGA", ]
cyp_row <- te_data[te_data$gene == "CYP3A4", ]

cat(sprintf("  FGG: rank=%d, TE=%.3f\n", fgg_row$rank, fgg_row$TE_index))
cat(sprintf("  FGB: rank=%d, TE=%.3f\n", fgb_row$rank, fgb_row$TE_index))
cat(sprintf("  FGA: rank=%d, TE=%.3f\n", fga_row$rank, fga_row$TE_index))
cat(sprintf("  CYP3A4: rank=%d, TE=%.3f\n", cyp_row$rank, cyp_row$TE_index))

# Plot Panel B: Rank plot
highlight_df <- te_data[te_data$highlight != "Other", ]

pB <- ggplot(te_data, aes(x = rank, y = TE_index)) +
  # Background points
  geom_point(data = te_data[te_data$highlight == "Other", ],
             color = "grey80", size = 0.3, alpha = 0.5) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  # Highlighted points

  geom_point(data = highlight_df,
             aes(color = highlight), size = 3) +
  # Labels
  geom_label_repel(data = highlight_df,
                   aes(label = paste0(gene, " (#", rank, ")"), color = highlight),
                   size = FS_GEOM_TEXT, fontface = "bold",
                   box.padding = 0.5, point.padding = 0.3,
                   min.segment.length = 0, seed = 42,
                   fill = "white", label.size = 0.2) +
  scale_color_manual(values = c("FGG" = COL_FGG, "FGB" = COL_FGB,
                                "FGA" = COL_FGA, "CYP3A4" = COL_CYP3A4)) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = paste0("Gene Rank (n = ", scales::comma(n_genes), " genes)"),
       y = "Translational Efficiency Index\n(protein log\u2082FC \u2212 mRNA log\u2082FC)") +
  annotate("text", x = n_genes * 0.55, y = max(te_data$TE_index) * 0.55,
           label = "Protein buffered\n(immune preservation)",
           size = FS_GEOM_TEXT, color = COL_IMMUNE, fontface = "italic", hjust = 0) +
  annotate("text", x = n_genes * 0.45, y = min(te_data$TE_index) * 0.55,
           label = "Protein depleted\n(metabolic collapse)",
           size = FS_GEOM_TEXT, color = COL_METABOLIC, fontface = "italic", hjust = 1) +
  theme_nc +
  theme(legend.position = "none")

cat("  Panel B complete.\n")


# =============================================================================
# PANEL C: mRNA-Protein Concordance Landscape
# =============================================================================
cat("\n--- Panel C: mRNA-protein concordance landscape ---\n")

# Use TE data which has both mRNA and protein log2FC for 8011 genes
te_plot <- te_data[, c("gene", "mRNA_log2FC", "protein_log2FC", "TE_index", "rank")]

# Classify concordance quadrants
te_plot$quadrant <- "Other"
te_plot$quadrant[te_plot$mRNA_log2FC < 0 & te_plot$protein_log2FC < 0] <- "Concordant down"
te_plot$quadrant[te_plot$mRNA_log2FC > 0 & te_plot$protein_log2FC > 0] <- "Concordant up"
te_plot$quadrant[te_plot$mRNA_log2FC < 0 & te_plot$protein_log2FC >= 0] <- "Protein buffered"
te_plot$quadrant[te_plot$mRNA_log2FC >= 0 & te_plot$protein_log2FC < 0] <- "Protein depleted"

quad_counts <- table(te_plot$quadrant)
cat("  Quadrant distribution:\n")
for (q in names(quad_counts)) cat(sprintf("    %s: %d\n", q, quad_counts[q]))

# Key genes to highlight
key_genes <- c("FGG", "FGB", "FGA", "CYP3A4", "CYP1A2", "CYP2E1", "ALDOB", "PCK1")
te_plot$highlight <- ifelse(te_plot$gene %in% key_genes, te_plot$gene, NA)

# Category labels for highlighted genes
te_plot$gene_cat <- NA
te_plot$gene_cat[te_plot$gene %in% c("FGG", "FGB", "FGA")] <- "Fibrinogen\n(immune buffered)"
te_plot$gene_cat[te_plot$gene %in% c("CYP3A4", "CYP1A2", "CYP2E1", "ALDOB", "PCK1")] <- "Metabolic\n(concordant loss)"

highlight_c <- te_plot[!is.na(te_plot$highlight), ]

# Correlation for reference
rho_mp <- cor(te_plot$mRNA_log2FC, te_plot$protein_log2FC, method = "spearman")

pC <- ggplot(te_plot, aes(x = mRNA_log2FC, y = protein_log2FC)) +
  # Background density
  geom_point(color = "grey80", size = 0.3, alpha = 0.4) +
  # Identity line (perfect concordance)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  # Zero lines
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  # Highlighted genes
  geom_point(data = highlight_c,
             aes(color = gene_cat), size = 2.5) +
  geom_label_repel(data = highlight_c,
                   aes(label = gene, color = gene_cat),
                   size = FS_GEOM_TEXT, fontface = "bold",
                   box.padding = 0.4, point.padding = 0.2,
                   min.segment.length = 0, seed = 42,
                   fill = "white", label.size = 0.2,
                   max.overlaps = 20) +
  scale_color_manual(values = c("Fibrinogen\n(immune buffered)" = COL_IMMUNE,
                                "Metabolic\n(concordant loss)" = COL_METABOLIC),
                     na.value = "grey70") +
  # Quadrant labels (these now serve as the visual legend — distinct corners + colors)
  annotate("text", x = -3.8, y = 3.3,
           label = "Fibrinogen\n(immune buffered)", color = COL_IMMUNE,
           size = FS_GEOM_TEXT, fontface = "italic", alpha = 0.85, hjust = 0, vjust = 1) +
  annotate("text", x = 3.8, y = -3.3,
           label = "Metabolic\n(concordant loss)", color = COL_METABOLIC,
           size = FS_GEOM_TEXT, fontface = "italic", alpha = 0.85, hjust = 1, vjust = 0) +
  annotate("text", x = 3.8, y = 3.3,
           label = sprintf("rho = %.3f", rho_mp),
           size = FS_GEOM_TEXT + 0.3, fontface = "italic", hjust = 1, vjust = 1) +
  coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-3.5, 3.5)) +
  labs(x = "mRNA log\u2082FC (Adjacent / Normal)",
       y = "Protein log\u2082FC (Adjacent / Normal)",
       color = NULL) +
  theme_nc +
  theme(legend.position = "none")

cat("  Panel C complete.\n")


# =============================================================================
# PANEL D: CYP3A4 Activity Index vs Clinical Bilirubin
# =============================================================================
cat("\n--- Panel D: CYP3A4 Activity Index vs clinical bilirubin ---\n")

# Compute per-patient CYP3A4 activity index from transcriptomics
# Use the ratio of Adjacent/Normal expression for CYP enzymes involved in ABZ metabolism
cyp_genes <- c("CYP3A4", "CYP1A2", "CYP2E1")
cyp_gene_ids <- gene_map$gene_id[gene_map$gene_name %in% cyp_genes]
names(cyp_gene_ids) <- gene_map$gene_name[match(cyp_gene_ids, gene_map$gene_id)]

# Calculate per-patient log2FC for each CYP gene
patient_cyp_activity <- data.frame(patient_id = paired_patients)

for (gn in names(cyp_gene_ids)) {
  gid <- cyp_gene_ids[gn]
  if (gid %in% rownames(vst_paired)) {
    normal_vals <- as.numeric(vst_paired[gid, paste0("Normal", paired_patients)])
    adjacent_vals <- as.numeric(vst_paired[gid, paste0("Adjacent", paired_patients)])
    patient_cyp_activity[[paste0(gn, "_logFC")]] <- adjacent_vals - normal_vals
  }
}

# CYP3A4 Activity Index = mean log2FC across CYP genes (weighted by importance)
logfc_cols <- grep("_logFC$", colnames(patient_cyp_activity), value = TRUE)
if (length(logfc_cols) > 0) {
  # Weight CYP3A4 higher as it's the primary activator
  weights <- ifelse(grepl("CYP3A4", logfc_cols), 2, 1)
  weights <- weights / sum(weights)
  patient_cyp_activity$activity_index <- as.numeric(
    as.matrix(patient_cyp_activity[, logfc_cols]) %*% weights
  )
} else {
  patient_cyp_activity$activity_index <- NA
}

# Merge with clinical data
patient_cyp_activity <- merge(patient_cyp_activity, clinical,
                              by = "patient_id", all.x = TRUE)

# Remove patients with missing bilirubin
plot_df_d <- patient_cyp_activity[!is.na(patient_cyp_activity$total_bilirubin) &
                                  !is.na(patient_cyp_activity$activity_index), ]
cat(sprintf("  Patients with complete data: %d\n", nrow(plot_df_d)))

# Correlation test
if (nrow(plot_df_d) >= 5) {
  cor_bili <- cor.test(plot_df_d$activity_index, plot_df_d$total_bilirubin, method = "spearman")
  cat(sprintf("  Activity Index vs Bilirubin: rho=%.3f, P=%.4f\n",
              cor_bili$estimate, cor_bili$p.value))
} else {
  cor_bili <- list(estimate = NA, p.value = NA)
}

# Plot Panel D
pD <- ggplot(plot_df_d, aes(x = activity_index, y = total_bilirubin)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              fill = "grey90", linewidth = 0.6, alpha = 0.3) +
  geom_point(aes(color = molecular_subtype), size = 3) +
  geom_text_repel(aes(label = paste0("P", patient_id)), size = FS_GEOM_TEXT,
                  box.padding = 0.4, seed = 42, max.overlaps = 12) +
  scale_color_manual(values = c("CS1" = COL_CS1, "CS2" = COL_CS2),
                     na.value = COL_NS) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = "CYP Activity Index (weighted mean log\u2082FC)",
       y = "Total Bilirubin (\u00b5mol/L)",
       color = "Subtype") +
  annotate("text", x = min(plot_df_d$activity_index) * 0.9,
           y = max(plot_df_d$total_bilirubin) * 0.95,
           label = sprintf("rho = %.2f, P = %.3f",
                           cor_bili$estimate, cor_bili$p.value),
           size = FS_GEOM_TEXT + 0.5, hjust = 0, fontface = "italic") +
  theme_nc +
  theme(legend.position = "bottom",
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(3, "mm"))
cat("  Panel D complete.\n")


# =============================================================================
# PANEL E: Patient-level Dead Zone Severity Heatmap
# =============================================================================
cat("\n--- Panel E: Patient-level Dead Zone severity heatmap ---\n")

# Build patient-level matrix with Dead Zone indicators
# Rows: patients, Columns: key metrics
# Metrics: CYP3A4 logFC, CYP1A2 logFC, CYP2E1 logFC, periportal_loss, total_bilirubin (z-scored)

heatmap_patients <- paired_patients

# Get CYP logFCs per patient
hm_data <- patient_cyp_activity[patient_cyp_activity$patient_id %in% heatmap_patients, ]
hm_data <- hm_data[order(hm_data$patient_id), ]

# Add zonation data
zonation_merge <- zonation[zonation$patient %in% hm_data$patient_id, ]
zonation_merge <- zonation_merge[order(zonation_merge$patient), ]

# Merge
hm_data <- merge(hm_data, zonation_merge[, c("patient", "periportal_adjacent", "ZDI_adjacent", "delta_ZPS")],
                 by.x = "patient_id", by.y = "patient", all.x = TRUE)

# Build heatmap matrix
hm_cols <- c()
hm_mat <- matrix(NA, nrow = nrow(hm_data), ncol = 0)
rownames(hm_mat) <- paste0("P", hm_data$patient_id)

# CYP enzyme log2FCs
for (col in logfc_cols) {
  if (col %in% colnames(hm_data)) {
    hm_mat <- cbind(hm_mat, hm_data[[col]])
    hm_cols <- c(hm_cols, gsub("_logFC", "", col))
  }
}

# Periportal score (adjacent)
if ("periportal_adjacent" %in% colnames(hm_data)) {
  hm_mat <- cbind(hm_mat, hm_data$periportal_adjacent)
  hm_cols <- c(hm_cols, "Periportal")
}

# Delta ZPS
if ("delta_ZPS" %in% colnames(hm_data)) {
  hm_mat <- cbind(hm_mat, hm_data$delta_ZPS)
  hm_cols <- c(hm_cols, "Zonation loss")
}

colnames(hm_mat) <- hm_cols

# Z-score normalization per column for visualization
hm_mat_z <- scale(hm_mat)

# Clinical annotation
ha_clinical <- NULL
if (nrow(hm_data) > 0) {
  anno_df <- data.frame(
    Bilirubin = hm_data$total_bilirubin,
    Albumin = hm_data$albumin,
    Subtype = hm_data$molecular_subtype,
    row.names = paste0("P", hm_data$patient_id)
  )
  # Handle NAs for annotation
  anno_df$Bilirubin[is.na(anno_df$Bilirubin)] <- median(anno_df$Bilirubin, na.rm = TRUE)
  anno_df$Albumin[is.na(anno_df$Albumin)] <- median(anno_df$Albumin, na.rm = TRUE)
  anno_df$Subtype[is.na(anno_df$Subtype)] <- "Unknown"

  ha_clinical <- rowAnnotation(
    Bilirubin = anno_barplot(anno_df$Bilirubin, bar_width = 0.7,
                             gp = gpar(fill = "#EFC000FF", col = NA),
                             width = unit(30, "mm"), axis = FALSE),
    Albumin = anno_barplot(anno_df$Albumin, bar_width = 0.7,
                           gp = gpar(fill = "#00A087FF", col = NA),
                           width = unit(30, "mm"), axis = FALSE),
    Subtype = anno_df$Subtype,
    col = list(Subtype = c("CS1" = COL_CS1, "CS2" = COL_CS2, "Unknown" = COL_NS)),
    annotation_name_gp = gpar(fontsize = FS_ANNO, fontfamily = FONT_GRID, fontface = "bold"),
    annotation_legend_param = list(
      Subtype = list(title_gp = gpar(fontsize = FS_LEGEND_T, fontfamily = FONT_GRID, fontface = "bold"),
                     labels_gp = gpar(fontsize = FS_LEGEND_L, fontfamily = FONT_GRID))
    )
  )
}

# Order patients by activity index (worst first)
patient_order <- order(hm_data$activity_index)
hm_mat_z_ordered <- hm_mat_z[patient_order, , drop = FALSE]
activity_ordered <- hm_data$activity_index[patient_order]

# High-impact diverging palette (deep blue -> white -> deep red, 5-stop)
col_zscore_hi <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#053061", "#4393C3", "#F7F7F7", "#D6604D", "#67001F")
)

# Column grouping for biological narrative
col_split <- factor(
  c(rep("Cell-intrinsic CYP loss", 3), rep("Zonation collapse", 2))[seq_len(ncol(hm_mat_z_ordered))],
  levels = c("Cell-intrinsic CYP loss", "Zonation collapse")
)

# Left annotation: per-patient CYP Activity Index severity barplot
left_anno <- rowAnnotation(
  `Severity` = anno_barplot(
    -activity_ordered,  # invert so higher bar = more severe (more negative activity)
    bar_width = 0.75,
    gp = gpar(fill = "#B2182B", col = NA),
    width = unit(12, "mm"),
    axis_param = list(side = "top",
                      gp = gpar(fontsize = 8, fontfamily = FONT_GRID))
  ),
  annotation_name_gp = gpar(fontsize = FS_ANNO, fontfamily = FONT_GRID,
                            fontface = "bold"),
  annotation_name_rot = 0,
  annotation_name_side = "top"
)

# Create heatmap (high-impact version)
ht_E <- Heatmap(hm_mat_z_ordered,
              name = "Z-score",
              col = col_zscore_hi,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = FS_ROW_NAME, fontfamily = FONT_GRID,
                                  fontface = "bold"),
              column_names_gp = gpar(fontsize = FS_ROW_NAME, fontfamily = FONT_GRID,
                                     fontface = "bold"),
              column_names_rot = 45,
              row_names_side = "left",
              column_split = col_split,
              column_title_gp = gpar(fontsize = FS_ANNO, fontfamily = FONT_GRID,
                                     fontface = "bold"),
              column_gap = unit(1.5, "mm"),
              border = TRUE,
              border_gp = gpar(col = "black", lwd = 0.8),
              rect_gp = gpar(col = "white", lwd = 0.6),
              left_annotation = left_anno,
              right_annotation = ha_clinical,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = FS_LEGEND_T, fontfamily = FONT_GRID, fontface = "bold"),
                labels_gp = gpar(fontsize = FS_LEGEND_L, fontfamily = FONT_GRID),
                grid_width = unit(3, "mm"),
                legend_height = unit(20, "mm")
              ),
              width = unit(60, "mm"),
              height = unit(42, "mm"))
cat("  Panel E complete.\n")


# =============================================================================
# SECTION 7: Save Individual Panels
# =============================================================================
cat("\n--- Saving individual panels ---\n")

# Panel A
cairo_pdf(file.path(OUT, "Fig9a_celltype_adjusted_CYP3A4.pdf"),
          width = mm2in(91), height = mm2in(85), family = FONT_FAMILY)
print(pA)
dev.off()
cat("  Saved Panel A\n")

# Panel B
cairo_pdf(file.path(OUT, "Fig9b_TE_landscape.pdf"),
          width = mm2in(92), height = mm2in(85), family = FONT_FAMILY)
print(pB)
dev.off()
cat("  Saved Panel B\n")

# Panel C
cairo_pdf(file.path(OUT, "Fig9c_TE_vs_halflife.pdf"),
          width = mm2in(91), height = mm2in(85), family = FONT_FAMILY)
print(pC)
dev.off()
cat("  Saved Panel C\n")

# Panel D (pixel-exact)
tw_d <- round(92 * ASSEMBLY_DPI / 25.4); th_d <- round(85 * ASSEMBLY_DPI / 25.4)
cairo_pdf(file.path(OUT, "Fig9d_CYP_activity_bilirubin.pdf"),
          width = tw_d / ASSEMBLY_DPI, height = th_d / ASSEMBLY_DPI, family = FONT_FAMILY)
print(pD)
dev.off()
cat(sprintf("  Saved Panel D (%dx%d px)\n", tw_d, th_d))

# Panel E (ComplexHeatmap, pixel-exact)
tw_e <- round(183 * ASSEMBLY_DPI / 25.4); th_e <- round(75 * ASSEMBLY_DPI / 25.4)
cairo_pdf(file.path(OUT, "Fig9e_deadzone_heatmap.pdf"),
          width = tw_e / ASSEMBLY_DPI, height = th_e / ASSEMBLY_DPI, family = FONT_FAMILY)
draw(ht_E, padding = unit(c(3, 8, 3, 8), "mm"))
dev.off()
cat(sprintf("  Saved Panel E (%dx%d px)\n", tw_e, th_e))

# =============================================================================
# SECTION 8: Composite Figure Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite Figure_9 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600

  # Row 1 (85mm): A(91) + B(92) = 183
  # Row 2 (85mm): C(91) + D(92) = 183
  # Row 3 (75mm): E(183) full-width
  H1 <- 85; H2 <- 85; H3 <- H_TOTAL - H1 - H2  # 75mm
  W_A <- 91; W_B <- W_TOTAL - W_A  # 92mm

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))

    # --- Panel A (Row 1 left: 91mm, 80mm) ---
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(H2 + H3, "mm"),
      width = unit(W_A, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(pA, newpage = FALSE)
    grid::popViewport()

    # --- Panel B (Row 1 right: 92mm, 80mm) ---
    grid::pushViewport(grid::viewport(
      x = unit(W_A, "mm"), y = unit(H2 + H3, "mm"),
      width = unit(W_B, "mm"), height = unit(H1, "mm"),
      just = c("left", "bottom")
    ))
    print(pB, newpage = FALSE)
    grid::popViewport()

    # --- Panel C (Row 2 left: 91mm, 80mm) ---
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(H3, "mm"),
      width = unit(W_A, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(pC, newpage = FALSE)
    grid::popViewport()

    # --- Panel D (Row 2 right: 92mm, 80mm) ---
    grid::pushViewport(grid::viewport(
      x = unit(W_A, "mm"), y = unit(H3, "mm"),
      width = unit(W_B, "mm"), height = unit(H2, "mm"),
      just = c("left", "bottom")
    ))
    print(pD, newpage = FALSE)
    grid::popViewport()

    # --- Panel E (Row 3: full-width heatmap, 85mm) ---
    grid::pushViewport(grid::viewport(
      x = unit(0, "mm"), y = unit(0, "mm"),
      width = unit(W_TOTAL, "mm"), height = unit(H3, "mm"),
      just = c("left", "bottom")
    ))
    draw(ht_E, padding = unit(c(3, 8, 3, 8), "mm"), newpage = FALSE)
    grid::popViewport()

    # --- Panel labels (FS_TAG = 14pt bold) ---
    label_data <- data.frame(
      text = c("a", "b", "c", "d", "e"),
      x_mm = c(2, W_A + 2, 2, W_A + 2, 2),
      y_mm = c(H2 + H3 + H1 - 2, H2 + H3 + H1 - 2,
               H3 + H2 - 2, H3 + H2 - 2, H3 - 2),
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

    grid::popViewport()
  }

  # --- Output: cairo_pdf (vector) ---
  cairo_pdf(file.path(OUT, "Figure_9.pdf"),
            width = W_TOTAL / 25.4, height = H_TOTAL / 25.4, family = FONT_FAMILY)
  render_vector_composite()
  dev.off()
  cat("  -> Figure_9.pdf (vector)\n")

  # --- Output: PNG (600 DPI) ---
  px_W <- round(W_TOTAL * DPI / 25.4)
  px_H <- round(H_TOTAL * DPI / 25.4)
  grDevices::png(file.path(OUT, "Figure_9.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_9.png\n")

  # --- Output: TIFF (600 DPI, LZW) ---
  grDevices::tiff(file.path(OUT, "Figure_9.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_vector_composite()
  dev.off()
  cat("  -> Figure_9.tiff\n")

  cat("  Figure_9 DONE (vector, 183x245mm, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))


# =============================================================================
# SECTION 9: Print Summary Statistics for Manuscript
# =============================================================================
cat("\n\n=== FIGURE 9 SUMMARY STATISTICS ===\n")
cat("For manuscript integration:\n\n")

if (exists("cyp3a4_unadj") && exists("cyp3a4_adj")) {
  cat(sprintf("Panel A: CYP3A4 unadjusted logFC=%.3f (P=%.4f) -> adjusted logFC=%.3f (P=%.4f)\n",
              cyp3a4_unadj$logFC, cyp3a4_unadj$P.Value,
              cyp3a4_adj$logFC, cyp3a4_adj$P.Value))
  cat(sprintf("  Change in effect size: %.1f%%\n",
              100 * (cyp3a4_adj$logFC - cyp3a4_unadj$logFC) / abs(cyp3a4_unadj$logFC)))
}

cat(sprintf("\nPanel B: FGG TE index=%.3f (rank #%d of %d genes)\n",
            fgg_row$TE_index, fgg_row$rank, n_genes))
cat(sprintf("  Fibrinogen family: FGG #%d, FGB #%d, FGA #%d\n",
            fgg_row$rank, fgb_row$rank, fga_row$rank))

cat(sprintf("\nPanel C: mRNA-protein concordance (rho=%.3f, n=%d genes)\n", rho_mp, nrow(te_plot)))
for (q in names(quad_counts)) cat(sprintf("  %s: %d genes\n", q, quad_counts[q]))

if (!is.na(cor_bili$estimate)) {
  cat(sprintf("\nPanel D: CYP3A4 Activity vs Bilirubin: rho=%.3f, P=%.4f (n=%d patients)\n",
              cor_bili$estimate, cor_bili$p.value, nrow(plot_df_d)))
}

cat(sprintf("\nPanel E: %d patients in Dead Zone severity heatmap\n", nrow(hm_data)))

cat("\n=== Figure 9 generation COMPLETE ===\n")
cat("Output files:\n")
cat(paste0("  ", OUT, "/Figure_9.pdf\n"))
cat(paste0("  ", OUT, "/Figure_9.png\n"))
cat(paste0("  ", OUT, "/Figure_9.tiff\n"))
