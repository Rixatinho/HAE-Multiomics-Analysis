#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_05_standalone.R
# Supplementary Figure 5: RBP Post-transcriptional Buffering +
#                          Immune Checkpoint Landscape
# HAE (Hepatic Alveolar Echinococcosis) Multi-omics Study
# Target: Nature Communications - vector PDF + PNG/TIFF (600 DPI) composite
# =============================================================================
# 8 panels (four-row layout, 183x245mm canvas):
#   A = RBP mRNA vs protein scatter (72 RBPs detected in both layers,
#       32 / 44.4% discordant; OASL/OAS1/G3BP2 labelled at mag > 0.45) — 91x60
#   B = Pathway discordance comparison after biological reclassification:
#       11/11 (100%) immune Hallmark vs 0/8 (0%) metabolic — 92x60
#   C = Discordant Hallmark NES horizontal bar (18 pathways, all TC<0 / PR>0;
#       max split = IFN-gamma response TC=-2.54 / PR=+2.07) — 91x65
#   D = Top 15 discordant RBPs ranked by disc_magnitude (OAS family dominates:
#       OASL/OAS1/OAS3/OAS2 at ranks 1, 2, 5, 8) — 92x65
#   E = Immune checkpoint expression (top 12 of 22 detected molecules;
#       VISTA/VSIR Protein log2FC = +0.521, P = 0.033 marked by *) — 91x60
#   F = Immunophenoscore (IPS) Normal vs peri-lesional boxplot
#       (Wilcoxon P = 0.347, NS) — 92x60
#   G = IPS sub-component box (MHC P=0.93, Effector P=0.18,
#       Checkpoint P=0.59, Suppressor P=0.80; all NS) — 91x60
#   H = ssGSEA infiltration box for 6 signatures (Antigen presentation,
#       Cytotoxic T, Myeloid suppression, T-cell activation, T-cell exhaustion,
#       Treg; all Wilcoxon P >= 0.32, NS) — 92x60
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_05_standalone.R
# =============================================================================

cat("=== Supplementary Figure 5: RBP Immune Buffering + Immune Checkpoint Landscape ===\n")
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
  library(ggsci)
  library(ggpubr)
  library(ggrepel)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)
cat("  Arial registered in PDF/PostScript font databases\n")

# =============================================================================
# SECTION 2: Project Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
RBP_DIR  <- file.path(RES, "rbp_immune_buffering")
CKPT_DIR <- file.path(RES, "immune_checkpoint_landscape")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_05")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 3: Font & Dimension Constants
# =============================================================================
FONT_FAMILY <- "Arial"
MM_PER_INCH <- 25.4

W_SINGLE <- 89; W_DOUBLE <- 183
H_STD <- 75; H_TALL <- 90

ASSEMBLY_DPI <- 600   # pixel-exact rendering constant

mm2in <- function(mm) mm / MM_PER_INCH
FS_TITLE <- 10; FS_AXIS_TITLE <- 9; FS_AXIS_TEXT <- 8
FS_LEGEND_T <- 8; FS_LEGEND_L <- 8; FS_TAG <- 12

# =============================================================================
# SECTION 4: Color Palette
# =============================================================================
COL_CONCORDANT <- "#868686FF"
COL_DISCORDANT <- "#CD534CFF"
COL_TC  <- "#0073C2FF"
COL_PR  <- "#CD534CFF"
COL_IMMUNE   <- "#CD534CFF"
COL_METABOLIC <- "#0073C2FF"
COL_NORMAL   <- "#7AA6DCFF"
COL_ADJACENT <- "#CD534CFF"
COL_SIG      <- "#C0392B"
PAL_CAT <- pal_jco("default")(10)

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
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5,
                              family = FONT_FAMILY,
                              margin = margin(b = 1.5, unit = "mm")),
    plot.subtitle = element_text(size = 7.5, color = "grey40", hjust = 0.5,
                                 family = FONT_FAMILY,
                                 margin = margin(b = 1.5, unit = "mm")),
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
  # Curated display map: preserves biological acronyms (TNF, IL, KRAS, G2/M, PI3K/AKT/mTOR,
  # E2F, UV, p53, IFN-gamma/alpha) for Hallmark pathways visible in Panel C.
  display_map <- c(
    HALLMARK_COMPLEMENT                  = "Complement",
    HALLMARK_INTERFERON_GAMMA_RESPONSE   = "Interferon-\u03b3 response",
    HALLMARK_ALLOGRAFT_REJECTION         = "Allograft rejection",
    HALLMARK_INTERFERON_ALPHA_RESPONSE   = "Interferon-\u03b1 response",
    HALLMARK_INFLAMMATORY_RESPONSE       = "Inflammatory response",
    HALLMARK_COAGULATION                 = "Coagulation",
    HALLMARK_APOPTOSIS                   = "Apoptosis",
    HALLMARK_TNFA_SIGNALING_VIA_NFKB     = "TNF-\u03b1 signaling via NF-\u03baB",
    HALLMARK_ANDROGEN_RESPONSE           = "Androgen response",
    HALLMARK_IL6_JAK_STAT3_SIGNALING     = "IL-6/JAK/STAT3 signaling",
    HALLMARK_PROTEIN_SECRETION           = "Protein secretion",
    HALLMARK_KRAS_SIGNALING_UP           = "KRAS signaling up",
    HALLMARK_G2M_CHECKPOINT              = "G2/M checkpoint",
    HALLMARK_IL2_STAT5_SIGNALING         = "IL-2/STAT5 signaling",
    HALLMARK_PI3K_AKT_MTOR_SIGNALING     = "PI3K/AKT/mTOR signaling",
    HALLMARK_E2F_TARGETS                 = "E2F targets",
    HALLMARK_UV_RESPONSE_UP              = "UV response up",
    HALLMARK_P53_PATHWAY                 = "p53 pathway"
  )
  out <- ifelse(x %in% names(display_map), display_map[x], NA_character_)
  # Generic fallback for any pathway not in map: strip prefix, replace underscores,
  # then convert to sentence case (first letter upper, rest lower) to avoid mangling acronyms.
  needs_fallback <- is.na(out)
  if (any(needs_fallback)) {
    fb <- gsub("^HALLMARK_|^REACTOME_|^KEGG_|^GOBP_|^GO_|^WP_", "", x[needs_fallback])
    fb <- gsub("_", " ", fb)
    fb <- tolower(fb)
    fb <- paste0(toupper(substr(fb, 1, 1)), substr(fb, 2, nchar(fb)))
    out[needs_fallback] <- fb
  }
  unname(out)
}

# =============================================================================
# SECTION 7: Data Loading
# =============================================================================
cat("  Loading data...\n")

# --- RBP data ---
rbp <- read.csv(file.path(RBP_DIR, "rbp_expression_summary.csv"), stringsAsFactors = FALSE)
rbp_both <- rbp %>% filter(detected_mRNA == "True" & detected_protein == "True")
cat("    RBPs detected in both layers:", nrow(rbp_both), "\n")
cat("    Discordant:", sum(rbp_both$is_discordant == "True"), "\n")

# --- Pathway discordance ---
pathways <- read.csv(file.path(RBP_DIR, "immune_pathway_discordance_rates.csv"), stringsAsFactors = FALSE)
# Convert concordant column to logical (CSV stores as "True"/"False" strings)
pathways$concordant_lgl <- pathways$concordant == "True" | pathways$concordant == TRUE

# Comprehensive reclassification of pathway_type by biological function
# (the source CSV labels only 5 of 30 pathways as Immune/Metabolic, leaving
# 14 functionally-immune/metabolic Hallmarks under "Other"; reclassify here
# for adequate denominators in Panel B)
.is_immune <- function(name) {
  grepl("INTERFERON|COMPLEMENT|ALLOGRAFT|INFLAMMATORY|COAGULATION|APOPTOSIS|TNFA|IL6_JAK|IL2_STAT5|TCR_|_TCR|JAK_STAT|IMMUNE",
        name, ignore.case = TRUE)
}
.is_metabolic <- function(name) {
  grepl("METABOLISM|XENOBIOTIC|BILE_ACID|FATTY_ACID|HEME|PEROXISOME|ADIPOGENESIS|OXIDATIVE_PHOSPHORYLATION|CHOLESTEROL|GLYCOLYSIS|MTORC1|GLUCONEOGENESIS",
        name, ignore.case = TRUE)
}
pathways$pathway_type <- ifelse(.is_immune(pathways$pathway), "Immune",
                          ifelse(.is_metabolic(pathways$pathway), "Metabolic", "Other"))

# --- Checkpoint expression ---
ckpt <- read.csv(file.path(CKPT_DIR, "checkpoint_expression.csv"), stringsAsFactors = FALSE)

# --- Immunophenoscore ---
ips <- read.csv(file.path(CKPT_DIR, "immunophenoscore.csv"), stringsAsFactors = FALSE)

# --- Immune cell infiltration scores ---
score_files <- c("Antigen_presentation_scores.csv", "Cytotoxic_T_scores.csv",
                 "Myeloid_suppression_scores.csv", "T-cell_activation_scores.csv",
                 "T-cell_exhaustion_scores.csv", "Treg_signature_scores.csv")
score_names <- c("Antigen\npresentation", "Cytotoxic T", "Myeloid\nsuppression",
                 "T-cell\nactivation", "T-cell\nexhaustion", "Treg\nsignature")
infiltration <- do.call(rbind, lapply(seq_along(score_files), function(i) {
  df <- read.csv(file.path(CKPT_DIR, score_files[i]), stringsAsFactors = FALSE)
  df$Signature <- score_names[i]
  df
}))

cat("  Data loading complete.\n\n")

# =============================================================================
# PANEL A: RBP mRNA vs Protein Scatter
# =============================================================================
cat("  Panel A: RBP mRNA vs protein scatter...\n")

rbp_both <- rbp_both %>%
  mutate(
    status = ifelse(is_discordant == "True", "Discordant", "Concordant"),
    label_gene = ifelse(is_discordant == "True" & disc_magnitude > 0.45,
                        gene, NA_character_)
  )

pA <- ggplot(rbp_both, aes(x = mRNA_logFC, y = protein_logFC, color = status)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              linewidth = 0.3, color = "grey75") +
  geom_point(size = 1.7, alpha = 0.85, stroke = 0) +
  ggrepel::geom_label_repel(
    aes(label = label_gene, fill = status),
    size = 2.7, fontface = "bold.italic", color = "white",
    label.padding = unit(0.7, "mm"), label.r = unit(0.5, "mm"),
    label.size = 0.15,
    box.padding = unit(0.6, "mm"), point.padding = unit(0.4, "mm"),
    segment.size = 0.3, segment.color = "grey35",
    min.segment.length = 0, force = 2, force_pull = 0.5,
    max.overlaps = 30, show.legend = FALSE, family = FONT_FAMILY) +
  scale_color_manual(values = c("Concordant" = COL_CONCORDANT,
                                "Discordant" = COL_DISCORDANT),
                     name = NULL) +
  scale_fill_manual(values = c("Concordant" = COL_CONCORDANT,
                               "Discordant" = COL_DISCORDANT),
                    guide = "none") +
  labs(x = "mRNA log\u2082 FC", y = "Protein log\u2082 FC",
       title = "RBP mRNA\u2013protein discordance",
       subtitle = "72 RBPs; 32 (44.4%) discordant") +
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(2, "mm"),
        legend.key.width  = unit(3, "mm"))

save_pdf("Supp05a_RBP_scatter.pdf", 91, 60, print(pA))

# =============================================================================
# PANEL B: Pathway Discordance Comparison
# =============================================================================
cat("  Panel B: Pathway discordance comparison...\n")

path_summary <- pathways %>%
  filter(pathway_type %in% c("Immune", "Metabolic")) %>%
  group_by(pathway_type) %>%
  summarise(
    n_total = n(),
    n_discordant = sum(!concordant_lgl),
    pct_discordant = n_discordant / n_total * 100,
    .groups = "drop"
  )

pB <- ggplot(path_summary, aes(x = pathway_type, y = pct_discordant,
                               fill = pathway_type)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.3, alpha = 0.9) +
  geom_text(aes(label = sprintf("%d / %d", n_discordant, n_total),
                y = pct_discordant + 5),
            size = 2.7, family = FONT_FAMILY, fontface = "bold") +
  scale_fill_manual(values = c("Immune" = COL_IMMUNE, "Metabolic" = COL_METABOLIC)) +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Discordance rate",
       title = "Pathway discordance: immune vs metabolic",
       subtitle = "TC\u2193 / PR\u2191 pattern") +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank())

save_pdf("Supp05b_pathway_discordance.pdf", 92, 60, print(pB))

# =============================================================================
# PANEL C: Discordant Hallmark NES Barplot
# =============================================================================
cat("  Panel C: Discordant Hallmark NES barplot...\n")

disc_paths <- pathways %>%
  filter(!concordant_lgl, database == "Hallmark") %>%
  arrange(TC_NES) %>%
  mutate(pathway_label = clean_pathway(pathway))

disc_long <- disc_paths %>%
  select(pathway_label, TC_NES, PR_NES) %>%
  pivot_longer(cols = c(TC_NES, PR_NES), names_to = "Layer", values_to = "NES") %>%
  mutate(Layer = ifelse(Layer == "TC_NES", "Transcriptomics", "Proteomics"))

disc_long$pathway_label <- factor(disc_long$pathway_label,
                                  levels = rev(disc_paths$pathway_label))

pC <- ggplot(disc_long, aes(x = NES, y = pathway_label, fill = Layer)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("Transcriptomics" = COL_TC, "Proteomics" = COL_PR)) +
  labs(x = "Normalised Enrichment Score (NES)", y = NULL,
       title = "Discordant Hallmark pathways",
       subtitle = "TC\u2193 vs PR\u2191 opposing enrichment") +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 8))

save_pdf("Supp05c_hallmark_NES.pdf", 91, 65, print(pC))

# =============================================================================
# PANEL D: Top 12 Discordant RBPs
# =============================================================================
cat("  Panel D: Top 12 discordant RBPs...\n")

top15 <- rbp_both %>%
  filter(is_discordant == "True") %>%
  arrange(desc(disc_magnitude)) %>%
  head(15)

top15_long <- top15 %>%
  select(gene, mRNA_logFC, protein_logFC) %>%
  pivot_longer(cols = c(mRNA_logFC, protein_logFC), names_to = "Layer", values_to = "logFC") %>%
  mutate(Layer = ifelse(Layer == "mRNA_logFC", "mRNA", "Protein"))

top15_long$gene <- factor(top15_long$gene, levels = rev(top15$gene))

pD <- ggplot(top15_long, aes(x = logFC, y = gene, fill = Layer)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("mRNA" = COL_TC, "Protein" = COL_PR)) +
  labs(x = "log\u2082FC (Adjacent/Normal)", y = NULL,
       title = "Top 15 discordant RBPs",
       subtitle = "Ranked by discordance magnitude") +
  theme(legend.position = "top")

save_pdf("Supp05d_top_RBPs.pdf", 92, 65, print(pD))

# =============================================================================
# PANEL E: Checkpoint Expression
# =============================================================================
cat("  Panel E: Checkpoint expression...\n")

ckpt_plot <- ckpt %>%
  filter(!is.na(RNA_logFC) | !is.na(Protein_logFC)) %>%
  mutate(
    sig_label = case_when(
      Protein_sig == "True" ~ "Protein sig",
      RNA_sig == "True" ~ "RNA sig",
      TRUE ~ "NS"
    ),
    max_abs = pmax(abs(RNA_logFC), abs(Protein_logFC), na.rm = TRUE)
  )

n_total_ckpt <- nrow(ckpt_plot)
ckpt_top <- ckpt_plot %>% arrange(desc(max_abs)) %>% head(12)

ckpt_long <- ckpt_top %>%
  select(Checkpoint, RNA_logFC, Protein_logFC, sig_label) %>%
  pivot_longer(cols = c(RNA_logFC, Protein_logFC), names_to = "Layer", values_to = "logFC") %>%
  mutate(Layer = ifelse(Layer == "RNA_logFC", "RNA", "Protein")) %>%
  filter(!is.na(logFC))

# Order by RNA logFC
ckpt_order <- ckpt_top %>% arrange(RNA_logFC) %>% pull(Checkpoint)
ckpt_long$Checkpoint <- factor(ckpt_long$Checkpoint, levels = rev(ckpt_order))

pE <- ggplot(ckpt_long, aes(x = logFC, y = Checkpoint, fill = Layer)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65,
           color = "black", linewidth = 0.18, alpha = 0.9) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "black") +
  geom_point(data = ckpt_long %>% filter(sig_label != "NS"),
             aes(x = logFC), shape = 8, size = 1.6, stroke = 0.4, color = COL_SIG,
             position = position_dodge(width = 0.7), show.legend = FALSE) +
  scale_fill_manual(values = c("RNA" = COL_TC, "Protein" = COL_PR), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = "log\u2082 FC (Adjacent / Normal)", y = NULL,
       title = "Immune checkpoint expression",
       subtitle = sprintf("Top 12 of %d molecules; * P < 0.05", n_total_ckpt)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(2, "mm"),
        legend.key.width  = unit(3, "mm"),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank())

save_pdf("Supp05e_checkpoint_expr.pdf", 91, 60, print(pE))

# =============================================================================
# PANEL F: Immunophenoscore Comparison
# =============================================================================
cat("  Panel F: Immunophenoscore comparison...\n")

wt_ips <- wilcox.test(ips$IPS_total[ips$Group == "Normal"],
                      ips$IPS_total[ips$Group == "Adjacent"])
p_label <- ifelse(wt_ips$p.value < 0.05,
                  sprintf("P = %.3f", wt_ips$p.value),
                  sprintf("P = %.2f (NS)", wt_ips$p.value))

pF <- ggplot(ips, aes(x = Group, y = IPS_total, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.4,
               color = "grey20", alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.75, shape = 21,
              color = "grey20", stroke = 0.25,
              aes(fill = Group), show.legend = FALSE) +
  scale_fill_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT)) +
  annotate("text", x = 1.5, y = max(ips$IPS_total) * 1.08,
           label = p_label, size = 2.7, family = FONT_FAMILY, fontface = "bold") +
  labs(x = NULL, y = "Immunophenoscore (IPS)",
       title = "IPS: Normal vs peri-lesional",
       subtitle = "Wilcoxon rank-sum test") +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank())

save_pdf("Supp05f_immunophenoscore.pdf", 92, 60, print(pF))

# =============================================================================
# PANEL G: IPS Component Analysis
# =============================================================================
cat("  Panel G: IPS component analysis...\n")

ips_comp <- ips %>%
  select(Sample, Group, MHC_molecules, Effector_cells, Checkpoints, Suppressor_cells) %>%
  pivot_longer(cols = c(MHC_molecules, Effector_cells, Checkpoints, Suppressor_cells),
               names_to = "Component", values_to = "Score") %>%
  mutate(Component = gsub("_", " ", Component))

comp_order <- c("MHC molecules", "Effector cells", "Checkpoints", "Suppressor cells")
ips_comp$Component <- factor(ips_comp$Component, levels = comp_order)

# Headroom for P-value labels (above max whisker)
g_max <- max(ips_comp$Score, na.rm = TRUE)
g_min <- min(ips_comp$Score, na.rm = TRUE)
g_pad <- (g_max - g_min) * 0.18

pG <- ggplot(ips_comp, aes(x = Component, y = Score, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.4,
               color = "grey20", alpha = 0.9,
               position = position_dodge(width = 0.72)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.72),
             size = 0.8, alpha = 0.65, shape = 21, color = "grey25",
             stroke = 0.2, aes(fill = Group), show.legend = FALSE) +
  scale_fill_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT), name = NULL) +
  scale_x_discrete(labels = c("MHC", "Effector", "Checkpoint", "Suppressor")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  stat_compare_means(aes(group = Group), method = "wilcox.test",
                     label = "p.format",
                     label.y = g_max + g_pad * 0.55,
                     size = 2.4, family = FONT_FAMILY, color = "grey30") +
  labs(x = NULL, y = "Component score",
       title = "IPS component breakdown",
       subtitle = "4 immunogenicity components (Wilcoxon P)") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(2, "mm"),
        legend.key.width  = unit(3, "mm"),
        axis.text.x = element_text(size = 7.5),
        panel.grid.major.x = element_blank())

save_pdf("Supp05g_IPS_components.pdf", 91, 60, print(pG))

# =============================================================================
# PANEL H: Immune Cell Infiltration Scores
# =============================================================================
cat("  Panel H: Immune cell infiltration scores...\n")

sig_order <- c("Antigen\npresentation", "Cytotoxic T", "Myeloid\nsuppression",
               "T-cell\nactivation", "T-cell\nexhaustion", "Treg\nsignature")
infiltration$Signature <- factor(infiltration$Signature, levels = sig_order)

# Short tick labels for x-axis readability
sig_short <- c("Antigen\npresentation" = "Antigen pres.",
               "Cytotoxic T"           = "Cytotoxic T",
               "Myeloid\nsuppression"  = "Myeloid supp.",
               "T-cell\nactivation"    = "T activation",
               "T-cell\nexhaustion"    = "T exhaustion",
               "Treg\nsignature"       = "Treg")

# Headroom for P-value labels (above max whisker)
h_max <- max(infiltration$Score, na.rm = TRUE)
h_min <- min(infiltration$Score, na.rm = TRUE)
h_pad <- (h_max - h_min) * 0.18

pH <- ggplot(infiltration, aes(x = Signature, y = Score, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.4,
               color = "grey20", alpha = 0.9,
               position = position_dodge(width = 0.72)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.72),
             size = 0.8, alpha = 0.65, shape = 21, color = "grey25",
             stroke = 0.2, aes(fill = Group), show.legend = FALSE) +
  scale_fill_manual(values = c("Normal" = COL_NORMAL, "Adjacent" = COL_ADJACENT), name = NULL) +
  scale_x_discrete(labels = sig_short) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  stat_compare_means(aes(group = Group), method = "wilcox.test",
                     label = "p.format",
                     label.y = h_max + h_pad * 0.55,
                     size = 2.4, family = FONT_FAMILY, color = "grey30") +
  labs(x = NULL, y = "Infiltration score (ssGSEA)",
       title = "Immune cell infiltration",
       subtitle = "6 signatures, Normal vs Adjacent (Wilcoxon P)") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(2, "mm"),
        legend.key.width  = unit(3, "mm"),
        axis.text.x = element_text(size = 7, angle = 25, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank())

save_pdf("Supp05h_infiltration.pdf", 92, 60, print(pH))

# =============================================================================
# =============================================================================
# Composite Assembly (Pure Vector Grid Viewport)
# =============================================================================
cat("\n--- Assembling composite SuppFig_05 (vector, 183x245mm, 600DPI) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245; DPI <- 600
  H1 <- 60; H2 <- 65; H3 <- 60; H4 <- H_TOTAL - H1 - H2 - H3
  W_L <- 91; W_R <- W_TOTAL - W_L

  render_vector_composite <- function() {
    grid::pushViewport(grid::viewport(
      width = unit(W_TOTAL, "mm"), height = unit(H_TOTAL, "mm"),
      xscale = c(0, W_TOTAL), yscale = c(0, H_TOTAL)
    ))
    # Panel A (Row1 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(pA, newpage=FALSE); grid::popViewport()
    # Panel B (Row1 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H2+H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H1,"mm"), just=c("left","bottom")))
    print(pB, newpage=FALSE); grid::popViewport()
    # Panel C (Row2 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(pC, newpage=FALSE); grid::popViewport()
    # Panel D (Row2 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H3+H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H2,"mm"), just=c("left","bottom")))
    print(pD, newpage=FALSE); grid::popViewport()
    # Panel E (Row3 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(H4,"mm"),
      width=unit(W_L,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(pE, newpage=FALSE); grid::popViewport()
    # Panel F (Row3 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(H4,"mm"),
      width=unit(W_R,"mm"), height=unit(H3,"mm"), just=c("left","bottom")))
    print(pF, newpage=FALSE); grid::popViewport()
    # Panel G (Row4 left)
    grid::pushViewport(grid::viewport(x=unit(0,"mm"), y=unit(0,"mm"),
      width=unit(W_L,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(pG, newpage=FALSE); grid::popViewport()
    # Panel H (Row4 right)
    grid::pushViewport(grid::viewport(x=unit(W_L,"mm"), y=unit(0,"mm"),
      width=unit(W_R,"mm"), height=unit(H4,"mm"), just=c("left","bottom")))
    print(pH, newpage=FALSE); grid::popViewport()
    # Labels
    label_data <- data.frame(
      text=c("a","b","c","d","e","f","g","h"),
      x_mm=c(2,W_L+2,2,W_L+2,2,W_L+2,2,W_L+2),
      y_mm=c(H2+H3+H4+H1-2, H2+H3+H4+H1-2, H3+H4+H2-2, H3+H4+H2-2,
             H4+H3-2, H4+H3-2, H4-2, H4-2), stringsAsFactors=FALSE)
    for(i in seq_len(nrow(label_data))) {
      grid::grid.text(label=label_data$text[i],
        x=unit(label_data$x_mm[i],"mm"), y=unit(label_data$y_mm[i],"mm"),
        just=c("left","top"),
        gp=grid::gpar(fontsize=FS_TAG, fontface="bold", fontfamily=FONT_FAMILY))
    }
    grid::popViewport()
  }

  cairo_pdf(file.path(OUT,"SuppFig_05.pdf"),
            width=W_TOTAL/25.4, height=H_TOTAL/25.4, family=FONT_FAMILY)
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_05.pdf (vector)\n")

  px_W <- round(W_TOTAL*DPI/25.4); px_H <- round(H_TOTAL*DPI/25.4)
  grDevices::png(file.path(OUT,"SuppFig_05.png"),
                 width=px_W, height=px_H, res=DPI, type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_05.png\n")

  grDevices::tiff(file.path(OUT,"SuppFig_05.tiff"),
                  width=px_W, height=px_H, res=DPI, compression="lzw", type="cairo")
  render_vector_composite(); dev.off()
  cat("  -> SuppFig_05.tiff\n")

  cat("  SuppFig_05 DONE (vector, AI-editable)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 5 rendering complete ===\n")
