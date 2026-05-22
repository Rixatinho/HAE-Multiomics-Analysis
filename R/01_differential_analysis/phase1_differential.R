#!/usr/bin/env Rscript
# ============================================================================
# Phase 1: Differential Expression Analysis
# HAE Multi-omics Integration Study
# ============================================================================
# Strategy:
#   Transcriptomics: limma-voom (paired) — better power with small n
#   Proteomics:      limma (paired)
#   Metabolomics:    limma (paired) + Wilcoxon (robustness check)
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

# ---- Paths ----
PROJECT   <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
DATA_ROOT <- "/Volumes/阿热-实验数据备份/工作/1.肝包虫/1.泡型包虫/1.2024-11-泡球蚴"
DATA_DIR  <- file.path(PROJECT, "analysis/data/processed")
OUT_DIR   <- file.path(PROJECT, "analysis/results/phase1_diff")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

col_group   <- c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")
col_volcano <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey70")

# Dual-threshold strategy for small-sample multi-omics:
#   Strict: padj < 0.05, |log2FC| > 1    (for reporting in text)
#   Relaxed: P < 0.05, |log2FC| > 0.585  (for heatmaps & pathway input)
# GSEA (Phase 2) uses full ranked list — no cutoff needed
FC_CUT       <- 1.0
FC_CUT_RELAX <- 0.585   # ~1.5-fold
PADJ_CUT     <- 0.05
P_CUT_RELAX  <- 0.05    # nominal P for exploratory screen

cat("========================================\n")
cat("Phase 1: Differential Expression Analysis\n")
cat("========================================\n\n")

extract_id <- function(x) as.integer(gsub("Normal|Adjacent", "", x))

# ============================================================================
# 1. TRANSCRIPTOMICS — limma-voom (paired design)
# ============================================================================
cat(">>> 1. Transcriptomics Differential Analysis (limma-voom, paired)\n")

# Read raw counts for limma-voom
tc_file <- file.path(DATA_ROOT,
  "1.转录组学测序/03.Result_X101SC24101695-Z02-F001_homo_sapiens/Result_X101SC24101695-Z02-F001_homo_sapiens/4.Quant/1.Count/gene_count.xls")
tc_raw <- read.delim(tc_file, check.names = FALSE)
rownames(tc_raw) <- tc_raw$gene_id

# Extract sample columns
sample_cols <- grep("^(Normal|Adjacent)", colnames(tc_raw), value = TRUE)
tc_counts <- tc_raw[, sample_cols]

# Gene name map
gene_name_map <- setNames(tc_raw$gene_name, tc_raw$gene_id)

# Keep protein-coding genes
if ("gene_biotype" %in% colnames(tc_raw)) {
  pc_genes <- rownames(tc_raw)[tc_raw$gene_biotype == "protein_coding"]
  tc_counts <- tc_counts[intersect(pc_genes, rownames(tc_counts)), ]
}

# Build paired sample set
normal_tc  <- sort(grep("^Normal", colnames(tc_counts), value = TRUE))
adjacent_tc <- sort(grep("^Adjacent", colnames(tc_counts), value = TRUE))
paired_ids <- intersect(extract_id(normal_tc), extract_id(adjacent_tc))
paired_samples <- c(paste0("Normal", paired_ids), paste0("Adjacent", paired_ids))
tc_paired <- tc_counts[, paired_samples]

cat(sprintf("  Paired patients: %d (%s)\n", length(paired_ids),
            paste(paired_ids, collapse = ", ")))

# Sample info
si <- data.frame(
  sample = colnames(tc_paired),
  group  = factor(ifelse(grepl("^Normal", colnames(tc_paired)), "Normal", "Adjacent"),
                  levels = c("Normal", "Adjacent")),
  patient = factor(extract_id(colnames(tc_paired))),
  stringsAsFactors = FALSE
)
rownames(si) <- si$sample

# DGEList + filtering
dge <- DGEList(counts = tc_paired, group = si$group)
keep <- filterByExpr(dge, group = si$group, min.count = 5)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
cat(sprintf("  After filtering: %d genes\n", nrow(dge)))

# Design: paired with duplicateCorrelation for better power
design_tc <- model.matrix(~ 0 + group, data = si)
colnames(design_tc) <- gsub("^group", "", colnames(design_tc))

# voom transformation
v <- voom(dge, design_tc, plot = FALSE)

# Use duplicateCorrelation for paired samples (more powerful than full blocking)
corfit <- duplicateCorrelation(v, design_tc, block = si$patient)
cat(sprintf("  Intra-patient correlation: %.3f\n", corfit$consensus.correlation))

# Refit voom with correlation
v <- voom(dge, design_tc, block = si$patient,
          correlation = corfit$consensus.correlation, plot = FALSE)
corfit <- duplicateCorrelation(v, design_tc, block = si$patient)

# lmFit with blocking
fit_tc <- lmFit(v, design_tc, block = si$patient,
                correlation = corfit$consensus.correlation)
contrast_tc <- makeContrasts(Adjacent - Normal, levels = design_tc)
fit2_tc <- contrasts.fit(fit_tc, contrast_tc)
fit2_tc <- eBayes(fit2_tc)

res_tc <- topTable(fit2_tc, number = Inf, sort.by = "P")
res_tc$gene_id <- rownames(res_tc)
res_tc$gene_name <- gene_name_map[res_tc$gene_id]

# Strict significance (BH-adjusted — for manuscript text)
res_tc$sig_strict <- "NS"
res_tc$sig_strict[res_tc$adj.P.Val < PADJ_CUT & res_tc$logFC > FC_CUT] <- "Up"
res_tc$sig_strict[res_tc$adj.P.Val < PADJ_CUT & res_tc$logFC < -FC_CUT] <- "Down"

# Relaxed significance (nominal P — for heatmaps/pathway input)
res_tc$significance <- "NS"
res_tc$significance[res_tc$P.Value < P_CUT_RELAX & res_tc$logFC > FC_CUT_RELAX] <- "Up"
res_tc$significance[res_tc$P.Value < P_CUT_RELAX & res_tc$logFC < -FC_CUT_RELAX] <- "Down"

n_strict_tc <- sum(res_tc$sig_strict != "NS")
n_up_tc   <- sum(res_tc$significance == "Up")
n_down_tc <- sum(res_tc$significance == "Down")
cat(sprintf("  Strict  (|log2FC|>%.1f, padj<%.2f): %d DEGs\n",
            FC_CUT, PADJ_CUT, n_strict_tc))
cat(sprintf("  Relaxed (|log2FC|>%.3f, P<%.2f): %d up + %d down = %d DEGs\n",
            FC_CUT_RELAX, P_CUT_RELAX, n_up_tc, n_down_tc, n_up_tc + n_down_tc))

# Save
write.csv(res_tc, file.path(OUT_DIR, "DEGs_Adjacent_vs_Normal.csv"), row.names = FALSE)
degs <- res_tc[res_tc$significance != "NS", ]
write.csv(degs, file.path(OUT_DIR, "DEGs_significant_relaxed.csv"), row.names = FALSE)

# Save VST matrix (from edgeR logCPM for heatmap use)
tc_logcpm <- cpm(dge, log = TRUE, prior.count = 1)
write.csv(tc_logcpm, file.path(DATA_DIR, "transcriptomics_logcpm_paired.csv"))

# Volcano plot (use nominal P on y-axis for visibility with small n)
top_genes_tc <- rbind(
  head(degs[degs$significance == "Up", ], 15),
  head(degs[degs$significance == "Down", ], 15)
)

p_volcano_tc <- ggplot(res_tc, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = significance), size = 0.8, alpha = 0.6) +
  geom_text_repel(
    data = top_genes_tc,
    aes(label = gene_name),
    size = 2.8, max.overlaps = 25, fontface = "italic",
    segment.size = 0.3, segment.color = "grey50"
  ) +
  scale_color_manual(values = col_volcano) +
  geom_vline(xintercept = c(-FC_CUT_RELAX, FC_CUT_RELAX), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(P_CUT_RELAX), linetype = "dashed", color = "grey40") +
  labs(
    title = "Transcriptomics: Adjacent vs Normal (limma-voom, paired)",
    subtitle = sprintf("%d Up / %d Down (|log2FC|>%.2f, P<%.2f)",
                       n_up_tc, n_down_tc, FC_CUT_RELAX, P_CUT_RELAX),
    x = "log2(Fold Change)", y = "-log10(nominal P-value)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "volcano_transcriptomics.pdf"), p_volcano_tc,
       width = 7, height = 6)
cat("  Volcano plot saved.\n\n")


# ============================================================================
# 2. PROTEOMICS — limma (paired design)
# ============================================================================
cat(">>> 2. Proteomics Differential Analysis (limma, paired)\n")

pr_norm <- as.matrix(read.csv(file.path(DATA_DIR, "proteomics_log2_norm.csv"),
                               row.names = 1, check.names = FALSE))

# Read annotation from raw data directly
pr_raw <- read.delim(file.path(DATA_ROOT,
  "2.蛋白质组学/Result-X101SC24101695_Z01_J001_B1_43/3.ProtExprQuantification/all_sample.xls"),
  check.names = FALSE)
gene_map_pr <- setNames(pr_raw$Gene, pr_raw$Protein)

# Sample info
si_pr <- data.frame(
  sample = colnames(pr_norm),
  group  = factor(ifelse(grepl("^Normal", colnames(pr_norm)), "Normal", "Adjacent"),
                  levels = c("Normal", "Adjacent")),
  patient = factor(extract_id(colnames(pr_norm))),
  stringsAsFactors = FALSE
)

# Paired limma with duplicateCorrelation
# NOTE: For non-voom data (PR/MT already log2-normalized), single iteration is sufficient
# as there are no observation weights to update between iterations
design_pr <- model.matrix(~ 0 + group, data = si_pr)
colnames(design_pr) <- gsub("^group", "", colnames(design_pr))

corfit_pr <- duplicateCorrelation(pr_norm, design_pr, block = si_pr$patient)
cat(sprintf("  Intra-patient correlation: %.3f\n", corfit_pr$consensus.correlation))

fit_pr <- lmFit(pr_norm, design_pr, block = si_pr$patient,
                correlation = corfit_pr$consensus.correlation)
contrast_pr <- makeContrasts(Adjacent - Normal, levels = design_pr)
fit2_pr <- contrasts.fit(fit_pr, contrast_pr)
fit2_pr <- eBayes(fit2_pr)

res_pr <- topTable(fit2_pr, number = Inf, sort.by = "P")
res_pr$Protein <- rownames(res_pr)
res_pr$gene_name <- gene_map_pr[res_pr$Protein]

# Dual thresholds
res_pr$sig_strict <- "NS"
res_pr$sig_strict[res_pr$adj.P.Val < PADJ_CUT & res_pr$logFC > FC_CUT] <- "Up"
res_pr$sig_strict[res_pr$adj.P.Val < PADJ_CUT & res_pr$logFC < -FC_CUT] <- "Down"

res_pr$significance <- "NS"
res_pr$significance[res_pr$P.Value < P_CUT_RELAX & res_pr$logFC > FC_CUT_RELAX] <- "Up"
res_pr$significance[res_pr$P.Value < P_CUT_RELAX & res_pr$logFC < -FC_CUT_RELAX] <- "Down"

n_strict_pr <- sum(res_pr$sig_strict != "NS")
n_up_pr   <- sum(res_pr$significance == "Up")
n_down_pr <- sum(res_pr$significance == "Down")
cat(sprintf("  Total tested: %d proteins\n", nrow(res_pr)))
cat(sprintf("  Strict  (|log2FC|>%.1f, padj<%.2f): %d DEPs\n",
            FC_CUT, PADJ_CUT, n_strict_pr))
cat(sprintf("  Relaxed (|log2FC|>%.3f, P<%.2f): %d up + %d down = %d DEPs\n",
            FC_CUT_RELAX, P_CUT_RELAX, n_up_pr, n_down_pr, n_up_pr + n_down_pr))

# Save
write.csv(res_pr, file.path(OUT_DIR, "DEPs_Adjacent_vs_Normal.csv"), row.names = FALSE)
deps <- res_pr[res_pr$significance != "NS", ]
write.csv(deps, file.path(OUT_DIR, "DEPs_significant.csv"), row.names = FALSE)

# Volcano
top_genes_pr <- rbind(
  head(deps[deps$significance == "Up", ], 15),
  head(deps[deps$significance == "Down", ], 15)
)

p_volcano_pr <- ggplot(res_pr, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = significance), size = 0.8, alpha = 0.6) +
  geom_text_repel(
    data = top_genes_pr,
    aes(label = gene_name),
    size = 2.8, max.overlaps = 25, fontface = "italic",
    segment.size = 0.3, segment.color = "grey50"
  ) +
  scale_color_manual(values = col_volcano) +
  geom_vline(xintercept = c(-FC_CUT_RELAX, FC_CUT_RELAX), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(P_CUT_RELAX), linetype = "dashed", color = "grey40") +
  labs(
    title = "Proteomics: Adjacent vs Normal (limma, paired)",
    subtitle = sprintf("%d Up / %d Down (|log2FC|>%.2f, P<%.2f)",
                       n_up_pr, n_down_pr, FC_CUT_RELAX, P_CUT_RELAX),
    x = "log2(Fold Change)", y = "-log10(nominal P-value)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "volcano_proteomics.pdf"), p_volcano_pr,
       width = 7, height = 6)
cat("  Volcano plot saved.\n\n")


# ============================================================================
# 3. METABOLOMICS — limma (paired) + Wilcoxon
# ============================================================================
cat(">>> 3. Metabolomics Differential Analysis (limma, paired)\n")

met_log2 <- as.matrix(read.csv(file.path(DATA_DIR, "metabolomics_log2_merged.csv"),
                                row.names = 1, check.names = FALSE))
met_anno <- read.csv(file.path(DATA_DIR, "metabolomics_annotation.csv"),
                      row.names = 1, check.names = FALSE)

si_met <- data.frame(
  sample = colnames(met_log2),
  group  = factor(ifelse(grepl("^Normal", colnames(met_log2)), "Normal", "Adjacent"),
                  levels = c("Normal", "Adjacent")),
  patient = factor(extract_id(colnames(met_log2))),
  stringsAsFactors = FALSE
)

# Paired limma with duplicateCorrelation
# NOTE: For non-voom data (MT already log2-normalized), single iteration is sufficient
# as there are no observation weights to update between iterations
design_met <- model.matrix(~ 0 + group, data = si_met)
colnames(design_met) <- gsub("^group", "", colnames(design_met))

corfit_met <- duplicateCorrelation(met_log2, design_met, block = si_met$patient)
cat(sprintf("  Intra-patient correlation: %.3f\n", corfit_met$consensus.correlation))

fit_met <- lmFit(met_log2, design_met, block = si_met$patient,
                 correlation = corfit_met$consensus.correlation)
contrast_met <- makeContrasts(Adjacent - Normal, levels = design_met)
fit2_met <- contrasts.fit(fit_met, contrast_met)
fit2_met <- eBayes(fit2_met)

res_met <- topTable(fit2_met, number = Inf, sort.by = "P")
res_met$Compound_ID <- rownames(res_met)

# Map metabolite names
if ("Name" %in% colnames(met_anno)) {
  name_map_met <- setNames(met_anno$Name, rownames(met_anno))
  res_met$metabolite_name <- name_map_met[res_met$Compound_ID]
}

# Paired Wilcoxon as robustness check
normal_met  <- sort(grep("^Normal", colnames(met_log2), value = TRUE))
adjacent_met <- sort(grep("^Adjacent", colnames(met_log2), value = TRUE))
paired_ids_met <- intersect(extract_id(normal_met), extract_id(adjacent_met))
n_cols_met <- paste0("Normal", paired_ids_met)
a_cols_met <- paste0("Adjacent", paired_ids_met)

wilcox_p <- sapply(rownames(met_log2), function(m) {
  tryCatch({
    wilcox.test(as.numeric(met_log2[m, a_cols_met]),
                as.numeric(met_log2[m, n_cols_met]),
                paired = TRUE)$p.value
  }, error = function(e) NA)
})
wilcox_padj <- p.adjust(wilcox_p, method = "BH")
res_met$wilcox_pval <- wilcox_p[res_met$Compound_ID]
res_met$wilcox_padj <- wilcox_padj[res_met$Compound_ID]

# Significance (dual threshold for metabolites)
# Strict: |log2FC|>1.0 (consistent with TC/PR), padj<0.05
# Relaxed: |log2FC|>0.585 (~1.5-fold), P<0.05 (for heatmaps & pathway input)
MET_FC_CUT_STRICT <- FC_CUT  # 1.0 - same as TC/PR strict threshold
MET_FC_CUT_RELAX  <- 0.585   # ~1.5-fold for relaxed/exploratory

res_met$sig_strict <- "NS"
res_met$sig_strict[res_met$adj.P.Val < PADJ_CUT & res_met$logFC > MET_FC_CUT_STRICT] <- "Up"
res_met$sig_strict[res_met$adj.P.Val < PADJ_CUT & res_met$logFC < -MET_FC_CUT_STRICT] <- "Down"

res_met$significance <- "NS"
res_met$significance[res_met$P.Value < P_CUT_RELAX & res_met$logFC > MET_FC_CUT_RELAX] <- "Up"
res_met$significance[res_met$P.Value < P_CUT_RELAX & res_met$logFC < -MET_FC_CUT_RELAX] <- "Down"

n_strict_met <- sum(res_met$sig_strict != "NS")
n_up_met   <- sum(res_met$significance == "Up")
n_down_met <- sum(res_met$significance == "Down")
cat(sprintf("  Total tested: %d metabolites\n", nrow(res_met)))
cat(sprintf("  Strict  (|log2FC|>%.3f, padj<%.2f): %d DEMs\n",
            MET_FC_CUT_STRICT, PADJ_CUT, n_strict_met))
cat(sprintf("  Relaxed (|log2FC|>%.3f, P<%.2f): %d up + %d down = %d DEMs\n",
            MET_FC_CUT_RELAX, P_CUT_RELAX, n_up_met, n_down_met, n_up_met + n_down_met))

# Save
write.csv(res_met, file.path(OUT_DIR, "DEMs_Adjacent_vs_Normal.csv"), row.names = FALSE)
dems <- res_met[res_met$significance != "NS", ]
write.csv(dems, file.path(OUT_DIR, "DEMs_significant.csv"), row.names = FALSE)

# Volcano
top_mets <- rbind(
  head(dems[dems$significance == "Up", ], 15),
  head(dems[dems$significance == "Down", ], 15)
)

p_volcano_met <- ggplot(res_met, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = significance), size = 0.8, alpha = 0.6) +
  geom_text_repel(
    data = top_mets,
    aes(label = metabolite_name),
    size = 2.5, max.overlaps = 25,
    segment.size = 0.3, segment.color = "grey50"
  ) +
  scale_color_manual(values = col_volcano) +
  geom_vline(xintercept = c(-MET_FC_CUT_RELAX, MET_FC_CUT_RELAX), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(P_CUT_RELAX), linetype = "dashed", color = "grey40") +
  labs(
    title = "Metabolomics: Adjacent vs Normal (limma, paired)",
    subtitle = sprintf("%d Up / %d Down (|log2FC|>%.2f, P<%.2f)",
                       n_up_met, n_down_met, MET_FC_CUT_RELAX, P_CUT_RELAX),
    x = "log2(Fold Change)", y = "-log10(nominal P-value)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "volcano_metabolomics.pdf"), p_volcano_met,
       width = 8, height = 6)
cat("  Volcano plot saved.\n\n")


# ============================================================================
# 4. TOP DIFFERENTIAL MOLECULES HEATMAPS (using relaxed DEGs/DEPs/DEMs)
# ============================================================================
cat(">>> 4. Heatmaps of Top Differential Molecules\n")

ann_colors <- list(Group = col_group)

# -- Transcriptomics top DEGs --
if (nrow(degs) > 0) {
  n_top <- min(30, nrow(degs))
  top_degs <- head(degs[order(degs$P.Value), ], n_top)

  hm_tc <- tc_logcpm[top_degs$gene_id, ]
  hm_tc <- t(scale(t(hm_tc)))
  # Handle NA/duplicate gene names
  rn <- top_degs$gene_name
  rn[is.na(rn)] <- top_degs$gene_id[is.na(rn)]
  rn <- make.unique(rn)
  rownames(hm_tc) <- rn

  ann_col <- data.frame(
    Group = ifelse(grepl("^Normal", colnames(hm_tc)), "Normal", "Adjacent"),
    row.names = colnames(hm_tc)
  )

  pdf(file.path(FIG_DIR, "heatmap_top_DEGs.pdf"), width = 10, height = 8)
  pheatmap(hm_tc, annotation_col = ann_col, annotation_colors = ann_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_cols = TRUE, cluster_rows = TRUE,
           fontsize_row = 8, fontsize_col = 7,
           main = sprintf("Top %d DEGs (Z-score)", n_top))
  dev.off()
  cat(sprintf("  Top %d DEGs heatmap saved.\n", n_top))
} else {
  cat("  No significant DEGs for heatmap.\n")
}

# -- Proteomics top DEPs --
if (nrow(deps) > 0) {
  n_top <- min(30, nrow(deps))
  top_deps <- head(deps[order(deps$P.Value), ], n_top)

  hm_pr <- pr_norm[top_deps$Protein, ]
  hm_pr <- t(scale(t(hm_pr)))
  rn_pr <- top_deps$gene_name
  rn_pr[is.na(rn_pr)] <- top_deps$Protein[is.na(rn_pr)]
  rn_pr <- make.unique(rn_pr)
  rownames(hm_pr) <- rn_pr

  ann_col_pr <- data.frame(
    Group = ifelse(grepl("^Normal", colnames(hm_pr)), "Normal", "Adjacent"),
    row.names = colnames(hm_pr)
  )

  pdf(file.path(FIG_DIR, "heatmap_top_DEPs.pdf"), width = 10, height = 8)
  pheatmap(hm_pr, annotation_col = ann_col_pr, annotation_colors = ann_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_cols = TRUE, cluster_rows = TRUE,
           fontsize_row = 8, fontsize_col = 7,
           main = sprintf("Top %d DEPs (Z-score)", n_top))
  dev.off()
  cat(sprintf("  Top %d DEPs heatmap saved.\n", n_top))
} else {
  cat("  No significant DEPs for heatmap.\n")
}

# -- Metabolomics top DEMs --
if (nrow(dems) > 0) {
  n_top <- min(30, nrow(dems))
  top_dems <- head(dems[order(dems$P.Value), ], n_top)

  hm_met <- met_log2[top_dems$Compound_ID, ]
  hm_met <- t(scale(t(hm_met)))
  labels_met <- top_dems$metabolite_name
  labels_met[is.na(labels_met)] <- top_dems$Compound_ID[is.na(labels_met)]
  labels_met <- substr(labels_met, 1, 45)
  labels_met <- make.unique(labels_met)
  rownames(hm_met) <- labels_met

  ann_col_met <- data.frame(
    Group = ifelse(grepl("^Normal", colnames(hm_met)), "Normal", "Adjacent"),
    row.names = colnames(hm_met)
  )

  pdf(file.path(FIG_DIR, "heatmap_top_DEMs.pdf"), width = 10, height = 8)
  pheatmap(hm_met, annotation_col = ann_col_met, annotation_colors = ann_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_cols = TRUE, cluster_rows = TRUE,
           fontsize_row = 7, fontsize_col = 7,
           main = sprintf("Top %d DEMs (Z-score)", n_top))
  dev.off()
  cat(sprintf("  Top %d DEMs heatmap saved.\n", n_top))
} else {
  cat("  No significant DEMs for heatmap.\n")
}


# ============================================================================
# 5. CROSS-OMICS OVERLAP
# ============================================================================
cat("\n>>> 5. Cross-omics Overlap Analysis\n")

if (nrow(degs) > 0 && nrow(deps) > 0) {
  common_genes <- intersect(degs$gene_name, deps$gene_name)
  common_genes <- common_genes[!is.na(common_genes)]
  cat(sprintf("  DEGs-DEPs shared gene names: %d\n", length(common_genes)))

  if (length(common_genes) > 0) {
    degs_sub <- degs[match(common_genes, degs$gene_name), ]
    deps_sub <- deps[match(common_genes, deps$gene_name), ]
    concordant <- sum(sign(degs_sub$logFC) == sign(deps_sub$logFC))
    cat(sprintf("  Concordant direction: %d / %d (%.1f%%)\n",
                concordant, length(common_genes),
                100 * concordant / length(common_genes)))

    # Save overlap table
    overlap_df <- data.frame(
      gene_name = common_genes,
      tc_logFC = degs_sub$logFC,
      tc_padj = degs_sub$adj.P.Val,
      pr_logFC = deps_sub$logFC,
      pr_padj = deps_sub$adj.P.Val,
      concordant = sign(degs_sub$logFC) == sign(deps_sub$logFC)
    )
    write.csv(overlap_df, file.path(OUT_DIR, "DEGs_DEPs_overlap.csv"), row.names = FALSE)
  }
}


# ============================================================================
# 6. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 1 SUMMARY\n")
cat("========================================\n")
cat("NOTE: With n=12-14 paired samples, BH-adjusted P-values are\n")
cat("conservative. We report both strict and relaxed thresholds.\n")
cat("GSEA (Phase 2) uses full ranked lists and is the primary approach.\n\n")
cat(sprintf("Transcriptomics (limma-voom, paired):\n"))
cat(sprintf("  Strict  (|log2FC|>%.1f, padj<%.2f): %d\n", FC_CUT, PADJ_CUT, n_strict_tc))
cat(sprintf("  Relaxed (|log2FC|>%.3f, P<%.2f): %d up + %d down = %d\n",
            FC_CUT_RELAX, P_CUT_RELAX, n_up_tc, n_down_tc, n_up_tc + n_down_tc))
cat(sprintf("\nProteomics (limma, paired):\n"))
cat(sprintf("  Strict  (|log2FC|>%.1f, padj<%.2f): %d\n", FC_CUT, PADJ_CUT, n_strict_pr))
cat(sprintf("  Relaxed (|log2FC|>%.3f, P<%.2f): %d up + %d down = %d\n",
            FC_CUT_RELAX, P_CUT_RELAX, n_up_pr, n_down_pr, n_up_pr + n_down_pr))
cat(sprintf("\nMetabolomics (limma + Wilcoxon, paired):\n"))
cat(sprintf("  Strict  (|log2FC|>%.3f, padj<%.2f): %d\n", MET_FC_CUT_STRICT, PADJ_CUT, n_strict_met))
cat(sprintf("  Relaxed (|log2FC|>%.3f, P<%.2f): %d up + %d down = %d\n",
            MET_FC_CUT_RELAX, P_CUT_RELAX, n_up_met, n_down_met, n_up_met + n_down_met))
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 1 COMPLETE.\n")
