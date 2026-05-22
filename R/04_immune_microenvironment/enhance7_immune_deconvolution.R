#!/usr/bin/env Rscript
# ============================================================================
# Enhancement 7: Immune Deconvolution & Microenvironment Analysis
# ============================================================================
# Estimate immune cell composition using ssGSEA with curated immune signatures
# and ESTIMATE-like stromal/immune scoring.

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
E1_DIR   <- file.path(PROJECT, "analysis/results/enhancement1_paired_subtyping")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhancement7_immune_deconvolution")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Enhancement 7: Immune Deconvolution\n")
cat("========================================\n\n")

# ---- 0. Load data ----
cat(">>> 0. Loading data\n")
degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)

tc_vst_raw <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                       check.names = FALSE, row.names = 1)
id2sym_tc <- setNames(degs$gene_name, degs$gene_id)
id2sym_tc <- id2sym_tc[!is.na(id2sym_tc) & id2sym_tc != ""]
tc_sym <- id2sym_tc[rownames(tc_vst_raw)]
keep <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst_raw[keep, ])
rownames(tc_mat) <- tc_sym[keep]

# Groups & subtypes
group_vec <- ifelse(grepl("^Normal", colnames(tc_mat)), "Normal", "Adjacent")
names(group_vec) <- colnames(tc_mat)

paired_sub <- read.csv(file.path(E1_DIR, "paired_subtype_K2.csv"), stringsAsFactors = FALSE)
subtype_map <- setNames(paired_sub$subtype, paired_sub$patient_id)
sample_subtype <- rep(NA, ncol(tc_mat))
names(sample_subtype) <- colnames(tc_mat)
for (s in colnames(tc_mat)) {
  pid <- gsub("^(Normal|Adjacent)", "", s)
  if (pid %in% names(subtype_map)) sample_subtype[s] <- subtype_map[pid]
}

cat(sprintf("  TC matrix: %d genes x %d samples\n", nrow(tc_mat), ncol(tc_mat)))

# ============================================================================
# 1. CURATED IMMUNE CELL SIGNATURES
# ============================================================================
cat("\n>>> 1. Building immune cell gene signatures\n")

# Comprehensive immune signatures curated from:
# - Bindea et al. 2013 (Immunity) - immune cell metagenes
# - Charoentong et al. 2017 (Cell Reports) - immunophenotype signatures
# - Davoli et al. 2017 (Science) - immune cell markers
# - Newman et al. 2015 (Nature Methods) - CIBERSORTx LM22
immune_sigs <- list(
  # Innate immune cells
  Macrophages_M1 = c("CD80", "CD86", "NOS2", "IL12A", "IL12B", "IL23A",
                      "TNF", "IL1B", "IL6", "CXCL9", "CXCL10", "CXCL11",
                      "CCL5", "IRF5", "SOCS1", "IDO1"),
  Macrophages_M2 = c("CD163", "MRC1", "CD209", "MSR1", "MARCO",
                      "IL10", "TGFB1", "CCL18", "CCL22", "ARG1",
                      "VEGFA", "MMP9", "IRF4", "STAT6"),
  Monocytes = c("CD14", "CSF1R", "FCGR1A", "LYZ", "S100A8", "S100A9",
                "VCAN", "FCN1", "CD163", "THBS1"),
  Dendritic_cells = c("CD1A", "CD1C", "CD1E", "FCER1A", "CLEC10A",
                       "CD83", "LAMP3", "CCR7", "IDO1", "CD80", "CD86",
                       "HLA-DRA", "HLA-DRB1", "BATF3"),
  Neutrophils = c("CEACAM8", "FUT4", "ITGAM", "CXCR1", "CXCR2",
                   "FCGR3B", "CSF3R", "S100A8", "S100A9", "S100A12",
                   "ELANE", "MPO", "MMP8"),
  NK_cells = c("NCAM1", "NKG7", "KLRD1", "KLRB1", "KLRC1", "KLRK1",
               "GNLY", "PRF1", "GZMA", "GZMB", "GZMH",
               "NCR1", "NCR3", "EOMES", "TBX21"),
  Mast_cells = c("KIT", "CPA3", "TPSAB1", "TPSB2", "MS4A2", "FCER1A",
                  "FCER1G", "HDC", "HPGDS", "IL4", "IL13"),
  Eosinophils = c("CCR3", "SIGLEC8", "IL5RA", "EPX", "RNASE2",
                   "CLC", "PRG2", "BMK", "ALOX15"),

  # Adaptive immune cells
  CD8_T_cells = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1",
                   "IFNG", "EOMES", "TBX21", "CXCR3", "NKG7",
                   "KLRG1", "FASLG"),
  CD4_T_helper = c("CD4", "IL2", "IL2RA", "ICOS", "CD28", "TNFRSF4",
                    "BCL6", "MAF", "BATF"),
  Th1 = c("TBX21", "IFNG", "TNF", "IL2", "CXCR3", "CCR5",
           "STAT1", "STAT4", "IL12RB1", "IL12RB2"),
  Th2 = c("GATA3", "IL4", "IL5", "IL13", "IL4R", "CCR4",
           "STAT6", "IL25", "IL33"),
  Th17 = c("RORC", "IL17A", "IL17F", "IL22", "IL23R", "CCR6",
            "STAT3", "IL6R", "IL21"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18",
            "IL10", "TGFB1", "ENTPD1", "NT5E", "TIGIT",
            "LAG3", "HAVCR2"),
  B_cells = c("CD19", "CD79A", "CD79B", "MS4A1", "CD22",
               "PAX5", "BLK", "BANK1", "CD72", "FCRL5"),
  Plasma_cells = c("SDC1", "XBP1", "IRF4", "PRDM1", "MZB1",
                    "JCHAIN", "IGHA1", "IGHG1", "IGHM"),

  # Functional signatures
  Cytotoxicity = c("PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "GNLY",
                    "NKG7", "FASLG", "IFNG"),
  Exhaustion = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
                  "TOX", "TOX2", "BTLA", "CD244", "ENTPD1"),
  Immune_checkpoint = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3",
                         "HAVCR2", "TIGIT", "BTLA", "VISTA", "IDO1"),
  Angiogenesis = c("VEGFA", "VEGFB", "VEGFC", "KDR", "FLT1", "FLT4",
                    "ANGPT1", "ANGPT2", "TEK", "PECAM1", "CDH5", "MCAM"),
  Fibrosis = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1",
                "ACTA2", "TAGLN", "TGFB1", "TGFB2", "CTGF",
                "LOXL2", "LOX", "TIMP1", "MMP2"),

  # ESTIMATE-like signatures
  Stromal_score = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL5A1",
                     "COL6A1", "FN1", "SPARC", "POSTN", "BGN", "THY1",
                     "ACTA2", "PDGFRA", "PDGFRB", "FAP", "VIM"),
  Immune_score = c("PTPRC", "CD2", "CD3D", "CD3E", "CD4", "CD8A",
                    "CD19", "CD68", "CD14", "LYZ", "NKG7", "GNLY",
                    "CCL5", "CXCL9", "CXCL10", "HLA-DRA", "HLA-DRB1")
)

# Filter to genes present in our data
for (nm in names(immune_sigs)) {
  immune_sigs[[nm]] <- intersect(immune_sigs[[nm]], rownames(tc_mat))
}
# Remove signatures with too few genes
sig_sizes <- sapply(immune_sigs, length)
immune_sigs <- immune_sigs[sig_sizes >= 3]

cat(sprintf("  %d immune signatures (min 3 genes each)\n", length(immune_sigs)))
for (nm in names(immune_sigs)) {
  cat(sprintf("    %s: %d genes\n", nm, length(immune_sigs[[nm]])))
}

# ============================================================================
# 2. ssGSEA SCORING
# ============================================================================
cat("\n>>> 2. Running ssGSEA immune cell scoring\n")

# Convert to GeneSetCollection
gs_list <- lapply(names(immune_sigs), function(nm) {
  GeneSet(immune_sigs[[nm]], setName = nm)
})
gs_collection <- GeneSetCollection(gs_list)

# Run ssGSEA using gsva() with method="ssgsea"
ssgsea_params <- ssgseaParam(tc_mat, gs_collection, normalize = TRUE)
ssgsea_scores <- gsva(ssgsea_params, verbose = FALSE)

cat(sprintf("  ssGSEA scores: %d signatures x %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))

write.csv(ssgsea_scores, file.path(OUT_DIR, "ssGSEA_immune_scores.csv"))

# ============================================================================
# 3. STATISTICAL COMPARISON: ADJACENT vs NORMAL
# ============================================================================
cat("\n>>> 3. Comparing immune cell scores: Adjacent vs Normal\n")

adj_samps <- colnames(ssgsea_scores)[grepl("^Adjacent", colnames(ssgsea_scores))]
nor_samps <- colnames(ssgsea_scores)[grepl("^Normal", colnames(ssgsea_scores))]

immune_diff <- data.frame()
for (sig in rownames(ssgsea_scores)) {
  tryCatch({
    # Build paired vectors for proper paired Wilcoxon test
    extract_id <- function(x) as.integer(gsub("Normal|Adjacent", "", x))
    adj_ids <- extract_id(adj_samps)
    nor_ids <- extract_id(nor_samps)
    paired_ids <- intersect(adj_ids, nor_ids)
    adj_p <- paste0("Adjacent", paired_ids)
    nor_p <- paste0("Normal", paired_ids)
    wt <- wilcox.test(ssgsea_scores[sig, adj_p], ssgsea_scores[sig, nor_p], paired = TRUE)
    immune_diff <- rbind(immune_diff, data.frame(
      signature = sig,
      mean_Adjacent = mean(ssgsea_scores[sig, adj_samps]),
      mean_Normal = mean(ssgsea_scores[sig, nor_samps]),
      diff = mean(ssgsea_scores[sig, adj_samps]) - mean(ssgsea_scores[sig, nor_samps]),
      pvalue = wt$p.value,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}
immune_diff$padj <- p.adjust(immune_diff$pvalue, method = "BH")
immune_diff <- immune_diff[order(immune_diff$pvalue), ]

cat(sprintf("  Significant (P<0.05): %d / %d\n",
            sum(immune_diff$pvalue < 0.05), nrow(immune_diff)))
cat(sprintf("  Significant (padj<0.05): %d\n", sum(immune_diff$padj < 0.05)))

write.csv(immune_diff, file.path(OUT_DIR, "immune_diff_Adjacent_vs_Normal.csv"), row.names = FALSE)

cat("\n  Immune signature differences:\n")
for (i in 1:nrow(immune_diff)) {
  sig <- ifelse(immune_diff$pvalue[i] < 0.05, "*", "")
  dir <- ifelse(immune_diff$diff[i] > 0, "UP", "DOWN")
  cat(sprintf("    %s: diff=%.3f (%s), P=%.4f %s\n",
              immune_diff$signature[i], immune_diff$diff[i], dir,
              immune_diff$pvalue[i], sig))
}

# ============================================================================
# 4. SUBTYPE COMPARISON: DS1 vs DS2
# ============================================================================
cat("\n>>> 4. Comparing immune scores by subtype (DS1 vs DS2)\n")

ds1_adj <- intersect(adj_samps, names(sample_subtype)[sample_subtype == "DS1" & !is.na(sample_subtype)])
ds2_adj <- intersect(adj_samps, names(sample_subtype)[sample_subtype == "DS2" & !is.na(sample_subtype)])

cat(sprintf("  DS1 Adjacent: %d, DS2 Adjacent: %d\n", length(ds1_adj), length(ds2_adj)))

immune_subtype <- data.frame()
if (length(ds1_adj) >= 3 && length(ds2_adj) >= 3) {
  for (sig in rownames(ssgsea_scores)) {
    tryCatch({
      wt <- wilcox.test(ssgsea_scores[sig, ds1_adj], ssgsea_scores[sig, ds2_adj])
      immune_subtype <- rbind(immune_subtype, data.frame(
        signature = sig,
        mean_DS1 = mean(ssgsea_scores[sig, ds1_adj]),
        mean_DS2 = mean(ssgsea_scores[sig, ds2_adj]),
        diff = mean(ssgsea_scores[sig, ds1_adj]) - mean(ssgsea_scores[sig, ds2_adj]),
        pvalue = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  immune_subtype$padj <- p.adjust(immune_subtype$pvalue, method = "BH")
  immune_subtype <- immune_subtype[order(immune_subtype$pvalue), ]

  cat(sprintf("  Significant (P<0.05): %d / %d\n",
              sum(immune_subtype$pvalue < 0.05), nrow(immune_subtype)))

  write.csv(immune_subtype, file.path(OUT_DIR, "immune_subtype_DS1_vs_DS2.csv"), row.names = FALSE)

  cat("\n  Top subtype-differential signatures:\n")
  for (i in 1:min(15, nrow(immune_subtype))) {
    sig <- ifelse(immune_subtype$pvalue[i] < 0.05, "*", "")
    cat(sprintf("    %s: DS1=%.3f, DS2=%.3f, P=%.4f %s\n",
                immune_subtype$signature[i], immune_subtype$mean_DS1[i],
                immune_subtype$mean_DS2[i], immune_subtype$pvalue[i], sig))
  }
}

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================
cat("\n>>> 5. Generating visualizations\n")

# ---- 5a. Immune score heatmap ----
samp_order <- c(sort(nor_samps), sort(adj_samps))
samp_order <- intersect(samp_order, colnames(ssgsea_scores))
hm_data <- ssgsea_scores[, samp_order]

anno_col <- data.frame(
  Group = ifelse(grepl("^Normal", samp_order), "Normal", "Adjacent")
)
rownames(anno_col) <- samp_order
for (s in samp_order) {
  if (!is.na(sample_subtype[s])) anno_col[s, "Subtype"] <- sample_subtype[s]
}

pdf(file.path(FIG_DIR, "immune_score_heatmap.pdf"), width = 14, height = 12)
pheatmap(hm_data,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         annotation_col = anno_col,
         annotation_colors = list(
           Group = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35"),
           Subtype = c("DS1" = "#3C5488", "DS2" = "#F39B7F")
         ),
         cluster_cols = FALSE, cluster_rows = TRUE,
         show_colnames = TRUE, fontsize_row = 8, fontsize_col = 7,
         scale = "row",
         main = "Immune Cell & Microenvironment Scores (ssGSEA)")
dev.off()
cat("  Immune heatmap saved.\n")

# ---- 5b. Boxplots for significant signatures ----
sig_sigs <- immune_diff$signature[immune_diff$pvalue < 0.1]
if (length(sig_sigs) > 0) {
  box_data <- data.frame()
  for (sig_name in sig_sigs) {
    box_data <- rbind(box_data, data.frame(
      signature = sig_name,
      score = c(ssgsea_scores[sig_name, nor_samps], ssgsea_scores[sig_name, adj_samps]),
      group = c(rep("Normal", length(nor_samps)), rep("Adjacent", length(adj_samps))),
      stringsAsFactors = FALSE
    ))
  }
  box_data$group <- factor(box_data$group, levels = c("Normal", "Adjacent"))

  p_box <- ggplot(box_data, aes(x = signature, y = score, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    scale_fill_manual(values = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")) +
    coord_flip() +
    labs(title = "Immune Signatures: Adjacent vs Normal (P<0.1)",
         x = "", y = "ssGSEA Score", fill = "Group") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, "immune_boxplot.pdf"), p_box,
         width = 10, height = max(4, length(sig_sigs) * 0.5 + 2))
  cat("  Immune boxplot saved.\n")
}

# ---- 5c. Immune-clinical correlation ----
# Correlate immune scores with clinical variables
clin_file <- file.path(E1_DIR, "paired_clinical_with_subtypes_K2.csv")
if (file.exists(clin_file)) {
  clin <- read.csv(clin_file, check.names = FALSE)
  cat("  Computing immune-clinical correlations...\n")

  # Numeric clinical variables
  num_vars <- c("age", "NLR", "PLR", "CRP", "ALT", "AST", "total_bilirubin",
                "direct_bilirubin", "albumin", "WBC", "PNM_sum")
  num_vars <- intersect(num_vars, colnames(clin))

  if (length(num_vars) > 0) {
    # Map immune scores to patients (Adjacent samples)
    immune_clin_cor <- data.frame()
    for (sig_name in rownames(ssgsea_scores)) {
      for (cv in num_vars) {
        vals_immune <- c()
        vals_clin <- c()
        for (s in adj_samps) {
          pid <- gsub("^Adjacent", "", s)
          clin_row <- which(clin$patient_id == pid | rownames(clin) == pid)
          if (length(clin_row) > 0) {
            clin_val <- as.numeric(clin[clin_row[1], cv])
            if (!is.na(clin_val)) {
              vals_immune <- c(vals_immune, ssgsea_scores[sig_name, s])
              vals_clin <- c(vals_clin, clin_val)
            }
          }
        }
        if (length(vals_immune) >= 5) {
          ct <- cor.test(vals_immune, vals_clin, method = "spearman", exact = FALSE)
          immune_clin_cor <- rbind(immune_clin_cor, data.frame(
            signature = sig_name, clinical_var = cv,
            rho = ct$estimate, pvalue = ct$p.value,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    if (nrow(immune_clin_cor) > 0) {
      immune_clin_cor$padj <- p.adjust(immune_clin_cor$pvalue, method = "BH")
      immune_clin_cor <- immune_clin_cor[order(immune_clin_cor$pvalue), ]
      write.csv(immune_clin_cor, file.path(OUT_DIR, "immune_clinical_correlation.csv"),
                row.names = FALSE)

      sig_corr <- immune_clin_cor[immune_clin_cor$pvalue < 0.05, ]
      cat(sprintf("  Significant immune-clinical correlations: %d\n", nrow(sig_corr)))
      if (nrow(sig_corr) > 0) {
        cat("  Top immune-clinical correlations:\n")
        for (i in 1:min(15, nrow(sig_corr))) {
          cat(sprintf("    %s ~ %s: rho=%.2f, P=%.4f\n",
                      sig_corr$signature[i], sig_corr$clinical_var[i],
                      sig_corr$rho[i], sig_corr$pvalue[i]))
        }
      }

      # Correlation heatmap
      cor_wide <- immune_clin_cor %>%
        select(signature, clinical_var, rho) %>%
        pivot_wider(names_from = clinical_var, values_from = rho) %>%
        as.data.frame()
      rownames(cor_wide) <- cor_wide$signature
      cor_wide$signature <- NULL
      cor_wide <- as.matrix(cor_wide)

      # Filter to interesting rows
      row_max <- apply(abs(cor_wide), 1, max, na.rm = TRUE)
      cor_wide <- cor_wide[row_max > 0.3, , drop = FALSE]

      if (nrow(cor_wide) >= 2 && ncol(cor_wide) >= 2) {
        cor_wide[is.na(cor_wide)] <- 0
        pdf(file.path(FIG_DIR, "immune_clinical_correlation.pdf"), width = 10, height = 10)
        pheatmap(cor_wide,
                 color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                 breaks = seq(-1, 1, length.out = 101),
                 cluster_rows = TRUE, cluster_cols = TRUE,
                 fontsize_row = 8, fontsize_col = 10,
                 main = "Immune Signature ~ Clinical Variable Correlation (Spearman)")
        dev.off()
        cat("  Immune-clinical correlation heatmap saved.\n")
      }
    }
  }
}

# ============================================================================
# 6. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Enhancement 7 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Immune signatures scored: %d (ssGSEA)\n", nrow(ssgsea_scores)))
cat(sprintf("Adjacent vs Normal significant (P<0.05): %d\n",
            sum(immune_diff$pvalue < 0.05)))
if (nrow(immune_subtype) > 0) {
  cat(sprintf("DS1 vs DS2 significant (P<0.05): %d\n",
              sum(immune_subtype$pvalue < 0.05)))
}
if (exists("sig_corr")) {
  cat(sprintf("Immune-clinical correlations (P<0.05): %d\n", nrow(sig_corr)))
}
cat(sprintf("Output: %s\n", OUT_DIR))
cat("Enhancement 7 COMPLETE.\n")
