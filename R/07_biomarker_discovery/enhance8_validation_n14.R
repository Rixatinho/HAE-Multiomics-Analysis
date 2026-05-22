#!/usr/bin/env Rscript
# ============================================================================
# Enhancement 8: Cell Type Deconvolution + n=14 Expansion
# HAE Multi-omics Integration Study
# ============================================================================
# Part 1: Liver cell type deconvolution (ssGSEA-based) + decoupling analysis
# Part 2: n=14 expansion (PR+MB MOFA2, subtyping, clinical correlation)
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(MOFA2)
  library(ConsensusClusterPlus)
  library(RColorBrewer)
})

PROJECT <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
DATA_DIR <- file.path(PROJECT, "analysis/data/processed")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhancement8_validation")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

col_group <- c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")
extract_id <- function(x) as.integer(gsub("Normal|Adjacent", "", x))

# Load GMT file for Hallmark pathways (pre-downloaded)
read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  gene_sets <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    gene_sets[[parts[1]]] <- parts[-(1:2)]
  }
  gene_sets
}

GMT_FILE <- "/tmp/h.all.v2024.1.Hs.symbols.gmt"
if (!file.exists(GMT_FILE)) {
  stop("Hallmark GMT file not found. Please download first.")
}
hallmark_sets <- read_gmt(GMT_FILE)
cat(sprintf("Loaded %d Hallmark gene sets from GMT\n", length(hallmark_sets)))

# Rank-based ssGSEA implementation
compute_ssgsea <- function(expr_mat, gene_sets) {
  ranks <- apply(expr_mat, 2, function(x) rank(-x))
  n_genes <- nrow(ranks)
  scores <- matrix(NA, nrow = length(gene_sets), ncol = ncol(ranks))
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(ranks)
  for (gs_name in names(gene_sets)) {
    idx <- which(rownames(ranks) %in% gene_sets[[gs_name]])
    if (length(idx) < 3) next
    for (j in seq_len(ncol(ranks))) {
      scores[gs_name, j] <- -(mean(ranks[idx, j]) - mean(ranks[-idx, j]))
    }
  }
  scores
}

cat("================================================================\n")
cat("Enhancement 8: Deconvolution + n=14 Expansion\n")
cat("================================================================\n\n")

# ============================================================================
# PART 1: Cell Type Deconvolution + Decoupling Mechanism Analysis
# ============================================================================
cat(">>> PART 1: Liver Cell Type Deconvolution\n")

tc_vst <- as.matrix(read.csv(file.path(DATA_DIR, "transcriptomics_vst_paired.csv"),
                              row.names = 1, check.names = FALSE))

# Map gene IDs to symbols (use project-internal DEG results first, fallback to external)
degs_file <- file.path(PROJECT, "analysis/results/phase1_diff/DEGs_Adjacent_vs_Normal.csv")
DATA_ROOT <- "/Volumes/阿热-实验数据备份/工作/1.肝包虫/1.泡型包虫/1.2024-11-泡球蚴"
tc_file <- file.path(DATA_ROOT,
  "1.转录组学测序/03.Result_X101SC24101695-Z02-F001_homo_sapiens/Result_X101SC24101695-Z02-F001_homo_sapiens/4.Quant/1.Count/gene_count.xls")

if (file.exists(degs_file)) {
  # Preferred: use project-internal DEG results for ID->symbol mapping
  deg_info <- read.csv(degs_file, check.names = FALSE)
  id_to_symbol <- setNames(deg_info$gene_name, deg_info$gene_id)
  cat("  Gene ID mapping source: phase1_diff DEGs (project-internal)\n")
} else if (file.exists(tc_file)) {
  tc_raw_info <- read.delim(tc_file, check.names = FALSE)
  id_to_symbol <- setNames(tc_raw_info$gene_name, tc_raw_info$gene_id)
  cat("  Gene ID mapping source: external drive (raw count file)\n")
} else {
  id_to_symbol <- setNames(rownames(tc_vst), rownames(tc_vst))
  cat("  WARNING: No gene mapping file found, using row names as-is\n")
}

# Map to gene symbols and deduplicate
tc_symbol <- tc_vst
mapped <- id_to_symbol[rownames(tc_symbol)]
valid <- !is.na(mapped) & mapped != "" & mapped != "---"
tc_symbol <- tc_symbol[valid, ]
rownames(tc_symbol) <- mapped[valid]
rv <- apply(tc_symbol, 1, var)
tc_symbol <- tc_symbol[order(rownames(tc_symbol), -rv), ]
tc_symbol <- tc_symbol[!duplicated(rownames(tc_symbol)), ]
cat(sprintf("  Expression matrix: %d genes x %d samples\n",
            nrow(tc_symbol), ncol(tc_symbol)))

# --- Liver cell type signatures ---
liver_sigs <- list(
  Hepatocytes = c("ALB","APOB","APOA1","APOC3","SERPINA1","TTR","TF","HP",
                   "FGA","FGB","FGG","CYP3A4","CYP2E1","CYP1A2","CYP2C9",
                   "UGT1A1","ADH1B","ALDOB","PCK1","G6PC","TAT","HPD","HAO1"),
  Cholangiocytes = c("KRT19","KRT7","EPCAM","SOX9","CFTR","TFF1","MUC1",
                      "FXYD2","DEFB1","ANXA4"),
  LSECs = c("CLEC4G","CLEC4M","STAB2","FCGR2B","LYVE1","CLEC1B","PECAM1",
             "ERG","FLT1","DNASE1L3"),
  Kupffer_cells = c("CD68","CD163","MARCO","TIMD4","VSIG4","C1QA","C1QB",
                     "C1QC","FOLR2","SLC40A1","GPNMB"),
  NK_NKT = c("NKG7","GNLY","GZMB","GZMA","GZMK","KLRD1","KLRB1","NCAM1",
              "CD7","PRF1"),
  T_cells = c("CD3D","CD3E","CD3G","CD2","CD8A","CD8B","CD4","TRAC","IL7R",
               "LEF1","TCF7"),
  B_Plasma = c("CD79A","CD79B","MS4A1","CD19","JCHAIN","MZB1","IGKC",
                "IGHG1","XBP1","PRDM1"),
  Stellate_Fibro = c("ACTA2","COL1A1","COL1A2","COL3A1","PDGFRB","DES",
                      "RGS5","TAGLN","DCN","BGN","LUM","COL6A1","COL6A2",
                      "FBLN1","VCAN","FAP"),
  Macrophages_M1 = c("TNF","IL1B","IL6","CCL2","CCL3","CCL4","CXCL10",
                      "IDO1","NOS2","CD86"),
  Macrophages_M2 = c("CD163","MRC1","MSR1","CD209","CLEC10A","CCL18",
                      "CCL22","TGFB1","IL10","ARG1"),
  Neutrophils = c("S100A8","S100A9","S100A12","FCGR3B","CXCR2","CSF3R",
                   "CEACAM8","MMP9"),
  Endothelial = c("PECAM1","CDH5","VWF","ENG","KDR","FLT1","EMCN","PLVAP")
)

# Run deconvolution
deconv <- compute_ssgsea(tc_symbol, liver_sigs)
deconv <- deconv[!apply(is.na(deconv), 1, all), ]
cat(sprintf("  Cell types quantified: %d\n", nrow(deconv)))
write.csv(deconv, file.path(OUT_DIR, "cell_type_ssgsea_scores.csv"))

# Paired differential (n=12)
normal_cols <- grep("^Normal", colnames(deconv), value = TRUE)
adj_cols <- grep("^Adjacent", colnames(deconv), value = TRUE)
paired_ids <- intersect(extract_id(normal_cols), extract_id(adj_cols))

deconv_diff <- data.frame(cell_type = rownames(deconv), mean_Adj = NA,
                          mean_Nor = NA, diff = NA, pvalue = NA)
for (i in seq_len(nrow(deconv))) {
  ct <- rownames(deconv)[i]
  n_v <- deconv[ct, paste0("Normal", paired_ids)]
  a_v <- deconv[ct, paste0("Adjacent", paired_ids)]
  deconv_diff$mean_Adj[i] <- mean(a_v)
  deconv_diff$mean_Nor[i] <- mean(n_v)
  deconv_diff$diff[i] <- mean(a_v) - mean(n_v)
  deconv_diff$pvalue[i] <- wilcox.test(a_v, n_v, paired = TRUE)$p.value
}
deconv_diff$padj <- p.adjust(deconv_diff$pvalue, method = "BH")
deconv_diff <- deconv_diff[order(deconv_diff$pvalue), ]
write.csv(deconv_diff, file.path(OUT_DIR, "cell_type_differential.csv"),
          row.names = FALSE)
cat("  Cell type differential:\n")
print(deconv_diff[, c("cell_type", "diff", "pvalue", "padj")])

# Heatmap
deconv_sc <- t(scale(t(deconv)))
ann_df <- data.frame(Group = ifelse(grepl("^Normal", colnames(deconv)),
                                     "Normal", "Adjacent"),
                     row.names = colnames(deconv))

pdf(file.path(FIG_DIR, "cell_type_deconvolution_heatmap.pdf"), width = 10, height = 7)
print(pheatmap(deconv_sc,
  annotation_col = ann_df,
  annotation_colors = list(Group = col_group),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 9,
  main = "Liver Cell Type Deconvolution (ssGSEA, n=12 paired)"
))
dev.off()
cat("  Cell type heatmap saved.\n")

# --- Boxplot of key cell types ---
deconv_long <- data.frame()
for (ct in rownames(deconv)) {
  for (s in colnames(deconv)) {
    deconv_long <- rbind(deconv_long, data.frame(
      cell_type = ct, sample = s, score = deconv[ct, s],
      group = ifelse(grepl("^Normal", s), "Normal", "Adjacent"),
      patient = extract_id(s)
    ))
  }
}

top_cts <- deconv_diff$cell_type[1:min(6, nrow(deconv_diff))]
p_box <- ggplot(deconv_long[deconv_long$cell_type %in% top_cts, ],
                aes(group, score, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_line(aes(group = patient), alpha = 0.3, linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = col_group) +
  labs(title = "Cell Type Abundance: Adjacent vs Normal (paired)",
       y = "ssGSEA Score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "cell_type_boxplots.pdf"), p_box, width = 9, height = 6)
cat("  Cell type boxplots saved.\n")

# --- Decoupling mechanism: correlate cell types with discordant pathways ---
cat("\n  Decoupling analysis...\n")

# Compute per-sample Hallmark ssGSEA scores
hallmark_scores <- compute_ssgsea(tc_symbol, hallmark_sets)

# Load our cross-omics results to identify discordant pathways
our_hall <- read.csv(file.path(PROJECT,
  "analysis/results/phase2_enrichment/cross_omics_convergent_Hallmark.csv"))
discordant <- our_hall[our_hall$concordant == FALSE &
                        our_hall$TC_padj < 0.05 &
                        our_hall$PR_padj < 0.05, ]
cat(sprintf("  Discordant pathways: %d\n", nrow(discordant)))

# Correlate cell type scores with discordant pathway scores
cor_results <- data.frame()
disc_pws <- discordant$pathway
disc_pws <- disc_pws[disc_pws %in% rownames(hallmark_scores)]

for (ct in rownames(deconv)) {
  for (pw in disc_pws) {
    cs <- intersect(colnames(deconv), colnames(hallmark_scores))
    ct_v <- as.numeric(deconv[ct, cs])
    pw_v <- as.numeric(hallmark_scores[pw, cs])
    if (sd(ct_v) > 0 & sd(pw_v) > 0) {
      ct_res <- cor.test(ct_v, pw_v, method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        cell_type = ct, pathway = gsub("^HALLMARK_", "", pw),
        rho = ct_res$estimate, pvalue = ct_res$p.value
      ))
    }
  }
}
cor_results$padj <- p.adjust(cor_results$pvalue, method = "BH")
cor_results <- cor_results[order(cor_results$pvalue), ]
write.csv(cor_results, file.path(OUT_DIR, "celltype_pathway_correlation.csv"),
          row.names = FALSE)

n_sig <- sum(cor_results$pvalue < 0.05)
cat(sprintf("  Cell type-pathway correlations: %d tested, %d P<0.05\n",
            nrow(cor_results), n_sig))
if (n_sig > 0) {
  cat("  Top correlations (P<0.05):\n")
  print(head(cor_results[cor_results$pvalue < 0.05, ], 10))
}

# Heatmap of correlations
if (nrow(cor_results) > 5) {
  rho_wide <- reshape(cor_results[, c("cell_type", "pathway", "rho")],
                      idvar = "cell_type", timevar = "pathway",
                      direction = "wide")
  rownames(rho_wide) <- rho_wide$cell_type
  rho_wide$cell_type <- NULL
  colnames(rho_wide) <- gsub("^rho\\.", "", colnames(rho_wide))
  rho_mat <- as.matrix(rho_wide)
  rho_mat[is.na(rho_mat)] <- 0

  pval_wide <- reshape(cor_results[, c("cell_type", "pathway", "pvalue")],
                       idvar = "cell_type", timevar = "pathway",
                       direction = "wide")
  rownames(pval_wide) <- pval_wide$cell_type
  pval_wide$cell_type <- NULL
  pval_mat <- as.matrix(pval_wide)

  sig_lab <- matrix("", nrow = nrow(rho_mat), ncol = ncol(rho_mat))
  sig_lab[pval_mat < 0.001] <- "***"
  sig_lab[pval_mat >= 0.001 & pval_mat < 0.01] <- "**"
  sig_lab[pval_mat >= 0.01 & pval_mat < 0.05] <- "*"

  pdf(file.path(FIG_DIR, "celltype_discordant_pathway_heatmap.pdf"),
      width = 12, height = 6)
  print(pheatmap(rho_mat,
    display_numbers = sig_lab,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = TRUE, cluster_cols = TRUE,
    fontsize = 9, fontsize_number = 11,
    main = "Cell Type vs Discordant Pathway Correlation (Spearman)"
  ))
  dev.off()
  cat("  Correlation heatmap saved.\n")
}

# --- Hepatocyte fraction vs discordant pathway scatter ---
hep_scores <- deconv["Hepatocytes", ]
ifng_pw <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
if (ifng_pw %in% rownames(hallmark_scores)) {
  cs <- intersect(names(hep_scores), colnames(hallmark_scores))
  scat_df <- data.frame(
    Hepatocyte_score = as.numeric(hep_scores[cs]),
    IFNg_score = as.numeric(hallmark_scores[ifng_pw, cs]),
    group = ifelse(grepl("^Normal", cs), "Normal", "Adjacent"),
    sample = cs
  )
  ct <- cor.test(scat_df$Hepatocyte_score, scat_df$IFNg_score, method = "spearman")

  p_scat <- ggplot(scat_df, aes(Hepatocyte_score, IFNg_score)) +
    geom_point(aes(color = group), size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8) +
    scale_color_manual(values = col_group) +
    labs(
      title = sprintf("Hepatocyte Fraction vs IFN-gamma (rho=%.2f, P=%.3f)",
                      ct$estimate, ct$p.value),
      x = "Hepatocyte ssGSEA Score", y = "IFN-gamma Pathway Score"
    ) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")

  ggsave(file.path(FIG_DIR, "hepatocyte_vs_IFNg_scatter.pdf"),
         p_scat, width = 6, height = 5)
  cat("  Hepatocyte vs IFN-gamma scatter saved.\n")
}

cat("  PART 1 COMPLETE.\n\n")


# ============================================================================
# PART 2: n=14 Expansion
# ============================================================================
cat(">>> PART 2: n=14 Expansion Analyses\n")

pr_norm <- as.matrix(read.csv(file.path(DATA_DIR, "proteomics_log2_norm.csv"),
                               row.names = 1, check.names = FALSE))
met_log2 <- as.matrix(read.csv(file.path(DATA_DIR, "metabolomics_log2_merged.csv"),
                                row.names = 1, check.names = FALSE))
cat(sprintf("  PR: %d x %d, MB: %d x %d\n",
            nrow(pr_norm), ncol(pr_norm), nrow(met_log2), ncol(met_log2)))

# --- 2A: MOFA2 (PR+MB, n=14) ---
cat("\n  2A: MOFA2 (PR+MB, n=14)...\n")

tryCatch({
  pr_adj <- pr_norm[, paste0("Adjacent", 1:14)]
  mb_adj <- met_log2[, paste0("Adjacent", 1:14)]

  # Top variable features
  pr_top <- names(sort(apply(pr_adj, 1, var), decreasing = TRUE))[1:2000]
  mb_top <- names(sort(apply(mb_adj, 1, var), decreasing = TRUE))[1:1000]

  pr_m <- pr_adj[pr_top, ]
  mb_m <- mb_adj[mb_top, ]
  colnames(pr_m) <- colnames(mb_m) <- paste0("P", 1:14)

  mofa_obj <- create_mofa(list(Proteomics = pr_m, Metabolomics = mb_m))
  model_opts <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- 5
  train_opts <- get_default_training_options(mofa_obj)
  train_opts$convergence_mode <- "slow"
  train_opts$seed <- 42
  train_opts$verbose <- FALSE

  mofa_obj <- prepare_mofa(mofa_obj,
    data_options = get_default_data_options(mofa_obj),
    model_options = model_opts,
    training_options = train_opts
  )

  mofa_out <- file.path(OUT_DIR, "MOFA2_n14_PRMB.hdf5")
  mofa_trained <- run_mofa(mofa_obj, use_basilisk = FALSE, outfile = mofa_out)

  factors <- get_factors(mofa_trained)$group1
  var_exp <- get_variance_explained(mofa_trained)
  cat(sprintf("    Factors: %d x %d\n", nrow(factors), ncol(factors)))
  cat("    Variance explained:\n")
  print(var_exp$r2_per_factor$group1)

  write.csv(factors, file.path(OUT_DIR, "MOFA2_n14_factors.csv"))

  # Clinical correlation (all 14 patients)
  clinical <- read.csv(file.path(PROJECT,
    "analysis/results/phase6_subtyping/clinical_parsed.csv"))

  num_vars <- c("age","BMI","disease_duration","lesion_size","ALT","AST",
                "ALP","GGT","total_bilirubin","albumin","WBC",
                "neutrophil_pct","lymphocyte_pct","CRP","IL6")

  fc_pval <- matrix(NA, ncol(factors), length(num_vars),
                    dimnames = list(colnames(factors), num_vars))
  fc_rho <- fc_pval

  for (cv in num_vars) {
    if (!(cv %in% colnames(clinical))) next
    for (fi in seq_len(ncol(factors))) {
      pid <- as.integer(gsub("^P", "", rownames(factors)))
      cval <- clinical[[cv]][match(pid, clinical$patient_id)]
      ok <- !is.na(cval) & !is.na(factors[, fi])
      if (sum(ok) >= 5) {
        ct <- cor.test(factors[ok, fi], cval[ok], method = "spearman")
        fc_pval[fi, cv] <- ct$p.value
        fc_rho[fi, cv] <- ct$estimate
      }
    }
  }

  write.csv(fc_pval, file.path(OUT_DIR, "MOFA2_n14_factor_clinical_pval.csv"))
  write.csv(fc_rho, file.path(OUT_DIR, "MOFA2_n14_factor_clinical_rho.csv"))

  cat("    Significant associations (P<0.05):\n")
  sig <- which(fc_pval < 0.05, arr.ind = TRUE)
  if (nrow(sig) > 0) {
    for (k in seq_len(nrow(sig))) {
      cat(sprintf("      %s vs %s: rho=%.3f, P=%.4f\n",
                  rownames(fc_pval)[sig[k,1]], colnames(fc_pval)[sig[k,2]],
                  fc_rho[sig[k,1], sig[k,2]], fc_pval[sig[k,1], sig[k,2]]))
    }
  } else { cat("      None\n") }

  # Factor-clinical heatmap
  rho_plot <- fc_rho
  rho_plot[is.na(rho_plot)] <- 0
  sig_lab2 <- matrix("", nrow(rho_plot), ncol(rho_plot))
  sig_lab2[fc_pval < 0.001] <- "***"
  sig_lab2[fc_pval >= 0.001 & fc_pval < 0.01] <- "**"
  sig_lab2[fc_pval >= 0.01 & fc_pval < 0.05] <- "*"

  pdf(file.path(FIG_DIR, "MOFA2_n14_factor_clinical_heatmap.pdf"), width = 10, height = 5)
  print(pheatmap(rho_plot,
    display_numbers = sig_lab2,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE, cluster_cols = TRUE,
    fontsize = 9, fontsize_number = 11,
    main = "MOFA2 Factors vs Clinical (Spearman, n=14)"
  ))
  dev.off()
  cat("    Factor-clinical heatmap saved.\n")

  # Variance explained barplot
  r2_df <- data.frame()
  r2_mat <- var_exp$r2_per_factor$group1
  for (vi in seq_len(ncol(r2_mat))) {
    for (fi in seq_len(nrow(r2_mat))) {
      r2_df <- rbind(r2_df, data.frame(
        Factor = rownames(r2_mat)[fi], View = colnames(r2_mat)[vi],
        R2 = r2_mat[fi, vi]
      ))
    }
  }
  p_var <- ggplot(r2_df, aes(Factor, R2, fill = View)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "MOFA2 Variance Explained (PR+MB, n=14)",
         y = "R2 (%)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, "MOFA2_n14_variance.pdf"), p_var, width = 7, height = 5)
  cat("    Variance plot saved.\n")

}, error = function(e) {
  cat(sprintf("    MOFA2 error: %s\n", e$message))
})


# --- 2B: Consensus subtyping (PR+MB, n=14) ---
cat("\n  2B: Consensus Subtyping (PR+MB, n=14)...\n")

tryCatch({
  pr_adj14 <- pr_norm[, paste0("Adjacent", 1:14)]
  mb_adj14 <- met_log2[, paste0("Adjacent", 1:14)]

  pr_sel <- names(sort(apply(pr_adj14, 1, var), decreasing = TRUE))[1:500]
  mb_sel <- names(sort(apply(mb_adj14, 1, var), decreasing = TRUE))[1:300]

  pr_sc <- t(scale(t(pr_adj14[pr_sel, ])))
  mb_sc <- t(scale(t(mb_adj14[mb_sel, ])))
  concat14 <- rbind(pr_sc, mb_sc)
  colnames(concat14) <- paste0("P", 1:14)
  concat14 <- concat14[complete.cases(concat14), ]
  cat(sprintf("    Matrix: %d features x %d samples\n",
              nrow(concat14), ncol(concat14)))

  graphics.off()
  cc_dir <- file.path(OUT_DIR, "consensus_n14")
  dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)

  cc <- ConsensusClusterPlus(as.matrix(concat14), maxK = 4, reps = 1000,
    pItem = 0.8, pFeature = 0.8, clusterAlg = "hc", distance = "pearson",
    seed = 42, plot = "pdf", title = cc_dir)

  k2 <- cc[[2]]$consensusClass
  cat(sprintf("    K=2: %s\n", paste(table(k2), collapse = " / ")))

  sub_n14 <- data.frame(patient_id = 1:14, subtype_K2 = paste0("CS", k2))
  write.csv(sub_n14, file.path(OUT_DIR, "subtype_n14_K2.csv"), row.names = FALSE)

  # Compare with n=12
  sub_n12 <- read.csv(file.path(PROJECT,
    "analysis/results/phase6_subtyping/subtype_K2.csv"))
  merged <- merge(sub_n14, sub_n12, by = "patient_id", suffixes = c("_n14", "_n12"),
                  all.x = TRUE)
  write.csv(merged, file.path(OUT_DIR, "subtype_n14_vs_n12.csv"), row.names = FALSE)
  cat("    n=14 vs n=12 comparison:\n")
  print(merged)

  # Clinical association (n=14)
  clinical <- read.csv(file.path(PROJECT,
    "analysis/results/phase6_subtyping/clinical_parsed.csv"))
  clin_sub <- merge(clinical, sub_n14, by = "patient_id")

  num_vars <- c("age","BMI","disease_duration","lesion_size","ALT","AST",
                "ALP","GGT","total_bilirubin","albumin","WBC",
                "neutrophil_pct","lymphocyte_pct","CRP","IL6")

  clin_assoc <- data.frame()
  for (cv in num_vars) {
    if (!(cv %in% colnames(clin_sub))) next
    g1 <- clin_sub[[cv]][clin_sub$subtype_K2 == "CS1"]
    g2 <- clin_sub[[cv]][clin_sub$subtype_K2 == "CS2"]
    g1 <- g1[!is.na(g1)]; g2 <- g2[!is.na(g2)]
    if (length(g1) >= 2 & length(g2) >= 2) {
      wt <- wilcox.test(g1, g2)
      clin_assoc <- rbind(clin_assoc, data.frame(
        variable = cv, mean_CS1 = mean(g1), mean_CS2 = mean(g2),
        pvalue = wt$p.value))
    }
  }
  if (nrow(clin_assoc) > 0) {
    clin_assoc$padj <- p.adjust(clin_assoc$pvalue, method = "BH")
    clin_assoc <- clin_assoc[order(clin_assoc$pvalue), ]
    write.csv(clin_assoc, file.path(OUT_DIR, "clinical_association_n14_K2.csv"),
              row.names = FALSE)
    cat("    Clinical associations (P<0.1):\n")
    print(clin_assoc[clin_assoc$pvalue < 0.1, ])
  }

  # Consensus heatmap
  graphics.off()
  cm <- cc[[2]]$consensusMatrix
  colnames(cm) <- rownames(cm) <- paste0("P", 1:14)

  pdf(file.path(FIG_DIR, "consensus_heatmap_n14_K2.pdf"), width = 8, height = 7)
  print(pheatmap(cm,
    annotation_col = data.frame(Subtype = paste0("CS", k2),
                                row.names = paste0("P", 1:14)),
    annotation_colors = list(Subtype = c(CS1="#E64B35", CS2="#4DBBD5")),
    color = colorRampPalette(c("white", "#2166AC"))(100),
    main = "Consensus Matrix K=2 (PR+MB, n=14)", fontsize = 10
  ))
  dev.off()
  cat("    Consensus heatmap saved.\n")

}, error = function(e) {
  cat(sprintf("    Consensus error: %s\n", e$message))
})


# --- 2C: n=12 vs n=14 DEP comparison ---
cat("\n  2C: n=12 vs n=14 DEP comparison...\n")

tc_ids <- c(1:9, 11, 13, 14)

pr_12 <- pr_norm[, c(paste0("Normal", tc_ids), paste0("Adjacent", tc_ids))]
si_12 <- data.frame(sample = colnames(pr_12),
  group = factor(ifelse(grepl("^Normal", colnames(pr_12)), "Normal", "Adjacent"),
                 levels = c("Normal", "Adjacent")),
  patient = factor(extract_id(colnames(pr_12))))

design_12 <- model.matrix(~ 0 + group, data = si_12)
colnames(design_12) <- c("Normal", "Adjacent")
corfit_12 <- duplicateCorrelation(pr_12, design_12, block = si_12$patient)
fit_12 <- lmFit(pr_12, design_12, block = si_12$patient,
                correlation = corfit_12$consensus.correlation)
fit2_12 <- eBayes(contrasts.fit(fit_12, makeContrasts(Adjacent - Normal,
                                                       levels = design_12)))
res_12 <- topTable(fit2_12, number = Inf)

res_14 <- read.csv(file.path(PROJECT,
  "analysis/results/phase1_diff/DEPs_Adjacent_vs_Normal.csv"))

comp <- data.frame(
  Metric = c("DEP_padj005", "DEP_relaxed_P005_FC0585", "DEP_padj005_FC1"),
  n12 = c(sum(res_12$adj.P.Val < 0.05),
          sum(res_12$P.Value < 0.05 & abs(res_12$logFC) > 0.585),
          sum(res_12$adj.P.Val < 0.05 & abs(res_12$logFC) > 1)),
  n14 = c(sum(res_14$adj.P.Val < 0.05, na.rm = TRUE),
          sum(res_14$P.Value < 0.05 & abs(res_14$logFC) > 0.585, na.rm = TRUE),
          sum(res_14$adj.P.Val < 0.05 & abs(res_14$logFC) > 1, na.rm = TRUE))
)
comp$gain <- comp$n14 - comp$n12
comp$gain_pct <- round(100 * comp$gain / pmax(comp$n12, 1), 1)
write.csv(comp, file.path(OUT_DIR, "n12_vs_n14_DEP_comparison.csv"),
          row.names = FALSE)
cat("  DEP comparison:\n")
print(comp)


# --- 2D: Convergent gene clinical correlation (n=14) ---
cat("\n  2D: Convergent Gene-Clinical Correlation (n=14)...\n")

pr_adj_all <- pr_norm[, paste0("Adjacent", 1:14)]
pr_anno <- read.csv(file.path(DATA_DIR, "proteomics_annotation.csv"),
                     row.names = 1, check.names = FALSE)
gene_map <- setNames(pr_anno$Gene, pr_anno$Protein)

convergent <- read.csv(file.path(PROJECT,
  "analysis/results/optimization_figures/cross_omics_convergent_genes.csv"))
top_genes <- convergent$gene[convergent$evidence_count >= 3]

clinical <- read.csv(file.path(PROJECT,
  "analysis/results/phase6_subtyping/clinical_parsed.csv"))

# Find convergent genes in PR
pr_genes <- gene_map[rownames(pr_adj_all)]
conv_idx <- which(pr_genes %in% top_genes)

if (length(conv_idx) > 0) {
  conv_expr <- pr_adj_all[conv_idx, ]
  rownames(conv_expr) <- pr_genes[conv_idx]

  gc_cor <- data.frame()
  for (gi in seq_len(nrow(conv_expr))) {
    gn <- rownames(conv_expr)[gi]
    for (cv in num_vars) {
      if (!(cv %in% colnames(clinical))) next
      cvals <- clinical[[cv]][match(1:14, clinical$patient_id)]
      gvals <- as.numeric(conv_expr[gi, ])
      ok <- !is.na(cvals) & !is.na(gvals)
      if (sum(ok) >= 5) {
        ct <- cor.test(gvals[ok], cvals[ok], method = "spearman")
        gc_cor <- rbind(gc_cor, data.frame(gene = gn, clinical_var = cv,
                                            rho = ct$estimate, pvalue = ct$p.value))
      }
    }
  }

  if (nrow(gc_cor) > 0) {
    gc_cor$padj <- p.adjust(gc_cor$pvalue, method = "BH")
    gc_cor <- gc_cor[order(gc_cor$pvalue), ]
    write.csv(gc_cor, file.path(OUT_DIR, "convergent_gene_clinical_n14.csv"),
              row.names = FALSE)
    cat(sprintf("    Tested: %d, P<0.05: %d, padj<0.05: %d\n",
                nrow(gc_cor), sum(gc_cor$pvalue < 0.05), sum(gc_cor$padj < 0.05)))
    cat("    Top 15 (P<0.05):\n")
    print(head(gc_cor[gc_cor$pvalue < 0.05, ], 15))

    # Heatmap
    sig_g <- unique(gc_cor$gene[gc_cor$pvalue < 0.05])
    if (length(sig_g) >= 3) {
      rho_w <- reshape(gc_cor[gc_cor$gene %in% sig_g, c("gene","clinical_var","rho")],
                       idvar = "gene", timevar = "clinical_var", direction = "wide")
      rownames(rho_w) <- rho_w$gene; rho_w$gene <- NULL
      colnames(rho_w) <- gsub("^rho\\.", "", colnames(rho_w))
      rho_m <- as.matrix(rho_w); rho_m[is.na(rho_m)] <- 0

      plot_n <- min(20, nrow(rho_m))
      pdf(file.path(FIG_DIR, "convergent_gene_clinical_heatmap_n14.pdf"),
          width = 10, height = max(5, plot_n * 0.4))
      print(pheatmap(rho_m[1:plot_n, , drop = FALSE],
        color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        breaks = seq(-1, 1, length.out = 101),
        cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 9,
        main = "Convergent Genes vs Clinical (Spearman rho, n=14)"
      ))
      dev.off()
      cat("    Gene-clinical heatmap saved.\n")
    }
  }
}

cat("  PART 2 COMPLETE.\n\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("================================================================\n")
cat("Enhancement 8 COMPLETE\n")
cat("================================================================\n")
cat("Part 1: Cell Type Deconvolution\n")
cat("  - 12 liver cell types quantified (ssGSEA)\n")
cat("  - Cell type vs discordant pathway correlation\n")
cat("  - Hepatocyte fraction vs IFN-gamma pathway\n")
cat("Part 2: n=14 Expansion\n")
cat("  - MOFA2 (PR+MB, n=14) + clinical correlation\n")
cat("  - Consensus subtyping K=2 (n=14)\n")
cat("  - n=12 vs n=14 DEP statistical comparison\n")
cat("  - Convergent gene-clinical correlation (n=14)\n")
cat(sprintf("Output: %s\n", OUT_DIR))
cat("================================================================\n")
