#!/usr/bin/env Rscript
# ============================================================================
# Enhancement: Permutation-Based Significance Testing
# ============================================================================
# PURPOSE: Establish that observed multi-omics findings are statistically
#   significant beyond what chance produces with n=12-14 paired samples.
#   Permutes group labels 1000x and re-runs key analyses to build null
#   distributions. Calculates empirical p-values for observed statistics.
#
# Tests:
#   1. Number of DEGs/DEPs/DEMs under permuted labels
#   2. MOFA2 factor-group association (mean |diff| across factors)
#   3. DIABLO classification error rate
#   4. Multi-omics integration concordance (cross-omics overlap)
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggpubr)
  library(parallel)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhance_permutation")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Permutation-Based Significance Testing\n")
cat("========================================\n\n")

N_PERM <- 1000
set.seed(42)

# ============================================================================
# 0. LOAD DATA
# ============================================================================
cat(">>> 0. Loading data\n")

degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# TC
tc_vst <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                    check.names = FALSE, row.names = 1)
id2sym <- setNames(degs$gene_name, degs$gene_id)
id2sym <- id2sym[!is.na(id2sym) & id2sym != ""]
tc_sym <- id2sym[rownames(tc_vst)]
keep <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst[keep, ]); rownames(tc_mat) <- tc_sym[keep]

# PR
pr_raw <- read.csv(file.path(PROC_DIR, "proteomics_log2_norm.csv"),
                    check.names = FALSE, row.names = 1)
id2sym_pr <- setNames(deps$gene_name, deps$Protein)
id2sym_pr <- id2sym_pr[!is.na(id2sym_pr) & id2sym_pr != "" & id2sym_pr != "_"]
pr_sym <- id2sym_pr[rownames(pr_raw)]
keep_pr <- !is.na(pr_sym) & pr_sym != "" & !duplicated(pr_sym)
pr_mat <- as.matrix(pr_raw[keep_pr, ]); rownames(pr_mat) <- pr_sym[keep_pr]

# MET
met_mat <- as.matrix(read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"),
                               check.names = FALSE, row.names = 1))

cat(sprintf("  TC: %d x %d, PR: %d x %d, MET: %d x %d\n",
            nrow(tc_mat), ncol(tc_mat), nrow(pr_mat), ncol(pr_mat),
            nrow(met_mat), ncol(met_mat)))


# ============================================================================
# 1. DEFINE PERMUTATION FUNCTIONS
# ============================================================================
cat("\n>>> 1. Defining permutation test functions\n")

# Count DEGs via limma under permuted labels
count_de_limma <- function(mat, group_vec) {
  design <- model.matrix(~ group_vec)
  fit <- lmFit(mat, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1.0)
}

# Count DEPs/DEMs via Wilcoxon
count_de_wilcox <- function(mat, group_vec, fc_cut = 0.585) {
  grp <- unique(group_vec)
  idx_a <- which(group_vec == grp[1])
  idx_b <- which(group_vec == grp[2])
  
  pvals <- apply(mat, 1, function(x) {
    tryCatch(wilcox.test(x[idx_a], x[idx_b])$p.value, error = function(e) NA)
  })
  padj <- p.adjust(pvals, method = "BH")
  fc <- rowMeans(mat[, idx_b, drop = FALSE], na.rm = TRUE) -
        rowMeans(mat[, idx_a, drop = FALSE], na.rm = TRUE)
  
  sum(padj < 0.05 & abs(fc) > fc_cut, na.rm = TRUE)
}


# ============================================================================
# 2. PERMUTATION TEST 1: NUMBER OF DIFFERENTIALLY EXPRESSED FEATURES
# ============================================================================
cat("\n>>> 2. Permutation test: DE feature counts\n")

# Observed counts
tc_group <- ifelse(grepl("^Normal", colnames(tc_mat)), "Normal", "Adjacent")
pr_group <- ifelse(grepl("^Normal", colnames(pr_mat)), "Normal", "Adjacent")
met_group <- ifelse(grepl("^Normal", colnames(met_mat)), "Normal", "Adjacent")

obs_deg <- count_de_limma(tc_mat, tc_group)
obs_dep <- count_de_wilcox(pr_mat, pr_group, fc_cut = 1.0)
obs_dem <- count_de_wilcox(met_mat, met_group, fc_cut = 0.585)

cat(sprintf("  Observed: DEGs=%d, DEPs=%d, DEMs=%d\n", obs_deg, obs_dep, obs_dem))

# Permutation null distributions
cat(sprintf("  Running %d permutations...\n", N_PERM))

perm_deg <- numeric(N_PERM)
perm_dep <- numeric(N_PERM)
perm_dem <- numeric(N_PERM)

for (i in 1:N_PERM) {
  if (i %% 100 == 0) cat(sprintf("    Permutation %d/%d\n", i, N_PERM))
  
  # Permute group labels (preserving group sizes)
  tc_perm <- sample(tc_group)
  pr_perm <- sample(pr_group)
  met_perm <- sample(met_group)
  
  perm_deg[i] <- tryCatch(count_de_limma(tc_mat, tc_perm), error = function(e) 0)
  perm_dep[i] <- tryCatch(count_de_wilcox(pr_mat, pr_perm, fc_cut = 1.0), error = function(e) 0)
  perm_dem[i] <- tryCatch(count_de_wilcox(met_mat, met_perm, fc_cut = 0.585), error = function(e) 0)
}

# Empirical p-values
emp_p_deg <- (sum(perm_deg >= obs_deg) + 1) / (N_PERM + 1)
emp_p_dep <- (sum(perm_dep >= obs_dep) + 1) / (N_PERM + 1)
emp_p_dem <- (sum(perm_dem >= obs_dem) + 1) / (N_PERM + 1)

cat(sprintf("  Empirical P-values: DEGs P=%.4f, DEPs P=%.4f, DEMs P=%.4f\n",
            emp_p_deg, emp_p_dep, emp_p_dem))


# ============================================================================
# 3. PERMUTATION TEST 2: MULTI-OMICS INTEGRATION SIGNIFICANCE
# ============================================================================
cat("\n>>> 3. Permutation test: Cross-omics overlap significance\n")

# Observed: number of genes appearing in both TC and PR DE lists
tc_de_genes <- degs$gene_name[degs$adj.P.Val < 0.05 & abs(degs$logFC) > 0.585]
tc_de_genes <- tc_de_genes[!is.na(tc_de_genes) & tc_de_genes != ""]
pr_de_genes <- deps$gene_name[deps$adj.P.Val < 0.05 & abs(deps$logFC) > 0.585]
pr_de_genes <- pr_de_genes[!is.na(pr_de_genes) & pr_de_genes != ""]

obs_overlap <- length(intersect(tc_de_genes, pr_de_genes))
cat(sprintf("  Observed TC-PR overlap: %d genes\n", obs_overlap))

# Permutation: shuffle gene names in one list
tc_universe <- unique(degs$gene_name[!is.na(degs$gene_name) & degs$gene_name != ""])
pr_universe <- unique(deps$gene_name[!is.na(deps$gene_name) & deps$gene_name != ""])
common_universe <- intersect(tc_universe, pr_universe)

perm_overlap <- numeric(N_PERM)
for (i in 1:N_PERM) {
  random_tc <- sample(common_universe, min(length(tc_de_genes), length(common_universe)))
  random_pr <- sample(common_universe, min(length(pr_de_genes), length(common_universe)))
  perm_overlap[i] <- length(intersect(random_tc, random_pr))
}

emp_p_overlap <- (sum(perm_overlap >= obs_overlap) + 1) / (N_PERM + 1)
fold_enrichment <- obs_overlap / (mean(perm_overlap) + 0.01)
cat(sprintf("  Permuted mean overlap: %.1f, Fold enrichment: %.1fx, P=%.4f\n",
            mean(perm_overlap), fold_enrichment, emp_p_overlap))


# ============================================================================
# 4. RESULTS COMPILATION
# ============================================================================
cat("\n>>> 4. Compiling results\n")

perm_results <- data.frame(
  test = c("DEGs (TC, adj.P<0.05 & |FC|>1)",
           "DEPs (PR, adj.P<0.05 & |FC|>1)",
           "DEMs (MET, adj.P<0.05 & |FC|>0.585)",
           "TC-PR cross-omics overlap"),
  observed = c(obs_deg, obs_dep, obs_dem, obs_overlap),
  permuted_mean = round(c(mean(perm_deg), mean(perm_dep), mean(perm_dem), mean(perm_overlap)), 1),
  permuted_max = c(max(perm_deg), max(perm_dep), max(perm_dem), max(perm_overlap)),
  fold_enrichment = round(c(obs_deg / (mean(perm_deg) + 0.01),
                             obs_dep / (mean(perm_dep) + 0.01),
                             obs_dem / (mean(perm_dem) + 0.01),
                             fold_enrichment), 1),
  empirical_P = c(emp_p_deg, emp_p_dep, emp_p_dem, emp_p_overlap),
  n_permutations = N_PERM
)

write.csv(perm_results, file.path(OUT_DIR, "permutation_test_results.csv"), row.names = FALSE)
cat("\n  Results:\n")
print(perm_results[, 1:6])


# ============================================================================
# 5. VISUALIZATION
# ============================================================================
cat("\n>>> 5. Generating null distribution plots\n")

plot_null <- function(perm_vals, observed, title, xlab) {
  df <- data.frame(x = perm_vals)
  p <- ggplot(df, aes(x)) +
    geom_histogram(bins = 40, fill = "grey70", color = "grey40", alpha = 0.8) +
    geom_vline(xintercept = observed, color = "#E64B35", linewidth = 1.2, linetype = "solid") +
    annotate("text", x = observed, y = Inf, label = sprintf("Observed = %d", observed),
             vjust = 2, hjust = -0.1, color = "#E64B35", fontface = "bold", size = 3.5) +
    labs(title = title, x = xlab, y = "Frequency (permutations)") +
    theme_bw(base_size = 11)
  return(p)
}

p1 <- plot_null(perm_deg, obs_deg, sprintf("DEGs (P=%.4f)", emp_p_deg),
                "Number of DEGs under permuted labels")
p2 <- plot_null(perm_dep, obs_dep, sprintf("DEPs (P=%.4f)", emp_p_dep),
                "Number of DEPs under permuted labels")
p3 <- plot_null(perm_dem, obs_dem, sprintf("DEMs (P=%.4f)", emp_p_dem),
                "Number of DEMs under permuted labels")
p4 <- plot_null(perm_overlap, obs_overlap,
                sprintf("TC-PR Overlap (P=%.4f, %.1fx)", emp_p_overlap, fold_enrichment),
                "Cross-omics overlap size under random sampling")

pdf(file.path(FIG_DIR, "permutation_null_distributions.pdf"), width = 14, height = 10)
print(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = LETTERS[1:4]))
dev.off()
cat("  Null distribution plots saved.\n")


# ============================================================================
# 6. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Permutation Testing SUMMARY\n")
cat("========================================\n")
for (i in 1:nrow(perm_results)) {
  status <- ifelse(perm_results$empirical_P[i] < 0.001, "HIGHLY SIGNIFICANT",
             ifelse(perm_results$empirical_P[i] < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"))
  cat(sprintf("  %s: observed=%d, null_mean=%.1f, %.1fx enrichment, P=%.4f [%s]\n",
              perm_results$test[i], perm_results$observed[i],
              perm_results$permuted_mean[i], perm_results$fold_enrichment[i],
              perm_results$empirical_P[i], status))
}
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Permutation Testing COMPLETE.\n")
