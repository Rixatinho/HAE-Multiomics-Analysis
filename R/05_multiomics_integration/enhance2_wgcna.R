#!/usr/bin/env Rscript
# ============================================================================
# Enhancement 2: WGCNA Cross-omics Co-expression Network
# ============================================================================
# Key improvement: Identify co-expressed gene/protein modules, relate to
# clinical traits and molecular subtypes, find hub genes for each module.
# ============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

allowWGCNAThreads(4)

PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
E1_DIR   <- file.path(PROJECT, "analysis/results/enhancement1_paired_subtyping")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhancement2_wgcna")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Enhancement 2: WGCNA Network Analysis\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA (all samples for WGCNA — more power)
# ============================================================================
cat(">>> 0. Loading data\n")

degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# Transcriptomics
tc_vst_raw <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                       check.names = FALSE, row.names = 1)
id2sym_tc <- setNames(degs$gene_name, degs$gene_id)
id2sym_tc <- id2sym_tc[!is.na(id2sym_tc) & id2sym_tc != ""]
tc_sym <- id2sym_tc[rownames(tc_vst_raw)]
keep_tc <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst_raw[keep_tc, ])
rownames(tc_mat) <- tc_sym[keep_tc]

# Proteomics
pr_raw <- read.csv(file.path(PROC_DIR, "proteomics_log2_norm.csv"),
                   check.names = FALSE, row.names = 1)
id2sym_pr <- setNames(deps$gene_name, deps$Protein)
id2sym_pr <- id2sym_pr[!is.na(id2sym_pr) & id2sym_pr != "" & id2sym_pr != "_"]
pr_sym <- id2sym_pr[rownames(pr_raw)]
keep_pr <- !is.na(pr_sym) & pr_sym != "" & !duplicated(pr_sym)
pr_mat <- as.matrix(pr_raw[keep_pr, ])
rownames(pr_mat) <- pr_sym[keep_pr]

# Common samples and genes
shared_samples <- intersect(colnames(tc_mat), colnames(pr_mat))
shared_genes <- intersect(rownames(tc_mat), rownames(pr_mat))

cat(sprintf("  Shared samples: %d, Shared genes (TC+PR): %d\n",
            length(shared_samples), length(shared_genes)))

# Group info
group_vec <- ifelse(grepl("^Normal", shared_samples), 0, 1)  # 0=Normal, 1=Adjacent
names(group_vec) <- shared_samples

# Load paired subtypes
paired_sub <- read.csv(file.path(E1_DIR, "paired_subtype_K2.csv"), stringsAsFactors = FALSE)
subtype_map <- setNames(ifelse(paired_sub$subtype == "DS1", 1, 2), paired_sub$patient_id)

# Map subtype to samples (Adjacent only)
sample_subtype <- rep(NA, length(shared_samples))
names(sample_subtype) <- shared_samples
for (s in shared_samples) {
  pid <- gsub("^(Normal|Adjacent)", "", s)
  if (pid %in% names(subtype_map)) sample_subtype[s] <- subtype_map[pid]
}


# ============================================================================
# 1. TRANSCRIPTOMICS WGCNA
# ============================================================================
cat("\n>>> 1. Transcriptomics WGCNA\n")

# Select top variable genes
tc_vars <- apply(tc_mat[, shared_samples], 1, var, na.rm = TRUE)
tc_top_genes <- names(head(sort(tc_vars, decreasing = TRUE), 5000))
tc_expr <- t(tc_mat[tc_top_genes, shared_samples])  # samples x genes

# Check for good genes/samples
gsg <- goodSamplesGenes(tc_expr, verbose = 0)
if (!gsg$allOK) {
  tc_expr <- tc_expr[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("  Removed bad genes/samples. Remaining: %d x %d\n", nrow(tc_expr), ncol(tc_expr)))
}
cat(sprintf("  TC WGCNA input: %d samples x %d genes\n", nrow(tc_expr), ncol(tc_expr)))

# Pick soft threshold
cat("  Picking soft threshold...\n")
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(tc_expr, powerVector = powers, verbose = 0, networkType = "signed")

# Find threshold where R^2 > 0.8
r2_vals <- sft$fitIndices$SFT.R.sq * sign(sft$fitIndices$slope)
best_power <- sft$powerEstimate
if (is.na(best_power) || best_power < 4) best_power <- 6
cat(sprintf("  Best power: %d (R2=%.2f)\n", best_power, sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == best_power]))

# Save power selection plot
pdf(file.path(FIG_DIR, "TC_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq * sign(sft$fitIndices$slope),
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n", main = "Scale Independence (TC)")
text(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq * sign(sft$fitIndices$slope),
     labels = powers, col = "red", cex = 0.9)
abline(h = 0.8, col = "blue", lty = 2)

plot(sft$fitIndices$Power, sft$fitIndices$mean.k,
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity (TC)")
text(sft$fitIndices$Power, sft$fitIndices$mean.k, labels = powers, col = "red", cex = 0.9)
dev.off()

# Build network
cat("  Building TC network (may take a few minutes)...\n")
tc_net <- blockwiseModules(tc_expr,
                            power = best_power,
                            networkType = "signed",
                            TOMType = "signed",
                            minModuleSize = 30,
                            reassignThreshold = 0,
                            mergeCutHeight = 0.25,
                            numericLabels = TRUE,
                            saveTOMs = FALSE,
                            verbose = 0)

tc_module_colors <- labels2colors(tc_net$colors)
cat(sprintf("  TC modules: %d (excluding grey)\n", length(unique(tc_module_colors)) - 1))
cat(sprintf("  Module sizes: %s\n",
            paste(sprintf("%s=%d", names(table(tc_module_colors)),
                          table(tc_module_colors)), collapse = ", ")))

# Save module assignments
tc_mod_df <- data.frame(
  gene = colnames(tc_expr),
  module_num = tc_net$colors,
  module_color = tc_module_colors,
  stringsAsFactors = FALSE
)
write.csv(tc_mod_df, file.path(OUT_DIR, "TC_module_assignments.csv"), row.names = FALSE)


# ============================================================================
# 2. MODULE-TRAIT CORRELATION
# ============================================================================
cat("\n>>> 2. Module-trait correlations\n")

# Module eigengenes
tc_MEs <- tc_net$MEs
colnames(tc_MEs) <- gsub("^ME", "", colnames(tc_MEs))

# Trait matrix
trait_mat <- data.frame(
  Group = group_vec[rownames(tc_expr)],
  Subtype = sample_subtype[rownames(tc_expr)]
)

# Load clinical data for Adjacent samples
clin_file <- file.path(E1_DIR, "paired_clinical_with_subtypes_K2.csv")
if (file.exists(clin_file)) {
  clin <- read.csv(clin_file, stringsAsFactors = FALSE)
  
  for (s in rownames(tc_expr)) {
    pid <- gsub("^(Normal|Adjacent)", "", s)
    idx <- match(pid, clin$patient_id)
    if (!is.na(idx)) {
      trait_mat[s, "Age"] <- as.numeric(clin$age[idx])
      trait_mat[s, "BMI"] <- as.numeric(clin$BMI[idx])
      trait_mat[s, "Duration"] <- as.numeric(clin$disease_duration[idx])
      trait_mat[s, "LesionSize"] <- as.numeric(clin$lesion_size[idx])
      trait_mat[s, "ALT"] <- as.numeric(clin$ALT[idx])
      trait_mat[s, "AST"] <- as.numeric(clin$AST[idx])
      trait_mat[s, "ALP"] <- as.numeric(clin$ALP[idx])
      trait_mat[s, "GGT"] <- as.numeric(clin$GGT[idx])
      trait_mat[s, "Bilirubin"] <- as.numeric(clin$total_bilirubin[idx])
      trait_mat[s, "Albumin"] <- as.numeric(clin$albumin[idx])
      trait_mat[s, "CRP"] <- as.numeric(clin$CRP[idx])
      trait_mat[s, "WBC"] <- as.numeric(clin$WBC[idx])
      trait_mat[s, "BileInvasion"] <- as.numeric(clin$bile_duct_invasion[idx])
    }
  }
}

trait_mat_num <- as.data.frame(lapply(trait_mat, function(x) as.numeric(as.character(x))))
rownames(trait_mat_num) <- rownames(trait_mat)

# Correlation
cor_result <- cor(tc_MEs, trait_mat_num, use = "pairwise.complete.obs")
p_result <- corPvalueStudent(cor_result, nrow(tc_expr))

# Save
write.csv(cor_result, file.path(OUT_DIR, "TC_module_trait_cor.csv"))
write.csv(p_result, file.path(OUT_DIR, "TC_module_trait_pval.csv"))

# Print significant module-trait associations
cat("  Significant module-trait associations (P<0.05):\n")
for (i in 1:nrow(cor_result)) {
  for (j in 1:ncol(cor_result)) {
    if (!is.na(p_result[i, j]) && p_result[i, j] < 0.05) {
      cat(sprintf("    ME%s ~ %s: r=%.3f, P=%.4f\n",
                  rownames(cor_result)[i], colnames(cor_result)[j],
                  cor_result[i, j], p_result[i, j]))
    }
  }
}

# Module-trait heatmap
# Create text matrix for display
text_mat <- paste(signif(cor_result, 2), "\n(",
                  signif(p_result, 1), ")", sep = "")
dim(text_mat) <- dim(cor_result)

# Remove grey module
keep_mods <- !grepl("^0$", rownames(cor_result))
cor_plot <- cor_result[keep_mods, , drop = FALSE]
p_plot <- p_result[keep_mods, , drop = FALSE]

# Rename modules
rownames(cor_plot) <- paste0("ME", labels2colors(as.numeric(rownames(cor_plot))))
rownames(p_plot) <- rownames(cor_plot)

pdf(file.path(FIG_DIR, "TC_module_trait_heatmap.pdf"), width = 12, height = max(6, nrow(cor_plot) * 0.5 + 2))
pheatmap(cor_plot,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = matrix(ifelse(p_plot < 0.05, sprintf("%.2f*", cor_plot),
                                          sprintf("%.2f", cor_plot)),
                                   nrow = nrow(cor_plot)),
         number_color = "black",
         fontsize_number = 7,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         main = "Module-Trait Relationships (TC WGCNA)",
         fontsize_row = 9,
         fontsize_col = 9)
dev.off()
cat("  Module-trait heatmap saved.\n")


# ============================================================================
# 3. HUB GENE IDENTIFICATION
# ============================================================================
cat("\n>>> 3. Hub gene identification\n")

# Identify group-associated modules
group_cor <- cor_result[, "Group", drop = TRUE]
sig_modules_group <- names(which(p_result[, "Group"] < 0.05 & abs(group_cor) > 0.3))

# Identify subtype-associated modules
if ("Subtype" %in% colnames(cor_result)) {
  sub_cor <- cor_result[, "Subtype", drop = TRUE]
  sig_modules_subtype <- names(which(!is.na(p_result[, "Subtype"]) &
                                       p_result[, "Subtype"] < 0.05 &
                                       abs(sub_cor) > 0.3))
} else {
  sig_modules_subtype <- character(0)
}

key_modules <- unique(c(sig_modules_group, sig_modules_subtype))
cat(sprintf("  Key modules (group/subtype-associated): %d\n", length(key_modules)))

hub_genes_all <- data.frame()

for (mod_num in key_modules) {
  mod_color <- labels2colors(as.numeric(mod_num))
  mod_genes <- colnames(tc_expr)[tc_net$colors == as.numeric(mod_num)]
  
  if (length(mod_genes) < 5) next
  
  # Module membership (kME)
  kME <- cor(tc_expr[, mod_genes], tc_MEs[, mod_num, drop = FALSE], use = "pairwise.complete.obs")
  
  # Gene significance for Group
  gs_group <- cor(tc_expr[, mod_genes], trait_mat_num[, "Group", drop = FALSE],
                  use = "pairwise.complete.obs")
  
  hub_df <- data.frame(
    gene = mod_genes,
    module = mod_color,
    module_num = as.numeric(mod_num),
    kME = kME[, 1],
    GS_group = gs_group[, 1],
    stringsAsFactors = FALSE
  )
  hub_df$abs_kME <- abs(hub_df$kME)
  hub_df$abs_GS <- abs(hub_df$GS_group)
  hub_df <- hub_df[order(-hub_df$abs_kME), ]
  
  # Top 10 hub genes (high kME and high GS)
  hub_df$hub_score <- hub_df$abs_kME * hub_df$abs_GS
  hub_df <- hub_df[order(-hub_df$hub_score), ]
  
  cat(sprintf("  Module %s (%d genes): top hubs = %s\n",
              mod_color, length(mod_genes),
              paste(head(hub_df$gene, 5), collapse = ", ")))
  
  hub_genes_all <- rbind(hub_genes_all, head(hub_df, 20))
}

write.csv(hub_genes_all, file.path(OUT_DIR, "TC_hub_genes.csv"), row.names = FALSE)
cat(sprintf("  Total hub genes identified: %d\n", nrow(hub_genes_all)))

# ---- Hub gene scatter plot (kME vs GS) ----
if (nrow(hub_genes_all) >= 5) {
  p_hub <- ggplot(hub_genes_all, aes(abs_kME, abs_GS, color = module)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_text_repel(data = hub_genes_all %>% group_by(module) %>% slice_max(hub_score, n = 3),
                    aes(label = gene), size = 3, max.overlaps = 20) +
    labs(title = "Hub Genes: Module Membership vs Gene Significance",
         x = "|kME| (Module Membership)",
         y = "|GS| (Gene Significance for Group)") +
    theme_bw(base_size = 12)
  ggsave(file.path(FIG_DIR, "TC_hub_gene_scatter.pdf"), p_hub, width = 10, height = 7)
  cat("  Hub gene scatter plot saved.\n")
}


# ============================================================================
# 4. PROTEOMICS WGCNA (ABBREVIATED)
# ============================================================================
cat("\n>>> 4. Proteomics WGCNA\n")

pr_vars <- apply(pr_mat[, shared_samples], 1, var, na.rm = TRUE)
pr_top <- names(head(sort(pr_vars, decreasing = TRUE), 3000))
pr_expr <- t(pr_mat[pr_top, shared_samples])

gsg_pr <- goodSamplesGenes(pr_expr, verbose = 0)
if (!gsg_pr$allOK) pr_expr <- pr_expr[gsg_pr$goodSamples, gsg_pr$goodGenes]
cat(sprintf("  PR WGCNA input: %d samples x %d proteins\n", nrow(pr_expr), ncol(pr_expr)))

sft_pr <- pickSoftThreshold(pr_expr, powerVector = powers, verbose = 0, networkType = "signed")
pr_power <- sft_pr$powerEstimate
if (is.na(pr_power) || pr_power < 4) pr_power <- 6
cat(sprintf("  PR best power: %d\n", pr_power))

cat("  Building PR network...\n")
pr_net <- blockwiseModules(pr_expr,
                            power = pr_power,
                            networkType = "signed",
                            TOMType = "signed",
                            minModuleSize = 20,
                            mergeCutHeight = 0.25,
                            numericLabels = TRUE,
                            saveTOMs = FALSE,
                            verbose = 0)

pr_module_colors <- labels2colors(pr_net$colors)
cat(sprintf("  PR modules: %d (excluding grey)\n", length(unique(pr_module_colors)) - 1))

# Save
pr_mod_df <- data.frame(gene = colnames(pr_expr), module_num = pr_net$colors,
                          module_color = pr_module_colors, stringsAsFactors = FALSE)
write.csv(pr_mod_df, file.path(OUT_DIR, "PR_module_assignments.csv"), row.names = FALSE)

# Module-trait correlation
pr_MEs <- pr_net$MEs
colnames(pr_MEs) <- gsub("^ME", "", colnames(pr_MEs))

cor_pr <- cor(pr_MEs, trait_mat_num[rownames(pr_expr), ], use = "pairwise.complete.obs")
p_pr <- corPvalueStudent(cor_pr, nrow(pr_expr))

write.csv(cor_pr, file.path(OUT_DIR, "PR_module_trait_cor.csv"))

cat("  PR significant module-trait associations:\n")
for (i in 1:nrow(cor_pr)) {
  for (j in 1:ncol(cor_pr)) {
    if (!is.na(p_pr[i, j]) && p_pr[i, j] < 0.05 && abs(cor_pr[i, j]) > 0.3) {
      cat(sprintf("    ME%s ~ %s: r=%.3f, P=%.4f\n",
                  labels2colors(as.numeric(rownames(cor_pr)[i])),
                  colnames(cor_pr)[j], cor_pr[i, j], p_pr[i, j]))
    }
  }
}


# ============================================================================
# 5. CROSS-OMICS MODULE PRESERVATION
# ============================================================================
cat("\n>>> 5. Cross-omics module overlap\n")

# Find overlap between TC and PR modules for shared genes
tc_mods <- setNames(tc_mod_df$module_color, tc_mod_df$gene)
pr_mods <- setNames(pr_mod_df$module_color, pr_mod_df$gene)
common <- intersect(names(tc_mods), names(pr_mods))
cat(sprintf("  Shared genes between TC and PR WGCNA: %d\n", length(common)))

# Initialize sig_overlap before conditional block to avoid reference error
sig_overlap <- data.frame()

if (length(common) > 50) {
  # Contingency table of module assignments
  ct <- table(TC = tc_mods[common], PR = pr_mods[common])
  
  # Fisher's test for enrichment of TC modules in PR modules
  overlap_results <- data.frame()
  for (tc_m in rownames(ct)) {
    if (tc_m == "grey") next
    for (pr_m in colnames(ct)) {
      if (pr_m == "grey") next
      n_overlap <- ct[tc_m, pr_m]
      if (n_overlap < 3) next
      
      # 2x2 table
      a <- n_overlap
      b <- sum(ct[tc_m, ]) - a
      c <- sum(ct[, pr_m]) - a
      d <- length(common) - a - b - c
      
      ft <- fisher.test(matrix(c(a, b, c, d), 2, 2))
      overlap_results <- rbind(overlap_results, data.frame(
        TC_module = tc_m, PR_module = pr_m,
        overlap = n_overlap, pvalue = ft$p.value, OR = ft$estimate,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  overlap_results$padj <- p.adjust(overlap_results$pvalue, method = "BH")
  overlap_results <- overlap_results[order(overlap_results$pvalue), ]
  write.csv(overlap_results, file.path(OUT_DIR, "cross_omics_module_overlap.csv"), row.names = FALSE)
  
  sig_overlap <- overlap_results[overlap_results$padj < 0.05, ]
  cat(sprintf("  Significant cross-omics module overlaps: %d\n", nrow(sig_overlap)))
  if (nrow(sig_overlap) > 0) {
    for (i in 1:min(10, nrow(sig_overlap))) {
      cat(sprintf("    TC_%s ~ PR_%s: %d genes, OR=%.1f, padj=%.4f\n",
                  sig_overlap$TC_module[i], sig_overlap$PR_module[i],
                  sig_overlap$overlap[i], sig_overlap$OR[i], sig_overlap$padj[i]))
    }
  }
}


# ============================================================================
# 6. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Enhancement 2 SUMMARY\n")
cat("========================================\n")
cat(sprintf("TC WGCNA: %d modules, power=%d, %d genes\n",
            length(unique(tc_module_colors)) - 1, best_power, ncol(tc_expr)))
cat(sprintf("PR WGCNA: %d modules, power=%d, %d proteins\n",
            length(unique(pr_module_colors)) - 1, pr_power, ncol(pr_expr)))
cat(sprintf("Hub genes: %d from %d key modules\n",
            nrow(hub_genes_all), length(key_modules)))
cat(sprintf("Cross-omics: %d shared genes, %d sig overlaps\n",
            length(common), nrow(sig_overlap)))
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Enhancement 2 COMPLETE.\n")
