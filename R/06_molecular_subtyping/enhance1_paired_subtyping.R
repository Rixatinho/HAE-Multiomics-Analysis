#!/usr/bin/env Rscript
# ============================================================================
# Enhancement 1: Paired-Difference Molecular Subtyping
# ============================================================================
# Key improvement: Instead of clustering raw Adjacent tissue profiles,
# compute per-patient log2(Adjacent/Normal) difference vectors across all
# 3 omics, then perform consensus clustering on these "disease perturbation"
# profiles. This removes inter-individual baseline variation and captures
# the magnitude/direction of disease-induced molecular changes.
# ============================================================================

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(cluster)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(RColorBrewer)
  library(dplyr)
  library(readxl)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
P6_DIR   <- file.path(PROJECT, "analysis/results/phase6_subtyping")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhancement1_paired_subtyping")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Enhancement 1: Paired-Difference Subtyping\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD & COMPUTE PAIRED DIFFERENCES
# ============================================================================
cat(">>> 0. Computing paired differences\n")

# Load ID mapping
degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# ---- Transcriptomics ----
tc_vst_raw <- read.csv(file.path(PROC_DIR, "transcriptomics_vst_paired.csv"),
                       check.names = FALSE, row.names = 1)
id2sym_tc <- setNames(degs$gene_name, degs$gene_id)
id2sym_tc <- id2sym_tc[!is.na(id2sym_tc) & id2sym_tc != ""]
tc_sym <- id2sym_tc[rownames(tc_vst_raw)]
keep_tc <- !is.na(tc_sym) & tc_sym != "" & !duplicated(tc_sym)
tc_mat <- as.matrix(tc_vst_raw[keep_tc, ])
rownames(tc_mat) <- tc_sym[keep_tc]

# ---- Proteomics ----
pr_raw <- read.csv(file.path(PROC_DIR, "proteomics_log2_norm.csv"),
                   check.names = FALSE, row.names = 1)
id2sym_pr <- setNames(deps$gene_name, deps$Protein)
id2sym_pr <- id2sym_pr[!is.na(id2sym_pr) & id2sym_pr != "" & id2sym_pr != "_"]
pr_sym <- id2sym_pr[rownames(pr_raw)]
keep_pr <- !is.na(pr_sym) & pr_sym != "" & !duplicated(pr_sym)
pr_mat <- as.matrix(pr_raw[keep_pr, ])
rownames(pr_mat) <- pr_sym[keep_pr]

# ---- Metabolomics ----
met_mat <- as.matrix(read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"),
                               check.names = FALSE, row.names = 1))

# ---- Identify paired samples ----
shared_all <- Reduce(intersect, list(colnames(tc_mat), colnames(pr_mat), colnames(met_mat)))
adj_samples <- grep("^Adjacent", shared_all, value = TRUE)
nor_samples <- grep("^Normal", shared_all, value = TRUE)

# Patient IDs
adj_pids <- gsub("^Adjacent", "", adj_samples)
nor_pids <- gsub("^Normal", "", nor_samples)
paired_pids <- intersect(adj_pids, nor_pids)
cat(sprintf("  Paired patients: %d\n", length(paired_pids)))

# ---- Compute per-patient differences (Adjacent - Normal) ----
# For VST/log2 data, difference = log2(fold change)
compute_paired_diff <- function(mat, pids) {
  diff_mat <- matrix(NA, nrow = nrow(mat), ncol = length(pids))
  rownames(diff_mat) <- rownames(mat)
  colnames(diff_mat) <- paste0("P", pids)
  
  for (i in seq_along(pids)) {
    adj_col <- paste0("Adjacent", pids[i])
    nor_col <- paste0("Normal", pids[i])
    if (adj_col %in% colnames(mat) && nor_col %in% colnames(mat)) {
      diff_mat[, i] <- mat[, adj_col] - mat[, nor_col]
    }
  }
  
  # Remove rows with all NA
  diff_mat <- diff_mat[rowSums(is.na(diff_mat)) < ncol(diff_mat), ]
  return(diff_mat)
}

tc_diff <- compute_paired_diff(tc_mat, paired_pids)
pr_diff <- compute_paired_diff(pr_mat, paired_pids)
met_diff <- compute_paired_diff(met_mat, paired_pids)

cat(sprintf("  TC diff: %d genes x %d patients\n", nrow(tc_diff), ncol(tc_diff)))
cat(sprintf("  PR diff: %d proteins x %d patients\n", nrow(pr_diff), ncol(pr_diff)))
cat(sprintf("  MET diff: %d metabolites x %d patients\n", nrow(met_diff), ncol(met_diff)))


# ============================================================================
# 1. FEATURE SELECTION & CONCATENATION
# ============================================================================
cat("\n>>> 1. Feature selection for paired-diff clustering\n")

select_top_var <- function(mat, n) {
  mat <- mat[complete.cases(mat), ]
  vars <- apply(mat, 1, var, na.rm = TRUE)
  top <- head(order(vars, decreasing = TRUE), min(n, nrow(mat)))
  mat[top, ]
}

tc_sel <- select_top_var(tc_diff, 1000)
pr_sel <- select_top_var(pr_diff, 1000)
met_sel <- select_top_var(met_diff, 300)

# Scale
tc_scaled <- t(scale(t(tc_sel)))
pr_scaled <- t(scale(t(pr_sel)))
met_scaled <- t(scale(t(met_sel)))

# Remove NA rows from scaling
tc_scaled <- tc_scaled[complete.cases(tc_scaled), ]
pr_scaled <- pr_scaled[complete.cases(pr_scaled), ]
met_scaled <- met_scaled[complete.cases(met_scaled), ]

# Prefix and concatenate
rownames(tc_scaled) <- paste0("TC_", rownames(tc_scaled))
rownames(pr_scaled) <- paste0("PR_", rownames(pr_scaled))
rownames(met_scaled) <- paste0("MET_", rownames(met_scaled))

concat_diff <- rbind(tc_scaled, pr_scaled, met_scaled)
cat(sprintf("  Concatenated diff matrix: %d features x %d patients\n",
            nrow(concat_diff), ncol(concat_diff)))


# ============================================================================
# 2. CONSENSUS CLUSTERING ON PAIRED DIFFERENCES
# ============================================================================
cat("\n>>> 2. Consensus clustering on paired differences\n")

cc_dir <- file.path(FIG_DIR, "consensus_paired")
dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)

cc_paired <- ConsensusClusterPlus(
  d = concat_diff,
  maxK = 5,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "pearson",
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2",
  seed = 42,
  title = cc_dir,
  plot = "pdf"
)

icl <- calcICL(cc_paired, title = cc_dir, plot = "pdf")

# ---- Metrics ----
compute_PAC <- function(cc_result, maxK) {
  pac <- numeric(maxK - 1)
  for (k in 2:maxK) {
    M <- cc_result[[k]]$consensusMatrix
    m_lower <- M[lower.tri(M)]
    pac[k - 1] <- mean(m_lower > 0.1 & m_lower < 0.9)
  }
  names(pac) <- paste0("K=", 2:maxK)
  return(pac)
}

compute_sil <- function(cc_result, maxK) {
  sil <- numeric(maxK - 1)
  for (k in 2:maxK) {
    M <- cc_result[[k]]$consensusMatrix
    d <- as.dist(1 - M)
    cl <- cc_result[[k]]$consensusClass
    if (length(unique(cl)) >= 2) {
      s <- silhouette(cl, d)
      sil[k - 1] <- mean(s[, "sil_width"])
    } else sil[k - 1] <- NA
  }
  names(sil) <- paste0("K=", 2:maxK)
  return(sil)
}

pac <- compute_PAC(cc_paired, 5)
sil <- compute_sil(cc_paired, 5)

cat("  PAC:", paste(sprintf("%s=%.3f", names(pac), pac), collapse = ", "), "\n")
cat("  Silhouette:", paste(sprintf("%s=%.3f", names(sil), sil), collapse = ", "), "\n")

best_k_pac <- which.min(pac) + 1
best_k_sil <- which.max(sil) + 1
cat(sprintf("  Best K by PAC: %d, by silhouette: %d\n", best_k_pac, best_k_sil))

# K selection plot
pdf(file.path(FIG_DIR, "paired_K_selection.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(2:5, pac, type = "b", pch = 16, col = "#3C5488", lwd = 2,
     xlab = "K", ylab = "PAC", main = "PAC (Paired-Difference)")
plot(2:5, sil, type = "b", pch = 16, col = "#E64B35", lwd = 2,
     xlab = "K", ylab = "Mean Silhouette Width", main = "Silhouette (Paired-Difference)")
dev.off()

# Save K metrics
write.csv(data.frame(K = 2:5, PAC = pac, Silhouette = sil),
          file.path(OUT_DIR, "paired_K_metrics.csv"), row.names = FALSE)


# ============================================================================
# 3. ASSIGN SUBTYPES FOR K=2 AND K=3
# ============================================================================
cat("\n>>> 3. Subtype assignment\n")

for (k in 2:min(3, 5)) {
  cl <- cc_paired[[k]]$consensusClass
  labels <- paste0("DS", cl)  # DS = Disease Subtype
  names(labels) <- names(cl)
  
  cat(sprintf("  K=%d: %s\n", k,
              paste(sprintf("DS%d=%d", 1:k, table(factor(cl, levels = 1:k))), collapse = ", ")))
  
  # Save assignment
  df <- data.frame(
    patient = names(cl),
    patient_id = gsub("^P", "", names(cl)),
    subtype = labels,
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(OUT_DIR, sprintf("paired_subtype_K%d.csv", k)), row.names = FALSE)
  
  # Consensus heatmap
  cons_mat <- cc_paired[[k]]$consensusMatrix
  rownames(cons_mat) <- colnames(cons_mat) <- names(cl)
  
  anno_df <- data.frame(Subtype = labels)
  rownames(anno_df) <- names(cl)
  subtype_colors <- c("DS1" = "#3C5488", "DS2" = "#E64B35", "DS3" = "#00A087")[1:k]
  
  pdf(file.path(FIG_DIR, sprintf("paired_consensus_K%d.pdf", k)), width = 8, height = 7)
  print(pheatmap(cons_mat,
           color = colorRampPalette(c("white", "#3C5488"))(100),
           annotation_col = anno_df, annotation_row = anno_df,
           annotation_colors = list(Subtype = subtype_colors),
           clustering_method = "ward.D2",
           main = sprintf("Paired-Difference Consensus (K=%d, 1000 reps)", k),
           fontsize = 10))
  dev.off()
  
  # Silhouette
  d <- as.dist(1 - cons_mat)
  s <- silhouette(cl, d)
  pdf(file.path(FIG_DIR, sprintf("paired_silhouette_K%d.pdf", k)), width = 8, height = 6)
  plot(s, col = subtype_colors[1:k],
       main = sprintf("Silhouette (K=%d, mean=%.3f)", k, mean(s[, "sil_width"])),
       border = NA)
  dev.off()
}


# ============================================================================
# 4. CLINICAL ASSOCIATION (KEY IMPROVEMENT)
# ============================================================================
cat("\n>>> 4. Clinical association analysis\n")

# Parse clinical data
clin_raw <- read_excel(file.path(PROJECT, "14 例肝泡型棘球蚴病患者临床信息.xlsx"), col_names = FALSE)
col_names <- as.character(clin_raw[3, ])
clin_data <- clin_raw[4:17, ]
colnames(clin_data) <- col_names

clinical <- data.frame(
  patient_id = as.character(clin_data[["编号"]]),
  age = as.numeric(clin_data[["年龄"]]),
  sex = as.character(clin_data[["性别"]]),
  ethnicity = as.character(clin_data[["民族"]]),
  endemic_area = as.character(clin_data[["居住地是否为疫区"]]),
  BMI = as.numeric(clin_data[["BMI"]]),
  disease_duration = as.numeric(clin_data[["确诊时间（年）"]]),
  surgery_history = as.character(clin_data[["是否接受过手术"]]),
  ABZ_treatment = as.character(clin_data[["阿苯达唑用药史"]]),
  lesion_size = as.numeric(clin_data[["大小（cm2）"]]),
  lesion_count = as.numeric(clin_data[["数目"]]),
  lesion_location = as.character(clin_data[["位置"]]),
  stringsAsFactors = FALSE
)

# PNM & PIVM staging
pnm_cols <- which(col_names %in% c("病灶", "临近器官", "转移病灶"))
if (length(pnm_cols) >= 3) {
  clinical$PNM_P <- as.character(clin_data[[pnm_cols[1]]])
  clinical$PNM_N <- as.character(clin_data[[pnm_cols[2]]])
  clinical$PNM_M <- as.character(clin_data[[pnm_cols[3]]])
}

pivm_cols <- which(col_names %in% c("病灶", "侵犯胆道", "血管侵犯", "转移病灶"))
if (length(pivm_cols) >= 4) {
  clinical$PIVM_stage <- as.character(clin_data[[pivm_cols[2] - 1]])
  clinical$bile_duct_invasion <- as.character(clin_data[[pivm_cols[2]]])
  clinical$vascular_invasion <- as.character(clin_data[[pivm_cols[3]]])
  clinical$PIVM_M <- as.character(clin_data[[pivm_cols[4]]])
}

# Lab values
clinical$WBC <- as.numeric(clin_data[["白细胞"]])
clinical$neutrophil_pct <- as.numeric(clin_data[["中性粒细胞（%）"]])
clinical$lymphocyte_pct <- as.numeric(clin_data[["淋巴细胞（%）"]])
clinical$ALT <- as.numeric(clin_data[["ALT(U/L)"]])
clinical$AST <- as.numeric(clin_data[["AST(U/L)"]])
clinical$ALP <- as.numeric(clin_data[["ALP(U)"]])
clinical$GGT <- as.numeric(clin_data[["GGT(U/L)"]])
clinical$total_bilirubin <- as.numeric(clin_data[["总胆红素（umol/L）"]])
clinical$albumin <- as.numeric(clin_data[["白蛋白(g/L)"]])

crp_raw <- as.character(clin_data[["CRP（mg/L)"]])
clinical$CRP <- suppressWarnings(as.numeric(gsub("[<>]", "", crp_raw)))
il6_raw <- as.character(clin_data[["IL-6（pg/ml）"]])
clinical$IL6 <- suppressWarnings(as.numeric(gsub("[<>]", "", il6_raw)))

clinical$glucose <- as.numeric(clin_data[["葡萄糖(mmol/L)"]])
clinical$lactate <- as.numeric(clin_data[["乳酸(mmol/L)"]])
clinical$cholesterol <- as.numeric(clin_data[["胆固醇（mmol/L）"]])

# Derived clinical variables
clinical$NLR <- clinical$neutrophil_pct / clinical$lymphocyte_pct  # Neutrophil-lymphocyte ratio
clinical$PNM_sum <- suppressWarnings(as.numeric(clinical$PNM_P)) +
                    suppressWarnings(as.numeric(clinical$PNM_N)) +
                    suppressWarnings(as.numeric(clinical$PNM_M))

# Categorize: advanced vs early
clinical$advanced_PNM <- ifelse(!is.na(clinical$PNM_sum) & clinical$PNM_sum >= 3, "Advanced", "Early")
clinical$high_inflammation <- ifelse(!is.na(clinical$CRP) & clinical$CRP > 10, "High", "Low")
clinical$liver_damage <- ifelse(!is.na(clinical$ALT) & clinical$ALT > 40, "Elevated", "Normal")

cat(sprintf("  Clinical data: %d patients x %d variables\n", nrow(clinical), ncol(clinical)))

# ---- Test associations for K=2 and K=3 ----
run_association <- function(clin_df, subtype_col) {
  results <- data.frame()
  st <- clin_df[[subtype_col]]
  
  if (length(unique(st[!is.na(st)])) < 2) return(results)
  
  # Continuous variables
  cont_vars <- c("age", "BMI", "disease_duration", "lesion_size", "lesion_count",
                  "WBC", "neutrophil_pct", "lymphocyte_pct", "NLR",
                  "ALT", "AST", "ALP", "GGT", "total_bilirubin", "albumin",
                  "CRP", "IL6", "glucose", "lactate", "cholesterol", "PNM_sum")
  
  for (v in cont_vars) {
    if (!v %in% colnames(clin_df)) next
    vals <- as.numeric(clin_df[[v]])
    if (sum(!is.na(vals)) < 4) next
    
    tryCatch({
      if (length(unique(st)) == 2) {
        test <- wilcox.test(vals ~ st)
        test_name <- "Wilcoxon"
      } else {
        test <- kruskal.test(vals ~ factor(st))
        test_name <- "Kruskal-Wallis"
      }
      
      means <- tapply(vals, st, mean, na.rm = TRUE)
      medians <- tapply(vals, st, median, na.rm = TRUE)
      
      results <- rbind(results, data.frame(
        variable = v, type = "continuous", test = test_name,
        pvalue = test$p.value,
        subtype_means = paste(sprintf("%s=%.2f", names(means), means), collapse="; "),
        subtype_medians = paste(sprintf("%s=%.2f", names(medians), medians), collapse="; "),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  
  # Categorical variables
  cat_vars <- c("sex", "ethnicity", "endemic_area", "surgery_history", "ABZ_treatment",
                "PNM_P", "PNM_N", "PNM_M", "bile_duct_invasion", "PIVM_stage", "PIVM_M",
                "advanced_PNM", "high_inflammation", "liver_damage")
  
  for (v in cat_vars) {
    if (!v %in% colnames(clin_df)) next
    vals <- clin_df[[v]]
    if (sum(!is.na(vals)) < 4 || length(unique(vals[!is.na(vals)])) < 2) next
    
    tryCatch({
      tab <- table(st, vals)
      ft <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
      
      results <- rbind(results, data.frame(
        variable = v, type = "categorical", test = "Fisher",
        pvalue = ft$p.value,
        subtype_means = "", subtype_medians = "",
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  
  results <- results[order(results$pvalue), ]
  return(results)
}

# Merge clinical with subtypes
for (k in 2:3) {
  sub_df <- read.csv(file.path(OUT_DIR, sprintf("paired_subtype_K%d.csv", k)), stringsAsFactors = FALSE)
  clin_merged <- merge(clinical, sub_df[, c("patient_id", "subtype")],
                       by = "patient_id", all.x = TRUE)
  clin_typed <- clin_merged[!is.na(clin_merged$subtype), ]
  
  cat(sprintf("\n  === K=%d Clinical Associations ===\n", k))
  cat(sprintf("  Patients with subtype: %d\n", nrow(clin_typed)))
  
  assoc <- run_association(clin_typed, "subtype")
  
  sig <- assoc[assoc$pvalue < 0.05, ]
  trend <- assoc[assoc$pvalue < 0.1 & assoc$pvalue >= 0.05, ]
  
  cat(sprintf("  Total tests: %d\n", nrow(assoc)))
  cat(sprintf("  Significant (P<0.05): %d\n", nrow(sig)))
  cat(sprintf("  Trends (P<0.1): %d\n", nrow(trend)))
  
  if (nrow(sig) > 0) {
    cat("  Significant associations:\n")
    for (i in 1:nrow(sig)) {
      cat(sprintf("    * %s: P=%.4f [%s] %s\n",
                  sig$variable[i], sig$pvalue[i], sig$type[i], sig$subtype_means[i]))
    }
  }
  
  if (nrow(trend) > 0) {
    cat("  Trend associations:\n")
    for (i in 1:nrow(trend)) {
      cat(sprintf("    ~ %s: P=%.4f [%s] %s\n",
                  trend$variable[i], trend$pvalue[i], trend$type[i], trend$subtype_means[i]))
    }
  }
  
  write.csv(assoc, file.path(OUT_DIR, sprintf("paired_clinical_assoc_K%d.csv", k)), row.names = FALSE)
  write.csv(clin_typed, file.path(OUT_DIR, sprintf("paired_clinical_with_subtypes_K%d.csv", k)), row.names = FALSE)
}


# ============================================================================
# 5. PAIRED-DIFFERENCE LANDSCAPE VISUALIZATION
# ============================================================================
cat("\n>>> 5. Visualization\n")

# ---- PCA on paired differences ----
pca_res <- prcomp(t(concat_diff), scale. = FALSE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$patient <- rownames(pca_df)
var_exp <- summary(pca_res)$importance[2, 1:2] * 100

# Add K=2 subtypes
cl2 <- cc_paired[[2]]$consensusClass
pca_df$subtype_K2 <- paste0("DS", cl2[pca_df$patient])

# Add K=3 if available
if (length(cc_paired) >= 3) {
  cl3 <- cc_paired[[3]]$consensusClass
  pca_df$subtype_K3 <- paste0("DS", cl3[pca_df$patient])
}

p1 <- ggplot(pca_df, aes(PC1, PC2, color = subtype_K2)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = patient), size = 3) +
  scale_color_manual(values = c("DS1" = "#3C5488", "DS2" = "#E64B35")) +
  labs(title = "Paired-Difference PCA (K=2)",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

p2 <- ggplot(pca_df, aes(PC1, PC2, color = subtype_K3)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = patient), size = 3) +
  scale_color_manual(values = c("DS1" = "#3C5488", "DS2" = "#E64B35", "DS3" = "#00A087")) +
  labs(title = "Paired-Difference PCA (K=3)",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

pdf(file.path(FIG_DIR, "paired_PCA.pdf"), width = 14, height = 6)
print(ggarrange(p1, p2, ncol = 2))
dev.off()
cat("  PCA plot saved.\n")

# ---- Multi-omics perturbation profile heatmap ----
# Select top discriminating features between K=2 subtypes
cl2 <- cc_paired[[2]]$consensusClass
st2 <- paste0("DS", cl2)
names(st2) <- names(cl2)

# Wilcoxon test on each feature
pvals_diff <- apply(concat_diff, 1, function(x) {
  tryCatch(wilcox.test(x[st2 == "DS1"], x[st2 == "DS2"])$p.value, error = function(e) NA)
})

fcs_diff <- rowMeans(concat_diff[, st2 == "DS1", drop = FALSE]) -
            rowMeans(concat_diff[, st2 == "DS2", drop = FALSE])

top_disc <- names(sort(pvals_diff))[1:min(60, sum(!is.na(pvals_diff)))]

# Order patients by subtype
samp_order <- c(
  sort(names(st2)[st2 == "DS1"]),
  sort(names(st2)[st2 == "DS2"])
)

mat_plot <- concat_diff[top_disc, samp_order]
mat_plot[mat_plot > 3] <- 3
mat_plot[mat_plot < -3] <- -3

# Omics annotation
feat_omics <- ifelse(grepl("^TC_", top_disc), "Transcriptomics",
               ifelse(grepl("^PR_", top_disc), "Proteomics", "Metabolomics"))
anno_row <- data.frame(Omics = feat_omics)
rownames(anno_row) <- top_disc

anno_col <- data.frame(Subtype = st2[samp_order])
rownames(anno_col) <- samp_order

pdf(file.path(FIG_DIR, "paired_perturbation_heatmap.pdf"), width = 10, height = 14)
pheatmap(mat_plot,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(
           Subtype = c("DS1" = "#3C5488", "DS2" = "#E64B35"),
           Omics = c("Transcriptomics" = "#4DBBD5", "Proteomics" = "#E64B35", "Metabolomics" = "#00A087")
         ),
         cluster_cols = FALSE, cluster_rows = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6, fontsize_col = 9,
         main = "Disease Perturbation Profiles (Adjacent - Normal)",
         breaks = seq(-3, 3, length.out = 101))
dev.off()
cat("  Perturbation heatmap saved.\n")

# ---- Clinical boxplots for best K ----
best_k <- ifelse(best_k_pac <= 3, best_k_pac, 2)
sub_df <- read.csv(file.path(OUT_DIR, sprintf("paired_subtype_K%d.csv", best_k)), stringsAsFactors = FALSE)
clin_merged <- merge(clinical, sub_df[, c("patient_id", "subtype")], by = "patient_id", all.x = TRUE)
clin_typed <- clin_merged[!is.na(clin_merged$subtype), ]

cont_vars <- c("age", "BMI", "lesion_size", "disease_duration",
               "ALT", "AST", "ALP", "GGT", "WBC", "NLR",
               "CRP", "albumin")
cont_vars <- intersect(cont_vars, colnames(clin_typed))

plot_list <- list()
for (v in cont_vars) {
  df_p <- clin_typed[, c(v, "subtype")]
  df_p[[v]] <- as.numeric(df_p[[v]])
  df_p <- df_p[!is.na(df_p[[v]]), ]
  if (nrow(df_p) < 4) next
  
  p <- ggplot(df_p, aes(x = subtype, y = .data[[v]], fill = subtype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("DS1" = "#3C5488", "DS2" = "#E64B35", "DS3" = "#00A087")) +
    labs(title = v, x = "", y = v) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
  
  tryCatch({
    if (length(unique(df_p$subtype)) == 2) {
      pval <- wilcox.test(as.numeric(df_p[[v]]) ~ df_p$subtype)$p.value
    } else {
      pval <- kruskal.test(as.numeric(df_p[[v]]) ~ factor(df_p$subtype))$p.value
    }
    p <- p + annotate("text", x = 1.5, y = max(df_p[[v]], na.rm = TRUE),
                      label = sprintf("P=%.3f", pval), size = 3)
  }, error = function(e) NULL)
  
  plot_list[[v]] <- p
}

if (length(plot_list) >= 4) {
  pdf(file.path(FIG_DIR, sprintf("paired_clinical_boxplots_K%d.pdf", best_k)), width = 16, height = 12)
  print(ggarrange(plotlist = plot_list, ncol = 4, nrow = 3))
  dev.off()
  cat("  Clinical boxplots saved.\n")
}


# ============================================================================
# 6. COMPARE ORIGINAL vs PAIRED SUBTYPING
# ============================================================================
cat("\n>>> 6. Comparing original vs paired-diff subtypes\n")

# Load original subtypes
orig_sub <- read.csv(file.path(P6_DIR, "subtype_K2.csv"), stringsAsFactors = FALSE)
orig_map <- setNames(orig_sub$subtype, orig_sub$patient_id)

# Paired subtypes
paired_sub <- read.csv(file.path(OUT_DIR, "paired_subtype_K2.csv"), stringsAsFactors = FALSE)
paired_map <- setNames(paired_sub$subtype, paired_sub$patient_id)

# Compare
shared_patients <- intersect(names(orig_map), names(paired_map))
if (length(shared_patients) > 0) {
  comp_df <- data.frame(
    patient = shared_patients,
    original = orig_map[shared_patients],
    paired = paired_map[shared_patients],
    stringsAsFactors = FALSE
  )
  comp_df$concordant <- comp_df$original == gsub("DS", "CS", comp_df$paired) |
                         comp_df$original == gsub("DS(\\d)", "CS\\1", comp_df$paired)
  
  # ARI between the two methods
  if (requireNamespace("mclust", quietly = TRUE)) {
    library(mclust)
    ari <- adjustedRandIndex(
      as.integer(gsub("CS", "", orig_map[shared_patients])),
      as.integer(gsub("DS", "", paired_map[shared_patients]))
    )
    cat(sprintf("  ARI (Original vs Paired): %.3f\n", ari))
  }
  
  cat("  Patient-level comparison:\n")
  print(table(Original = orig_map[shared_patients], Paired = paired_map[shared_patients]))
  
  write.csv(comp_df, file.path(OUT_DIR, "original_vs_paired_comparison.csv"), row.names = FALSE)
}


# ============================================================================
# 7. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Enhancement 1 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Paired patients: %d\n", length(paired_pids)))
cat(sprintf("Diff matrix: %d features x %d patients\n", nrow(concat_diff), ncol(concat_diff)))
cat(sprintf("\nOptimal K: PAC=%d, Silhouette=%d\n", best_k_pac, best_k_sil))

for (k in 2:3) {
  cl <- cc_paired[[k]]$consensusClass
  cat(sprintf("\nK=%d: %s\n", k,
              paste(sprintf("DS%d=%d", 1:k, table(factor(cl, levels = 1:k))), collapse=", ")))
  
  assoc <- read.csv(file.path(OUT_DIR, sprintf("paired_clinical_assoc_K%d.csv", k)), stringsAsFactors = FALSE)
  cat(sprintf("  Clinical associations: %d sig (P<0.05), %d trends (P<0.1)\n",
              sum(assoc$pvalue < 0.05), sum(assoc$pvalue < 0.1)))
}

cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Enhancement 1 COMPLETE.\n")
