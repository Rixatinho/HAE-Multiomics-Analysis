#!/usr/bin/env Rscript
# ============================================================================
# Enhancement 3: TF Activity (DoRothEA via decoupleR) + Pathway Activity (PROGENy)
# ============================================================================

suppressPackageStartupMessages({
  library(decoupleR)
  library(OmnipathR)
  library(progeny)
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
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhancement3_tf_pathway")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Enhancement 3: TF & Pathway Activity\n")
cat("========================================\n\n")

# ---- Load data ----
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

# Groups
group_vec <- ifelse(grepl("^Normal", colnames(tc_mat)), "Normal", "Adjacent")
names(group_vec) <- colnames(tc_mat)

# Subtypes
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
# 1. TF ACTIVITY INFERENCE (DoRothEA via decoupleR)
# ============================================================================
cat("\n>>> 1. TF activity inference (DoRothEA)\n")

# Get DoRothEA regulons
cat("  Loading DoRothEA regulons...\n")
dorothea_net <- tryCatch({
  get_dorothea(organism = "human", levels = c("A", "B", "C"))
}, error = function(e) {
  cat("  Error loading DoRothEA:", e$message, "\n")
  cat("  Trying alternative...\n")
  get_collectri(organism = "human")
})

cat(sprintf("  Regulon network: %d interactions, %d TFs\n",
            nrow(dorothea_net), length(unique(dorothea_net$source))))

# Run ULM (Univariate Linear Model) for TF activity
cat("  Running TF activity inference (ULM)...\n")
tf_activity <- run_ulm(mat = tc_mat, net = dorothea_net, .source = "source",
                        .target = "target", .mor = "mor", minsize = 5)

# Pivot to matrix
tf_mat <- tf_activity %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  as.data.frame()

rownames(tf_mat) <- tf_mat$source
tf_mat$source <- NULL
tf_mat <- as.matrix(tf_mat)

cat(sprintf("  TF activity matrix: %d TFs x %d samples\n", nrow(tf_mat), ncol(tf_mat)))

# ---- Test TF activity: Adjacent vs Normal (paired Wilcoxon) ----
adj_samps <- colnames(tf_mat)[grepl("^Adjacent", colnames(tf_mat))]
nor_samps <- colnames(tf_mat)[grepl("^Normal", colnames(tf_mat))]

# Build paired sample vectors
extract_id <- function(x) as.integer(gsub("Normal|Adjacent", "", x))
adj_ids <- extract_id(adj_samps)
nor_ids <- extract_id(nor_samps)
paired_ids <- intersect(adj_ids, nor_ids)
adj_paired <- paste0("Adjacent", paired_ids)
nor_paired <- paste0("Normal", paired_ids)

tf_diff <- data.frame()
for (tf in rownames(tf_mat)) {
  tryCatch({
    wt <- wilcox.test(tf_mat[tf, adj_paired], tf_mat[tf, nor_paired], paired = TRUE)
    fc <- mean(tf_mat[tf, adj_samps]) - mean(tf_mat[tf, nor_samps])
    tf_diff <- rbind(tf_diff, data.frame(
      TF = tf, mean_Adjacent = mean(tf_mat[tf, adj_samps]),
      mean_Normal = mean(tf_mat[tf, nor_samps]),
      diff = fc, pvalue = wt$p.value, stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}
tf_diff$padj <- p.adjust(tf_diff$pvalue, method = "BH")
tf_diff <- tf_diff[order(tf_diff$pvalue), ]

cat(sprintf("  TFs sig (P<0.05): %d, (padj<0.05): %d\n",
            sum(tf_diff$pvalue < 0.05), sum(tf_diff$padj < 0.05)))

write.csv(tf_diff, file.path(OUT_DIR, "TF_activity_diff.csv"), row.names = FALSE)

# ---- Test TF: DS1 vs DS2 ----
ds1_samps <- intersect(adj_samps, names(sample_subtype)[sample_subtype == "DS1"])
ds2_samps <- intersect(adj_samps, names(sample_subtype)[sample_subtype == "DS2"])

tf_subtype <- data.frame()
for (tf in rownames(tf_mat)) {
  if (length(ds1_samps) < 2 || length(ds2_samps) < 2) next
  tryCatch({
    wt <- wilcox.test(tf_mat[tf, ds1_samps], tf_mat[tf, ds2_samps])
    tf_subtype <- rbind(tf_subtype, data.frame(
      TF = tf, mean_DS1 = mean(tf_mat[tf, ds1_samps]),
      mean_DS2 = mean(tf_mat[tf, ds2_samps]),
      diff = mean(tf_mat[tf, ds1_samps]) - mean(tf_mat[tf, ds2_samps]),
      pvalue = wt$p.value, stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}
tf_subtype$padj <- p.adjust(tf_subtype$pvalue, method = "BH")
tf_subtype <- tf_subtype[order(tf_subtype$pvalue), ]

cat(sprintf("  TFs diff between subtypes (P<0.05): %d\n", sum(tf_subtype$pvalue < 0.05)))
write.csv(tf_subtype, file.path(OUT_DIR, "TF_activity_subtype.csv"), row.names = FALSE)

# Top TFs
cat("\n  Top TFs (Adjacent vs Normal):\n")
for (i in 1:min(10, nrow(tf_diff))) {
  cat(sprintf("    %s: diff=%.2f, P=%.4f\n", tf_diff$TF[i], tf_diff$diff[i], tf_diff$pvalue[i]))
}

# ---- TF activity heatmap ----
top_tfs <- head(tf_diff$TF, 30)
tf_plot <- tf_mat[top_tfs, ]

samp_order <- c(sort(nor_samps), sort(adj_samps))
samp_order <- intersect(samp_order, colnames(tf_plot))
tf_plot <- tf_plot[, samp_order]

anno_col <- data.frame(
  Group = ifelse(grepl("^Normal", samp_order), "Normal", "Adjacent")
)
rownames(anno_col) <- samp_order

# Add subtype for Adjacent
for (s in samp_order) {
  if (!is.na(sample_subtype[s])) anno_col[s, "Subtype"] <- sample_subtype[s]
}

pdf(file.path(FIG_DIR, "TF_activity_heatmap.pdf"), width = 12, height = 10)
pheatmap(tf_plot,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         annotation_col = anno_col,
         annotation_colors = list(
           Group = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35"),
           Subtype = c("DS1" = "#3C5488", "DS2" = "#F39B7F")
         ),
         cluster_cols = FALSE, cluster_rows = TRUE,
         show_colnames = TRUE, fontsize_row = 8, fontsize_col = 7,
         main = "Top 30 Transcription Factor Activities (DoRothEA/ULM)",
         scale = "row")
dev.off()
cat("  TF heatmap saved.\n")


# ============================================================================
# 2. PATHWAY ACTIVITY (PROGENy - direct progeny package)
# ============================================================================
cat("\n>>> 2. Pathway activity inference (PROGENy)\n")

# Use progeny() function directly (avoids OmnipathR network + decoupleR run_mlm issues)
cat("  Running PROGENy pathway scoring (top 500 genes, 10000 permutations)...\n")
pw_scores <- progeny(tc_mat, scale = FALSE, organism = "Human", top = 500,
                     perm = 10000, verbose = TRUE)

# progeny() returns samples x pathways; transpose to pathways x samples
pw_mat <- t(pw_scores)
# Fix pathway names (JAK.STAT -> JAK-STAT)
rownames(pw_mat) <- gsub("\\.", "-", rownames(pw_mat))

# Remove pathways with zero variance (constant across all samples)
pw_var <- apply(pw_mat, 1, var, na.rm = TRUE)
pw_mat <- pw_mat[pw_var > 1e-10, , drop = FALSE]

cat(sprintf("  Pathway activity matrix: %d pathways x %d samples\n", nrow(pw_mat), ncol(pw_mat)))

# ---- Test pathway activity: Adjacent vs Normal (paired Wilcoxon) ----
pw_diff <- data.frame()
for (pw in rownames(pw_mat)) {
  tryCatch({
    wt <- wilcox.test(pw_mat[pw, adj_paired], pw_mat[pw, nor_paired], paired = TRUE)
    fc <- mean(pw_mat[pw, adj_samps]) - mean(pw_mat[pw, nor_samps])
    pw_diff <- rbind(pw_diff, data.frame(
      pathway = pw, mean_Adjacent = mean(pw_mat[pw, adj_samps]),
      mean_Normal = mean(pw_mat[pw, nor_samps]),
      diff = fc, pvalue = wt$p.value, stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}
pw_diff$padj <- p.adjust(pw_diff$pvalue, method = "BH")
pw_diff <- pw_diff[order(pw_diff$pvalue), ]

cat(sprintf("  Pathways sig (P<0.05): %d / %d\n", sum(pw_diff$pvalue < 0.05, na.rm = TRUE), nrow(pw_diff)))
write.csv(pw_diff, file.path(OUT_DIR, "PROGENy_pathway_diff.csv"), row.names = FALSE)

cat("\n  PROGENy pathway differences:\n")
for (i in 1:nrow(pw_diff)) {
  sig_mark <- ifelse(pw_diff$pvalue[i] < 0.05, "*", "")
  cat(sprintf("    %s: diff=%.2f, P=%.4f %s\n",
              pw_diff$pathway[i], pw_diff$diff[i], pw_diff$pvalue[i], sig_mark))
}

# ---- Test by subtype ----
pw_subtype <- data.frame()
for (pw in rownames(pw_mat)) {
  if (length(ds1_samps) < 2 || length(ds2_samps) < 2) next
  tryCatch({
    wt <- wilcox.test(pw_mat[pw, ds1_samps], pw_mat[pw, ds2_samps])
    pw_subtype <- rbind(pw_subtype, data.frame(
      pathway = pw, mean_DS1 = mean(pw_mat[pw, ds1_samps]),
      mean_DS2 = mean(pw_mat[pw, ds2_samps]),
      diff = mean(pw_mat[pw, ds1_samps]) - mean(pw_mat[pw, ds2_samps]),
      pvalue = wt$p.value, stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}
pw_subtype$padj <- p.adjust(pw_subtype$pvalue, method = "BH")
pw_subtype <- pw_subtype[order(pw_subtype$pvalue), ]
write.csv(pw_subtype, file.path(OUT_DIR, "PROGENy_pathway_subtype.csv"), row.names = FALSE)

cat(sprintf("  Pathways diff between subtypes (P<0.05): %d\n", sum(pw_subtype$pvalue < 0.05)))

# ---- Pathway activity heatmap ----
pw_plot <- pw_mat[, samp_order]

pdf(file.path(FIG_DIR, "PROGENy_heatmap.pdf"), width = 12, height = 6)
pheatmap(pw_plot,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         annotation_col = anno_col,
         annotation_colors = list(
           Group = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35"),
           Subtype = c("DS1" = "#3C5488", "DS2" = "#F39B7F")
         ),
         cluster_cols = FALSE, cluster_rows = TRUE,
         show_colnames = TRUE, fontsize_row = 10, fontsize_col = 7,
         main = "PROGENy Pathway Activities",
         scale = "row")
dev.off()
cat("  PROGENy heatmap saved.\n")

# ---- Pathway boxplots ----
pw_long <- data.frame()
for (pw in rownames(pw_mat)) {
  pw_long <- rbind(pw_long, data.frame(
    pathway = pw,
    score = c(pw_mat[pw, nor_samps], pw_mat[pw, adj_samps]),
    group = c(rep("Normal", length(nor_samps)), rep("Adjacent", length(adj_samps))),
    stringsAsFactors = FALSE
  ))
}

p_pw <- ggplot(pw_long, aes(x = pathway, y = score, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35")) +
  coord_flip() +
  labs(title = "PROGENy Pathway Activities: Adjacent vs Normal",
       x = "", y = "Pathway Activity Score", fill = "Group") +
  theme_bw(base_size = 10)

ggsave(file.path(FIG_DIR, "PROGENy_boxplot.pdf"), p_pw, width = 10, height = 8)
cat("  PROGENy boxplot saved.\n")


# ============================================================================
# 3. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Enhancement 3 SUMMARY\n")
cat("========================================\n")
cat(sprintf("TF activity: %d TFs scored, %d sig (P<0.05), %d sig subtypes\n",
            nrow(tf_diff), sum(tf_diff$pvalue < 0.05), sum(tf_subtype$pvalue < 0.05)))
cat(sprintf("PROGENy: %d pathways scored, %d sig (P<0.05), %d sig subtypes\n",
            nrow(pw_diff), sum(pw_diff$pvalue < 0.05, na.rm = TRUE), sum(pw_subtype$pvalue < 0.05, na.rm = TRUE)))
cat(sprintf("Output: %s\n", OUT_DIR))
cat("Enhancement 3 COMPLETE.\n")
