#!/usr/bin/env Rscript
# ============================================================================
# Phase 3: Metabolic Reprogramming Deep-Dive
# HAE Multi-omics Integration Study
# ============================================================================
# Objectives:
#   1. Metabolite class enrichment (over-representation of classes in DEMs)
#   2. KEGG metabolic pathway enrichment for metabolomics
#   3. Lipid metabolism profiling (subclass-level)
#   4. Amino acid & energy metabolism analysis
#   5. Cross-omics metabolic enzyme-metabolite integration
#   6. Metabolic network visualization
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
  library(dplyr)
  library(tidyr)
  library(VennDiagram)
  library(grid)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
OUT_DIR  <- file.path(PROJECT, "analysis/results/phase3_metabolic")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 3: Metabolic Reprogramming\n")
cat("========================================\n\n")

# ============================================================================
# 0. LOAD DATA
# ============================================================================
cat(">>> 0. Loading data\n")

# Metabolomics differential results
dems_all <- read.csv(file.path(DIFF_DIR, "DEMs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# Metabolomics annotation (rich: ClassI/II/III, KEGG_ID, HMDB_ID, etc.)
annot <- read.csv(file.path(PROC_DIR, "metabolomics_annotation.csv"), check.names = FALSE)

# Merge annotation with differential results
dems <- merge(dems_all, annot, by = "Compound_ID", all.x = TRUE)
cat(sprintf("  Total metabolites: %d\n", nrow(dems)))

# Define significance using BH-adjusted P-values (corrected from P.Value)
# With 3299 metabolites, unadjusted P would yield ~165 false positives at alpha=0.05
dems$sig <- dems$adj.P.Val < 0.05 & abs(dems$logFC) > 0.585
dems$direction <- ifelse(dems$logFC > 0, "Up", "Down")
dems$direction[!dems$sig] <- "NS"

n_up   <- sum(dems$direction == "Up")
n_down <- sum(dems$direction == "Down")
cat(sprintf("  Significant DEMs: %d up, %d down, %d total\n", n_up, n_down, n_up + n_down))

# Load transcriptomics and proteomics for metabolic enzyme integration
degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)


# ============================================================================
# 1. METABOLITE CLASS ENRICHMENT (Fisher's exact test)
# ============================================================================
cat("\n>>> 1. Metabolite Class Enrichment Analysis\n")

run_class_enrichment <- function(dems, class_col, label, min_class_size = 5) {
  # Clean class column
  dems[[class_col]][dems[[class_col]] == "-" | dems[[class_col]] == "" | is.na(dems[[class_col]])] <- "Unclassified"
  
  all_classes <- unique(dems[[class_col]])
  all_classes <- all_classes[all_classes != "Unclassified"]
  
  results <- data.frame()
  n_total <- nrow(dems)
  n_sig <- sum(dems$sig)
  
  for (cls in all_classes) {
    in_class <- dems[[class_col]] == cls
    n_class <- sum(in_class)
    if (n_class < min_class_size) next
    
    n_sig_in_class <- sum(dems$sig & in_class)
    n_sig_not_class <- n_sig - n_sig_in_class
    n_nosig_in_class <- n_class - n_sig_in_class
    n_nosig_not_class <- n_total - n_sig - n_nosig_in_class
    
    mat <- matrix(c(n_sig_in_class, n_sig_not_class,
                     n_nosig_in_class, n_nosig_not_class),
                   nrow = 2)
    test <- fisher.test(mat, alternative = "two.sided")
    
    # Direction: mean logFC of sig metabolites in this class
    mean_fc <- mean(dems$logFC[in_class & dems$sig], na.rm = TRUE)
    n_up_cls <- sum(dems$direction == "Up" & in_class)
    n_down_cls <- sum(dems$direction == "Down" & in_class)
    
    results <- rbind(results, data.frame(
      class = cls,
      n_total_in_class = n_class,
      n_sig_in_class = n_sig_in_class,
      expected = round(n_class * n_sig / n_total, 1),
      fold_enrichment = round((n_sig_in_class / n_class) / (n_sig / n_total), 2),
      n_up = n_up_cls,
      n_down = n_down_cls,
      mean_logFC_sig = round(mean_fc, 3),
      pvalue = test$p.value,
      OR = round(test$estimate, 2),
      stringsAsFactors = FALSE
    ))
  }
  
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(results$pvalue), ]
  
  cat(sprintf("  %s: %d classes tested, %d enriched (padj<0.1)\n",
              label, nrow(results), sum(results$padj < 0.1)))
  return(results)
}

# ClassI enrichment
enrich_c1 <- run_class_enrichment(dems, "ClassI", "ClassI")
write.csv(enrich_c1, file.path(OUT_DIR, "metabolite_ClassI_enrichment.csv"), row.names = FALSE)

# ClassII enrichment
enrich_c2 <- run_class_enrichment(dems, "ClassII", "ClassII", min_class_size = 3)
write.csv(enrich_c2, file.path(OUT_DIR, "metabolite_ClassII_enrichment.csv"), row.names = FALSE)

# Visualization: ClassI enrichment barplot
enrich_c1_plot <- enrich_c1[enrich_c1$n_sig_in_class >= 3, ]
enrich_c1_plot <- head(enrich_c1_plot[order(-enrich_c1_plot$n_sig_in_class), ], 15)
enrich_c1_plot$class_short <- substr(enrich_c1_plot$class, 1, 40)
enrich_c1_plot$sig_label <- paste0(enrich_c1_plot$n_up, " up / ", enrich_c1_plot$n_down, " down")

p_class <- ggplot(enrich_c1_plot, aes(n_sig_in_class, reorder(class_short, n_sig_in_class))) +
  geom_col(aes(fill = fold_enrichment), width = 0.7) +
  geom_text(aes(label = sig_label), hjust = -0.1, size = 3) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "grey90", high = "#E64B35", midpoint = 1) +
  labs(title = "Metabolite Class Enrichment in DEMs",
       x = "Number of significant metabolites", y = NULL,
       fill = "Fold\nEnrichment") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 9)) +
  xlim(0, max(enrich_c1_plot$n_sig_in_class) * 1.4)
ggsave(file.path(FIG_DIR, "metabolite_ClassI_enrichment.pdf"), p_class, width = 10, height = 6)
cat("  ClassI enrichment plot saved.\n")


# ============================================================================
# 2. KEGG METABOLIC PATHWAY ENRICHMENT
# ============================================================================
cat("\n>>> 2. KEGG Metabolic Pathway Enrichment\n")

# Parse KEGG_MapID (semicolon-separated)
dems$kegg_maps <- dems$KEGG_MapID
dems$has_kegg <- !is.na(dems$kegg_maps) & dems$kegg_maps != "-" & dems$kegg_maps != ""

cat(sprintf("  Metabolites with KEGG pathway mapping: %d / %d\n",
            sum(dems$has_kegg), nrow(dems)))

# Build pathway-metabolite mapping
pw_met <- data.frame()
for (i in which(dems$has_kegg)) {
  maps <- trimws(unlist(strsplit(dems$kegg_maps[i], ";")))
  if (length(maps) > 0) {
    pw_met <- rbind(pw_met, data.frame(
      pathway = maps,
      Compound_ID = dems$Compound_ID[i],
      sig = dems$sig[i],
      logFC = dems$logFC[i],
      stringsAsFactors = FALSE
    ))
  }
}

# Fisher's exact for each KEGG pathway
pw_counts <- pw_met %>%
  group_by(pathway) %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig),
    n_up = sum(sig & logFC > 0),
    n_down = sum(sig & logFC < 0),
    mean_logFC = mean(logFC[sig], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_total >= 3) %>%
  arrange(desc(n_sig))

n_total_mapped <- sum(dems$has_kegg)
n_sig_mapped <- sum(dems$has_kegg & dems$sig)

pw_counts$pvalue <- sapply(1:nrow(pw_counts), function(i) {
  mat <- matrix(c(
    pw_counts$n_sig[i],
    n_sig_mapped - pw_counts$n_sig[i],
    pw_counts$n_total[i] - pw_counts$n_sig[i],
    n_total_mapped - n_sig_mapped - (pw_counts$n_total[i] - pw_counts$n_sig[i])
  ), nrow = 2)
  mat[mat < 0] <- 0
  fisher.test(mat, alternative = "greater")$p.value
})

pw_counts$padj <- p.adjust(pw_counts$pvalue, method = "BH")

# Add pathway names via KEGG API
cat("  Fetching KEGG pathway names...\n")
tryCatch({
  kegg_info <- clusterProfiler::download_KEGG("hsa")
  pw_names <- setNames(kegg_info$KEGGPATHID2NAME[,2], 
                        paste0("map", gsub("hsa", "", kegg_info$KEGGPATHID2NAME[,1])))
  pw_counts$pathway_name <- pw_names[pw_counts$pathway]
  pw_counts$pathway_name[is.na(pw_counts$pathway_name)] <- pw_counts$pathway[is.na(pw_counts$pathway_name)]
}, error = function(e) {
  cat("  KEGG name fetch failed, using pathway IDs\n")
  pw_counts$pathway_name <<- pw_counts$pathway
})

pw_counts <- pw_counts[order(pw_counts$pvalue), ]
write.csv(pw_counts, file.path(OUT_DIR, "KEGG_metabolic_pathway_enrichment.csv"), row.names = FALSE)

cat(sprintf("  KEGG pathways tested: %d\n", nrow(pw_counts)))
cat(sprintf("  Significant (padj<0.1): %d\n", sum(pw_counts$padj < 0.1, na.rm = TRUE)))

# Top pathways plot
pw_plot <- head(pw_counts[pw_counts$n_sig >= 2, ], 20)
if (nrow(pw_plot) > 0) {
  pw_plot$name_short <- substr(pw_plot$pathway_name, 1, 45)
  
  p_kegg <- ggplot(pw_plot, aes(-log10(pvalue), reorder(name_short, -log10(pvalue)))) +
    geom_point(aes(size = n_sig, color = n_up / (n_up + n_down)), alpha = 0.8) +
    scale_color_gradient2(low = "#4DBBD5", mid = "grey80", high = "#E64B35",
                          midpoint = 0.5, limits = c(0, 1)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(title = "KEGG Metabolic Pathway Enrichment (Metabolomics)",
         x = "-log10(p-value)", y = NULL,
         size = "DEMs in\npathway", color = "Fraction\nUpregulated") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(FIG_DIR, "KEGG_metabolic_pathway_enrichment.pdf"), p_kegg, width = 10, height = 7)
  cat("  KEGG pathway enrichment plot saved.\n")
}


# ============================================================================
# 3. LIPID METABOLISM PROFILING
# ============================================================================
cat("\n>>> 3. Lipid Metabolism Profiling\n")

lipids <- dems[grepl("Lipid|lipid", dems$ClassI), ]
cat(sprintf("  Total lipids: %d, significant: %d\n", nrow(lipids), sum(lipids$sig)))

# Lipid subclass distribution
lipid_sub <- lipids %>%
  group_by(ClassII) %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig),
    n_up = sum(direction == "Up"),
    n_down = sum(direction == "Down"),
    mean_logFC = mean(logFC, na.rm = TRUE),
    mean_logFC_sig = mean(logFC[sig], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_total >= 3) %>%
  arrange(desc(n_sig))

write.csv(lipid_sub, file.path(OUT_DIR, "lipid_subclass_summary.csv"), row.names = FALSE)

# Lipid subclass barplot (stacked up/down)
lipid_bar <- lipids[lipids$sig, ] %>%
  group_by(ClassII, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(!is.na(ClassII), ClassII != "-")

if (nrow(lipid_bar) > 0) {
  p_lipid <- ggplot(lipid_bar, aes(n, reorder(ClassII, n), fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5")) +
    labs(title = "Differentially Abundant Lipid Subclasses",
         x = "Number of significant metabolites", y = NULL, fill = "Direction") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG_DIR, "lipid_subclass_barplot.pdf"), p_lipid, width = 8, height = 5)
  cat("  Lipid subclass barplot saved.\n")
}

# Fatty Acyl detailed breakdown (ClassIII)
fa <- lipids[grepl("Fatty Acyl", lipids$ClassII), ]
if (nrow(fa) > 0) {
  fa_sub <- fa %>%
    group_by(ClassIII) %>%
    summarise(
      n_total = n(), n_sig = sum(sig),
      n_up = sum(direction == "Up"), n_down = sum(direction == "Down"),
      .groups = "drop"
    ) %>%
    filter(n_total >= 2) %>%
    arrange(desc(n_sig))
  write.csv(fa_sub, file.path(OUT_DIR, "fatty_acyl_ClassIII_summary.csv"), row.names = FALSE)
  cat(sprintf("  Fatty Acyls: %d total, %d significant, %d subclasses\n",
              nrow(fa), sum(fa$sig), nrow(fa_sub)))
}


# ============================================================================
# 4. AMINO ACID METABOLISM
# ============================================================================
cat("\n>>> 4. Amino Acid Metabolism Analysis\n")

aa <- dems[grepl("Amino acid|amino acid|Peptide|peptide", dems$ClassII, ignore.case = TRUE) |
           grepl("amino acid", dems$ClassI, ignore.case = TRUE) |
           grepl("Organic acid", dems$ClassI, ignore.case = TRUE), ]
cat(sprintf("  Amino acids & organic acids: %d total, %d significant\n",
            nrow(aa), sum(aa$sig)))

aa_sig <- aa[aa$sig, c("Compound_ID", "Name", "ClassII", "ClassIII", "logFC", "P.Value", "KEGG_ID")]
aa_sig <- aa_sig[order(aa_sig$P.Value), ]
write.csv(aa_sig, file.path(OUT_DIR, "amino_acid_organic_acid_DEMs.csv"), row.names = FALSE)


# ============================================================================
# 5. CROSS-OMICS METABOLIC ENZYME-METABOLITE INTEGRATION
# ============================================================================
cat("\n>>> 5. Cross-omics Metabolic Enzyme Integration\n")

# Define key metabolic gene sets (manually curated for HAE relevance)
metabolic_genes <- list(
  # Lipid metabolism
  fatty_acid_synthesis = c("FASN", "ACACA", "ACACB", "SCD", "ELOVL1", "ELOVL2", 
                            "ELOVL5", "ELOVL6", "FADS1", "FADS2"),
  fatty_acid_oxidation = c("CPT1A", "CPT1B", "CPT2", "ACADL", "ACADM", "ACADS",
                            "HADHA", "HADHB", "ACOX1", "ACOX2"),
  cholesterol_metabolism = c("HMGCR", "HMGCS1", "HMGCS2", "MVK", "FDPS", "SQLE",
                              "LSS", "CYP51A1", "DHCR7", "DHCR24"),
  bile_acid_metabolism = c("CYP7A1", "CYP8B1", "CYP27A1", "CYP7B1", "SLC10A1",
                            "SLC10A2", "ABCB11", "NR1H4", "BAAT", "SLC27A5"),
  # Amino acid metabolism
  glutamine_metabolism = c("GLS", "GLS2", "GLUL", "GLUD1", "GLUD2", "GOT1", "GOT2",
                            "GPT", "GPT2", "SLC1A5"),
  tryptophan_metabolism = c("TDO2", "IDO1", "IDO2", "KMO", "KYNU", "HAAO",
                             "AFMID", "TPH1", "TPH2", "DDC"),
  # Energy metabolism
  glycolysis = c("HK1", "HK2", "GPI", "PFKL", "PFKM", "ALDOA", "ALDOB",
                  "GAPDH", "PGK1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB"),
  tca_cycle = c("CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "OGDH",
                 "SUCLA2", "SDHA", "SDHB", "FH", "MDH1", "MDH2"),
  oxidative_phosphorylation = c("NDUFA1", "NDUFB1", "NDUFS1", "SDHA", "UQCRC1",
                                 "COX5A", "COX7A1", "ATP5F1A", "ATP5F1B", "ATP5MC1"),
  # Immune-metabolic
  arginine_metabolism = c("ARG1", "ARG2", "NOS1", "NOS2", "NOS3", "ASS1", "ASL",
                           "OTC", "CPS1", "ODC1")
)

# Extract metabolic enzyme expression from transcriptomics and proteomics
extract_enzyme_data <- function(de_results, gene_sets, label) {
  result <- data.frame()
  for (pathway in names(gene_sets)) {
    genes <- gene_sets[[pathway]]
    matched <- de_results[de_results$gene_name %in% genes, ]
    if (nrow(matched) > 0) {
      matched$pathway <- pathway
      matched$omics <- label
      result <- rbind(result, matched[, c("gene_name", "logFC", "t", "P.Value", 
                                            "adj.P.Val", "pathway", "omics")])
    }
  }
  return(result)
}

enzyme_tc <- extract_enzyme_data(degs, metabolic_genes, "Transcriptomics")
enzyme_pr <- extract_enzyme_data(deps, metabolic_genes, "Proteomics")
enzyme_all <- rbind(enzyme_tc, enzyme_pr)

cat(sprintf("  Metabolic enzymes found: TC=%d, PR=%d\n", nrow(enzyme_tc), nrow(enzyme_pr)))

write.csv(enzyme_all, file.path(OUT_DIR, "metabolic_enzyme_expression.csv"), row.names = FALSE)

# Heatmap of metabolic enzyme expression across pathways
if (nrow(enzyme_tc) > 0 && nrow(enzyme_pr) > 0) {
  # Pivot to matrix: genes x (TC_logFC, PR_logFC)
  tc_fc <- setNames(enzyme_tc$logFC, enzyme_tc$gene_name)
  pr_fc <- setNames(enzyme_pr$logFC, enzyme_pr$gene_name)
  
  shared_genes <- intersect(names(tc_fc), names(pr_fc))
  if (length(shared_genes) > 0) {
    mat <- cbind(
      TC_logFC = tc_fc[shared_genes],
      PR_logFC = pr_fc[shared_genes]
    )
    colnames(mat) <- c("mRNA logFC", "Protein logFC")
    
    # Annotation: pathway membership
    gene_pw <- enzyme_tc[enzyme_tc$gene_name %in% shared_genes, c("gene_name", "pathway")]
    gene_pw <- gene_pw[!duplicated(gene_pw$gene_name), ]
    rownames(gene_pw) <- gene_pw$gene_name
    
    ann_row <- data.frame(
      Pathway = gsub("_", " ", gene_pw[shared_genes, "pathway"]),
      row.names = shared_genes
    )
    
    # Significance markers
    tc_pval <- setNames(enzyme_tc$P.Value, enzyme_tc$gene_name)
    pr_pval <- setNames(enzyme_pr$P.Value, enzyme_pr$gene_name)
    
    # Limit to max 60 genes for readability
    if (nrow(mat) > 60) {
      # Prioritize significant genes
      combined_p <- pmin(tc_pval[shared_genes], pr_pval[shared_genes], na.rm = TRUE)
      top_idx <- order(combined_p)[1:60]
      mat <- mat[top_idx, ]
      ann_row <- ann_row[rownames(mat), , drop = FALSE]
    }
    
    pdf(file.path(FIG_DIR, "heatmap_metabolic_enzymes_cross_omics.pdf"),
        width = 8, height = max(8, nrow(mat) * 0.25 + 3))
    pheatmap(
      mat,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      annotation_row = ann_row,
      fontsize_row = 7,
      main = "Metabolic Enzyme Expression (Adjacent vs Normal)",
      breaks = seq(-max(abs(mat), na.rm=TRUE), max(abs(mat), na.rm=TRUE), length.out = 101)
    )
    dev.off()
    cat(sprintf("  Metabolic enzyme heatmap saved (%d genes).\n", nrow(mat)))
  }
}


# ============================================================================
# 6. PATHWAY-LEVEL METABOLIC SUMMARY
# ============================================================================
cat("\n>>> 6. Pathway-level Metabolic Summary\n")

# For each metabolic pathway, summarize enzyme + metabolite changes
pathway_summary <- data.frame()
for (pw in names(metabolic_genes)) {
  genes <- metabolic_genes[[pw]]
  
  # Transcriptomics
  tc_match <- degs[degs$gene_name %in% genes, ]
  tc_sig <- sum(tc_match$P.Value < 0.05 & abs(tc_match$logFC) > 0.585)
  tc_mean_fc <- mean(tc_match$logFC, na.rm = TRUE)
  
  # Proteomics
  pr_match <- deps[deps$gene_name %in% genes, ]
  pr_sig <- sum(pr_match$P.Value < 0.05 & abs(pr_match$logFC) > 0.585)
  pr_mean_fc <- mean(pr_match$logFC, na.rm = TRUE)
  
  pathway_summary <- rbind(pathway_summary, data.frame(
    pathway = pw,
    n_genes_defined = length(genes),
    TC_detected = nrow(tc_match),
    TC_sig = tc_sig,
    TC_mean_logFC = round(tc_mean_fc, 3),
    PR_detected = nrow(pr_match),
    PR_sig = pr_sig,
    PR_mean_logFC = round(pr_mean_fc, 3),
    stringsAsFactors = FALSE
  ))
}

write.csv(pathway_summary, file.path(OUT_DIR, "metabolic_pathway_summary.csv"), row.names = FALSE)
cat("  Pathway summary:\n")
for (i in 1:nrow(pathway_summary)) {
  cat(sprintf("    %-30s  TC: %d/%d sig (mean FC=%.2f)  PR: %d/%d sig (mean FC=%.2f)\n",
              pathway_summary$pathway[i],
              pathway_summary$TC_sig[i], pathway_summary$TC_detected[i],
              pathway_summary$TC_mean_logFC[i],
              pathway_summary$PR_sig[i], pathway_summary$PR_detected[i],
              pathway_summary$PR_mean_logFC[i]))
}

# Bubble plot: pathway activity across omics
pw_long <- pathway_summary %>%
  select(pathway, TC_sig, TC_mean_logFC, PR_sig, PR_mean_logFC) %>%
  pivot_longer(cols = c(TC_sig, PR_sig), names_to = "omics_sig", values_to = "n_sig") %>%
  pivot_longer(cols = c(TC_mean_logFC, PR_mean_logFC), names_to = "omics_fc", values_to = "mean_logFC") %>%
  filter((grepl("TC", omics_sig) & grepl("TC", omics_fc)) |
         (grepl("PR", omics_sig) & grepl("PR", omics_fc))) %>%
  mutate(omics = ifelse(grepl("TC", omics_sig), "Transcriptomics", "Proteomics"))

pw_long$pathway_label <- gsub("_", " ", pw_long$pathway)
pw_long$pathway_label <- tools::toTitleCase(pw_long$pathway_label)

p_bubble <- ggplot(pw_long, aes(omics, reorder(pathway_label, mean_logFC))) +
  geom_point(aes(size = n_sig + 0.5, color = mean_logFC), alpha = 0.8) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Metabolic Pathway Activity (Adjacent vs Normal)",
       x = NULL, y = NULL, size = "Sig.\ngenes", color = "Mean\nlogFC") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10, angle = 0))
ggsave(file.path(FIG_DIR, "metabolic_pathway_bubble.pdf"), p_bubble, width = 8, height = 7)
cat("  Pathway bubble plot saved.\n")


# ============================================================================
# 7. TOP DIFFERENTIAL METABOLITES HEATMAP
# ============================================================================
cat("\n>>> 7. Top Metabolite Heatmap\n")

# Load normalized metabolomics matrix
met_mat <- read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"), 
                    check.names = FALSE, row.names = 1)

# Get sample info
sample_corr <- read.csv(file.path(PROC_DIR, "sample_correspondence.csv"), check.names = FALSE)

# Top 40 DEMs by P-value
top_dems <- dems[dems$sig, ]
top_dems <- head(top_dems[order(top_dems$P.Value), ], 40)

# Extract expression matrix
met_ids <- top_dems$Compound_ID
met_sub <- met_mat[rownames(met_mat) %in% met_ids, ]

if (nrow(met_sub) > 2) {
  # Use metabolite names for row labels
  name_map <- setNames(top_dems$Name, top_dems$Compound_ID)
  rownames(met_sub) <- make.unique(substr(
    ifelse(name_map[rownames(met_sub)] != "" & !is.na(name_map[rownames(met_sub)]),
           name_map[rownames(met_sub)], rownames(met_sub)),
    1, 40))
  
  # Column annotation
  col_ann <- data.frame(
    Group = ifelse(grepl("_N$|_N[0-9]|Normal|normal", colnames(met_sub)), "Normal", "Adjacent"),
    row.names = colnames(met_sub)
  )
  # Try to assign groups based on sample naming
  col_ann$Group <- ifelse(grepl("_T$|_T[0-9]|Adjacent|adjacent|_P$|_P[0-9]", colnames(met_sub)), 
                           "Adjacent", "Normal")
  
  ann_colors <- list(Group = c("Normal" = "#4DBBD5", "Adjacent" = "#E64B35"))
  
  # Scale by row
  met_scaled <- t(scale(t(met_sub)))
  
  pdf(file.path(FIG_DIR, "heatmap_top40_DEMs.pdf"),
      width = 10, height = max(8, nrow(met_scaled) * 0.3 + 3))
  pheatmap(
    met_scaled,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_col = col_ann,
    annotation_colors = ann_colors,
    fontsize_row = 7,
    show_colnames = FALSE,
    main = "Top 40 Differentially Abundant Metabolites"
  )
  dev.off()
  cat("  Top 40 DEMs heatmap saved.\n")
}


# ============================================================================
# 8. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 3 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total metabolites: %d\n", nrow(dems)))
cat(sprintf("Significant DEMs: %d (up=%d, down=%d)\n", n_up + n_down, n_up, n_down))
cat(sprintf("\nClassI enrichment: %d classes tested\n", nrow(enrich_c1)))
cat(sprintf("  Top class: %s (%d sig DEMs)\n",
            enrich_c1$class[1], enrich_c1$n_sig_in_class[1]))
cat(sprintf("\nKEGG metabolic pathways: %d tested, %d sig (padj<0.1)\n",
            nrow(pw_counts), sum(pw_counts$padj < 0.1, na.rm = TRUE)))
cat(sprintf("\nLipid metabolism: %d lipids, %d significant\n",
            nrow(lipids), sum(lipids$sig)))
cat(sprintf("  Fatty Acyls: %d sig\n", sum(fa$sig)))
cat(sprintf("\nMetabolic enzymes detected across omics: %d shared\n",
            length(intersect(enzyme_tc$gene_name, enzyme_pr$gene_name))))
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 3 COMPLETE.\n")
