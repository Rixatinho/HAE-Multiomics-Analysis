#!/usr/bin/env Rscript
# =============================================================================
# HAE Multi-omics: Enhanced Cell Deconvolution & Cell Communication Analysis
# Module M3: High-resolution cellular microenvironment characterization
# =============================================================================

suppressPackageStartupMessages({
  library(TOAST)         # Cell type deconvolution
  library(Seurat)        # Single-cell analysis
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
  library(ggalluvial)
  library(igraph)
  library(Matrix)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# =============================================================================
# Setup and Data Loading
# =============================================================================
cat("=" %+% paste(rep("=", 79), collapse="") %+% "\n")
cat("HAE M3: Cell Deconvolution & Cell Communication Analysis\n")
cat("=" %+% paste(rep("=", 79), collapse="") %+% "\n")
`%+%` <- function(a, b) paste0(a, b)

base_dir <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
data_dir <- file.path(base_dir, "analysis/data/processed")
results_dir <- file.path(base_dir, "analysis/results/enhancement_deconvolution")
ref_dir <- file.path(base_dir, "analysis/data/reference/liver_sc_atlas")
prev_results <- file.path(base_dir, "analysis/results/enhancement7_immune_deconvolution")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Load bulk RNA-seq data
cat("\n[1] Loading bulk RNA-seq data...\n")
bulk_counts <- read.csv(file.path(data_dir, "transcriptomics_norm_counts.csv"), row.names = 1)
bulk_vst <- read.csv(file.path(data_dir, "transcriptomics_vst_paired.csv"), row.names = 1)
cat("  - Dimensions:", nrow(bulk_counts), "genes x", ncol(bulk_counts), "samples\n")

# Clean gene names (remove ENSEMBL suffix)
clean_gene_names <- function(genes) {
  sapply(strsplit(genes, "_"), function(x) {
    if(length(x) > 1) {
      # Try to get gene symbol via biomaRt or use fallback
      x[1]  # Return ENSEMBL ID for now
    } else {
      x[1]
    }
  })
}

# Load subtype information
subtype_file <- file.path(base_dir, "analysis/results/enhancement1_paired_subtyping/paired_subtype_K2.csv")
subtypes <- read.csv(subtype_file)
cat("  - Loaded subtype info:", nrow(subtypes), "patients\n")

# Define sample groups
normal_samples <- grep("^Normal", colnames(bulk_counts), value = TRUE)
adjacent_samples <- grep("^Adjacent", colnames(bulk_counts), value = TRUE)
cat("  - Normal samples:", length(normal_samples), "\n")
cat("  - Adjacent samples:", length(adjacent_samples), "\n")

# =============================================================================
# Part A: Build Liver Cell Type Reference Signature Matrix
# =============================================================================
cat("\n[2] Building liver cell type reference signatures...\n")

# Define marker genes for liver cell types (based on MacParland 2018 & literature)
liver_markers <- list(
  # Parenchymal cells
  Hepatocyte = c("ALB", "APOA1", "APOB", "APOC3", "TTR", "TF", "HP", "SERPINA1", 
                 "CYP3A4", "CYP2E1", "CYP1A2", "FGB", "FGG", "AGT", "AHSG", "ORM1",
                 "C3", "HPD", "TAT", "HAL", "ASS1", "CPS1"),
  
  Cholangiocyte = c("KRT19", "KRT7", "EPCAM", "SOX9", "CFTR", "MUC1", "TFF1", 
                    "TACSTD2", "FXYD2", "SPP1", "CLDN4"),
  
  # Non-parenchymal cells
  Kupffer_cell = c("CD68", "MARCO", "VSIG4", "CD163", "TIMD4", "CLEC4F", "CD14",
                   "CSF1R", "MAFB", "SPI1", "IRF8", "FCGR1A", "FCGR3A", "MERTK"),
  
  Macrophage_infiltrating = c("S100A8", "S100A9", "LYZ", "VCAN", "FCN1", "CD14",
                               "ITGAM", "CCR2", "IL1B", "TNF", "CXCL8"),
  
  Stellate_cell = c("RGS5", "ACTA2", "PDGFRB", "COL1A1", "COL1A2", "COL3A1",
                    "LRAT", "TAGLN", "CYGB", "DES", "TIMP1", "VIM"),
  
  Endothelial = c("PECAM1", "VWF", "CDH5", "KDR", "FLT1", "CLEC4G", "CLEC4M",
                  "LYVE1", "STAB1", "STAB2", "FCN2", "FCN3", "PLVAP"),
  
  # Immune cells
  NK_cell = c("NCAM1", "NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "KLRD1", "KLRF1",
              "KLRB1", "KLRC1", "NCR1", "FCGR3A"),
  
  T_cell_CD8 = c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "GZMK", "EOMES", 
                  "TBX21", "IFNG", "CCL5"),
  
  T_cell_CD4 = c("CD4", "CD3D", "CD3E", "IL7R", "LEF1", "TCF7", "FOXP3", 
                  "IL2RA", "CTLA4", "ICOS"),
  
  B_cell = c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BANK1", "CD22", 
             "BLK", "FCER2"),
  
  Plasma_cell = c("JCHAIN", "MZB1", "SDC1", "TNFRSF17", "XBP1", "PRDM1",
                  "IGHA1", "IGHG1", "IGLC2"),
  
  # Dendritic cells
  DC_conventional = c("CLEC9A", "XCR1", "BATF3", "IRF8", "THBD", "CD1C", 
                       "FCER1A", "CLEC10A"),
  
  DC_plasmacytoid = c("LILRA4", "IRF7", "IL3RA", "CLEC4C", "TCF4", "GZMB"),
  
  Neutrophil = c("FCGR3B", "CSF3R", "CXCR1", "CXCR2", "S100A8", "S100A9",
                 "CEACAM8", "MMP8", "MMP9")
)

cat("  - Defined", length(liver_markers), "cell types with marker genes\n")

# =============================================================================
# Part B: Prepare Expression Data for Deconvolution
# =============================================================================
cat("\n[3] Preparing expression data for deconvolution...\n")

# Convert ENSEMBL IDs to gene symbols using annotation
# First, try to map ENSEMBL IDs to symbols
ensembl_to_symbol <- function(bulk_data) {
  # Extract gene IDs
  gene_ids <- rownames(bulk_data)
  
  # Try to use org.Hs.eg.db if available
  if(requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    library(org.Hs.eg.db)
    ensembl_ids <- gsub("_.*", "", gene_ids)
    symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, 
                      column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    symbols[is.na(symbols)] <- ensembl_ids[is.na(symbols)]
    return(symbols)
  } else {
    # Fallback: use ENSEMBL IDs directly
    return(gsub("_.*", "", gene_ids))
  }
}

# Map genes
gene_symbols <- ensembl_to_symbol(bulk_counts)
rownames(bulk_counts) <- make.unique(gene_symbols)
rownames(bulk_vst) <- make.unique(gene_symbols)

# =============================================================================
# Part C: Create Signature Matrix for TOAST Deconvolution
# =============================================================================
cat("\n[4] Creating signature matrix from marker genes...\n")

# Build signature matrix based on marker genes
create_signature_matrix <- function(expr_matrix, marker_list, method = "mean") {
  sig_matrix <- matrix(0, nrow = nrow(expr_matrix), ncol = length(marker_list),
                       dimnames = list(rownames(expr_matrix), names(marker_list)))
  
  for(ct in names(marker_list)) {
    markers <- marker_list[[ct]]
    found_markers <- intersect(markers, rownames(expr_matrix))
    
    if(length(found_markers) > 0) {
      if(method == "mean") {
        # Use mean expression as pseudo-signature
        sig_matrix[found_markers, ct] <- 1
      }
    }
    cat("    -", ct, ":", length(found_markers), "/", length(markers), "markers found\n")
  }
  
  return(sig_matrix)
}

# Create binary marker matrix
marker_matrix <- create_signature_matrix(bulk_counts, liver_markers)

# Calculate mean expression for each cell type marker set
sig_values <- sapply(names(liver_markers), function(ct) {
  markers <- liver_markers[[ct]]
  found <- intersect(markers, rownames(bulk_counts))
  if(length(found) > 0) {
    colMeans(log2(bulk_counts[found, , drop=FALSE] + 1))
  } else {
    rep(NA, ncol(bulk_counts))
  }
})

sig_values <- sig_values[, !apply(sig_values, 2, function(x) all(is.na(x)))]
cat("  - Signature values computed for", ncol(sig_values), "cell types\n")

# =============================================================================
# Part D: TOAST-based Cell Type Deconvolution
# =============================================================================
cat("\n[5] Running TOAST deconvolution...\n")

# Prepare bulk data for TOAST
bulk_for_deconv <- as.matrix(bulk_counts)
bulk_for_deconv[is.na(bulk_for_deconv)] <- 0
bulk_for_deconv <- bulk_for_deconv[rowSums(bulk_for_deconv) > 0, ]

# Build reference-free deconvolution using marker gene enrichment approach
# Since we don't have single-cell reference, use ssGSEA-like approach

# Function to calculate enrichment score (similar to ssGSEA)
calculate_cell_scores <- function(expr_mat, gene_sets, method = "gsva") {
  # Normalize expression
  expr_scaled <- t(scale(t(log2(expr_mat + 1))))
  
  scores <- sapply(names(gene_sets), function(gs_name) {
    genes <- gene_sets[[gs_name]]
    found <- intersect(genes, rownames(expr_scaled))
    
    if(length(found) >= 3) {
      # Calculate mean z-score of marker genes
      gene_scores <- colMeans(expr_scaled[found, , drop = FALSE], na.rm = TRUE)
      return(gene_scores)
    } else {
      return(rep(NA, ncol(expr_mat)))
    }
  })
  
  scores <- scores[, !apply(scores, 2, function(x) all(is.na(x)))]
  return(as.data.frame(scores))
}

# Calculate cell type scores
cell_scores <- calculate_cell_scores(bulk_for_deconv, liver_markers)
cat("  - Calculated scores for", ncol(cell_scores), "cell types\n")

# Convert scores to proportions using softmax-like transformation
scores_to_proportions <- function(scores) {
  # Shift to positive
  scores_pos <- scores - min(scores, na.rm = TRUE) + 0.01
  
  # Normalize rows to sum to 1
  props <- t(apply(scores_pos, 1, function(x) {
    if(all(is.na(x))) return(rep(NA, length(x)))
    x[is.na(x)] <- 0
    x / sum(x)
  }))
  
  colnames(props) <- colnames(scores)
  return(as.data.frame(props))
}

cell_proportions <- scores_to_proportions(cell_scores)

cat("\n  Cell type proportion summary:\n")
print(round(colMeans(cell_proportions, na.rm = TRUE) * 100, 2))

# Save deconvolution results
write.csv(cell_scores, file.path(results_dir, "cell_type_scores.csv"))
write.csv(cell_proportions, file.path(results_dir, "cell_type_proportions.csv"))

# =============================================================================
# Part E: Statistical Analysis - Adjacent vs Normal
# =============================================================================
cat("\n[6] Statistical analysis: Adjacent vs Normal comparison...\n")

# Prepare data for paired analysis
get_paired_data <- function(data, normal_cols, adjacent_cols) {
  # Extract patient IDs
  normal_ids <- gsub("Normal", "", normal_cols)
  adjacent_ids <- gsub("Adjacent", "", adjacent_cols)
  
  # Find paired samples
  paired_ids <- intersect(normal_ids, adjacent_ids)
  
  normal_paired <- paste0("Normal", paired_ids)
  adjacent_paired <- paste0("Adjacent", paired_ids)
  
  list(
    normal = data[normal_paired, , drop = FALSE],
    adjacent = data[adjacent_paired, , drop = FALSE],
    patient_ids = paired_ids
  )
}

paired_data <- get_paired_data(cell_proportions, normal_samples, adjacent_samples)
cat("  - Found", length(paired_data$patient_ids), "paired samples\n")

# Paired Wilcoxon test for each cell type
diff_results <- data.frame(
  cell_type = colnames(cell_proportions),
  mean_Adjacent = colMeans(paired_data$adjacent, na.rm = TRUE),
  mean_Normal = colMeans(paired_data$normal, na.rm = TRUE),
  stringsAsFactors = FALSE
)

diff_results$diff <- diff_results$mean_Adjacent - diff_results$mean_Normal
diff_results$fold_change <- diff_results$mean_Adjacent / diff_results$mean_Normal

# Paired Wilcoxon test
diff_results$pvalue <- sapply(colnames(cell_proportions), function(ct) {
  adj_vals <- paired_data$adjacent[, ct]
  norm_vals <- paired_data$normal[, ct]
  
  if(sum(!is.na(adj_vals) & !is.na(norm_vals)) >= 3) {
    tryCatch({
      wilcox.test(adj_vals, norm_vals, paired = TRUE)$p.value
    }, error = function(e) NA)
  } else {
    NA
  }
})

diff_results$padj <- p.adjust(diff_results$pvalue, method = "BH")
diff_results <- diff_results[order(diff_results$pvalue), ]

cat("\n  Adjacent vs Normal - Top differential cell types:\n")
print(head(diff_results, 10))

write.csv(diff_results, file.path(results_dir, "deconv_diff_Adjacent_vs_Normal.csv"), row.names = FALSE)

# =============================================================================
# Part F: Subtype Analysis (DS1 vs DS2)
# =============================================================================
cat("\n[7] Subtype analysis: DS1 vs DS2...\n")

# Merge cell proportions with subtype info
cell_props_adjacent <- cell_proportions[adjacent_samples, ]
cell_props_adjacent$sample <- rownames(cell_props_adjacent)
cell_props_adjacent$patient_id <- as.numeric(gsub("Adjacent", "", cell_props_adjacent$sample))

subtype_merged <- merge(cell_props_adjacent, subtypes[, c("patient_id", "subtype")], by = "patient_id")
rownames(subtype_merged) <- subtype_merged$sample

# Compare DS1 vs DS2
ds1_samples <- subtype_merged$sample[subtype_merged$subtype == "DS1"]
ds2_samples <- subtype_merged$sample[subtype_merged$subtype == "DS2"]

cat("  - DS1 samples:", length(ds1_samples), "\n")
cat("  - DS2 samples:", length(ds2_samples), "\n")

subtype_diff <- data.frame(
  cell_type = colnames(cell_proportions),
  mean_DS1 = colMeans(cell_proportions[ds1_samples, ], na.rm = TRUE),
  mean_DS2 = colMeans(cell_proportions[ds2_samples, ], na.rm = TRUE),
  stringsAsFactors = FALSE
)

subtype_diff$diff <- subtype_diff$mean_DS2 - subtype_diff$mean_DS1
subtype_diff$fold_change <- subtype_diff$mean_DS2 / subtype_diff$mean_DS1

# Wilcoxon test
subtype_diff$pvalue <- sapply(colnames(cell_proportions), function(ct) {
  ds1_vals <- cell_proportions[ds1_samples, ct]
  ds2_vals <- cell_proportions[ds2_samples, ct]
  
  if(length(ds1_vals) >= 2 & length(ds2_vals) >= 2) {
    tryCatch({
      wilcox.test(ds1_vals, ds2_vals)$p.value
    }, error = function(e) NA)
  } else {
    NA
  }
})

subtype_diff$padj <- p.adjust(subtype_diff$pvalue, method = "BH")
subtype_diff <- subtype_diff[order(subtype_diff$pvalue), ]

cat("\n  DS1 vs DS2 - Top differential cell types:\n")
print(head(subtype_diff, 10))

write.csv(subtype_diff, file.path(results_dir, "deconv_diff_DS1_vs_DS2.csv"), row.names = FALSE)

# =============================================================================
# Part G: Cell Communication Analysis
# =============================================================================
cat("\n[8] Cell communication analysis using ligand-receptor pairs...\n")

# Define ligand-receptor pairs (curated from CellTalkDB and CellPhoneDB)
ligand_receptor_pairs <- data.frame(
  ligand = c("TGFB1", "TGFB1", "TGFB1", "TNF", "TNF", "IL6", "IL1B", "IL1B",
             "CXCL12", "CCL2", "CCL5", "VEGFA", "VEGFA", "HGF", "PDGFB",
             "WNT5A", "NOTCH1", "DLL4", "JAG1", "FGF2", "EGF", "IGF1",
             "MIF", "SPP1", "ANGPT1", "ANGPT2", "CSF1", "IFNG", "IL10",
             "PGF", "CTGF", "BMP4", "BMP2", "SHH", "WNT3A", "CXCL8",
             "CCL3", "CCL4", "IL4", "IL13", "TGFb1", "PDGFA", "FGF1"),
  
  receptor = c("TGFBR1", "TGFBR2", "ACVRL1", "TNFRSF1A", "TNFRSF1B", "IL6R", "IL1R1", "IL1R2",
               "CXCR4", "CCR2", "CCR5", "KDR", "FLT1", "MET", "PDGFRB",
               "FZD5", "NOTCH1", "NOTCH1", "NOTCH1", "FGFR1", "EGFR", "IGF1R",
               "CD74", "CD44", "TIE1", "TIE2", "CSF1R", "IFNGR1", "IL10RA",
               "FLT1", "TGFBR2", "BMPR1A", "BMPR2", "PTCH1", "FZD1", "CXCR1",
               "CCR1", "CCR5", "IL4R", "IL13RA1", "ITGB1", "PDGFRA", "FGFR2"),
  
  pathway = c("TGFb", "TGFb", "TGFb", "TNF", "TNF", "IL6", "IL1", "IL1",
              "Chemokine", "Chemokine", "Chemokine", "VEGF", "VEGF", "HGF", "PDGF",
              "WNT", "NOTCH", "NOTCH", "NOTCH", "FGF", "EGF", "IGF",
              "Migration", "Adhesion", "Angiopoietin", "Angiopoietin", "CSF", "Interferon", "IL10",
              "VEGF", "TGFb", "BMP", "BMP", "Hedgehog", "WNT", "Chemokine",
              "Chemokine", "Chemokine", "Th2", "Th2", "ECM", "PDGF", "FGF"),
  stringsAsFactors = FALSE
)

# Remove duplicate pairs
ligand_receptor_pairs <- unique(ligand_receptor_pairs)
cat("  - Using", nrow(ligand_receptor_pairs), "ligand-receptor pairs\n")

# Cell type - ligand/receptor associations (which cell type expresses what)
cell_lr_associations <- list(
  Hepatocyte = list(
    ligands = c("ALB", "APOA1", "AGT", "HGF", "IGF1", "IGFBP1"),
    receptors = c("MET", "EGFR", "TGFBR1", "TGFBR2")
  ),
  Kupffer_cell = list(
    ligands = c("TNF", "IL1B", "IL6", "CXCL8", "CCL2", "TGFB1"),
    receptors = c("CSF1R", "CD163", "MARCO", "TLR4", "FCGR1A")
  ),
  Macrophage_infiltrating = list(
    ligands = c("IL1B", "TNF", "CXCL8", "CCL3", "CCL4", "MIF"),
    receptors = c("CCR2", "CSF1R", "IFNGR1", "IL10RA")
  ),
  Stellate_cell = list(
    ligands = c("TGFB1", "CTGF", "PDGFB", "COL1A1", "FGF2"),
    receptors = c("PDGFRB", "TGFBR1", "TGFBR2", "FGFR1", "ITGB1")
  ),
  Endothelial = list(
    ligands = c("VEGFA", "ANGPT1", "ANGPT2", "DLL4"),
    receptors = c("KDR", "FLT1", "TIE1", "TIE2", "NOTCH1")
  ),
  NK_cell = list(
    ligands = c("IFNG", "TNF", "GZMB", "PRF1"),
    receptors = c("KLRD1", "KLRC1", "NCR1", "IL2RB")
  ),
  T_cell_CD8 = list(
    ligands = c("IFNG", "TNF", "GZMB", "CCL5"),
    receptors = c("CD8A", "PDCD1", "CTLA4", "IL2RB")
  ),
  T_cell_CD4 = list(
    ligands = c("IL4", "IL13", "IL10", "IFNG", "IL17A"),
    receptors = c("CD4", "IL7R", "ICOS", "CTLA4")
  )
)

# Calculate cell communication scores
calculate_communication <- function(expr_data, cell_props, lr_pairs, condition) {
  # Get expression values for ligands and receptors
  comm_scores <- list()
  
  for(i in 1:nrow(lr_pairs)) {
    lig <- lr_pairs$ligand[i]
    rec <- lr_pairs$receptor[i]
    path <- lr_pairs$pathway[i]
    
    if(lig %in% rownames(expr_data) & rec %in% rownames(expr_data)) {
      # Calculate mean expression across samples
      lig_expr <- mean(as.numeric(expr_data[lig, ]), na.rm = TRUE)
      rec_expr <- mean(as.numeric(expr_data[rec, ]), na.rm = TRUE)
      
      # Communication score = sqrt(ligand * receptor)
      score <- sqrt(lig_expr * rec_expr)
      
      comm_scores[[paste(lig, rec, sep = "_")]] <- data.frame(
        ligand = lig,
        receptor = rec,
        pathway = path,
        ligand_expr = lig_expr,
        receptor_expr = rec_expr,
        comm_score = score,
        condition = condition,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, comm_scores)
}

# Calculate for Normal and Adjacent
bulk_log <- log2(bulk_counts + 1)
comm_normal <- calculate_communication(bulk_log[, normal_samples], cell_proportions, 
                                        ligand_receptor_pairs, "Normal")
comm_adjacent <- calculate_communication(bulk_log[, adjacent_samples], cell_proportions,
                                          ligand_receptor_pairs, "Adjacent")

# Combine results
comm_all <- rbind(comm_normal, comm_adjacent)
cat("  - Calculated communication scores for", nrow(comm_normal), "L-R pairs\n")

# Compare communication between conditions
comm_wide <- comm_all %>%
  dplyr::select(ligand, receptor, pathway, comm_score, condition) %>%
  pivot_wider(names_from = condition, values_from = comm_score)

comm_wide$diff <- comm_wide$Adjacent - comm_wide$Normal
comm_wide$log2FC <- log2((comm_wide$Adjacent + 0.01) / (comm_wide$Normal + 0.01))
comm_wide <- comm_wide[order(-abs(comm_wide$diff)), ]

cat("\n  Top differential L-R interactions:\n")
print(head(comm_wide[, c("ligand", "receptor", "pathway", "Normal", "Adjacent", "log2FC")], 15))

write.csv(comm_wide, file.path(results_dir, "cell_communication_diff.csv"), row.names = FALSE)
write.csv(comm_all, file.path(results_dir, "cell_communication_scores.csv"), row.names = FALSE)

# =============================================================================
# Part H: Build Cell-Cell Communication Network
# =============================================================================
cat("\n[9] Building cell-cell communication network...\n")

# Aggregate communication by pathway
pathway_comm <- comm_all %>%
  group_by(pathway, condition) %>%
  summarise(
    mean_score = mean(comm_score, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = condition, values_from = mean_score)

pathway_comm$diff <- pathway_comm$Adjacent - pathway_comm$Normal
pathway_comm <- pathway_comm[order(-abs(pathway_comm$diff)), ]

cat("\n  Pathway-level communication changes:\n")
print(pathway_comm)

write.csv(pathway_comm, file.path(results_dir, "pathway_communication_summary.csv"), row.names = FALSE)

# =============================================================================
# Part I: Compare with Previous Enhancement7 Results
# =============================================================================
cat("\n[10] Comparing with previous immune deconvolution results...\n")

# Load previous ssGSEA results
prev_ssgsea <- read.csv(file.path(prev_results, "ssGSEA_immune_scores.csv"), row.names = 1)
prev_ssgsea <- t(prev_ssgsea)

# Find common cell types for comparison
# Map current cell types to previous
comparison_map <- list(
  "NK_cell" = "NK_cells",
  "T_cell_CD8" = "CD8_T_cells",
  "T_cell_CD4" = "CD4_T_helper",
  "B_cell" = "B_cells",
  "Plasma_cell" = "Plasma_cells",
  "Kupffer_cell" = "Macrophages_M2",
  "Macrophage_infiltrating" = "Macrophages_M1",
  "DC_conventional" = "Dendritic_cells",
  "Neutrophil" = "Neutrophils"
)

# Calculate correlations where possible
correlation_results <- data.frame(
  current_celltype = character(),
  previous_celltype = character(),
  correlation = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

common_samples <- intersect(rownames(cell_scores), rownames(prev_ssgsea))

for(curr_ct in names(comparison_map)) {
  prev_ct <- comparison_map[[curr_ct]]
  
  if(curr_ct %in% colnames(cell_scores) & prev_ct %in% colnames(prev_ssgsea)) {
    curr_vals <- cell_scores[common_samples, curr_ct]
    prev_vals <- prev_ssgsea[common_samples, prev_ct]
    
    if(sum(!is.na(curr_vals) & !is.na(prev_vals)) >= 5) {
      cor_test <- cor.test(curr_vals, prev_vals, method = "spearman")
      
      correlation_results <- rbind(correlation_results, data.frame(
        current_celltype = curr_ct,
        previous_celltype = prev_ct,
        correlation = cor_test$estimate,
        pvalue = cor_test$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n  Correlation with previous ssGSEA results:\n")
print(correlation_results)

write.csv(correlation_results, file.path(results_dir, "validation_correlation_with_ssGSEA.csv"), row.names = FALSE)

# =============================================================================
# Part J: Visualizations
# =============================================================================
cat("\n[11] Generating visualizations...\n")

# Set figure theme
theme_publication <- theme_bw() +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

# 1. Stacked barplot - Cell type composition per sample
cat("  - Generating stacked barplot...\n")

props_long <- cell_proportions %>%
  mutate(sample = rownames(cell_proportions)) %>%
  pivot_longer(cols = -sample, names_to = "cell_type", values_to = "proportion") %>%
  mutate(condition = ifelse(grepl("Normal", sample), "Normal", "Adjacent"))

# Define color palette for cell types
n_celltypes <- length(unique(props_long$cell_type))
cell_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_celltypes)
names(cell_colors) <- unique(props_long$cell_type)

p_stacked <- ggplot(props_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cell_colors) +
  facet_wrap(~condition, scales = "free_x") +
  labs(x = "", y = "Cell Type Proportion", fill = "Cell Type",
       title = "Cell Type Composition in HAE Samples") +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(results_dir, "fig_stacked_barplot_celltype.pdf"), p_stacked, 
       width = 10, height = 6)

# 2. Boxplot - Adjacent vs Normal for each cell type
cat("  - Generating boxplot comparison...\n")

p_boxplot <- ggplot(props_long, aes(x = cell_type, y = proportion, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_point(aes(group = condition), position = position_dodge(width = 0.6), 
             size = 0.8, alpha = 0.6) +
  scale_fill_manual(values = c("Normal" = "#3B9AB2", "Adjacent" = "#E67E50")) +
  labs(x = "", y = "Cell Type Proportion", fill = "Condition",
       title = "Cell Type Proportions: Adjacent vs Normal") +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "fig_boxplot_Adjacent_vs_Normal.pdf"), p_boxplot,
       width = 10, height = 5)

# 3. Heatmap - Cell type proportions
cat("  - Generating proportion heatmap...\n")

# Prepare annotation
sample_anno <- data.frame(
  Condition = ifelse(grepl("Normal", rownames(cell_proportions)), "Normal", "Adjacent"),
  row.names = rownames(cell_proportions)
)

# Add subtype for Adjacent samples
sample_anno$Subtype <- NA
for(i in 1:nrow(subtypes)) {
  adj_sample <- paste0("Adjacent", subtypes$patient_id[i])
  if(adj_sample %in% rownames(sample_anno)) {
    sample_anno[adj_sample, "Subtype"] <- subtypes$subtype[i]
  }
}

anno_colors <- list(
  Condition = c("Normal" = "#3B9AB2", "Adjacent" = "#E67E50"),
  Subtype = c("DS1" = "#66C2A5", "DS2" = "#FC8D62")
)

# Scale proportions for heatmap
props_scaled <- t(scale(cell_proportions))

pdf(file.path(results_dir, "fig_heatmap_celltype_proportions.pdf"), width = 10, height = 6)
pheatmap(
  props_scaled,
  annotation_col = sample_anno,
  annotation_colors = anno_colors,
  clustering_method = "ward.D2",
  show_colnames = TRUE,
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 6,
  main = "Cell Type Proportions Heatmap",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)
dev.off()

# 4. Alluvial plot - Subtype vs cell composition
cat("  - Generating alluvial plot...\n")

# Prepare data for alluvial
alluvial_data <- subtype_merged %>%
  pivot_longer(
    cols = -c(patient_id, sample, subtype),
    names_to = "cell_type",
    values_to = "proportion"
  ) %>%
  group_by(subtype, cell_type) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop")

# Discretize proportions for alluvial
alluvial_data$prop_level <- cut(alluvial_data$mean_prop, 
                                 breaks = c(0, 0.05, 0.10, 0.15, 1),
                                 labels = c("Low", "Medium", "High", "Very High"))

p_alluvial <- ggplot(alluvial_data,
                     aes(axis1 = subtype, axis2 = cell_type, y = mean_prop)) +
  geom_alluvium(aes(fill = cell_type), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
  scale_fill_manual(values = cell_colors) +
  labs(title = "Cell Type Flow Between Subtypes", y = "Mean Proportion") +
  theme_publication +
  theme(legend.position = "none")

ggsave(file.path(results_dir, "fig_alluvial_subtype_celltype.pdf"), p_alluvial,
       width = 8, height = 6)

# 5. Cell communication network
cat("  - Generating communication network...\n")

# Build network from pathway communication
net_data <- pathway_comm %>%
  filter(!is.na(Adjacent) & !is.na(Normal)) %>%
  mutate(weight = abs(diff))

# Create simple pathway-pathway network based on shared genes
# For visualization, create a network of pathways
pathway_nodes <- unique(comm_all$pathway)
n_pathways <- length(pathway_nodes)

# Create edges based on co-expression patterns
edges <- data.frame()
for(i in 1:(n_pathways-1)) {
  for(j in (i+1):n_pathways) {
    p1 <- pathway_nodes[i]
    p2 <- pathway_nodes[j]
    
    # Get L-R pairs for each pathway
    pairs1 <- comm_wide[comm_wide$pathway == p1, ]
    pairs2 <- comm_wide[comm_wide$pathway == p2, ]
    
    # Calculate correlation of their log2FC
    if(nrow(pairs1) > 0 & nrow(pairs2) > 0) {
      mean_fc1 <- mean(pairs1$log2FC, na.rm = TRUE)
      mean_fc2 <- mean(pairs2$log2FC, na.rm = TRUE)
      
      edges <- rbind(edges, data.frame(
        from = p1, to = p2,
        weight = abs(mean_fc1) + abs(mean_fc2)
      ))
    }
  }
}

# Create igraph network
if(nrow(edges) > 0) {
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = pathway_nodes)
  
  # Calculate node attributes
  node_attrs <- pathway_comm %>%
    mutate(
      change = Adjacent - Normal,
      direction = ifelse(change > 0, "Up in Adjacent", "Down in Adjacent")
    )
  
  V(g)$change <- node_attrs$diff[match(V(g)$name, node_attrs$pathway)]
  V(g)$color <- ifelse(V(g)$change > 0, "#E67E50", "#3B9AB2")
  V(g)$size <- abs(V(g)$change) * 20 + 5
  
  pdf(file.path(results_dir, "fig_pathway_communication_network.pdf"), width = 8, height = 8)
  set.seed(42)
  layout <- layout_with_fr(g)
  plot(g, layout = layout,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       edge.width = E(g)$weight,
       main = "Pathway Communication Network\n(Red=Up in Adjacent, Blue=Down)")
  legend("bottomright", 
         legend = c("Up in Adjacent", "Down in Adjacent"),
         fill = c("#E67E50", "#3B9AB2"),
         bty = "n", cex = 0.8)
  dev.off()
}

# 6. L-R pair heatmap
cat("  - Generating L-R pair heatmap...\n")

# Top differential L-R pairs
top_lr <- head(comm_wide[order(-abs(comm_wide$log2FC)), ], 30)

lr_mat <- as.matrix(top_lr[, c("Normal", "Adjacent")])
rownames(lr_mat) <- paste(top_lr$ligand, top_lr$receptor, sep = " → ")

pdf(file.path(results_dir, "fig_heatmap_LR_pairs.pdf"), width = 6, height = 8)
pheatmap(
  lr_mat,
  cluster_cols = FALSE,
  scale = "row",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  fontsize = 8,
  fontsize_row = 7,
  main = "Top Differential Ligand-Receptor Pairs",
  labels_col = c("Normal", "Adjacent")
)
dev.off()

# 7. Correlation scatter plot with previous results
cat("  - Generating validation correlation plot...\n")

if(nrow(correlation_results) > 0) {
  # Create scatter plot for one representative cell type
  best_corr <- correlation_results[which.max(abs(correlation_results$correlation)), ]
  
  curr_ct <- best_corr$current_celltype
  prev_ct <- best_corr$previous_celltype
  
  scatter_data <- data.frame(
    current = cell_scores[common_samples, curr_ct],
    previous = prev_ssgsea[common_samples, prev_ct],
    sample = common_samples
  )
  scatter_data$condition <- ifelse(grepl("Normal", scatter_data$sample), "Normal", "Adjacent")
  
  p_scatter <- ggplot(scatter_data, aes(x = current, y = previous, color = condition)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "grey50") +
    scale_color_manual(values = c("Normal" = "#3B9AB2", "Adjacent" = "#E67E50")) +
    labs(
      x = paste("Current:", curr_ct, "(marker score)"),
      y = paste("Previous:", prev_ct, "(ssGSEA score)"),
      title = paste("Validation: r =", round(best_corr$correlation, 3))
    ) +
    theme_publication
  
  ggsave(file.path(results_dir, "fig_validation_correlation.pdf"), p_scatter,
         width = 6, height = 5)
}

# =============================================================================
# Summary Report
# =============================================================================
cat("\n" %+% paste(rep("=", 80), collapse="") %+% "\n")
cat("ANALYSIS SUMMARY\n")
cat(paste(rep("=", 80), collapse="") %+% "\n")

cat("\n[A] Cell Type Deconvolution Results:\n")
cat("  - Total cell types analyzed:", ncol(cell_proportions), "\n")
cat("  - Mean proportions (top 5):\n")
top_props <- sort(colMeans(cell_proportions, na.rm = TRUE), decreasing = TRUE)
for(i in 1:min(5, length(top_props))) {
  cat("    ", names(top_props)[i], ":", round(top_props[i] * 100, 1), "%\n")
}

cat("\n[B] Adjacent vs Normal Differential Cell Types:\n")
sig_diff <- diff_results[diff_results$pvalue < 0.1 & !is.na(diff_results$pvalue), ]
if(nrow(sig_diff) > 0) {
  for(i in 1:nrow(sig_diff)) {
    direction <- ifelse(sig_diff$diff[i] > 0, "↑", "↓")
    cat("  ", direction, sig_diff$cell_type[i], 
        "(p=", round(sig_diff$pvalue[i], 3), 
        ", FC=", round(sig_diff$fold_change[i], 2), ")\n")
  }
} else {
  cat("  - No cell types with p < 0.1\n")
  cat("  - Top changes (by effect size):\n")
  top_changes <- head(diff_results[order(-abs(diff_results$diff)), ], 3)
  for(i in 1:nrow(top_changes)) {
    direction <- ifelse(top_changes$diff[i] > 0, "↑", "↓")
    cat("    ", direction, top_changes$cell_type[i],
        "(diff=", round(top_changes$diff[i], 3), ")\n")
  }
}

cat("\n[C] Subtype (DS1 vs DS2) Differential Cell Types:\n")
sig_subtype <- subtype_diff[subtype_diff$pvalue < 0.1 & !is.na(subtype_diff$pvalue), ]
if(nrow(sig_subtype) > 0) {
  for(i in 1:min(5, nrow(sig_subtype))) {
    direction <- ifelse(sig_subtype$diff[i] > 0, "↑DS2", "↑DS1")
    cat("  ", direction, ":", sig_subtype$cell_type[i],
        "(p=", round(sig_subtype$pvalue[i], 3), ")\n")
  }
} else {
  cat("  - No cell types with p < 0.1\n")
}

cat("\n[D] Key Cell Communication Findings:\n")
cat("  - Top upregulated pathways in Adjacent:\n")
up_pathways <- pathway_comm[pathway_comm$diff > 0, ]
up_pathways <- up_pathways[order(-up_pathways$diff), ]
for(i in 1:min(3, nrow(up_pathways))) {
  cat("    ↑", up_pathways$pathway[i], "(+", round(up_pathways$diff[i], 2), ")\n")
}

cat("  - Top downregulated pathways in Adjacent:\n")
down_pathways <- pathway_comm[pathway_comm$diff < 0, ]
down_pathways <- down_pathways[order(down_pathways$diff), ]
for(i in 1:min(3, nrow(down_pathways))) {
  cat("    ↓", down_pathways$pathway[i], "(", round(down_pathways$diff[i], 2), ")\n")
}

cat("\n[E] Validation with Previous ssGSEA:\n")
if(nrow(correlation_results) > 0) {
  sig_corr <- correlation_results[correlation_results$pvalue < 0.05, ]
  cat("  - Significant correlations found:", nrow(sig_corr), "/", nrow(correlation_results), "\n")
  if(nrow(sig_corr) > 0) {
    for(i in 1:nrow(sig_corr)) {
      cat("    ", sig_corr$current_celltype[i], "vs", sig_corr$previous_celltype[i],
          ": r =", round(sig_corr$correlation[i], 2), "\n")
    }
  }
} else {
  cat("  - No correlation data available\n")
}

cat("\n[F] Output Files:\n")
cat("  Results directory:", results_dir, "\n")
output_files <- list.files(results_dir)
for(f in output_files) {
  cat("    -", f, "\n")
}

cat("\n" %+% paste(rep("=", 80), collapse="") %+% "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse="") %+% "\n")
