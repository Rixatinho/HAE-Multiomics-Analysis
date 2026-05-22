#!/usr/bin/env Rscript
# =============================================================================
# Fig1_multiomics_landscape.R -- Multi-omics Landscape Heatmap (Panel b)
# Figure 1b | 183 x 158 mm | ComplexHeatmap stacked (TC+PR+MB+Immune)
# =============================================================================

slot <- get_panel_slot("Fig1_multiomics_landscape.pdf")
W <- slot$w; H <- slot$h   # 183 x 158 mm

# --- Data loading ---
tc_mat <- get_tc_mat(); pr_mat <- get_pr_mat(); mb_mat <- get_mb_mat()
degs <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
deps <- read.csv(file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
dems <- read.csv(file.path(RES, "phase1_diff/DEMs_Adjacent_vs_Normal.csv"), stringsAsFactors = FALSE)
subtypes <- read.csv(file.path(RES, "enhancement1_paired_subtyping/paired_subtype_K2.csv"), stringsAsFactors = FALSE)
clinical <- read.csv(file.path(RES, "enhancement1_paired_subtyping/paired_clinical_with_subtypes_K2.csv"), stringsAsFactors = FALSE)
ssgsea <- read.csv(file.path(RES, "enhancement7_immune_deconvolution/ssGSEA_immune_scores.csv"), row.names = 1, check.names = FALSE)

# --- Common patients across 3 omics ---
get_pid <- function(s) as.integer(gsub("^(Adjacent|Normal)", "", s))
adj_tc <- grep("^Adjacent", colnames(tc_mat), value = TRUE)
adj_pr <- grep("^Adjacent", colnames(pr_mat), value = TRUE)
adj_mb <- grep("^Adjacent", colnames(mb_mat), value = TRUE)
common_pids <- sort(Reduce(intersect, list(get_pid(adj_tc), get_pid(adj_pr), get_pid(adj_mb))))
cat("  Common patients:", length(common_pids), "\n")
adj_tc_s <- paste0("Adjacent", common_pids)
adj_pr_s <- paste0("Adjacent", common_pids)
adj_mb_s <- paste0("Adjacent", common_pids)

# --- Subtype mapping ---
subtype_map <- setNames(subtypes$subtype, paste0("P", subtypes$patient_id))
patient_subtypes <- subtype_map[paste0("P", common_pids)]
patient_subtypes <- gsub("^DS", "CS", patient_subtypes)

# --- Layer 1: Top DEGs (8 up + 7 down = 15) ---
degs_s <- degs[order(degs$P.Value), ]
degs_s <- degs_s[!is.na(degs_s$gene_name) & degs_s$gene_name != "", ]
degs_s$gene_name <- ensg_to_symbol(degs_s$gene_name)
degs_s <- degs_s[!grepl("^ENSG", degs_s$gene_name), ]
top_degs <- rbind(head(degs_s[degs_s$logFC > 0, ], 8), head(degs_s[degs_s$logFC < 0, ], 7))
tc_ids <- top_degs$gene_id[top_degs$gene_id %in% rownames(tc_mat)]
tc_sub <- tc_mat[tc_ids, adj_tc_s, drop = FALSE]
rownames(tc_sub) <- top_degs$gene_name[match(rownames(tc_sub), top_degs$gene_id)]
tc_scaled <- scale_rows(tc_sub, 2)

# --- Layer 2: Top DEPs (5 up + 5 down = 10) ---
deps_s <- deps[order(deps$P.Value), ]
deps_s <- deps_s[!is.na(deps_s$gene_name) & deps_s$gene_name != "", ]
deps_s$gene_name <- ensg_to_symbol(deps_s$gene_name)
deps_s <- deps_s[!grepl("^ENSG", deps_s$gene_name), ]
top_deps <- rbind(head(deps_s[deps_s$logFC > 0, ], 5), head(deps_s[deps_s$logFC < 0, ], 5))
pr_ids <- top_deps$Protein[top_deps$Protein %in% rownames(pr_mat)]
pr_sub <- pr_mat[pr_ids, adj_pr_s, drop = FALSE]
rownames(pr_sub) <- top_deps$gene_name[match(rownames(pr_sub), top_deps$Protein)]
pr_scaled <- scale_rows(pr_sub, 2)

# --- Layer 3: Top DEMs (5 up + 5 down = 10) ---
dems_s <- dems[order(dems$P.Value), ]
dems_s <- dems_s[!is.na(dems_s$metabolite_name) & dems_s$metabolite_name != "", ]
top_dems <- rbind(head(dems_s[dems_s$logFC > 0, ], 5), head(dems_s[dems_s$logFC < 0, ], 5))
mb_ids <- top_dems$Compound_ID[top_dems$Compound_ID %in% rownames(mb_mat)]
mb_sub <- mb_mat[mb_ids, adj_mb_s, drop = FALSE]
mb_mapping <- data.frame(id = top_dems$Compound_ID, name = top_dems$metabolite_name, stringsAsFactors = FALSE)
mb_names <- format_feature_names(rownames(mb_sub), mb_mapping, type = "MT", max_len = 35)
rownames(mb_sub) <- make.unique(mb_names)
mb_scaled <- scale_rows(mb_sub, 2)

# --- Layer 4: Immune signatures (10 fixed) ---
immune_sigs <- c("Macrophages_M1", "Macrophages_M2", "Monocytes", "Dendritic_cells",
                 "NK_cells", "CD8_T_cells", "Treg", "Th1", "Th2", "Th17")
immune_avail <- immune_sigs[immune_sigs %in% rownames(ssgsea)]
imm_sub <- ssgsea[immune_avail, adj_tc_s, drop = FALSE]
rownames(imm_sub) <- clean_celltype(rownames(imm_sub))
imm_scaled <- scale_rows(imm_sub, 2)

cat("  Layers: TC=", nrow(tc_scaled), " PR=", nrow(pr_scaled),
    " MB=", nrow(mb_scaled), " IMM=", nrow(imm_scaled), "\n")

# --- Clinical annotation ---
ha_col <- get_nc_annotation(patient_subtypes, clinical, common_pids)
col_expr <- colorRamp2(c(-2, 0, 2), c(COL_DOWN, "white", COL_UP))

# --- Row height: fit all layers into ~120mm (158mm - ~38mm overhead) ---
total_rows <- nrow(tc_scaled) + nrow(pr_scaled) + nrow(mb_scaled) + nrow(imm_scaled)
RH <- min(3, floor(120 / total_rows * 10) / 10)  # auto-fit, max 3mm
cat("  Row height:", RH, "mm  | Total rows:", total_rows, "\n")

# --- Build heatmaps ---
ht_tc <- Heatmap(tc_scaled, name = "Gene\nExpr (z)", col = col_expr,
  row_names_gp = gp_row_names(7), show_column_names = FALSE,
  cluster_columns = TRUE, cluster_rows = TRUE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_title = "Transcriptomics", row_title_gp = gp_row_title(),
  top_annotation = ha_col, height = unit(nrow(tc_scaled) * RH, "mm"),
  heatmap_legend_param = std_legend_param())

ht_pr <- Heatmap(pr_scaled, name = "Protein\nExpr (z)", col = col_expr,
  row_names_gp = gp_row_names(5.5), show_column_names = FALSE,
  cluster_columns = FALSE, cluster_rows = TRUE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_title = "Proteomics", row_title_gp = gp_row_title(),
  height = unit(nrow(pr_scaled) * RH, "mm"),
  heatmap_legend_param = std_legend_param())

ht_mb <- Heatmap(mb_scaled, name = "Metab\nExpr (z)", col = col_expr,
  row_names_gp = gp_row_names(5), show_column_names = FALSE,
  cluster_columns = FALSE, cluster_rows = TRUE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_title = "Metabolomics", row_title_gp = gp_row_title(),
  height = unit(nrow(mb_scaled) * RH, "mm"),
  heatmap_legend_param = std_legend_param())

ht_imm <- Heatmap(imm_scaled, name = "Immune\nScore (z)", col = col_immune,
  row_names_gp = gp_row_names(7), show_column_names = TRUE,
  column_names_gp = gp_col_names(5.5, bold = FALSE),
  cluster_columns = FALSE, cluster_rows = TRUE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_title = "Immune", row_title_gp = gp_row_title(),
  height = unit(nrow(imm_scaled) * RH, "mm"),
  heatmap_legend_param = std_legend_param())

ht_list <- ht_tc %v% ht_pr %v% ht_mb %v% ht_imm

# --- Render at exact slot dimensions ---
s <- sp("Fig1_multiomics_landscape.pdf", W, H)
cairo_pdf(s$path, width = s$width, height = s$height, family = FONT_FAMILY)
draw(ht_list,
     column_title = "Multi-omics Molecular Landscape of HAE Adjacent Liver Tissue",
     column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_GRID),
     padding = unit(c(2, 2, 2, 10), "mm"))
dev.off()
cat("  -> Fig1_multiomics_landscape.pdf (", W, "x", H, "mm)\n")
