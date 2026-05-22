#!/usr/bin/env Rscript
# =============================================================================
# COMPILE SUPPLEMENTARY TABLES into single Excel workbook
# DEPRECATED: Early-stage prototype (ST1-ST18 only). Final submission uses
# 42 tables (ST1-ST42) in 07_supplementary/Supplementary_Tables.xlsx.
# Retained for historical reference. See REVISION_CHECKLIST.md.
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
})

BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
RES  <- file.path(BASE, "analysis/results")
OUT  <- file.path(BASE, "07_supplementary")

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", fontSize = 11, 
                            fgFill = "#4472C4", fontColour = "white",
                            halign = "center", border = "TopBottomLeftRight")
note_style <- createStyle(fontSize = 10, fontColour = "#666666", wrapText = TRUE)

safe_read <- function(path) {
  if (!file.exists(path)) {
    cat("  WARNING: missing", basename(path), "\n")
    return(data.frame(NOTE = paste("File not found:", basename(path))))
  }
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) data.frame(ERROR = conditionMessage(e)))
}

add_table <- function(wb, sheet_name, data, title = NULL) {
  addWorksheet(wb, sheet_name)
  start_row <- 1
  if (!is.null(title)) {
    writeData(wb, sheet_name, title, startRow = 1, startCol = 1)
    addStyle(wb, sheet_name, note_style, rows = 1, cols = 1)
    start_row <- 3
  }
  writeData(wb, sheet_name, data, startRow = start_row, headerStyle = header_style)
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
}

cat("================================================================\n")
cat("   COMPILING SUPPLEMENTARY TABLES\n")
cat("================================================================\n\n")

# ST1: Clinical characteristics
cat("--- ST1: Clinical characteristics ---\n")
clin_file <- file.path(BASE, "14 例肝泡型棘球蚴病患者临床信息.xlsx")
if (file.exists(clin_file)) {
  clin <- tryCatch(read.xlsx(clin_file, sheet = 1), error = function(e) NULL)
  if (!is.null(clin)) {
    add_table(wb, "ST1_Clinical", clin,
              "Supplementary Table 1. Clinical characteristics of 14 HAE patients")
    cat("  Added (", nrow(clin), " patients)\n")
  }
} else {
  addWorksheet(wb, "ST1_Clinical")
  writeData(wb, "ST1_Clinical", "Clinical data to be added from source file")
}

# ST2: DEGs
cat("--- ST2: DEGs ---\n")
degs_all <- safe_read(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"))
add_table(wb, "ST2_DEGs_all", degs_all,
          "Supplementary Table 2. All differentially expressed genes (transcriptomics)")
cat("  Added (", nrow(degs_all), " genes)\n")

# ST3: DEPs
cat("--- ST3: DEPs ---\n")
deps_all <- safe_read(file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"))
add_table(wb, "ST3_DEPs_all", deps_all,
          "Supplementary Table 3. All differentially expressed proteins (proteomics)")
cat("  Added (", nrow(deps_all), " proteins)\n")

# ST4: DEMs
cat("--- ST4: DEMs ---\n")
dems_all <- safe_read(file.path(RES, "phase1_diff/DEMs_Adjacent_vs_Normal.csv"))
add_table(wb, "ST4_DEMs_all", dems_all,
          "Supplementary Table 4. All differentially abundant metabolites (metabolomics)")
cat("  Added (", nrow(dems_all), " metabolites)\n")

# ST5: Strict DEGs/DEPs/DEMs
cat("--- ST5: Significant molecules ---\n")
degs_sig <- safe_read(file.path(RES, "phase1_diff/DEGs_significant.csv"))
deps_sig <- safe_read(file.path(RES, "phase1_diff/DEPs_significant.csv"))
dems_sig <- safe_read(file.path(RES, "phase1_diff/DEMs_significant.csv"))
add_table(wb, "ST5a_DEGs_sig", degs_sig, "ST5a. Significant DEGs (strict threshold)")
add_table(wb, "ST5b_DEPs_sig", deps_sig, "ST5b. Significant DEPs (strict threshold)")
add_table(wb, "ST5c_DEMs_sig", dems_sig, "ST5c. Significant DEMs (strict threshold)")
cat("  Added\n")

# ST6: Cross-omics convergent Hallmark pathways
cat("--- ST6: Convergent pathways ---\n")
hallmark <- safe_read(file.path(RES, "phase2_enrichment/cross_omics_convergent_Hallmark.csv"))
add_table(wb, "ST6_Hallmark", hallmark,
          "Supplementary Table 6. Cross-omics convergent Hallmark pathways")
cat("  Added (", nrow(hallmark), " pathways)\n")

# ST7: GSEA results (all databases)
cat("--- ST7: GSEA results ---\n")
gsea_tc_h <- safe_read(file.path(RES, "phase2_enrichment/GSEA_transcriptomics_Hallmark.csv"))
gsea_pr_h <- safe_read(file.path(RES, "phase2_enrichment/GSEA_proteomics_Hallmark.csv"))
add_table(wb, "ST7a_GSEA_TC_Hallmark", gsea_tc_h, "ST7a. GSEA Transcriptomics Hallmark")
add_table(wb, "ST7b_GSEA_PR_Hallmark", gsea_pr_h, "ST7b. GSEA Proteomics Hallmark")
cat("  Added\n")

# ST8: Metabolic pathway summary
cat("--- ST8: Metabolic pathways ---\n")
metab <- safe_read(file.path(RES, "phase3_metabolic/metabolic_pathway_summary.csv"))
add_table(wb, "ST8_Metabolic", metab,
          "Supplementary Table 8. Metabolic pathway gene-level analysis")
cat("  Added (", nrow(metab), " pathways)\n")

# ST9: WGCNA module overlap
cat("--- ST9: WGCNA overlap ---\n")
wgcna <- safe_read(file.path(RES, "enhancement2_wgcna/cross_omics_module_overlap.csv"))
add_table(wb, "ST9_WGCNA_overlap", wgcna,
          "Supplementary Table 9. Cross-omics WGCNA module overlap")
cat("  Added (", nrow(wgcna), " module pairs)\n")

# ST10: K selection metrics
cat("--- ST10: K selection ---\n")
ksel <- safe_read(file.path(RES, "phase6_subtyping/K_selection_metrics.csv"))
add_table(wb, "ST10_K_selection", ksel,
          "Supplementary Table 10. Consensus clustering metrics for K=2 to K=6")
cat("  Added\n")

# ST11: Clinical associations K=2
cat("--- ST11: Clinical associations ---\n")
clin_k2 <- safe_read(file.path(RES, "phase6_subtyping/clinical_association_K2.csv"))
add_table(wb, "ST11_Clinical_K2", clin_k2,
          "Supplementary Table 11. Clinical variable associations with molecular subtypes")
cat("  Added\n")

# ST12: Druggable targets
cat("--- ST12: Druggable targets ---\n")
drugs <- safe_read(file.path(RES, "phase7_characterization/druggable_targets_by_subtype.csv"))
if (nrow(drugs) == 1 && "NOTE" %in% names(drugs)) {
  drugs <- safe_read(file.path(RES, "enhancement6_drug_repurposing/drug_target_scored.csv"))
}
add_table(wb, "ST12_Druggable", drugs,
          "Supplementary Table 12. Druggable targets by subtype")
cat("  Added (", nrow(drugs), " targets)\n")

# ST13: Diagnostic panel
cat("--- ST13: Diagnostic panel ---\n")
panel <- safe_read(file.path(RES, "phase8_biomarker/multiomics_diagnostic_panel.csv"))
add_table(wb, "ST13_Diagnostic_panel", panel,
          "Supplementary Table 13. Multi-omics diagnostic panel (18 features)")
cat("  Added\n")

# ST14: Individual biomarkers
cat("--- ST14: Biomarkers ---\n")
biomark <- safe_read(file.path(RES, "phase8_biomarker/disease_biomarkers.csv"))
add_table(wb, "ST14_Biomarkers", biomark,
          "Supplementary Table 14. Individual biomarker performance")
cat("  Added (", nrow(biomark), " biomarkers)\n")

# ST15: Meta-analysis cohort summary
cat("--- ST15: Meta cohorts ---\n")
meta_sum <- safe_read(file.path(RES, "enhance_meta_analysis/effective_sample_summary.csv"))
add_table(wb, "ST15_Meta_cohorts", meta_sum,
          "Supplementary Table 15. Six-cohort meta-analysis details")
cat("  Added\n")

# ST16: Consensus DEGs from meta-analysis (top 500 to keep file manageable)
cat("--- ST16: Consensus DEGs ---\n")
cons_degs <- safe_read(file.path(RES, "enhance_meta_analysis/meta_analysis_consensus_DEGs.csv"))
add_table(wb, "ST16_Consensus_DEGs", cons_degs,
          "Supplementary Table 16. Consensus DEGs from 6-cohort meta-analysis (2344 genes)")
cat("  Added (", nrow(cons_degs), " genes)\n")

# ST17: External validation AUC
cat("--- ST17: Validation AUC ---\n")
val_auc <- safe_read(file.path(RES, "enhance_meta_analysis/signature_validation_AUC.csv"))
add_table(wb, "ST17_Validation_AUC", val_auc,
          "Supplementary Table 17. External validation AUC for diagnostic signature")
cat("  Added\n")

# ST18: mRNA-protein discordance
cat("--- ST18: Discordance ---\n")
discord <- safe_read(file.path(RES, "enhancement13_discordance/discordant_genes.csv"))
add_table(wb, "ST18_Discordance", discord,
          "Supplementary Table 18. mRNA-protein discordant genes")
cat("  Added (", nrow(discord), " genes)\n")

# Save workbook
out_file <- file.path(OUT, "Supplementary_Tables.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
cat("\n================================================================\n")
cat("   SUPPLEMENTARY TABLES COMPILED\n")
cat("   Output:", out_file, "\n")
cat("   Sheets:", length(names(wb)), "\n")
cat("================================================================\n")
out_file <- file.path(OUT, "Supplementary_Tables.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
cat("\n================================================================\n")
cat("   SUPPLEMENTARY TABLES COMPILED\n")
cat("   Output:", out_file, "\n")
cat("   Sheets:", length(names(wb)), "\n")
cat("================================================================\n")
