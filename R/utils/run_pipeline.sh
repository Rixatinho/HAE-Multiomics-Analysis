# NOTE (2026-09-27, release v1.2.1): the SuppFig_NN directory names below are PRE-RENUMBERING
# labels. All 32 supplementary figures were renumbered by first citation order in the
# 2026-09 revision; the authoritative numbering is the Supplementary Figure Legends block in
# the submitted Supplementary Information. Do not treat these folder names as current figure
# numbers. This script is retained for reproducibility of the underlying panel PDFs only.
#!/bin/bash
cd "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
LOG="pipeline_run.log"
echo "=== START $(date) ===" > "$LOG"

echo "Running nc_production_pipeline.R..." >> "$LOG"
conda run -n multiomics Rscript analysis/scripts/nc_production_pipeline.R >> "$LOG" 2>&1
echo "Exit code: $?" >> "$LOG"

echo "" >> "$LOG"
echo "Running enhance_metabolomics_enrichment.R..." >> "$LOG"
conda run -n multiomics Rscript analysis/scripts/enhance_metabolomics_enrichment.R >> "$LOG" 2>&1
echo "Exit code: $?" >> "$LOG"

echo "" >> "$LOG"
echo "Running enhance_subtype_robustness_v2.R..." >> "$LOG"
conda run -n multiomics Rscript analysis/scripts/enhance_subtype_robustness_v2.R >> "$LOG" 2>&1
echo "Exit code: $?" >> "$LOG"

echo "" >> "$LOG"
echo "Running enhance_biomarker_robustness.R..." >> "$LOG"
conda run -n multiomics Rscript analysis/scripts/enhance_biomarker_robustness.R >> "$LOG" 2>&1
echo "Exit code: $?" >> "$LOG"

echo "" >> "$LOG"
echo "Running enhance_mofa_annotation.R..." >> "$LOG"
conda run -n multiomics Rscript analysis/scripts/enhance_mofa_annotation.R >> "$LOG" 2>&1
echo "Exit code: $?" >> "$LOG"

echo "" >> "$LOG"
echo "=== ALL SCRIPTS DONE $(date) ===" >> "$LOG"

# Copy enhance outputs to per-figure subdirectories
echo "Copying enhance outputs..." >> "$LOG"

# SuppFig_19
mkdir -p "图片/SuppFig_19"
for src in analysis/results/metabolomics_enrichment/figures/*.pdf; do
  fname=$(basename "$src")
  case "$fname" in
    *ora*|*ORA*) cp "$src" "图片/SuppFig_19/Supp19a_metabolic_pathway_ora.pdf" ;;
    *chemical*) cp "$src" "图片/SuppFig_19/Supp19b_chemical_class_enrichment.pdf" ;;
    *msea*|*MSEA*) cp "$src" "图片/SuppFig_19/Supp19c_msea_dotplot.pdf" ;;
    *consistency*|*cross*) cp "$src" "图片/SuppFig_19/Supp19d_cross_omics_consistency.pdf" ;;
  esac
done

# SuppFig_20
mkdir -p "图片/SuppFig_20"
for src in analysis/results/subtype_robustness_v2/figures/*.pdf; do
  fname=$(basename "$src")
  case "$fname" in
    *k_selection*) cp "$src" "图片/SuppFig_20/Supp20a_k_selection_comprehensive.pdf" ;;
    *effect*) cp "$src" "图片/SuppFig_20/Supp20b_molecular_effect_sizes.pdf" ;;
  esac
done

# SuppFig_21 (only a and b, no c)
mkdir -p "图片/SuppFig_21"
for src in analysis/results/biomarker_robustness/*.pdf; do
  fname=$(basename "$src")
  case "$fname" in
    *loo*|*roc*) cp "$src" "图片/SuppFig_21/Supp21a_biomarker_loo_roc.pdf" ;;
    *model*|*comparison*) cp "$src" "图片/SuppFig_21/Supp21b_biomarker_model_comparison.pdf" ;;
  esac
done
# Remove old panel c if it exists
rm -f "图片/SuppFig_21/Supp21c_biomarker_gse124362_validation.pdf"

# SuppFig_22
mkdir -p "图片/SuppFig_22"
for src in analysis/results/mofa_annotation/figures/*.pdf; do
  fname=$(basename "$src")
  case "$fname" in
    *variance*) cp "$src" "图片/SuppFig_22/Supp22a_factor_variance.pdf" ;;
    *pathway*) cp "$src" "图片/SuppFig_22/Supp22b_factor_pathway.pdf" ;;
    *clinical*) cp "$src" "图片/SuppFig_22/Supp22c_factor_clinical.pdf" ;;
    *features*|*top*) cp "$src" "图片/SuppFig_22/Supp22d_factor_top_features.pdf" ;;
  esac
done

echo "Copy done" >> "$LOG"

# Run assembly
echo "Running assemble_figures_v3.py..." >> "$LOG"
conda run -n multiomics python analysis/scripts/assemble_figures_v3.py >> "$LOG" 2>&1
echo "Assembly exit code: $?" >> "$LOG"

# Convert to PNG (in-place within each subdirectory)
echo "Converting to PNG..." >> "$LOG"
conda run -n multiomics python analysis/scripts/convert_pdf_png.py >> "$LOG" 2>&1
echo "PNG conversion exit code: $?" >> "$LOG"

echo "=== COMPLETE $(date) ===" >> "$LOG"
