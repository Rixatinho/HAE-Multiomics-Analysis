# 09_causal_network_medicine

Analysis scripts for Results Section R8: Causal Inference and Network Medicine.

## Source Scripts (in 02_analysis/scripts/02_enhancements/)

| Script | Analysis |
|--------|----------|
| enhance_mr_analysis.R | Mendelian Randomization (initial) |
| enhance_mr_twosamplemr_FIXED.R | Two-sample MR |
| compute_mr_fstats.R | MR instrument strength (F-statistics) |
| enhance_network_medicine.R | Network medicine drug proximity |
| enhance6_drug_repurposing.R | Drug repurposing via target overlap |
| enhance_pan_liver_comparison.R | Pan-liver disease comparison |
| enhance_pan_liver_v2.R | Pan-liver disease comparison (v2, final) |

## Results Directories (in 02_analysis/results/)

- `enhancement18_mr/` - MR results (forest plots, funnel plots, sensitivity)
- `enhancement19_network_medicine/` - Drug proximity z-scores, candidate ranking
- `molecular_docking/` - AutoDock Vina docking outputs
- `gdsc_pharmacogenomic/` - GDSC cell line validation
- `enhancement20_pan_liver/` - Cross-disease overlap, forest plots

## Supplementary Figures

- **Fig. S13**: Molecular Docking + Pharmacogenomic Validation (10 panels)
- **Fig. S14**: MR/Drug Repurposing + Pan-Liver Disease Comparison (9 panels)
