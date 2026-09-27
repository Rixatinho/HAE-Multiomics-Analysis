# 09_causal_network_medicine

Analysis scripts for Results Section R8: Network Medicine and Drug Repurposing.

> **Note (2026-09-27, release v1.2.1).** The Mendelian randomisation scripts previously shipped under
> `R/07_biomarker_discovery/` (`enhance_mr_analysis.R`, `enhance_mr_twosamplemr.R`,
> `enhance_mr_real_gwas.R`, `enhance_mr_minimal_test.R`) have been removed from this repository.
> The corresponding analysis was excised from the submitted manuscript; causal inference now rests on
> within-patient mediation analysis (`enhancements/enhancement36_mediation.py`) and on Bayesian /
> physiologically based pharmacokinetic modelling (`enhancements/enhancement48_pbpk_cyp3a4.R`).

## Source Scripts (in 02_analysis/scripts/02_enhancements/)

| Script | Analysis |
|--------|----------|
| enhance_network_medicine.R | Network medicine drug proximity |
| enhance6_drug_repurposing.R | Drug repurposing via target overlap |
| enhance_pan_liver_comparison.R | Pan-liver disease comparison |
| enhance_pan_liver_v2.R | Pan-liver disease comparison (v2, final) |

## Results Directories (in 02_analysis/results/)

- `enhancement19_network_medicine/` - Drug proximity z-scores, candidate ranking
- `molecular_docking/` - AutoDock Vina docking outputs
- `gdsc_pharmacogenomic/` - GDSC cell line validation
- `enhancement20_pan_liver/` - Cross-disease overlap, forest plots

## Supplementary Figures (current numbering, 2026-09 revision)

- **Supplementary Fig. 24**: Molecular docking validation and pharmacogenomic analysis
- **Supplementary Fig. 26**: Drug repurposing and pan-liver disease comparison
