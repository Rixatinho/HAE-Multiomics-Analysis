# HAE-Multiomics-Analysis

Source code accompanying the manuscript:

> **Peri-lesional loss of hepatic CYP3A4 may underlie subtherapeutic albendazole bioactivation in human alveolar echinococcosis**

This repository contains the R analysis pipeline and figure-generation scripts used to produce all main and supplementary figures of the study. The pipeline integrates paired bulk transcriptomic, proteomic, and metabolomic data from 14 patients with hepatic alveolar echinococcosis (HAE), together with cell-type deconvolution, gene regulatory network inference, multi-omics integration, biomarker discovery, and quantitative pharmacological modelling.

---

## Repository layout

```
HAE-Multiomics-Analysis/
├── R/                              Main analysis pipeline
│   ├── 00_preprocessing/           Quality control, normalisation, batch correction
│   ├── 01_differential_analysis/   limma / DESeq2 / edgeR differential testing
│   ├── 02_pathway_enrichment/      GSEA, GO, KEGG, Reactome, ssGSEA
│   ├── 03_metabolic_reprogramming/ Metabolic pathway analysis, drug repurposing
│   ├── 04_immune_microenvironment/ BayesPrism, CellChat, TCR/BCR repertoire
│   ├── 05_multiomics_integration/  MOFA2, DIABLO, SNF, WGCNA, GRN, network medicine
│   ├── 06_molecular_subtyping/     Consensus clustering and characterisation
│   ├── 07_biomarker_discovery/     Elastic-net, ML ensemble, validation
│   ├── 08_visualization/           Shared visual primitives (panel factories)
│   ├── 09_causal_network_medicine/ Random walk with restart, drug-target docking I/O
│   └── utils/                      Themes, layouts, helper functions
├── figures/                        Per-figure standalone scripts
│   ├── main/                       Figure_1_standalone.R … Figure_6_standalone.R
│   └── supplementary/              SuppFig_01_standalone.R … SuppFig_16_standalone.R
├── enhancements/                   v1.1.0+v1.2.0 extension analyses (Supplementary
│                                   Figures 16–32, Supplementary Tables 44–51)
├── environment.yml                 Conda environment specification
├── install_packages.R              Bioconductor / CRAN installer
├── LICENSE
└── README.md
```

Each script under `figures/` is fully self-contained: `library()` calls, palette
constants, helper functions, data loading, panel construction, and final
composition all live in one file. Running any standalone script with the
`multiomics` conda environment reproduces the corresponding figure end to end.

---

## Environment

The complete software environment is captured in `environment.yml`. To recreate
it on a new machine:

```bash
conda env create -f environment.yml
conda activate multiomics
Rscript install_packages.R
```

Tested under macOS 14 (Apple Silicon) and Linux x86_64 with R 4.5.0 and
Python 3.12.

---

## Reproducing a single figure

```bash
conda activate multiomics
Rscript figures/main/Figure_4/Figure_4_standalone.R
Rscript figures/supplementary/SuppFig_07/SuppFig_07_standalone.R
```

Output files (`.pdf`, `.png`, `.tiff`) are written to the same directory as the
script.

## Running an analysis stage

The numbered modules under `R/` follow the order in which they appear in the
Methods section. They can be executed individually:

```bash
Rscript R/01_differential_analysis/phase1_differential.R
Rscript R/04_immune_microenvironment/phase4_immune.R
Rscript R/05_multiomics_integration/phase5_integration.R
```

Most modules expect the processed data matrices produced by
`R/00_preprocessing/phase0_preprocessing.R`. Please refer to the manuscript
Methods for the exact processing parameters.

---

## Enhancement analyses (v1.1.0 + v1.2.0)

The `enhancements/` directory contains the extension analyses added during
revision, which produce Supplementary Figures 16–32 and Supplementary
Tables 44–51 (numbering follows the revised 2026-09 submission):

| Script | Output |
|--------|--------|
| `enhancement35_rifampicin_rescue.py` | Bayesian rifampicin-rescue simulation (Supp. Fig. 27, 29–31, ST48, ST50) |
| `enhancement36_mediation.py` | Patient-level paired-delta mediation analysis (Supp. Fig. 32) |
| `enhancement37_lincs_cmap.R` / `enhancement37_lincs_cmap_enrichr.py` | LINCS L1000 / CMap drug-reversal scoring (ST45) |
| `enhancement38_liver_atlas_deconv.py` / `enhancement39_cross_deconv.R` | Cross-atlas cell-type deconvolution (Supp. Fig. 13, 16, ST21) |
| `enhancement41_clinical_anchor.py` | Clinical-severity anchoring heatmap (Supp. Fig. 7, ST27) |
| `enhancement42_dorothea_progeny.R` / `enhancement42b_progeny_isolated.R` | DoRothEA / PROGENy pathway activity (Supp. Fig. 18, ST26) |
| `enhancement43_lincs_crosslib.py` | LINCS cross-library validation (ST46) |
| `enhancement44_celltype_clinical.py` / `enhancement45_tf_clinical.py` / `enhancement46_ges_validation.py` | Cell-type / TF clinical correlations and GES validation |
| `enhancement47_gtex_genomewide.py` | GTEx v10 genome-wide liver validation, 262 donors (Supp. Fig. 22, ST36) |
| `enhancement48_pbpk_cyp3a4.R` | CYP3A4-mediated PBPK mechanistic model (Supp. Fig. 28, ST49) |
| `enhancement51_zonation_composition.py` | Zonation composition-adjustment sensitivity analysis (Supp. Fig. 30) |
| `enhancement52_human_subgroup_meta.py` | Human subgroup random-effects meta-analysis (Supp. Fig. 31, ST50) |

---

## Data availability

Raw sequencing reads, mass-spectrometry spectra, and clinical metadata are
deposited in publicly accessible repositories (accession numbers are listed in
the Data availability statement of the manuscript). The processed,
de-identified data tables required to rerun every analysis script are provided
as Supplementary Tables 1–51 alongside the published article.

---

## How to cite

If this code is useful for your work, please cite the original article. The
final citation will be added once the article is in press.

---

## License

Released under the MIT License. See `LICENSE` for details.

---

## Release notes

- **v1.2.1 (2026-09-27)**: completed the removal signposted in v1.2.0 — the five Mendelian randomisation scripts under `R/07_biomarker_discovery/` are deleted; `R/09_causal_network_medicine` re-scoped to network medicine and drug repurposing; the stale `Figure_7` / `Figure_8` (MR) and `SuppFig_13` (MR extended) layout blocks removed from `R/utils/figure_layouts.R`; `enhancement47` and `enhancement48` relabelled to their current figure numbers (22 and 28); `enhancement31`–`enhancement34` added (Bayesian MCMC pharmacokinetics, healthy-liver reference, evidence synthesis, real-world consistency).
- **v1.2.0 (2026-09-26)**: sync with the revised EBioMedicine submission — Mendelian
  randomisation analyses removed throughout; main figures consolidated from 7 to 6
  (former Figure 7 promoted to Figure 6, former Figure 6 moved to Supplementary
  Fig. 19); supplementary figure panel tags standardised to uppercase; new
  enhancement analyses 51 (zonation composition adjustment) and 52 (human subgroup
  meta-analysis); supplementary tables renumbered to ST1-51. Manuscript: 58 references.
- **v1.1.0**: sync with EBioMedicine submission (7 main + 27 supplementary figures,
  ST1-52).
