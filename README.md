# HAE-Multiomics-Analysis

Source code accompanying the manuscript:

> **A pharmacological dead zone underlies albendazole failure in hepatic alveolar echinococcosis**

This repository contains the R analysis pipeline and figure-generation scripts used to produce all main and supplementary figures of the study. The pipeline integrates paired bulk transcriptomic, proteomic, and metabolomic data from 14 patients with hepatic alveolar echinococcosis (HAE), together with cell-type deconvolution, gene regulatory network inference, multi-omics integration, biomarker discovery, Mendelian randomisation, and quantitative pharmacological modelling.

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
│   ├── 07_biomarker_discovery/     Elastic-net, ML ensemble, MR, validation
│   ├── 08_visualization/           Shared visual primitives (panel factories)
│   ├── 09_causal_network_medicine/ Random walk with restart, drug-target docking I/O
│   └── utils/                      Themes, layouts, helper functions
├── figures/                        Per-figure standalone scripts
│   ├── main/                       Figure_1_standalone.R … Figure_10_standalone.R
│   └── supplementary/              SuppFig_01_standalone.R … SuppFig_15_standalone.R
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

## Data availability

Raw sequencing reads, mass-spectrometry spectra, and clinical metadata are
deposited in publicly accessible repositories (accession numbers are listed in
the Data availability statement of the manuscript). The processed,
de-identified data tables required to rerun every analysis script are provided
as Supplementary Tables 1–41 alongside the published article.

---

## How to cite

If this code is useful for your work, please cite the original article. The
final citation will be added once the article is in press.

---

## License

Released under the MIT License. See `LICENSE` for details.
