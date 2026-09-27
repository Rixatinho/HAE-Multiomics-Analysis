#!/usr/bin/env Rscript
# =============================================================================
# figure_layouts.R -- Slot Dimensions for All 21 Figures (R2重组版)
# =============================================================================
# Pre-calculated exact panel sizes (mm) so assembly does ZERO scaling.
# R2重组版: 8 Main Figures + 13 Supp Figures
# NC/J Hep limits: max_w = 183mm, max_h = 247mm, gap = 3mm
# Updated: 2026-03-27 (R2重组后)
# =============================================================================

# --- NC page constants ---
NC_MAX_W <- 183   # double-column width (mm)
NC_MAX_H <- 247   # max page height (mm)
NC_GAP   <- 3     # inter-panel gap (mm)

# --- Slot dimension helpers ---
slot_2col <- function(h) {
  # 2-column grid: each panel = (183-3)/2 = 90mm wide
  w <- (NC_MAX_W - NC_GAP) / 2
  list(w = w, h = h)
}

slot_3col <- function(h) {
  # 3-column grid: each panel = (183 - 2*3)/3 = 59mm wide
  w <- (NC_MAX_W - 2 * NC_GAP) / 3
  list(w = w, h = h)
}

slot_full <- function(h) {
  # Full-width panel
  list(w = NC_MAX_W, h = h)
}

# =============================================================================
# LAYOUTS -- each figure: list of panel entries
#   name: output PDF filename (without path)
#   w, h: exact dimensions in mm
# =============================================================================

LAYOUTS <- list(

  # ---- MAIN FIGURES (8张) ----

  # Figure 1: Study Design and Multi-omics Profiling (6 panels: a-f)
  Figure_1 = list(
    rows = list(
      list(
        list(name = "Fig1a_workflow.pdf",               w = NC_MAX_W, h = 80)
      ),
      list(
        list(name = "Fig1b_volcano_TC.pdf",             w = 90, h = 85),
        list(name = "Fig1c_volcano_PR.pdf",             w = 90, h = 85)
      ),
      list(
        list(name = "Fig1d_volcano_MT.pdf",             w = 90, h = 85),
        list(name = "Fig1e_top_DEGs.pdf",               w = 90, h = 85)
      ),
      list(
        list(name = "Fig1f_Hallmark_NES_heatmap.pdf",   w = NC_MAX_W, h = 65)
      )
    )
  ),

  # Figure 2: Metabolic Reprogramming (6 panels: a-f)
  Figure_2 = list(
    rows = list(
      list(
        list(name = "Fig2a_metabolic_bubble.pdf",       w = 90, h = 85),
        list(name = "Fig2b_enzyme_heatmap.pdf",         w = 90, h = 85)
      ),
      list(
        list(name = "Fig2c_top_DEMs_heatmap.pdf",       w = 90, h = 100),
        list(name = "Fig2d_msea_dotplot.pdf",           w = 90, h = 100)
      ),
      list(
        list(name = "Fig2e_prot_metab_network.pdf",     w = 90, h = 80),
        list(name = "Fig2f_GSEA_metab.pdf",             w = 90, h = 80)
      )
    )
  ),

  # Figure 3: Immune Microenvironment (5 panels: a-e)
  # 来源: old Fig3 panels a,b,d,g,h -> new a,b,c,d,e
  Figure_3 = list(
    rows = list(
      list(
        list(name = "Fig3a_ssGSEA_heatmap.pdf",         w = 90, h = 90),
        list(name = "Fig3b_immune_barplot.pdf",         w = 90, h = 90)
      ),
      list(
        list(name = "Fig3c_immune_NES.pdf",             w = 90, h = 80),
        list(name = "Fig3d_checkpoint_heatmap.pdf",     w = 90, h = 80)
      ),
      list(
        list(name = "Fig3e_immune_metab_coupling.pdf",  w = NC_MAX_W, h = 80)
      )
    )
  ),

  # Figure 4: Cell Communication and Master Regulators (7 panels: a-g)
  # Panel relabeling 2026-03-29: a=BayesPrism, b=CellChat, c=LR heatmap
  Figure_4 = list(
    rows = list(
      list(
        list(name = "Fig4a_stacked_barplot.pdf",          w = 90, h = 80),
        list(name = "Fig4b_cell_communication.pdf",       w = 90, h = 80)
      ),
      list(
        list(name = "Fig4c_LR_heatmap.pdf",              w = 90, h = 80),
        list(name = "Fig4d_LR_bubble.pdf",               w = 90, h = 80)
      ),
      list(
        list(name = "Fig4e_diversity_boxplot.pdf",        w = 90, h = 80),
        list(name = "Fig4f_master_regulator_network.pdf", w = 90, h = 80)
      ),
      list(
        list(name = "Fig4g_tf_centrality_ranking.pdf",   w = NC_MAX_W, h = 70)
      )
    )
  ),

  # Figure 5: Integrative Multi-omics Analysis (7 panels: a-g)
  # Panel relabeling 2026-03-29: d=MOFA pathway, e=MOFA clinical, f=variance, g=WGCNA
  Figure_5 = list(
    rows = list(
      list(
        list(name = "Fig5a_DIABLO_circos.pdf",           w = 90, h = 90),
        list(name = "Fig5b_DIABLO_samples.pdf",         w = 90, h = 90)
      ),
      list(
        list(name = "Fig5c_MOFA2_variance.pdf",         w = 90, h = 80),
        list(name = "Fig5d_factor_pathway_heatmap.pdf", w = 90, h = 80)
      ),
      list(
        list(name = "Fig5e_factor_clinical.pdf",        w = 59, h = 75),
        list(name = "Fig5f_variance_contribution.pdf",  w = 59, h = 75),
        list(name = "Fig5g_WGCNA_TC_trait.pdf",         w = 59, h = 75)
      )
    )
  ),

  # Figure 6: Consensus Molecular Subtyping (7 panels: a-g)
  # 来源: old Fig5 a-g -> new a-g
  Figure_6 = list(
    rows = list(
      list(
        list(name = "Fig6a_consensus_K2.pdf",            w = 90, h = 80),
        list(name = "Fig6b_silhouette_K2.pdf",           w = 90, h = 80)
      ),
      list(
        list(name = "Fig6c_subtype_PCA.pdf",             w = 90, h = 80),
        list(name = "Fig6d_clinical_boxplots_K2.pdf",    w = 90, h = 80)
      ),
      list(
        list(name = "Fig6e_subtype_profile.pdf",         w = NC_MAX_W, h = 85)
      ),
      list(
        list(name = "Fig6f_druggable_volcano.pdf",       w = 90, h = 80),
        list(name = "Fig6g_comprehensive_heatmap.pdf",   w = 90, h = 80)
      )
    )
  ),

  # ---- SUPPLEMENTARY FIGURES (13张) ----

  # Supp Fig 1: Quality Control (10 panels: a-j)
  SuppFig_01 = list(
    rows = list(
      list(
        list(name = "Supp1a_cor_TC.pdf",                w = 59, h = 90),
        list(name = "Supp1b_cor_PR.pdf",                w = 59, h = 90),
        list(name = "Supp1c_cor_MT.pdf",                w = 59, h = 90)
      ),
      list(
        list(name = "Supp1d_PCA_TC.pdf",                w = 59, h = 90),
        list(name = "Supp1e_PCA_PR.pdf",                w = 59, h = 90),
        list(name = "Supp1f_PCA_MT.pdf",                w = 59, h = 90)
      ),
      list(
        list(name = "Supp1g_PCA_combined.pdf",          w = 90, h = 90)
      ),
      list(
        list(name = "Supp1h_top_DEGs.pdf",              w = 59, h = 90),
        list(name = "Supp1i_top_DEPs.pdf",              w = 59, h = 90),
        list(name = "Supp1j_top_DEMs.pdf",              w = 59, h = 90)
      )
    )
  ),

  # Supp Fig 2: Extended Differential Analysis (6 panels: a-f)
  SuppFig_02 = list(
    rows = list(
      list(
        list(name = "Supp2a_top_DEPs.pdf",               w = 90, h = 100),
        list(name = "Supp2b_top_DEMs.pdf",               w = 90, h = 100)
      ),
      list(
        list(name = "Supp2c_mRNA_protein_correlation.pdf",w = 90, h = 80),
        list(name = "Supp2d_Reactome_GSEA_TC.pdf",       w = 90, h = 80)
      ),
      list(
        list(name = "Supp2e_Reactome_GSEA_PR.pdf",       w = 90, h = 80),
        list(name = "Supp2f_concordance.pdf",            w = 90, h = 80)
      )
    )
  ),

  # Supp Fig 3: Post-Transcriptional and KEGG Analysis (8 panels: a-h)
  # 来源: old Supp3 a,b + old Supp3 d,e,f->c,d,e + old Supp4 a,b,c->f,g,h
  SuppFig_03 = list(
    rows = list(
      list(
        list(name = "Supp3a_top_DEPs.pdf",               w = 90, h = 90),
        list(name = "Supp3b_top_DEMs.pdf",               w = 90, h = 90)
      ),
      list(
        list(name = "Supp3c_mRNA_protein_correlation.pdf",w = 90, h = 80),
        list(name = "Supp3d_external_validation.pdf",    w = 90, h = 80)
      ),
      list(
        list(name = "Supp3e_TC_KEGG_enrichment.pdf",     w = 59, h = 75),
        list(name = "Supp3f_PR_KEGG_enrichment.pdf",     w = 59, h = 75),
        list(name = "Supp3g_MT_KEGG_enrichment.pdf",     w = 59, h = 75)
      ),
      list(
        list(name = "Supp3h_mRNA_prot_cor_density.pdf",  w = 90, h = 80)
      )
    )
  ),

  # Supp Fig 4: TF and PROGENy Analysis (3 panels: a-c)
  # 来源: old Fig1 g,h,i
  SuppFig_04 = list(
    rows = list(
      list(
        list(name = "Supp4a_DoRothEA_TF.pdf",            w = 90, h = 90),
        list(name = "Supp4b_PROGENy_heatmap.pdf",        w = 90, h = 90)
      ),
      list(
        list(name = "Supp4c_PROGENy_boxplot.pdf",        w = NC_MAX_W, h = 80)
      )
    )
  ),

  # Supp Fig 5: Network and Venn Analysis (10 panels: a-j)
  # 来源: old Supp5 a,b + old Supp6 a-f->c,d,e,f,g,h + old Fig2 g,h->i,j
  SuppFig_05 = list(
    rows = list(
      list(
        list(name = "Supp5a_degree_dist.pdf",            w = 90, h = 80),
        list(name = "Supp5b_mRNA_prot_top_scatter.pdf",  w = 90, h = 80)
      ),
      list(
        list(name = "Supp5c_MT_KEGG_enrichment.pdf",     w = 90, h = 80),
        list(name = "Supp5d_immune_infiltration_boxplot.pdf",w = 90, h = 80)
      ),
      list(
        list(name = "Supp5e_Tcell_functional_boxplot.pdf",w = 90, h = 80),
        list(name = "Supp5f_GSEA_fibrosis.pdf",          w = 90, h = 80)
      ),
      list(
        list(name = "Supp5g_TME_score_boxplot.pdf",      w = 59, h = 75),
        list(name = "Supp5h_concordance.pdf",            w = 59, h = 75),
        list(name = "Supp5i_pathway_heatmap.pdf",        w = 59, h = 75)
      ),
      list(
        list(name = "Supp5j_Venn_Hallmark.pdf",          w = 90, h = 80),
        list(name = "Supp5k_GSEA_complement.pdf",        w = 90, h = 80)
      )
    )
  ),

  # Supp Fig 6: Extended Immune Analysis (7 panels: a-g)
  # 来源: old Supp8 a,b,c,d,h,j,o
  SuppFig_06 = list(
    rows = list(
      list(
        list(name = "Supp6a_checkpoint_expression.pdf",  w = 90, h = 75),
        list(name = "Supp6b_immune_infiltration_boxplot.pdf",w = 90, h = 75)
      ),
      list(
        list(name = "Supp6c_Tcell_functional_boxplot.pdf",w = 90, h = 75),
        list(name = "Supp6d_PCD_heatmap.pdf",            w = 90, h = 75)
      ),
      list(
        list(name = "Supp6e_ICP_heatmap.pdf",            w = 90, h = 75),
        list(name = "Supp6f_clinical_correlation.pdf",   w = 90, h = 75)
      ),
      list(
        list(name = "Supp6g_immune_correlation.pdf",      w = NC_MAX_W, h = 75)
      )
    )
  ),

  # Supp Fig 7: Cell Communication Extended (12 panels: a-l)
  # 来源: old Supp7 a,b,c,d + old Supp8 e,f,g,i,k,l,m,n
  SuppFig_07 = list(
    rows = list(
      list(
        list(name = "Supp7a_drug_categories.pdf",        w = 90, h = 75),
        list(name = "Supp7b_LR_heatmap.pdf",             w = 90, h = 75)
      ),
      list(
        list(name = "Supp7c_Venn_Hallmark.pdf",          w = 90, h = 75),
        list(name = "Supp7d_Venn_KEGG.pdf",              w = 90, h = 75)
      ),
      list(
        list(name = "Supp7e_LR_differential_dotplot.pdf",w = 90, h = 75),
        list(name = "Supp7f_boxplot_comparison.pdf",     w = 90, h = 75)
      ),
      list(
        list(name = "Supp7g_heatmap_proportions.pdf",    w = 90, h = 75),
        list(name = "Supp7h_TME_score_boxplot.pdf",      w = 90, h = 75)
      ),
      list(
        list(name = "Supp7i_pathway_LR_network.pdf",     w = 90, h = 75),
        list(name = "Supp7j_diversity_boxplot.pdf",      w = 90, h = 75)
      ),
      list(
        list(name = "Supp7k_diversity_boxplot.pdf",      w = 90, h = 75),
        list(name = "Supp7l_immune_metab_coupling.pdf",  w = 90, h = 75)
      )
    )
  ),

  # Supp Fig 8: Extended Multi-omics Integration (11 panels: a-k)
  # 来源: old Supp9 a-i + old Supp11 g->j + old Fig4 j->k
  SuppFig_08 = list(
    rows = list(
      list(
        list(name = "Supp8a_top_features.pdf",           w = 90, h = 80),
        list(name = "Supp8b_rf_importance.pdf",          w = 90, h = 80)
      ),
      list(
        list(name = "Supp8c_snf_mofa_confusion.pdf",     w = 90, h = 80),
        list(name = "Supp8d_mofa_scatter.pdf",           w = 90, h = 80)
      ),
      list(
        list(name = "Supp8e_MOFA_clustering.pdf",        w = 90, h = 75),
        list(name = "Supp8f_factor_pathway_heatmap.pdf", w = 90, h = 75)
      ),
      list(
        list(name = "Supp8g_condition_comparison.pdf",   w = 59, h = 75),
        list(name = "Supp8h_condition_comparison.pdf",   w = 59, h = 75),
        list(name = "Supp8i_feature_venn.pdf",           w = 59, h = 75)
      ),
      list(
        list(name = "Supp8j_tf_centrality_ranking.pdf",  w = 90, h = 75),
        list(name = "Supp8k_tf_target_heatmap.pdf",      w = 90, h = 75)
      )
    )
  ),

  # Supp Fig 9: Extended Subtyping Analysis (8 panels: a-h)
  # 来源: old Supp10 a,d,e,f,g,h,i,j->a,b,c,d,e,f,g,h
  SuppFig_09 = list(
    rows = list(
      list(
        list(name = "Supp9a_PAC_vs_K.pdf",               w = 90, h = 75),
        list(name = "Supp9b_consensus_heatmap.pdf",      w = 90, h = 75)
      ),
      list(
        list(name = "Supp9c_silhouette_plot.pdf",        w = 90, h = 75),
        list(name = "Supp9d_pac_comparison.pdf",         w = 90, h = 75)
      ),
      list(
        list(name = "Supp9e_silhouette_stability.pdf",   w = 90, h = 75),
        list(name = "Supp9f_gap_statistic.pdf",          w = 90, h = 75)
      ),
      list(
        list(name = "Supp9g_bootstrap_jaccard.pdf",      w = 90, h = 75),
        list(name = "Supp9h_bootstrap_stability.pdf",    w = 90, h = 75)
      )
    )
  ),

  # Supp Fig 10: Drug and Clustering Extended (7 panels: a-g)
  # 来源: old Supp10 b,c,k,l,m->a,b,c,d,e + old Fig5 h,i->f,g
  SuppFig_10 = list(
    rows = list(
      list(
        list(name = "Supp10a_drug_target_network.pdf",   w = 90, h = 75),
        list(name = "Supp10b_druggable_boxplot.pdf",     w = 90, h = 75)
      ),
      list(
        list(name = "Supp10c_drug_strategy.pdf",         w = 90, h = 75),
        list(name = "Supp10d_family_heatmap.pdf",        w = 90, h = 75)
      ),
      list(
        list(name = "Supp10e_silhouette_K2.pdf",         w = 90, h = 75),
        list(name = "Supp10f_heatmap_proportions.pdf",   w = 90, h = 75)
      ),
      list(
        list(name = "Supp10g_alluvial.pdf",              w = NC_MAX_W, h = 75)
      )
    )
  ),

  # Supp Fig 11: WGCNA and Validation (7 panels: a-g)
  # 来源: old Supp11 a,b,c,d,i,j,k
  SuppFig_11 = list(
    rows = list(
      list(
        list(name = "Supp11a_WGCNA_hubs.pdf",            w = 90, h = 70),
        list(name = "Supp11b_loocv_roc.pdf",             w = 90, h = 70)
      ),
      list(
        list(name = "Supp11c_direction_heatmap.pdf",     w = 90, h = 70),
        list(name = "Supp11d_concordance.pdf",           w = 90, h = 70)
      ),
      list(
        list(name = "Supp11e_dataset_overview.pdf",      w = 90, h = 70),
        list(name = "Supp11f_cluster_validation.pdf",    w = 90, h = 70)
      ),
      list(
        list(name = "Supp11g_bootstrap_stability.pdf",   w = NC_MAX_W, h = 70)
      )
    )
  ),

  # Supp Fig 12: Model Performance (6 panels: a-f)
  # 来源: old Supp11 e,l,m,n,o,p->a,b,c,d,e,f
  SuppFig_12 = list(
    rows = list(
      list(
        list(name = "Supp12a_roc_curves.pdf",            w = 90, h = 70),
        list(name = "Supp12b_importance_heatmap.pdf",    w = 90, h = 70)
      ),
      list(
        list(name = "Supp12c_performance_comparison.pdf",w = 90, h = 70),
        list(name = "Supp12d_importance_heatmap.pdf",    w = 90, h = 70)
      ),
      list(
        list(name = "Supp12e_external_concordance.pdf",  w = 90, h = 70),
        list(name = "Supp12f_cross_dataset_heatmap.pdf", w = 90, h = 70)
      )
    )
  ),

)

# =============================================================================
# Helper: get slot dimensions for a specific panel
# =============================================================================
get_panel_slot <- function(panel_name) {
  for (fig_name in names(LAYOUTS)) {
    for (row in LAYOUTS[[fig_name]]$rows) {
      for (entry in row) {
        if (entry$name == panel_name) {
          return(list(figure = fig_name, w = entry$w, h = entry$h))
        }
      }
    }
  }
  return(NULL)  # not found
}

# =============================================================================
# Helper: list all panel names across all figures
# =============================================================================
list_all_panels <- function() {
  panels <- character()
  for (fig_name in names(LAYOUTS)) {
    for (row in LAYOUTS[[fig_name]]$rows) {
      for (entry in row) {
        panels <- c(panels, entry$name)
      }
    }
  }
  panels
}

cat("[figure_layouts.R]", length(LAYOUTS), "figures,",
    length(list_all_panels()), "panel slots registered.\n")
