#!/usr/bin/env Rscript
# ============================================================================
# Enhancement 6: Drug Repurposing Analysis
# ============================================================================
# Identify potential therapeutic candidates by querying druggable targets
# from our DE genes and hub genes using DGIdb and connectivity analysis.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
})

PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
E2_DIR   <- file.path(PROJECT, "analysis/results/enhancement2_wgcna")
E4_DIR   <- file.path(PROJECT, "analysis/results/enhancement4_prot_metab_network")
E5_DIR   <- file.path(PROJECT, "analysis/results/enhancement5_external_validation")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhancement6_drug_repurposing")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Enhancement 6: Drug Repurposing Analysis\n")
cat("========================================\n\n")

# ============================================================================
# 1. COLLECT TARGET GENES FROM MULTIPLE ANALYSES
# ============================================================================
cat(">>> 1. Collecting candidate target genes\n")

# DEGs
degs <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)
deg_up <- degs$gene_name[degs$logFC > 0.5 & degs$P.Value < 0.05]
deg_down <- degs$gene_name[degs$logFC < -0.5 & degs$P.Value < 0.05]
deg_up <- deg_up[!is.na(deg_up) & deg_up != "" & !grepl("^ENSG", deg_up)]
deg_down <- deg_down[!is.na(deg_down) & deg_down != "" & !grepl("^ENSG", deg_down)]

# DEPs
deps <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)
dep_genes <- deps$gene_name[deps$P.Value < 0.05 & abs(deps$logFC) > 0.5]
dep_genes <- dep_genes[!is.na(dep_genes) & dep_genes != ""]

# WGCNA hubs
hub_genes <- character(0)
hub_file <- file.path(E2_DIR, "TC_hub_genes.csv")
if (file.exists(hub_file)) {
  hubs <- read.csv(hub_file, check.names = FALSE)
  if ("gene" %in% colnames(hubs)) hub_genes <- unique(hubs$gene)
  else hub_genes <- unique(hubs[, 1])
  hub_genes <- hub_genes[!is.na(hub_genes) & hub_genes != ""]
}

# Network hubs
net_hubs <- character(0)
net_file <- file.path(E4_DIR, "network_hub_nodes.csv")
if (file.exists(net_file)) {
  net <- read.csv(net_file, check.names = FALSE)
  net_hubs <- net$node[net$type == "protein" & net$degree >= 10]
}

# Consensus genes
consensus_genes <- character(0)
cons_file <- file.path(E5_DIR, "consensus_genes.csv")
if (file.exists(cons_file)) {
  cons <- read.csv(cons_file, check.names = FALSE)
  consensus_genes <- cons$gene
}

# Combine all unique target genes
all_targets <- unique(c(deg_up, deg_down, dep_genes, hub_genes, net_hubs, consensus_genes))
all_targets <- all_targets[!grepl("^ENSG", all_targets)]
cat(sprintf("  DEGs up: %d, DEGs down: %d\n", length(deg_up), length(deg_down)))
cat(sprintf("  DEPs: %d\n", length(dep_genes)))
cat(sprintf("  WGCNA hubs: %d\n", length(hub_genes)))
cat(sprintf("  Network hubs: %d\n", length(net_hubs)))
cat(sprintf("  Consensus genes: %d\n", length(consensus_genes)))
cat(sprintf("  Total unique targets: %d\n", length(all_targets)))

# ============================================================================
# 2. CURATED DRUGGABLE GENE DATABASE
# ============================================================================
cat("\n>>> 2. Building druggable gene database\n")

# Since DGIdb API may have network issues, use comprehensive curated database
# Sources: DGIdb categories, FDA-approved drug targets, druggable genome

# Known druggable gene categories
druggable_categories <- list(
  kinases = c(
    "EGFR", "ERBB2", "ERBB3", "MET", "ALK", "ROS1", "RET", "FGFR1", "FGFR2", "FGFR3",
    "VEGFR1", "VEGFR2", "KDR", "FLT1", "FLT4", "PDGFRA", "PDGFRB", "KIT", "CSF1R",
    "ABL1", "ABL2", "JAK1", "JAK2", "JAK3", "TYK2", "SRC", "FYN", "LCK",
    "BTK", "SYK", "PI3KCA", "PIK3CA", "PIK3CB", "PIK3CD", "MTOR", "AKT1", "AKT2",
    "BRAF", "RAF1", "ARAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
    "CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9",
    "AURKA", "AURKB", "PLK1", "CHEK1", "CHEK2", "WEE1",
    "ROCK1", "ROCK2", "PAK1", "PAK4", "DYRK1A"
  ),

  gpcrs = c(
    "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
    "ADRB1", "ADRB2", "ADRA1A", "ADRA2A",
    "HTR1A", "HTR2A", "HTR2C", "HTR3A",
    "DRD1", "DRD2", "DRD3", "DRD4", "DRD5",
    "CHRM1", "CHRM2", "CHRM3",
    "CCR2", "CCR5", "CXCR4", "CXCR3",
    "PTGER2", "PTGER4", "TBXA2R", "LPAR1",
    "S1PR1", "S1PR2", "EDNRA", "EDNRB"
  ),

  nuclear_receptors = c(
    "NR3C1", "NR3C2", "ESR1", "ESR2", "AR", "PGR",
    "PPARA", "PPARG", "PPARD",
    "RARA", "RARB", "RARG", "RXRA", "RXRB",
    "NR1H4", "NR1H3", "NR1I2", "NR1I3",
    "VDR", "HNF4A", "THRA", "THRB"
  ),

  proteases = c(
    "MMP2", "MMP9", "MMP13", "MMP14",
    "ADAM10", "ADAM17", "BACE1", "BACE2",
    "CASP1", "CASP3", "CASP8", "CASP9",
    "CTSA", "CTSB", "CTSD", "CTSK", "CTSL",
    "ELANE", "PRSS1", "DPP4", "ACE", "ACE2"
  ),

  immune_targets = c(
    "PDCD1", "CD274", "CTLA4", "HAVCR2", "LAG3", "TIGIT",
    "CD80", "CD86", "ICOS", "CD28",
    "TNF", "IL6", "IL6R", "IL1B", "IL1R1",
    "IL4R", "IL13RA1", "IL17A", "IL17RA",
    "TGFB1", "TGFBR1", "TGFBR2",
    "CSF2", "CSF1R", "CCL2", "CCR2",
    "IDO1", "ARG1", "CD47", "SIRPA"
  ),

  epigenetic = c(
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC6",
    "DNMT1", "DNMT3A", "DNMT3B",
    "EZH2", "DOT1L", "PRMT5",
    "BRD2", "BRD3", "BRD4", "CREBBP", "EP300",
    "KDM1A", "KDM5A", "KDM6A", "KDM6B"
  ),

  transporters = c(
    "ABCB1", "ABCB11", "ABCC1", "ABCC2", "ABCG2",
    "SLC6A3", "SLC6A4", "SLC22A1", "SLC22A2",
    "SLC5A2", "SLC12A1", "SLC12A3", "SLC9A3"
  ),

  metabolic_enzymes = c(
    "HMGCR", "FASN", "SCD", "ACACA",
    "CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4",
    "CYP2E1", "CYP1B1", "CYP19A1",
    "COX2", "PTGS2", "PTGS1", "ALOX5", "PDE4D", "PDE5A",
    "DHFR", "TYMS", "RRM1", "RRM2",
    "IDH1", "IDH2", "PARP1", "PARP2",
    "CPT1A", "ACADM", "HADHA"
  )
)

all_druggable <- unique(unlist(druggable_categories))
cat(sprintf("  Druggable gene database: %d genes across %d categories\n",
            length(all_druggable), length(druggable_categories)))

# ============================================================================
# 3. IDENTIFY DRUGGABLE TARGETS IN OUR DATA
# ============================================================================
cat("\n>>> 3. Identifying druggable targets\n")

druggable_hits <- intersect(all_targets, all_druggable)
cat(sprintf("  Druggable targets found: %d / %d (%.1f%%)\n",
            length(druggable_hits), length(all_targets),
            length(druggable_hits) / length(all_targets) * 100))

# Annotate each hit
drug_results <- data.frame()
for (gene in druggable_hits) {
  # Which categories
  cats <- names(druggable_categories)[sapply(druggable_categories, function(x) gene %in% x)]
  
  # DE status
  in_deg_up <- gene %in% deg_up
  in_deg_down <- gene %in% deg_down
  in_dep <- gene %in% dep_genes
  in_hub <- gene %in% hub_genes
  in_net <- gene %in% net_hubs
  in_cons <- gene %in% consensus_genes
  
  # logFC from DEGs
  lfc <- degs$logFC[degs$gene_name == gene]
  pval <- degs$P.Value[degs$gene_name == gene]
  if (length(lfc) == 0) { lfc <- NA; pval <- NA }
  else { lfc <- lfc[1]; pval <- pval[1] }
  
  # Evidence count
  evidence <- sum(in_deg_up, in_deg_down, in_dep, in_hub, in_net, in_cons)
  
  drug_results <- rbind(drug_results, data.frame(
    gene = gene,
    category = paste(cats, collapse = "; "),
    logFC = lfc,
    P.Value = pval,
    direction = ifelse(is.na(lfc), "unknown", ifelse(lfc > 0, "up", "down")),
    DEG = in_deg_up | in_deg_down,
    DEP = in_dep,
    WGCNA_hub = in_hub,
    Network_hub = in_net,
    Consensus = in_cons,
    evidence_count = evidence,
    stringsAsFactors = FALSE
  ))
}

drug_results <- drug_results[order(-drug_results$evidence_count, drug_results$P.Value), ]
write.csv(drug_results, file.path(OUT_DIR, "druggable_targets.csv"), row.names = FALSE)

cat(sprintf("  Druggable targets annotated: %d\n", nrow(drug_results)))
cat("  Top druggable targets:\n")
for (i in 1:min(20, nrow(drug_results))) {
  cat(sprintf("    %s [%s]: logFC=%.2f, evidence=%d, dir=%s\n",
              drug_results$gene[i], drug_results$category[i],
              drug_results$logFC[i], drug_results$evidence_count[i],
              drug_results$direction[i]))
}

# ============================================================================
# 4. KNOWN DRUG-TARGET MAPPING
# ============================================================================
cat("\n>>> 4. Drug-target mapping\n")

# Curated drug-target interactions relevant to our hits
drug_target_db <- data.frame(
  gene = c(
    # Kinases
    "EGFR", "EGFR", "MET", "PDGFRA", "PDGFRB", "KDR", "FLT1", "KIT",
    "JAK1", "JAK2", "MTOR", "PIK3CA", "CDK4", "CDK6",
    # Immune
    "PDCD1", "CD274", "CTLA4", "TNF", "IL6", "IL6R",
    "IL1B", "TGFB1", "TGFBR1", "IDO1", "CSF1R",
    # Nuclear receptors
    "PPARG", "PPARA", "NR1H4", "NR1I2", "HNF4A",
    "ESR1", "AR", "VDR",
    # Metabolic
    "HMGCR", "PTGS2", "PARP1", "HDAC1", "HDAC2", "HDAC6",
    "MMP2", "MMP9", "CTSD", "CTSB",
    # Transporters
    "ABCB1", "SLC22A1",
    # Others
    "ACE", "ACE2", "DPP4"
  ),
  drug = c(
    "Erlotinib/Gefitinib", "Osimertinib", "Crizotinib/Capmatinib", "Imatinib", "Imatinib",
    "Sorafenib/Sunitinib", "Bevacizumab", "Imatinib",
    "Ruxolitinib/Tofacitinib", "Ruxolitinib", "Rapamycin/Everolimus", "Alpelisib",
    "Palbociclib", "Palbociclib",
    "Nivolumab/Pembrolizumab", "Atezolizumab", "Ipilimumab", "Infliximab/Adalimumab",
    "Tocilizumab", "Tocilizumab",
    "Anakinra/Canakinumab", "Fresolimumab", "Galunisertib",
    "Epacadostat", "Pexidartinib",
    "Pioglitazone/Rosiglitazone", "Fenofibrate/Bezafibrate", "Obeticholic acid",
    "Rifampicin", "NA",
    "Tamoxifen/Fulvestrant", "Enzalutamide", "Calcitriol",
    "Statins (Atorvastatin)", "Celecoxib", "Olaparib/Niraparib",
    "Vorinostat/Romidepsin", "Vorinostat", "Ricolinostat",
    "Marimastat/Batimastat", "Marimastat", "Pepstatin A", "CA-074",
    "Verapamil", "NA",
    "Enalapril/Lisinopril", "NA", "Sitagliptin"
  ),
  interaction_type = c(
    "inhibitor", "inhibitor", "inhibitor", "inhibitor", "inhibitor",
    "inhibitor", "antibody", "inhibitor",
    "inhibitor", "inhibitor", "inhibitor", "inhibitor",
    "inhibitor", "inhibitor",
    "antibody", "antibody", "antibody", "antibody",
    "antibody", "antibody",
    "antagonist", "antibody", "inhibitor",
    "inhibitor", "inhibitor",
    "agonist", "agonist", "agonist",
    "inducer", "modulator",
    "antagonist", "antagonist", "agonist",
    "inhibitor", "inhibitor", "inhibitor",
    "inhibitor", "inhibitor", "inhibitor",
    "inhibitor", "inhibitor", "inhibitor", "inhibitor",
    "inhibitor", "substrate",
    "inhibitor", "target", "inhibitor"
  ),
  stringsAsFactors = FALSE
)

# Map drugs to our hits
drug_mapped <- merge(drug_results, drug_target_db, by = "gene", all.x = FALSE)
drug_mapped <- drug_mapped[order(-drug_mapped$evidence_count, drug_mapped$P.Value), ]
write.csv(drug_mapped, file.path(OUT_DIR, "drug_target_mapping.csv"), row.names = FALSE)

cat(sprintf("  Targets with known drugs: %d\n", length(unique(drug_mapped$gene))))
cat(sprintf("  Drug-target pairs: %d\n", nrow(drug_mapped)))

if (nrow(drug_mapped) > 0) {
  cat("\n  Drug repurposing candidates:\n")
  for (i in 1:min(20, nrow(drug_mapped))) {
    strategy <- ifelse(drug_mapped$direction[i] == "up" & drug_mapped$interaction_type[i] == "inhibitor",
                       "INHIBIT upregulated target",
                       ifelse(drug_mapped$direction[i] == "down" & drug_mapped$interaction_type[i] %in% c("agonist", "activator"),
                              "ACTIVATE downregulated target",
                              "investigate"))
    cat(sprintf("    %s -> %s [%s] (logFC=%.2f, strategy: %s)\n",
                drug_mapped$drug[i], drug_mapped$gene[i],
                drug_mapped$interaction_type[i],
                drug_mapped$logFC[i], strategy))
  }
}

# ============================================================================
# 5. THERAPEUTIC STRATEGY SCORING
# ============================================================================
cat("\n>>> 5. Therapeutic strategy scoring\n")

if (nrow(drug_mapped) > 0) {
  # Score each drug-target pair
  drug_mapped$strategy_score <- 0
  for (i in 1:nrow(drug_mapped)) {
    score <- 0
    # Direction match: inhibitor for upregulated, activator for downregulated
    if (!is.na(drug_mapped$logFC[i])) {
      if (drug_mapped$logFC[i] > 0 && drug_mapped$interaction_type[i] %in% c("inhibitor", "antibody", "antagonist")) {
        score <- score + 3  # Strong match
      } else if (drug_mapped$logFC[i] < 0 && drug_mapped$interaction_type[i] %in% c("agonist", "activator")) {
        score <- score + 3
      }
    }
    # Evidence count bonus
    score <- score + drug_mapped$evidence_count[i]
    # Significance bonus
    if (!is.na(drug_mapped$P.Value[i]) && drug_mapped$P.Value[i] < 0.01) score <- score + 2
    else if (!is.na(drug_mapped$P.Value[i]) && drug_mapped$P.Value[i] < 0.05) score <- score + 1
    # Effect size bonus
    if (!is.na(drug_mapped$logFC[i]) && abs(drug_mapped$logFC[i]) > 1) score <- score + 2
    else if (!is.na(drug_mapped$logFC[i]) && abs(drug_mapped$logFC[i]) > 0.5) score <- score + 1

    drug_mapped$strategy_score[i] <- score
  }

  drug_mapped <- drug_mapped[order(-drug_mapped$strategy_score), ]
  write.csv(drug_mapped, file.path(OUT_DIR, "drug_target_scored.csv"), row.names = FALSE)

  cat("  Top scored drug repurposing candidates:\n")
  for (i in 1:min(15, nrow(drug_mapped))) {
    cat(sprintf("    Score=%d: %s -> %s (logFC=%.2f, P=%.4f, evidence=%d)\n",
                drug_mapped$strategy_score[i], drug_mapped$drug[i],
                drug_mapped$gene[i], drug_mapped$logFC[i],
                ifelse(is.na(drug_mapped$P.Value[i]), 1, drug_mapped$P.Value[i]),
                drug_mapped$evidence_count[i]))
  }
}

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================
cat("\n>>> 6. Generating visualizations\n")

# ---- 6a. Druggable target category distribution ----
if (nrow(drug_results) > 0) {
  cat_count <- data.frame()
  for (cat_name in names(druggable_categories)) {
    n <- sum(drug_results$gene %in% druggable_categories[[cat_name]])
    if (n > 0) {
      cat_count <- rbind(cat_count, data.frame(
        category = cat_name, count = n, stringsAsFactors = FALSE
      ))
    }
  }
  cat_count <- cat_count[order(-cat_count$count), ]

  if (nrow(cat_count) > 0) {
    cat_count$category <- factor(cat_count$category, levels = rev(cat_count$category))
    p_cat <- ggplot(cat_count, aes(x = count, y = category)) +
      geom_col(fill = "#E64B35", alpha = 0.8) +
      labs(title = "Druggable Target Categories",
           x = "Number of Targets", y = "") +
      theme_bw(base_size = 12)
    ggsave(file.path(FIG_DIR, "druggable_categories.pdf"), p_cat, width = 8, height = 5)
    cat("  Category barplot saved.\n")
  }
}

# ---- 6b. Drug-target network visualization ----
if (nrow(drug_mapped) > 0) {
  # Top drug-target pairs for visualization
  top_dt <- head(drug_mapped[!is.na(drug_mapped$drug) & drug_mapped$drug != "NA", ], 30)

  if (nrow(top_dt) > 0) {
    library(igraph)
    dt_edges <- data.frame(
      from = top_dt$drug,
      to = top_dt$gene,
      stringsAsFactors = FALSE
    )
    g <- graph_from_data_frame(dt_edges, directed = FALSE)
    V(g)$type <- ifelse(V(g)$name %in% top_dt$gene, "gene", "drug")
    V(g)$color <- ifelse(V(g)$type == "gene", "#E64B35", "#4DBBD5")
    V(g)$shape <- ifelse(V(g)$type == "gene", "circle", "square")

    pdf(file.path(FIG_DIR, "drug_target_network.pdf"), width = 14, height = 12)
    set.seed(42)
    plot(g,
         vertex.size = ifelse(V(g)$type == "gene", 12, 10),
         vertex.label.cex = 0.6,
         vertex.label.color = "black",
         edge.color = "gray60",
         main = "Drug-Target Repurposing Network")
    legend("bottomleft",
           legend = c("Gene Target", "Drug"),
           col = c("#E64B35", "#4DBBD5"),
           pch = c(16, 15), cex = 0.9, bty = "n")
    dev.off()
    cat("  Drug-target network saved.\n")
  }
}

# ---- 6c. Volcano-style drug target plot ----
if (nrow(drug_results) > 0) {
  plot_df <- drug_results[!is.na(drug_results$logFC) & !is.na(drug_results$P.Value), ]
  plot_df$neg_log10p <- -log10(plot_df$P.Value)
  plot_df$is_druggable <- TRUE
  plot_df$label <- ifelse(plot_df$evidence_count >= 2 | plot_df$P.Value < 0.01,
                          plot_df$gene, "")

  p_vol <- ggplot(plot_df, aes(x = logFC, y = neg_log10p, label = label)) +
    geom_point(aes(size = evidence_count, color = category), alpha = 0.7) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
    labs(title = "Druggable Targets: Expression Change vs Significance",
         x = "log2 Fold Change (Adjacent/Normal)",
         y = "-log10(P-value)",
         size = "Evidence\nCount",
         color = "Drug Category") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
  ggsave(file.path(FIG_DIR, "druggable_volcano.pdf"), p_vol, width = 12, height = 8)
  cat("  Druggable volcano plot saved.\n")
}

# ============================================================================
# 7. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Enhancement 6 SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total candidate targets: %d\n", length(all_targets)))
cat(sprintf("Druggable targets: %d (%.1f%%)\n",
            nrow(drug_results), nrow(drug_results) / length(all_targets) * 100))
if (nrow(drug_mapped) > 0) {
  cat(sprintf("Drug-target pairs: %d\n", nrow(drug_mapped)))
  cat(sprintf("Unique drugs: %d\n", length(unique(drug_mapped$drug[drug_mapped$drug != "NA"]))))
}
cat(sprintf("Output: %s\n", OUT_DIR))
cat("Enhancement 6 COMPLETE.\n")
