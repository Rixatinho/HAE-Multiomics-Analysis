#!/usr/bin/env Rscript
# ============================================================================
# Enhancement 4: Protein-Metabolite Correlation Network
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(igraph)
})

PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
PROC_DIR <- file.path(PROJECT, "analysis/data/processed")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
E1_DIR   <- file.path(PROJECT, "analysis/results/enhancement1_paired_subtyping")
OUT_DIR  <- file.path(PROJECT, "analysis/results/enhancement4_prot_metab_network")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Source nc_theme.R for helper functions (ensp_to_symbol, etc.)
NC_THEME <- file.path(PROJECT, "analysis/scripts/nc_theme.R")
if (file.exists(NC_THEME)) {
  source(NC_THEME)
  cat("  nc_theme.R loaded for ID mapping helpers\n")
}

cat("========================================\n")
cat("Enhancement 4: Protein-Metabolite Network\n")
cat("========================================\n\n")

# ---- 0. Load data ----
cat(">>> 0. Loading data\n")

# Proteomics
pr_mat <- as.matrix(read.csv(file.path(PROC_DIR, "proteomics_log2_norm.csv"),
                              check.names = FALSE, row.names = 1))
pr_diff <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

# Metabolomics
mb_mat <- as.matrix(read.csv(file.path(PROC_DIR, "metabolomics_log2_merged.csv"),
                              check.names = FALSE, row.names = 1))
mb_diff <- read.csv(file.path(DIFF_DIR, "DEMs_Adjacent_vs_Normal.csv"), check.names = FALSE)
mb_anno <- read.csv(file.path(PROC_DIR, "metabolomics_annotation.csv"), check.names = FALSE)

# Subtypes
paired_sub <- read.csv(file.path(E1_DIR, "paired_subtype_K2.csv"), stringsAsFactors = FALSE)
subtype_map <- setNames(paired_sub$subtype, paired_sub$patient_id)

# Align samples
common_samps <- intersect(colnames(pr_mat), colnames(mb_mat))
pr_mat <- pr_mat[, common_samps]
mb_mat <- mb_mat[, common_samps]

cat(sprintf("  Proteomics: %d proteins x %d samples\n", nrow(pr_mat), ncol(pr_mat)))
cat(sprintf("  Metabolomics: %d metabolites x %d samples\n", nrow(mb_mat), ncol(mb_mat)))

# ============================================================================
# 1. SELECT DIFFERENTIALLY EXPRESSED FEATURES
# ============================================================================
cat("\n>>> 1. Selecting DE features\n")

# DEPs: P < 0.05 & |logFC| > 0.5
sig_prots <- pr_diff$Protein[pr_diff$P.Value < 0.05 & abs(pr_diff$logFC) > 0.5]
sig_prots <- intersect(sig_prots, rownames(pr_mat))

# Build comprehensive ENSP → gene symbol mapping
prot2gene <- setNames(pr_diff$gene_name, pr_diff$Protein)
# Clean version without ENSP suffix for additional matching
prot2gene_clean <- setNames(pr_diff$gene_name, gsub("\\.[0-9]+$", "", pr_diff$Protein))

# Helper function to convert ENSP ID to gene symbol with fallback
ensp_to_gene <- function(ensp_id) {
  # Try direct match first
  sym <- prot2gene[ensp_id]
  if (!is.na(sym) && sym != "" && sym != "_") return(sym)
  # Try without version suffix
  ensp_clean <- gsub("\\.[0-9]+$", "", ensp_id)
  sym <- prot2gene_clean[ensp_clean]
  if (!is.na(sym) && sym != "" && sym != "_") return(sym)
  # If nc_theme.R provides ensp_to_symbol, use it
  if (exists("ensp_to_symbol")) {
    sym <- ensp_to_symbol(ensp_id)
    if (!is.na(sym) && sym != ensp_id && sym != "") return(sym)
  }
  # Return original ID if no mapping found
  return(ensp_id)
}

# DEMs: P < 0.05 & |logFC| > 0.5
sig_mets <- mb_diff$Compound_ID[mb_diff$P.Value < 0.05 & abs(mb_diff$logFC) > 0.5]
sig_mets <- intersect(sig_mets, rownames(mb_mat))

# Metabolite name map - also read from DEMs_significant.csv if available
met2name <- setNames(mb_anno$Name, mb_anno$Compound_ID)
met2class <- setNames(mb_anno$ClassI, mb_anno$Compound_ID)
met2kegg <- setNames(mb_anno$KEGG_ID, mb_anno$Compound_ID)

# Supplement with DEMs_significant.csv if it has better names
dems_sig_file <- file.path(DIFF_DIR, "DEMs_significant.csv")
if (file.exists(dems_sig_file)) {
  dems_sig <- read.csv(dems_sig_file, stringsAsFactors = FALSE)
  if ("metabolite_name" %in% colnames(dems_sig) && "Compound_ID" %in% colnames(dems_sig)) {
    extra_names <- setNames(dems_sig$metabolite_name, dems_sig$Compound_ID)
    for (cid in names(extra_names)) {
      if (is.na(met2name[cid]) || met2name[cid] == "" || met2name[cid] == cid) {
        met2name[cid] <- extra_names[cid]
      }
    }
  }
}

# Helper function to format metabolite names (truncate long names)
format_met_name <- function(met_id, max_len = 35) {
  nm <- met2name[met_id]
  if (is.na(nm) || nm == "" || nm == met_id) return(met_id)
  if (nchar(nm) > max_len) {
    return(paste0(substr(nm, 1, max_len - 3), "..."))
  }
  return(nm)
}

cat(sprintf("  Significant proteins: %d\n", length(sig_prots)))
cat(sprintf("  Significant metabolites: %d\n", length(sig_mets)))

# Cap features if too many (top by P-value)
MAX_PROTS <- 200
MAX_METS  <- 200
if (length(sig_prots) > MAX_PROTS) {
  top_pr <- pr_diff[pr_diff$Protein %in% sig_prots, ]
  top_pr <- top_pr[order(top_pr$P.Value), ]
  sig_prots <- head(top_pr$Protein, MAX_PROTS)
  cat(sprintf("  Capped to top %d proteins by P-value\n", MAX_PROTS))
}
if (length(sig_mets) > MAX_METS) {
  top_mb <- mb_diff[mb_diff$Compound_ID %in% sig_mets, ]
  top_mb <- top_mb[order(top_mb$P.Value), ]
  sig_mets <- head(top_mb$Compound_ID, MAX_METS)
  cat(sprintf("  Capped to top %d metabolites by P-value\n", MAX_METS))
}

pr_sub <- pr_mat[sig_prots, ]
mb_sub <- mb_mat[sig_mets, ]

# ============================================================================
# 2. COMPUTE SPEARMAN CORRELATION MATRIX
# ============================================================================
cat("\n>>> 2. Computing Spearman correlations (protein x metabolite)\n")

n_samp <- ncol(pr_sub)
cor_mat <- matrix(NA, nrow = length(sig_prots), ncol = length(sig_mets))
pval_mat <- matrix(NA, nrow = length(sig_prots), ncol = length(sig_mets))
rownames(cor_mat) <- rownames(pval_mat) <- sig_prots
colnames(cor_mat) <- colnames(pval_mat) <- sig_mets

for (i in seq_along(sig_prots)) {
  for (j in seq_along(sig_mets)) {
    ct <- cor.test(pr_sub[i, ], mb_sub[j, ], method = "spearman", exact = FALSE)
    cor_mat[i, j] <- ct$estimate
    pval_mat[i, j] <- ct$p.value
  }
}

# FDR correction across all tests
padj_vec <- p.adjust(as.vector(pval_mat), method = "BH")
padj_mat <- matrix(padj_vec, nrow = nrow(pval_mat), ncol = ncol(pval_mat))
rownames(padj_mat) <- rownames(cor_mat)
colnames(padj_mat) <- colnames(cor_mat)

cat(sprintf("  Total tests: %d\n", length(padj_vec)))
cat(sprintf("  Significant (P<0.05): %d\n", sum(pval_mat < 0.05, na.rm = TRUE)))
cat(sprintf("  Significant (padj<0.05): %d\n", sum(padj_mat < 0.05, na.rm = TRUE)))
cat(sprintf("  Strong (|rho|>0.6 & P<0.05): %d\n",
            sum(abs(cor_mat) > 0.6 & pval_mat < 0.05, na.rm = TRUE)))

# ============================================================================
# 3. BUILD EDGE LIST & NETWORK
# ============================================================================
cat("\n>>> 3. Building correlation network\n")

# Use |rho| > 0.6 & P < 0.05 as edge threshold
RHO_CUT <- 0.6
P_CUT   <- 0.05

edges <- data.frame()
for (i in seq_along(sig_prots)) {
  for (j in seq_along(sig_mets)) {
    if (!is.na(cor_mat[i, j]) && abs(cor_mat[i, j]) > RHO_CUT && pval_mat[i, j] < P_CUT) {
      # Use improved mapping functions
      prot_name <- ensp_to_gene(sig_prots[i])
      met_name_full <- met2name[sig_mets[j]]
      if (is.na(met_name_full) || met_name_full == "") met_name_full <- sig_mets[j]
      met_name_short <- format_met_name(sig_mets[j], 35)  # Short version for display
      
      edges <- rbind(edges, data.frame(
        protein_id = sig_prots[i],           # Original ENSP ID for reference
        protein = prot_name,                  # Gene symbol (for display)
        metabolite_id = sig_mets[j],          # Original Compound ID for reference
        metabolite = met_name_short,          # Metabolite name short (for display)
        metabolite_name_full = met_name_full, # Full metabolite name
        rho = cor_mat[i, j],
        pvalue = pval_mat[i, j],
        padj = padj_mat[i, j],
        direction = ifelse(cor_mat[i, j] > 0, "positive", "negative"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat(sprintf("  Network edges: %d (|rho|>%.1f, P<%.2f)\n", nrow(edges), RHO_CUT, P_CUT))

# If too few edges, relax threshold
if (nrow(edges) < 50) {
  RHO_CUT <- 0.5
  cat(sprintf("  Relaxing threshold to |rho|>%.1f...\n", RHO_CUT))
  edges <- data.frame()
  for (i in seq_along(sig_prots)) {
    for (j in seq_along(sig_mets)) {
      if (!is.na(cor_mat[i, j]) && abs(cor_mat[i, j]) > RHO_CUT && pval_mat[i, j] < P_CUT) {
        prot_name <- ensp_to_gene(sig_prots[i])
        met_name_full <- met2name[sig_mets[j]]
        if (is.na(met_name_full) || met_name_full == "") met_name_full <- sig_mets[j]
        met_name_short <- format_met_name(sig_mets[j], 35)
        
        edges <- rbind(edges, data.frame(
          protein_id = sig_prots[i],
          protein = prot_name,
          metabolite_id = sig_mets[j],
          metabolite = met_name_short,
          metabolite_name_full = met_name_full,
          rho = cor_mat[i, j],
          pvalue = pval_mat[i, j],
          padj = padj_mat[i, j],
          direction = ifelse(cor_mat[i, j] > 0, "positive", "negative"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  cat(sprintf("  Edges after relaxation: %d\n", nrow(edges)))
}

write.csv(edges, file.path(OUT_DIR, "prot_metab_edges.csv"), row.names = FALSE)

# Also save a legacy format for compatibility (protein_name, metabolite_name columns)
edges_legacy <- edges
colnames(edges_legacy)[colnames(edges_legacy) == "protein"] <- "protein_name"
colnames(edges_legacy)[colnames(edges_legacy) == "metabolite"] <- "metabolite_name"
colnames(edges_legacy)[colnames(edges_legacy) == "protein_id"] <- "protein"
colnames(edges_legacy)[colnames(edges_legacy) == "metabolite_id"] <- "metabolite"
# Reorder columns to match expected format
edges_legacy <- edges_legacy[, c("protein", "protein_name", "metabolite", "metabolite_name", 
                                  "metabolite_name_full", "rho", "pvalue", "padj", "direction")]
write.csv(edges_legacy, file.path(OUT_DIR, "prot_metab_edges_legacy.csv"), row.names = FALSE)

# Unique nodes (use gene symbols now)
prot_nodes <- unique(edges$protein)
met_nodes  <- unique(edges$metabolite)
cat(sprintf("  Proteins in network: %d\n", length(prot_nodes)))
cat(sprintf("  Metabolites in network: %d\n", length(met_nodes)))

# ============================================================================
# 4. NETWORK ANALYSIS (igraph)
# ============================================================================
cat("\n>>> 4. Network topology analysis\n")

if (nrow(edges) > 0) {
  # Build bipartite graph using readable names (gene symbols / metabolite names)
  g <- graph_from_data_frame(
    edges[, c("protein", "metabolite", "rho", "direction")],
    directed = FALSE
  )

  # Set node types
  V(g)$type <- ifelse(V(g)$name %in% edges$protein, "protein", "metabolite")
  V(g)$color <- ifelse(V(g)$type == "protein", "#E64B35", "#4DBBD5")
  V(g)$shape <- ifelse(V(g)$type == "protein", "circle", "square")
  E(g)$color <- ifelse(E(g)$direction == "positive", "#E64B3580", "#4DBBD580")
  E(g)$width <- abs(E(g)$rho) * 3

  # Degree, betweenness (with robust type conversion)
  deg <- igraph::degree(g)
  btw <- igraph::betweenness(g, normalized = TRUE)
  # Ensure numeric vectors (igraph 2.x returns proper numeric)
  deg <- as.numeric(deg)
  btw <- as.numeric(btw)
  node_names <- V(g)$name
  names(deg) <- node_names
  names(btw) <- node_names

  # Hub nodes (top degree)
  hub_df <- data.frame(
    node = V(g)$name,
    type = V(g)$type,
    degree = deg,
    betweenness = btw,
    stringsAsFactors = FALSE
  )
  hub_df <- hub_df[order(-hub_df$degree), ]
  write.csv(hub_df, file.path(OUT_DIR, "network_hub_nodes.csv"), row.names = FALSE)

  cat(sprintf("  Total nodes: %d, Total edges: %d\n", vcount(g), ecount(g)))
  cat(sprintf("  Network density: %.4f\n", edge_density(g)))
  cat(sprintf("  Components: %d\n", count_components(g)))
  cat(sprintf("  Mean degree: %.2f\n", mean(deg)))

  # Top hub proteins
  hub_prots <- hub_df[hub_df$type == "protein", ]
  cat("\n  Top hub proteins:\n")
  for (i in 1:min(10, nrow(hub_prots))) {
    cat(sprintf("    %s: degree=%d, betweenness=%.4f\n",
                hub_prots$node[i], hub_prots$degree[i], hub_prots$betweenness[i]))
  }

  # Top hub metabolites
  hub_mets <- hub_df[hub_df$type == "metabolite", ]
  cat("\n  Top hub metabolites:\n")
  for (i in 1:min(10, nrow(hub_mets))) {
    cat(sprintf("    %s: degree=%d, betweenness=%.4f\n",
                hub_mets$node[i], hub_mets$degree[i], hub_mets$betweenness[i]))
  }

  # ---- Community detection ----
  comms <- cluster_louvain(g)
  V(g)$community <- membership(comms)
  cat(sprintf("\n  Louvain communities: %d (modularity=%.3f)\n",
              length(comms), modularity(comms)))

  comm_df <- data.frame(
    node = V(g)$name,
    type = V(g)$type,
    community = V(g)$community,
    degree = deg,  # Use pre-computed degree vector
    stringsAsFactors = FALSE
  )
  write.csv(comm_df, file.path(OUT_DIR, "network_communities.csv"), row.names = FALSE)

  # Community composition
  for (c in sort(unique(comm_df$community))) {
    sub <- comm_df[comm_df$community == c, ]
    n_prot <- sum(sub$type == "protein")
    n_met <- sum(sub$type == "metabolite")
    cat(sprintf("    Community %d: %d proteins, %d metabolites\n", c, n_prot, n_met))
  }

  # ============================================================================
  # 5. VISUALIZATIONS
  # ============================================================================
  cat("\n>>> 5. Generating visualizations\n")

  # ---- 5a. Network plot ----
  pdf(file.path(FIG_DIR, "prot_metab_network.pdf"), width = 14, height = 12)
  set.seed(42)
  layout <- layout_with_fr(g)

  # Size by degree
  V(g)$size <- pmin(3 + deg * 1.5, 20)
  V(g)$label.cex <- ifelse(deg >= quantile(deg, 0.8), 0.6, 0.4)
  V(g)$label.color <- "black"

  plot(g, layout = layout,
       vertex.label = ifelse(deg >= quantile(deg, 0.75), V(g)$name, NA),
       edge.curved = 0.1,
       main = sprintf("Protein-Metabolite Correlation Network\n(%d nodes, %d edges, |rho|>%.1f, P<%.2f)",
                       vcount(g), ecount(g), RHO_CUT, P_CUT))
  legend("bottomleft",
         legend = c("Protein", "Metabolite", "Positive corr.", "Negative corr."),
         col = c("#E64B35", "#4DBBD5", "#E64B35", "#4DBBD5"),
         pch = c(16, 15, NA, NA), lty = c(NA, NA, 1, 1),
         pt.cex = 1.5, cex = 0.8, bty = "n")
  dev.off()
  cat("  Network plot saved.\n")

  # ---- 5b. Correlation heatmap (top features) ----
  # Select top proteins and metabolites by degree for heatmap
  top_n <- 30
  top_p <- head(hub_prots$node, min(top_n, nrow(hub_prots)))
  top_m <- head(hub_mets$node, min(top_n, nrow(hub_mets)))

  # Map back to IDs (names are now gene symbols, need to map to ENSP IDs for matrix lookup)
  name2prot <- setNames(edges$protein_id, edges$protein)
  name2met  <- setNames(edges$metabolite_id, edges$metabolite)

  hm_prots <- unique(name2prot[top_p])
  hm_mets  <- unique(name2met[top_m])
  hm_prots <- intersect(hm_prots, rownames(cor_mat))
  hm_mets  <- intersect(hm_mets, colnames(cor_mat))

  if (length(hm_prots) >= 3 && length(hm_mets) >= 3) {
    hm_data <- cor_mat[hm_prots, hm_mets, drop = FALSE]
    # Use gene symbols for rows (improved mapping)
    rn <- sapply(rownames(hm_data), ensp_to_gene)
    rownames(hm_data) <- make.unique(rn)  # Ensure unique names
    # Use metabolite names for columns (improved formatting)
    cn <- sapply(colnames(hm_data), function(x) format_met_name(x, 35))
    colnames(hm_data) <- make.unique(cn)  # Ensure unique names

    pdf(file.path(FIG_DIR, "prot_metab_correlation_heatmap.pdf"), width = 14, height = 10)
    pheatmap(hm_data,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             breaks = seq(-1, 1, length.out = 101),
             cluster_rows = TRUE, cluster_cols = TRUE,
             fontsize_row = 7, fontsize_col = 7,
             main = "Protein-Metabolite Spearman Correlation (Hub Features)")
    dev.off()
    cat("  Correlation heatmap saved.\n")
  }

  # ---- 5c. Degree distribution ----
  deg_df <- data.frame(node = names(deg), degree = deg,
                       type = V(g)$type, stringsAsFactors = FALSE)
  p_deg <- ggplot(deg_df, aes(x = degree, fill = type)) +
    geom_histogram(binwidth = 1, alpha = 0.7, position = "dodge") +
    scale_fill_manual(values = c("protein" = "#E64B35", "metabolite" = "#4DBBD5")) +
    labs(title = "Node Degree Distribution", x = "Degree", y = "Count", fill = "Type") +
    theme_bw(base_size = 12)
  ggsave(file.path(FIG_DIR, "degree_distribution.pdf"), p_deg, width = 8, height = 5)
  cat("  Degree distribution saved.\n")

} else {
  cat("  WARNING: No edges found. Skipping network analysis.\n")
}

# ============================================================================
# 6. METABOLITE CLASS ENRICHMENT IN NETWORK
# ============================================================================
cat("\n>>> 6. Metabolite class enrichment\n")

if (nrow(edges) > 0) {
  net_met_ids <- unique(edges$metabolite)
  all_met_ids <- rownames(mb_mat)

  # Get classes for network metabolites
  net_classes <- met2class[net_met_ids]
  net_classes <- net_classes[!is.na(net_classes) & net_classes != "" & net_classes != "-"]
  all_classes <- met2class[all_met_ids]
  all_classes <- all_classes[!is.na(all_classes) & all_classes != "" & all_classes != "-"]

  if (length(net_classes) > 0 && length(all_classes) > 0) {
    class_enrich <- data.frame()
    for (cls in unique(net_classes)) {
      n_in_net <- sum(net_classes == cls)
      n_in_bg <- sum(all_classes == cls)
      n_net <- length(net_classes)
      n_bg <- length(all_classes)

      ft <- fisher.test(matrix(c(n_in_net, n_net - n_in_net,
                                  n_in_bg - n_in_net, n_bg - n_net - n_in_bg + n_in_net),
                                nrow = 2), alternative = "greater")
      class_enrich <- rbind(class_enrich, data.frame(
        class = cls, in_network = n_in_net, in_background = n_in_bg,
        network_pct = n_in_net / n_net * 100,
        background_pct = n_in_bg / n_bg * 100,
        OR = ft$estimate, pvalue = ft$p.value, stringsAsFactors = FALSE
      ))
    }
    class_enrich$padj <- p.adjust(class_enrich$pvalue, method = "BH")
    class_enrich <- class_enrich[order(class_enrich$pvalue), ]
    write.csv(class_enrich, file.path(OUT_DIR, "metabolite_class_enrichment.csv"), row.names = FALSE)

    cat("  Metabolite class enrichment:\n")
    for (i in 1:nrow(class_enrich)) {
      sig <- ifelse(class_enrich$pvalue[i] < 0.05, "*", "")
      cat(sprintf("    %s: %d/%d (%.1f%%), OR=%.2f, P=%.4f %s\n",
                  class_enrich$class[i], class_enrich$in_network[i],
                  class_enrich$in_background[i], class_enrich$network_pct[i],
                  class_enrich$OR[i], class_enrich$pvalue[i], sig))
    }
  }
}

# ============================================================================
# 7. SUBTYPE-SPECIFIC CORRELATIONS
# ============================================================================
cat("\n>>> 7. Subtype-specific correlation analysis\n")

adj_samps <- grep("^Adjacent", common_samps, value = TRUE)
samp_sub <- sapply(adj_samps, function(s) {
  pid <- gsub("^Adjacent", "", s)
  if (pid %in% names(subtype_map)) subtype_map[pid] else NA
})
ds1_samps <- adj_samps[samp_sub == "DS1" & !is.na(samp_sub)]
ds2_samps <- adj_samps[samp_sub == "DS2" & !is.na(samp_sub)]

cat(sprintf("  DS1 Adjacent samples: %d, DS2 Adjacent samples: %d\n",
            length(ds1_samps), length(ds2_samps)))

if (length(ds1_samps) >= 4 && length(ds2_samps) >= 4 && nrow(edges) > 0) {
  # For top edges, compare correlation in DS1 vs DS2
  top_edges <- head(edges[order(edges$pvalue), ], min(100, nrow(edges)))

  subtype_cor <- data.frame()
  for (k in 1:nrow(top_edges)) {
    p_id <- top_edges$protein_id[k]   # ENSP ID for matrix lookup
    m_id <- top_edges$metabolite_id[k]  # Compound ID for matrix lookup
    p_name <- top_edges$protein[k]     # Gene symbol for display
    m_name <- top_edges$metabolite[k]  # Metabolite name for display
    m_name_full <- top_edges$metabolite_name_full[k]  # Full metabolite name
    
    if (p_id %in% rownames(pr_mat) && m_id %in% rownames(mb_mat)) {
      rho_ds1 <- cor(pr_mat[p_id, ds1_samps], mb_mat[m_id, ds1_samps], method = "spearman")
      rho_ds2 <- cor(pr_mat[p_id, ds2_samps], mb_mat[m_id, ds2_samps], method = "spearman")
      subtype_cor <- rbind(subtype_cor, data.frame(
        protein_id = p_id,
        protein = p_name,             # Gene symbol (for display)
        metabolite_id = m_id,
        metabolite = m_name,           # Short name (for display)
        metabolite_name_full = m_name_full,
        rho_all = top_edges$rho[k],
        rho_DS1 = rho_ds1, rho_DS2 = rho_ds2,
        delta_rho = abs(rho_ds1 - rho_ds2),
        stringsAsFactors = FALSE
      ))
    }
  }
  subtype_cor <- subtype_cor[order(-subtype_cor$delta_rho), ]
  write.csv(subtype_cor, file.path(OUT_DIR, "subtype_specific_correlations.csv"), row.names = FALSE)

  # Also save legacy format with protein_name / metabolite_name columns
  subtype_cor_legacy <- subtype_cor
  colnames(subtype_cor_legacy)[colnames(subtype_cor_legacy) == "protein"] <- "protein_name"
  colnames(subtype_cor_legacy)[colnames(subtype_cor_legacy) == "metabolite"] <- "metabolite_name"
  colnames(subtype_cor_legacy)[colnames(subtype_cor_legacy) == "protein_id"] <- "protein"
  colnames(subtype_cor_legacy)[colnames(subtype_cor_legacy) == "metabolite_id"] <- "metabolite"
  write.csv(subtype_cor_legacy, file.path(OUT_DIR, "subtype_specific_correlations_legacy.csv"), row.names = FALSE)

  cat("  Subtype-differential correlations (top 10 by delta_rho):\n")
  for (i in 1:min(10, nrow(subtype_cor))) {
    cat(sprintf("    %s ~ %s: DS1=%.2f, DS2=%.2f, delta=%.2f\n",
                subtype_cor$protein[i], subtype_cor$metabolite[i],
                subtype_cor$rho_DS1[i], subtype_cor$rho_DS2[i],
                subtype_cor$delta_rho[i]))
  }
} else {
  cat("  Skipping subtype-specific analysis (insufficient samples per subtype)\n")
}

# ============================================================================
# 8. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Enhancement 4 SUMMARY\n")
cat("========================================\n")
cat(sprintf("DE Proteins: %d, DE Metabolites: %d\n", length(sig_prots), length(sig_mets)))
cat(sprintf("Network edges: %d (|rho|>%.1f, P<%.2f)\n", nrow(edges), RHO_CUT, P_CUT))
if (nrow(edges) > 0) {
  cat(sprintf("Network nodes: %d proteins + %d metabolites\n",
              length(prot_nodes), length(met_nodes)))
  cat(sprintf("Louvain communities: %d\n", length(unique(comm_df$community))))
}
cat(sprintf("Output: %s\n", OUT_DIR))
cat("Enhancement 4 COMPLETE.\n")
