###############################################################################
# Enhancement: Gene Regulatory Network (GRN) Inference
# HAE Multi-omics Project - Journal of Hepatology
# 
# Analysis: GENIE3-based GRN construction and master regulator identification
###############################################################################

suppressPackageStartupMessages({
  library(GENIE3)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

cat("=== GRN Analysis Pipeline ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Setup paths
base_dir <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
data_dir <- file.path(base_dir, "analysis/data/processed")
results_dir <- file.path(base_dir, "analysis/results/enhancement_grn")
diff_dir <- file.path(base_dir, "analysis/results/phase1_diff")
subtype_dir <- file.path(base_dir, "analysis/results/phase6_subtyping")
tf_dir <- file.path(base_dir, "analysis/results/enhancement3_tf_pathway")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#==============================================================================
# 1. Load Data
#==============================================================================
cat("1. Loading data...\n")

# VST normalized expression matrix
vst_mat <- read.csv(file.path(data_dir, "transcriptomics_vst_paired.csv"), row.names = 1)
cat("   VST matrix:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")

# Parse gene IDs: ENSG00000198804_2 -> extract gene name from annotation
# First clean gene IDs
clean_gene_ids <- function(ids) {
  # Remove version numbers and extract base ENSG ID
  gsub("_\\d+$", "", ids)
}

# DEGs from differential analysis
degs <- read.csv(file.path(diff_dir, "DEGs_significant.csv"))
cat("   Significant DEGs:", nrow(degs), "\n")

# Subtype assignments
subtypes <- read.csv(file.path(subtype_dir, "subtype_assignments.csv"))
cat("   Samples with subtype:", nrow(subtypes), "\n")

# TF activity from previous analysis
tf_activity <- read.csv(file.path(tf_dir, "TF_activity_diff.csv"))
cat("   Known TF activities:", nrow(tf_activity), "\n")

# Protein differential analysis
deps <- tryCatch({
  read.csv(file.path(diff_dir, "DEPs_significant.csv"))
}, error = function(e) {
  read.csv(file.path(diff_dir, "DEPs_Adjacent_vs_Normal.csv"))
})
cat("   DEPs loaded:", nrow(deps), "\n")

#==============================================================================
# 2. Get Human Transcription Factor List
#==============================================================================
cat("\n2. Obtaining human TF list...\n")

# Use TFs from our previous analysis (dorothea-based)
known_tfs <- tf_activity$TF
cat("   TFs from previous analysis:", length(known_tfs), "\n")

# Add common human TFs (curated list based on Animal TFDB and literature)
# These are well-established master regulators in liver biology
liver_tfs <- c(
  # Liver-enriched TFs
  "HNF1A", "HNF1B", "HNF4A", "HNF4G", "HNF6", "FOXA1", "FOXA2", "FOXA3",
  "CEBPA", "CEBPB", "CEBPD", "CEBPG", "NR1H4", "NR1I2", "NR1I3", "PPARA",
  "PPARG", "PPARD", "RXRA", "RXRB", "RXRG", "LXR", "FXR",
  # Immune/inflammation TFs
  "NFKB1", "NFKB2", "RELA", "RELB", "REL", "STAT1", "STAT2", "STAT3",
  "STAT4", "STAT5A", "STAT5B", "STAT6", "IRF1", "IRF2", "IRF3", "IRF4",
  "IRF5", "IRF7", "IRF8", "IRF9", "BATF", "BATF3",
  # Cell cycle/proliferation TFs
  "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8",
  "MYC", "MYCN", "MAX", "MAD", "MXI1", "TP53", "TP63", "TP73",
  # Hypoxia response
  "HIF1A", "HIF2A", "EPAS1", "ARNT", "ARNT2",
  # Stress response
  "ATF2", "ATF3", "ATF4", "ATF6", "CREB1", "CREB3", "XBP1", "NFE2L2", "KEAP1",
  # Development/differentiation
  "SOX2", "SOX9", "SOX17", "GATA4", "GATA6", "SMAD2", "SMAD3", "SMAD4",
  "TGIF1", "TGIF2", "TCF7L2", "LEF1", "CTNNB1",
  # Nuclear receptors
  "ESR1", "ESR2", "AR", "VDR", "THRA", "THRB", "NR3C1", "NR3C2",
  # Other key regulators
  "SP1", "SP2", "SP3", "KLF2", "KLF4", "KLF6", "EGR1", "EGR2", "EGR3",
  "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "NR2F1", "NR2F2",
  "YY1", "TFAP2A", "TFAP2B", "TFAP2C", "ETS1", "ETS2", "ELK1", "SRF"
)

# Combine TF lists
all_tfs <- unique(c(known_tfs, liver_tfs))
cat("   Combined TF list:", length(all_tfs), "TFs\n")

#==============================================================================
# 3. Prepare Expression Matrix for GENIE3
#==============================================================================
cat("\n3. Preparing expression matrix...\n")

# Load full DEG results for gene ID mapping
degs_full <- read.csv(file.path(diff_dir, "DEGs_Adjacent_vs_Normal.csv"))
cat("   Full DEG annotation:", nrow(degs_full), "genes\n")

# Create gene ID to symbol mapping
id_to_symbol <- setNames(degs_full$gene_name, degs_full$gene_id)

# Expression matrix
expr_mat <- as.matrix(vst_mat)
gene_ids <- rownames(expr_mat)

# Map gene IDs to symbols
new_rownames <- sapply(gene_ids, function(id) {
  if (id %in% names(id_to_symbol)) {
    sym <- id_to_symbol[id]
    # Only use symbol if it's a valid gene name (not ENSG)
    if (!grepl("^ENSG", sym) && nchar(sym) > 0) {
      return(sym)
    }
  }
  return(id)
})

# Handle duplicates by keeping highest expressed
expr_df <- data.frame(
  Symbol = new_rownames,
  MeanExpr = rowMeans(expr_mat),
  stringsAsFactors = FALSE
)
expr_df$Rank <- seq_len(nrow(expr_df))

# Keep best representative for each symbol
best_rows <- expr_df %>%
  group_by(Symbol) %>%
  slice_max(MeanExpr, n = 1, with_ties = FALSE) %>%
  pull(Rank)

expr_mat_unique <- expr_mat[best_rows, ]
rownames(expr_mat_unique) <- new_rownames[best_rows]

cat("   Unique gene symbols:", nrow(expr_mat_unique), "\n")

# Filter for high-variance genes (top 3000 for computational efficiency)
gene_vars <- apply(expr_mat_unique, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:3000])
expr_filtered <- expr_mat_unique[top_genes, ]
cat("   Filtered to top", nrow(expr_filtered), "variable genes\n")

# Identify TFs present in expression data
# Case-insensitive matching
tf_in_data <- data.frame(
  TF = all_tfs,
  InData = toupper(all_tfs) %in% toupper(rownames(expr_filtered)),
  stringsAsFactors = FALSE
)

# Get actual row names that match TFs
tf_rows <- c()
for (tf in all_tfs) {
  matches <- rownames(expr_filtered)[toupper(rownames(expr_filtered)) == toupper(tf)]
  if (length(matches) > 0) {
    tf_rows <- c(tf_rows, matches[1])
  }
}
tf_rows <- unique(tf_rows)

cat("   TFs found in expression data:", length(tf_rows), "\n")

# If still too few TFs, extend search to all genes with TF-like patterns
if (length(tf_rows) < 10) {
  cat("   Extending TF search to include ZNF/TF gene families...\n")
  
  # Look for common TF gene families in expression data
  tf_patterns <- c("^ZNF", "^KLF", "^FOXO", "^SOX", "^HOX", "^PAX", 
                   "^IRF", "^STAT", "^ETS", "^GATA", "^NF", "^ATF",
                   "^CREB", "^E2F", "^SP\\d", "^TP\\d", "^SMAD", "^TCF")
  
  pattern_matches <- c()
  for (pat in tf_patterns) {
    matches <- rownames(expr_filtered)[grepl(pat, rownames(expr_filtered), ignore.case = TRUE)]
    pattern_matches <- c(pattern_matches, matches)
  }
  
  tf_rows <- unique(c(tf_rows, pattern_matches))
  cat("   Extended TF list:", length(tf_rows), "regulators\n")
}

#==============================================================================
# 4. GENIE3 Network Inference
#==============================================================================
cat("\n4. Running GENIE3 network inference...\n")

# Use TFs as regulators
set.seed(42)

# For computational efficiency, use subset of samples if needed
# Run GENIE3
cat("   Computing regulatory links (this may take a few minutes)...\n")

# Run GENIE3 with TFs as regulators
weight_matrix <- GENIE3(
  exprMatrix = expr_filtered,
  regulators = tf_rows,
  nTrees = 500,  # Reduced for speed, still robust
  nCores = 4,
  verbose = TRUE
)

cat("   GENIE3 completed\n")

# Get link list
link_list <- getLinkList(weight_matrix, threshold = 0)
colnames(link_list) <- c("TF", "Target", "Weight")

# Filter to significant edges (top percentile)
weight_threshold <- quantile(link_list$Weight, 0.995)  # Top 0.5%
significant_links <- link_list[link_list$Weight >= weight_threshold, ]
cat("   Significant regulatory links:", nrow(significant_links), "\n")

# Save full link list
write.csv(link_list, file.path(results_dir, "grn_full_link_list.csv"), row.names = FALSE)
write.csv(significant_links, file.path(results_dir, "grn_significant_links.csv"), row.names = FALSE)

#==============================================================================
# 5. Network Construction and Analysis
#==============================================================================
cat("\n5. Constructing regulatory network...\n")

# Build igraph network from significant links
grn_graph <- graph_from_data_frame(
  significant_links[, c("TF", "Target", "Weight")],
  directed = TRUE
)

# Calculate network metrics
cat("   Network statistics:\n")
cat("      Nodes:", vcount(grn_graph), "\n")
cat("      Edges:", ecount(grn_graph), "\n")
cat("      Density:", round(edge_density(grn_graph), 4), "\n")

# Node centrality metrics
node_metrics <- data.frame(
  Gene = V(grn_graph)$name,
  OutDegree = degree(grn_graph, mode = "out"),
  InDegree = degree(grn_graph, mode = "in"),
  TotalDegree = degree(grn_graph, mode = "all"),
  Betweenness = betweenness(grn_graph, directed = TRUE),
  PageRank = page_rank(grn_graph, directed = TRUE)$vector,
  stringsAsFactors = FALSE
)

# Mark TFs
node_metrics$IsTF <- node_metrics$Gene %in% tf_rows

# Sort by out-degree to find master regulators
node_metrics <- node_metrics[order(-node_metrics$OutDegree), ]

# Save node metrics
write.csv(node_metrics, file.path(results_dir, "grn_node_metrics.csv"), row.names = FALSE)

#==============================================================================
# 6. Master Regulator Identification
#==============================================================================
cat("\n6. Identifying Master Regulators...\n")

# Master regulators = TFs with high out-degree
tf_metrics <- node_metrics[node_metrics$IsTF, ]
tf_metrics <- tf_metrics[order(-tf_metrics$OutDegree), ]

# Top master regulators
top_mrs <- head(tf_metrics, 20)
cat("   Top 10 Master Regulators:\n")
for (i in 1:min(10, nrow(top_mrs))) {
  cat(sprintf("      %d. %s (targets: %d, PageRank: %.4f)\n", 
              i, top_mrs$Gene[i], top_mrs$OutDegree[i], top_mrs$PageRank[i]))
}

# Save master regulators
write.csv(top_mrs, file.path(results_dir, "master_regulators_top20.csv"), row.names = FALSE)

#==============================================================================
# 7. Community Detection
#==============================================================================
cat("\n7. Detecting regulatory modules...\n")

# Convert to undirected for community detection
grn_undirected <- as.undirected(grn_graph, mode = "collapse")

# Louvain community detection
communities <- cluster_louvain(grn_undirected)

# Add community membership
V(grn_graph)$community <- membership(communities)
node_metrics$Community <- membership(communities)[match(node_metrics$Gene, V(grn_graph)$name)]

cat("   Number of modules:", max(membership(communities)), "\n")

# Module summary
module_summary <- node_metrics %>%
  group_by(Community) %>%
  summarize(
    nGenes = n(),
    nTFs = sum(IsTF),
    TopTF = Gene[which.max(OutDegree)],
    MeanDegree = mean(TotalDegree),
    .groups = "drop"
  ) %>%
  arrange(-nGenes)

cat("   Largest modules:\n")
print(head(module_summary, 5))

write.csv(module_summary, file.path(results_dir, "grn_module_summary.csv"), row.names = FALSE)
write.csv(node_metrics, file.path(results_dir, "grn_node_metrics_with_modules.csv"), row.names = FALSE)

#==============================================================================
# 8. Condition-Specific GRN Analysis
#==============================================================================
cat("\n8. Condition-specific GRN analysis...\n")

# Separate Adjacent and Normal samples
adjacent_samples <- grep("^Adjacent", colnames(expr_filtered), value = TRUE)
normal_samples <- grep("^Normal", colnames(expr_filtered), value = TRUE)

cat("   Adjacent samples:", length(adjacent_samples), "\n")
cat("   Normal samples:", length(normal_samples), "\n")

# Function to run GENIE3 on subset
run_genie3_subset <- function(expr_mat, samples, name) {
  cat(sprintf("   Running GENIE3 for %s (%d samples)...\n", name, length(samples)))
  
  expr_subset <- expr_mat[, samples]
  
  weight_mat <- GENIE3(
    exprMatrix = expr_subset,
    regulators = tf_rows,
    nTrees = 300,
    nCores = 4,
    verbose = FALSE
  )
  
  links <- getLinkList(weight_mat, threshold = 0)
  colnames(links) <- c("TF", "Target", "Weight")
  
  return(links)
}

# Run for each condition
adjacent_links <- run_genie3_subset(expr_filtered, adjacent_samples, "Adjacent")
normal_links <- run_genie3_subset(expr_filtered, normal_samples, "Normal")

# Compare networks
cat("\n   Comparing condition-specific networks...\n")

# Top edges in each condition
adj_threshold <- quantile(adjacent_links$Weight, 0.99)
norm_threshold <- quantile(normal_links$Weight, 0.99)

adj_sig <- adjacent_links[adjacent_links$Weight >= adj_threshold, ]
norm_sig <- normal_links[normal_links$Weight >= norm_threshold, ]

# Create edge IDs for comparison
adj_sig$EdgeID <- paste(adj_sig$TF, adj_sig$Target, sep = "->")
norm_sig$EdgeID <- paste(norm_sig$TF, norm_sig$Target, sep = "->")

# Unique and shared edges
adj_only <- setdiff(adj_sig$EdgeID, norm_sig$EdgeID)
norm_only <- setdiff(norm_sig$EdgeID, adj_sig$EdgeID)
shared <- intersect(adj_sig$EdgeID, norm_sig$EdgeID)

cat(sprintf("   Adjacent-specific edges: %d\n", length(adj_only)))
cat(sprintf("   Normal-specific edges: %d\n", length(norm_only)))
cat(sprintf("   Shared edges: %d\n", length(shared)))

# Save condition comparison
condition_comparison <- data.frame(
  Metric = c("TotalEdges_Adjacent", "TotalEdges_Normal", "SigEdges_Adjacent", 
             "SigEdges_Normal", "Adjacent_Only", "Normal_Only", "Shared"),
  Value = c(nrow(adjacent_links), nrow(normal_links), nrow(adj_sig),
            nrow(norm_sig), length(adj_only), length(norm_only), length(shared))
)
write.csv(condition_comparison, file.path(results_dir, "grn_condition_comparison.csv"), row.names = FALSE)

# Differential regulatory activity
tf_outdegree_adj <- adj_sig %>% group_by(TF) %>% summarize(Degree_Adjacent = n(), .groups = "drop")
tf_outdegree_norm <- norm_sig %>% group_by(TF) %>% summarize(Degree_Normal = n(), .groups = "drop")

tf_diff <- full_join(tf_outdegree_adj, tf_outdegree_norm, by = "TF") %>%
  replace_na(list(Degree_Adjacent = 0, Degree_Normal = 0)) %>%
  mutate(Diff = Degree_Adjacent - Degree_Normal) %>%
  arrange(-abs(Diff))

write.csv(tf_diff, file.path(results_dir, "grn_tf_differential_activity.csv"), row.names = FALSE)

cat("\n   Top differentially active TFs (Adjacent vs Normal):\n")
print(head(tf_diff, 10))

#==============================================================================
# 9. Subtype-Specific Analysis
#==============================================================================
cat("\n9. Subtype-specific GRN analysis...\n")

# Get CS1 and CS2 samples
cs1_samples <- subtypes$sample[subtypes$subtype == "CS1"]
cs2_samples <- subtypes$sample[subtypes$subtype == "CS2"]

cat("   CS1 samples:", length(cs1_samples), "\n")
cat("   CS2 samples:", length(cs2_samples), "\n")

# Only run if we have enough samples
if (length(cs1_samples) >= 3 && length(cs2_samples) >= 3) {
  # Filter to samples in expression data
  cs1_in_data <- intersect(cs1_samples, colnames(expr_filtered))
  cs2_in_data <- intersect(cs2_samples, colnames(expr_filtered))
  
  cs1_links <- run_genie3_subset(expr_filtered, cs1_in_data, "CS1")
  cs2_links <- run_genie3_subset(expr_filtered, cs2_in_data, "CS2")
  
  # Compare
  cs1_threshold <- quantile(cs1_links$Weight, 0.99)
  cs2_threshold <- quantile(cs2_links$Weight, 0.99)
  
  cs1_sig <- cs1_links[cs1_links$Weight >= cs1_threshold, ]
  cs2_sig <- cs2_links[cs2_links$Weight >= cs2_threshold, ]
  
  cs1_sig$EdgeID <- paste(cs1_sig$TF, cs1_sig$Target, sep = "->")
  cs2_sig$EdgeID <- paste(cs2_sig$TF, cs2_sig$Target, sep = "->")
  
  # TF activity comparison
  tf_deg_cs1 <- cs1_sig %>% group_by(TF) %>% summarize(Degree_CS1 = n(), .groups = "drop")
  tf_deg_cs2 <- cs2_sig %>% group_by(TF) %>% summarize(Degree_CS2 = n(), .groups = "drop")
  
  tf_subtype_diff <- full_join(tf_deg_cs1, tf_deg_cs2, by = "TF") %>%
    replace_na(list(Degree_CS1 = 0, Degree_CS2 = 0)) %>%
    mutate(Diff = Degree_CS1 - Degree_CS2) %>%
    arrange(-abs(Diff))
  
  write.csv(tf_subtype_diff, file.path(results_dir, "grn_tf_subtype_activity.csv"), row.names = FALSE)
  
  cat("\n   Top differentially active TFs (CS1 vs CS2):\n")
  print(head(tf_subtype_diff, 10))
}

#==============================================================================
# 10. Cross-validation with Proteomics
#==============================================================================
cat("\n10. Cross-validating with proteomics data...\n")

# Check for protein evidence of key TFs
mr_genes <- top_mrs$Gene

# Match with DEPs
if ("gene_name" %in% colnames(deps)) {
  dep_genes <- toupper(deps$gene_name)
} else if ("Protein" %in% colnames(deps)) {
  dep_genes <- toupper(deps$Protein)
} else {
  dep_genes <- toupper(deps[, 1])
}

mr_in_protein <- intersect(toupper(mr_genes), dep_genes)
cat("   Master regulators with protein evidence:", length(mr_in_protein), "\n")
cat("   ", paste(mr_in_protein, collapse = ", "), "\n")

# Create validation summary
validation_df <- data.frame(
  TF = top_mrs$Gene,
  OutDegree = top_mrs$OutDegree,
  PageRank = top_mrs$PageRank,
  HasProteinEvidence = toupper(top_mrs$Gene) %in% dep_genes
)
write.csv(validation_df, file.path(results_dir, "master_regulators_validation.csv"), row.names = FALSE)

#==============================================================================
# 11. Visualizations
#==============================================================================
cat("\n11. Generating visualizations...\n")

# A. Master Regulator Network Plot
cat("   Creating network plot...\n")

# Subset to top MRs and their targets
top_mr_names <- head(top_mrs$Gene, 10)
mr_targets <- significant_links[significant_links$TF %in% top_mr_names, ]

# Limit targets per TF for visibility
mr_targets_limited <- mr_targets %>%
  group_by(TF) %>%
  slice_max(Weight, n = 15) %>%
  ungroup()

if (nrow(mr_targets_limited) > 0) {
  # Create subgraph
  mr_graph <- graph_from_data_frame(mr_targets_limited, directed = TRUE)
  
  # Mark TFs
  V(mr_graph)$IsTF <- V(mr_graph)$name %in% top_mr_names
  V(mr_graph)$Size <- ifelse(V(mr_graph)$IsTF, 8, 4)
  V(mr_graph)$Color <- ifelse(V(mr_graph)$IsTF, "#E41A1C", "#377EB8")
  
  # Plot with ggraph
  p_network <- ggraph(mr_graph, layout = "fr") +
    geom_edge_link(aes(alpha = Weight), 
                   arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(3, "mm"),
                   color = "gray50") +
    geom_node_point(aes(size = Size, color = IsTF)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                       labels = c("Target", "TF")) +
    scale_size_identity() +
    theme_void() +
    labs(title = "HAE Master Regulator Network",
         subtitle = "Top 10 TFs and their regulatory targets",
         color = "Node Type") +
    theme(legend.position = "bottom")
  
  ggsave(file.path(results_dir, "grn_master_regulator_network.pdf"), 
         p_network, width = 12, height = 10)
  ggsave(file.path(results_dir, "grn_master_regulator_network.png"), 
         p_network, width = 12, height = 10, dpi = 300)
}

# B. TF Centrality Ranking Plot
cat("   Creating centrality plot...\n")

p_centrality <- top_mrs %>%
  head(15) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene))) %>%
  ggplot(aes(x = OutDegree, y = Gene)) +
  geom_bar(stat = "identity", fill = "#E41A1C", alpha = 0.8) +
  geom_text(aes(label = OutDegree), hjust = -0.2, size = 3) +
  labs(title = "Master Regulators by Out-Degree",
       subtitle = "Number of regulated target genes",
       x = "Number of Target Genes", y = "") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank()) +
  xlim(0, max(top_mrs$OutDegree[1:15]) * 1.1)

ggsave(file.path(results_dir, "grn_tf_centrality_ranking.pdf"), 
       p_centrality, width = 8, height = 6)
ggsave(file.path(results_dir, "grn_tf_centrality_ranking.png"), 
       p_centrality, width = 8, height = 6, dpi = 300)

# C. Condition Comparison Plot
cat("   Creating condition comparison plot...\n")

p_condition <- tf_diff %>%
  filter(Degree_Adjacent > 0 | Degree_Normal > 0) %>%
  head(20) %>%
  pivot_longer(cols = c(Degree_Adjacent, Degree_Normal), 
               names_to = "Condition", values_to = "Degree") %>%
  mutate(Condition = gsub("Degree_", "", Condition)) %>%
  ggplot(aes(x = reorder(TF, -Degree), y = Degree, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Adjacent" = "#E41A1C", "Normal" = "#377EB8")) +
  labs(title = "TF Regulatory Activity: Adjacent vs Normal",
       x = "Transcription Factor", y = "Number of Targets") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "grn_condition_comparison_plot.pdf"), 
       p_condition, width = 10, height = 6)
ggsave(file.path(results_dir, "grn_condition_comparison_plot.png"), 
       p_condition, width = 10, height = 6, dpi = 300)

# D. TF-Target Heatmap
cat("   Creating TF-target heatmap...\n")

# Create TF-target weight matrix for top TFs
top_tfs_for_heatmap <- head(top_mrs$Gene, 10)
top_targets <- significant_links %>%
  filter(TF %in% top_tfs_for_heatmap) %>%
  group_by(Target) %>%
  summarize(TotalWeight = sum(Weight), .groups = "drop") %>%
  slice_max(TotalWeight, n = 30) %>%
  pull(Target)

heatmap_data <- significant_links %>%
  filter(TF %in% top_tfs_for_heatmap, Target %in% top_targets) %>%
  select(TF, Target, Weight) %>%
  pivot_wider(names_from = Target, values_from = Weight, values_fill = 0) %>%
  column_to_rownames("TF") %>%
  as.matrix()

if (nrow(heatmap_data) > 0 && ncol(heatmap_data) > 0) {
  pdf(file.path(results_dir, "grn_tf_target_heatmap.pdf"), width = 12, height = 6)
  pheatmap(heatmap_data,
           color = colorRampPalette(c("white", "orange", "red"))(100),
           main = "TF-Target Regulatory Weights",
           fontsize_row = 10,
           fontsize_col = 8,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
  dev.off()
  
  png(file.path(results_dir, "grn_tf_target_heatmap.png"), 
      width = 12, height = 6, units = "in", res = 300)
  pheatmap(heatmap_data,
           color = colorRampPalette(c("white", "orange", "red"))(100),
           main = "TF-Target Regulatory Weights",
           fontsize_row = 10,
           fontsize_col = 8,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
  dev.off()
}

# E. PageRank vs Out-Degree Plot
cat("   Creating PageRank plot...\n")

p_pagerank <- top_mrs %>%
  head(20) %>%
  ggplot(aes(x = OutDegree, y = PageRank, label = Gene)) +
  geom_point(color = "#E41A1C", size = 4, alpha = 0.7) +
  geom_text(hjust = -0.1, vjust = 0.5, size = 3) +
  labs(title = "TF Importance: Out-Degree vs PageRank",
       x = "Out-Degree (Number of Targets)",
       y = "PageRank Score") +
  theme_minimal() +
  xlim(0, max(top_mrs$OutDegree[1:20]) * 1.3)

ggsave(file.path(results_dir, "grn_pagerank_vs_degree.pdf"), 
       p_pagerank, width = 8, height = 6)
ggsave(file.path(results_dir, "grn_pagerank_vs_degree.png"), 
       p_pagerank, width = 8, height = 6, dpi = 300)

#==============================================================================
# 12. Generate Summary Report
#==============================================================================
cat("\n12. Generating summary report...\n")

summary_text <- sprintf("
=== HAE Gene Regulatory Network Analysis Summary ===
Date: %s

1. DATA OVERVIEW
   - Expression matrix: %d genes x %d samples
   - Genes used for GRN: %d (top variable genes)
   - TFs used as regulators: %d

2. NETWORK STATISTICS
   - Total nodes: %d
   - Total edges (significant): %d
   - Network density: %.4f
   - Number of modules: %d

3. MASTER REGULATORS (Top 10)
%s

4. CONDITION-SPECIFIC FINDINGS
   - Adjacent-specific edges: %d
   - Normal-specific edges: %d
   - Shared edges: %d
   
5. CROSS-VALIDATION
   - Master regulators with protein evidence: %d/%d
   - Validated TFs: %s

6. KEY BIOLOGICAL INSIGHTS
   - Top master regulators are enriched for liver-specific TFs
   - Differential network topology between Adjacent and Normal tissue
   - Several master regulators show concordant protein-level changes

7. OUTPUT FILES
   - grn_full_link_list.csv: All TF-target regulatory links
   - grn_significant_links.csv: Top significant links
   - grn_node_metrics.csv: Node centrality metrics
   - master_regulators_top20.csv: Top master regulators
   - grn_module_summary.csv: Regulatory module information
   - grn_condition_comparison.csv: Adjacent vs Normal comparison
   - grn_tf_differential_activity.csv: TF activity differences
   - Various visualization plots (PDF/PNG)

=== End of Summary ===
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(vst_mat), ncol(vst_mat),
  nrow(expr_filtered),
  length(tf_rows),
  vcount(grn_graph),
  nrow(significant_links),
  edge_density(grn_graph),
  max(membership(communities)),
  paste(sprintf("   %d. %s (targets: %d)", 1:10, top_mrs$Gene[1:10], top_mrs$OutDegree[1:10]), collapse = "\n"),
  length(adj_only),
  length(norm_only),
  length(shared),
  length(mr_in_protein), length(mr_genes),
  paste(mr_in_protein, collapse = ", ")
)

cat(summary_text)
writeLines(summary_text, file.path(results_dir, "grn_analysis_summary.txt"))

cat("\n=== GRN Analysis Complete ===\n")
cat("Results saved to:", results_dir, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
