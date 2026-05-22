#!/usr/bin/env Rscript
# ============================================================================
# Phase 2: Pathway Enrichment & Cross-omics Convergence
# HAE Multi-omics Integration Study
# ============================================================================
# Strategy:
#   - GSEA (fgsea) on full ranked gene/protein lists — PRIMARY approach
#   - GO-BP via clusterProfiler::gseGO (local org.Hs.eg.db)
#   - Cross-omics pathway convergence (Venn + hypergeometric)
#   - Gene sets: Hallmark, KEGG, Reactome (GMT), GO-BP (local db)
# ============================================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
  library(grid)
})

# ---- Paths ----
PROJECT  <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
DIFF_DIR <- file.path(PROJECT, "analysis/results/phase1_diff")
OUT_DIR  <- file.path(PROJECT, "analysis/results/phase2_enrichment")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 2: Pathway Enrichment & GSEA\n")
cat("========================================\n\n")

# ============================================================================
# 0. PREPARE GENE SETS
# ============================================================================
cat(">>> 0. Preparing gene sets\n")

# --- 0a. GMT files from Broad Institute (Hallmark, KEGG, Reactome) ---
GMT_DIR <- file.path(PROJECT, "analysis/data/gmt")
dir.create(GMT_DIR, recursive = TRUE, showWarnings = FALSE)

MSIGDB_VER <- "2024.1.Hs"
BASE_URL   <- paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/", MSIGDB_VER)

download_gmt <- function(filename) {
  local_path <- file.path(GMT_DIR, filename)
  if (!file.exists(local_path)) {
    url <- paste0(BASE_URL, "/", filename)
    cat(sprintf("  Downloading %s ...\n", filename))
    download.file(url, local_path, quiet = TRUE, method = "auto")
  } else {
    cat(sprintf("  Using cached %s\n", filename))
  }
  gmtPathways(local_path)
}

hallmark_list <- download_gmt(paste0("h.all.v", MSIGDB_VER, ".symbols.gmt"))
kegg_list     <- download_gmt(paste0("c2.cp.kegg_medicus.v", MSIGDB_VER, ".symbols.gmt"))
reactome_list <- download_gmt(paste0("c2.cp.reactome.v", MSIGDB_VER, ".symbols.gmt"))

cat(sprintf("  Hallmark: %d pathways\n", length(hallmark_list)))
cat(sprintf("  KEGG: %d pathways\n", length(kegg_list)))
cat(sprintf("  Reactome: %d pathways\n", length(reactome_list)))
cat("  GO-BP: will use clusterProfiler::gseGO (local org.Hs.eg.db)\n")


# ============================================================================
# HELPER: Symbol -> Entrez ranked list for clusterProfiler gseGO
# ============================================================================
symbol_to_entrez_ranked <- function(ranked_symbols) {
  gene_df <- bitr(names(ranked_symbols), fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_df <- gene_df[!duplicated(gene_df$SYMBOL), ]
  ranked_entrez <- ranked_symbols[gene_df$SYMBOL]
  names(ranked_entrez) <- gene_df$ENTREZID
  ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)
  return(ranked_entrez)
}


# ============================================================================
# 1. TRANSCRIPTOMICS GSEA
# ============================================================================
cat("\n>>> 1. Transcriptomics GSEA\n")

res_tc <- read.csv(file.path(DIFF_DIR, "DEGs_Adjacent_vs_Normal.csv"), check.names = FALSE)

res_tc_clean <- res_tc[!is.na(res_tc$gene_name) & res_tc$gene_name != "", ]
res_tc_clean <- res_tc_clean[!duplicated(res_tc_clean$gene_name), ]
ranked_tc <- setNames(res_tc_clean$t, res_tc_clean$gene_name)
ranked_tc <- sort(ranked_tc, decreasing = TRUE)
cat(sprintf("  Ranked gene list: %d genes\n", length(ranked_tc)))

# fgsea for Hallmark, KEGG, Reactome
run_gsea <- function(ranked, pathways, name, min_size = 15, max_size = 500) {
  set.seed(42)
  res <- fgsea(pathways = pathways, stats = ranked,
               minSize = min_size, maxSize = max_size, nPermSimple = 10000)
  res <- res[order(res$pval), ]
  cat(sprintf("  %s: %d significant (padj<0.05), %d at padj<0.25\n",
              name, sum(res$padj < 0.05), sum(res$padj < 0.25)))
  return(as.data.frame(res))
}

gsea_tc_hallmark <- run_gsea(ranked_tc, hallmark_list, "Hallmark")
gsea_tc_kegg     <- run_gsea(ranked_tc, kegg_list, "KEGG")
gsea_tc_reactome <- run_gsea(ranked_tc, reactome_list, "Reactome")

# GO-BP via clusterProfiler gseGO
cat("  Running GO-BP (gseGO with org.Hs.eg.db)...\n")
ranked_tc_entrez <- symbol_to_entrez_ranked(ranked_tc)
cat(sprintf("  Mapped to Entrez: %d genes\n", length(ranked_tc_entrez)))

gsea_tc_gobp_raw <- gseGO(geneList = ranked_tc_entrez, OrgDb = org.Hs.eg.db,
                           ont = "BP", minGSSize = 15, maxGSSize = 300,
                           pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE)
gsea_tc_gobp <- as.data.frame(gsea_tc_gobp_raw)
# Rename to match fgsea format
if (nrow(gsea_tc_gobp) > 0) {
  colnames(gsea_tc_gobp)[colnames(gsea_tc_gobp) == "ID"] <- "pathway"
  colnames(gsea_tc_gobp)[colnames(gsea_tc_gobp) == "Description"] <- "pathway_name"
  colnames(gsea_tc_gobp)[colnames(gsea_tc_gobp) == "enrichmentScore"] <- "ES"
  colnames(gsea_tc_gobp)[colnames(gsea_tc_gobp) == "NES"] <- "NES"
  colnames(gsea_tc_gobp)[colnames(gsea_tc_gobp) == "p.adjust"] <- "padj"
  colnames(gsea_tc_gobp)[colnames(gsea_tc_gobp) == "pvalue"] <- "pval"
}
cat(sprintf("  GO-BP: %d significant (padj<0.05), %d at padj<0.25\n",
            sum(gsea_tc_gobp$padj < 0.05, na.rm = TRUE),
            sum(gsea_tc_gobp$padj < 0.25, na.rm = TRUE)))

# Save results
save_gsea <- function(df, filepath) {
  df_save <- df
  if ("leadingEdge" %in% colnames(df_save)) {
    df_save$leadingEdge <- sapply(df_save$leadingEdge, function(x) paste(x, collapse = ";"))
  }
  if ("core_enrichment" %in% colnames(df_save)) {
    # already string from clusterProfiler
  }
  write.csv(df_save, filepath, row.names = FALSE)
}

save_gsea(gsea_tc_hallmark, file.path(OUT_DIR, "GSEA_transcriptomics_Hallmark.csv"))
save_gsea(gsea_tc_kegg, file.path(OUT_DIR, "GSEA_transcriptomics_KEGG.csv"))
save_gsea(gsea_tc_reactome, file.path(OUT_DIR, "GSEA_transcriptomics_Reactome.csv"))
save_gsea(gsea_tc_gobp, file.path(OUT_DIR, "GSEA_transcriptomics_GOBP.csv"))


# ============================================================================
# 2. PROTEOMICS GSEA
# ============================================================================
cat("\n>>> 2. Proteomics GSEA\n")

res_pr <- read.csv(file.path(DIFF_DIR, "DEPs_Adjacent_vs_Normal.csv"), check.names = FALSE)

res_pr_clean <- res_pr[!is.na(res_pr$gene_name) & res_pr$gene_name != "", ]
res_pr_clean <- res_pr_clean[!duplicated(res_pr_clean$gene_name), ]
ranked_pr <- setNames(res_pr_clean$t, res_pr_clean$gene_name)
ranked_pr <- sort(ranked_pr, decreasing = TRUE)
cat(sprintf("  Ranked protein list: %d proteins\n", length(ranked_pr)))

gsea_pr_hallmark <- run_gsea(ranked_pr, hallmark_list, "Hallmark")
gsea_pr_kegg     <- run_gsea(ranked_pr, kegg_list, "KEGG")
gsea_pr_reactome <- run_gsea(ranked_pr, reactome_list, "Reactome")

# GO-BP
cat("  Running GO-BP (gseGO with org.Hs.eg.db)...\n")
ranked_pr_entrez <- symbol_to_entrez_ranked(ranked_pr)
cat(sprintf("  Mapped to Entrez: %d proteins\n", length(ranked_pr_entrez)))

gsea_pr_gobp_raw <- gseGO(geneList = ranked_pr_entrez, OrgDb = org.Hs.eg.db,
                           ont = "BP", minGSSize = 15, maxGSSize = 300,
                           pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE)
gsea_pr_gobp <- as.data.frame(gsea_pr_gobp_raw)
if (nrow(gsea_pr_gobp) > 0) {
  colnames(gsea_pr_gobp)[colnames(gsea_pr_gobp) == "ID"] <- "pathway"
  colnames(gsea_pr_gobp)[colnames(gsea_pr_gobp) == "Description"] <- "pathway_name"
  colnames(gsea_pr_gobp)[colnames(gsea_pr_gobp) == "enrichmentScore"] <- "ES"
  colnames(gsea_pr_gobp)[colnames(gsea_pr_gobp) == "NES"] <- "NES"
  colnames(gsea_pr_gobp)[colnames(gsea_pr_gobp) == "p.adjust"] <- "padj"
  colnames(gsea_pr_gobp)[colnames(gsea_pr_gobp) == "pvalue"] <- "pval"
}
cat(sprintf("  GO-BP: %d significant (padj<0.05), %d at padj<0.25\n",
            sum(gsea_pr_gobp$padj < 0.05, na.rm = TRUE),
            sum(gsea_pr_gobp$padj < 0.25, na.rm = TRUE)))

save_gsea(gsea_pr_hallmark, file.path(OUT_DIR, "GSEA_proteomics_Hallmark.csv"))
save_gsea(gsea_pr_kegg, file.path(OUT_DIR, "GSEA_proteomics_KEGG.csv"))
save_gsea(gsea_pr_reactome, file.path(OUT_DIR, "GSEA_proteomics_Reactome.csv"))
save_gsea(gsea_pr_gobp, file.path(OUT_DIR, "GSEA_proteomics_GOBP.csv"))


# ============================================================================
# 3. GSEA VISUALIZATION
# ============================================================================
cat("\n>>> 3. GSEA Visualization\n")

plot_gsea_dotplot <- function(df, title, n_show = 20, prefix_strip = NULL) {
  if (nrow(df) == 0) { cat("  [skip] No results for:", title, "\n"); return(NULL) }
  df_sig <- df[df$padj < 0.25, ]
  if (nrow(df_sig) == 0) df_sig <- head(df[order(df$pval), ], n_show)
  df_plot <- head(df_sig[order(df_sig$pval), ], n_show)

  if (!is.null(prefix_strip)) {
    df_plot$pathway_short <- gsub(prefix_strip, "", df_plot$pathway)
  } else if ("pathway_name" %in% colnames(df_plot)) {
    df_plot$pathway_short <- df_plot$pathway_name
  } else {
    df_plot$pathway_short <- df_plot$pathway
  }
  df_plot$pathway_short <- gsub("_", " ", df_plot$pathway_short)
  df_plot$pathway_short <- substr(df_plot$pathway_short, 1, 55)
  df_plot$direction <- ifelse(df_plot$NES > 0, "Up in Adjacent", "Down in Adjacent")

  p <- ggplot(df_plot, aes(NES, reorder(pathway_short, NES))) +
    geom_point(aes(size = -log10(padj), color = direction), alpha = 0.8) +
    scale_color_manual(values = c("Up in Adjacent" = "#E64B35",
                                   "Down in Adjacent" = "#4DBBD5")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = title, x = "Normalized Enrichment Score (NES)",
         y = NULL, size = "-log10(padj)") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
  return(p)
}

# Hallmark dotplots
p1 <- plot_gsea_dotplot(gsea_tc_hallmark, "Transcriptomics GSEA - Hallmark",
                         prefix_strip = "^HALLMARK_")
if (!is.null(p1)) ggsave(file.path(FIG_DIR, "GSEA_TC_Hallmark_dotplot.pdf"), p1, width = 10, height = 7)

p2 <- plot_gsea_dotplot(gsea_pr_hallmark, "Proteomics GSEA - Hallmark",
                         prefix_strip = "^HALLMARK_")
if (!is.null(p2)) ggsave(file.path(FIG_DIR, "GSEA_PR_Hallmark_dotplot.pdf"), p2, width = 10, height = 7)

# KEGG dotplots
p3 <- plot_gsea_dotplot(gsea_tc_kegg, "Transcriptomics GSEA - KEGG",
                         prefix_strip = "^KEGG_MEDICUS_")
if (!is.null(p3)) ggsave(file.path(FIG_DIR, "GSEA_TC_KEGG_dotplot.pdf"), p3, width = 10, height = 7)

p4 <- plot_gsea_dotplot(gsea_pr_kegg, "Proteomics GSEA - KEGG",
                         prefix_strip = "^KEGG_MEDICUS_")
if (!is.null(p4)) ggsave(file.path(FIG_DIR, "GSEA_PR_KEGG_dotplot.pdf"), p4, width = 10, height = 7)

# Reactome dotplots
p5 <- plot_gsea_dotplot(gsea_tc_reactome, "Transcriptomics GSEA - Reactome",
                         prefix_strip = "^REACTOME_")
if (!is.null(p5)) ggsave(file.path(FIG_DIR, "GSEA_TC_Reactome_dotplot.pdf"), p5, width = 11, height = 7)

p6 <- plot_gsea_dotplot(gsea_pr_reactome, "Proteomics GSEA - Reactome",
                         prefix_strip = "^REACTOME_")
if (!is.null(p6)) ggsave(file.path(FIG_DIR, "GSEA_PR_Reactome_dotplot.pdf"), p6, width = 11, height = 7)

# GO-BP dotplots
p7 <- plot_gsea_dotplot(gsea_tc_gobp, "Transcriptomics GSEA - GO Biological Process")
if (!is.null(p7)) ggsave(file.path(FIG_DIR, "GSEA_TC_GOBP_dotplot.pdf"), p7, width = 11, height = 7)

p8 <- plot_gsea_dotplot(gsea_pr_gobp, "Proteomics GSEA - GO Biological Process")
if (!is.null(p8)) ggsave(file.path(FIG_DIR, "GSEA_PR_GOBP_dotplot.pdf"), p8, width = 11, height = 7)

cat("  GSEA dotplots saved.\n")


# ============================================================================
# 4. CROSS-OMICS PATHWAY CONVERGENCE
# ============================================================================
cat("\n>>> 4. Cross-omics Pathway Convergence\n")

# ---- 4a. KEGG convergence ----
sig_tc_kegg <- gsea_tc_kegg$pathway[gsea_tc_kegg$padj < 0.25]
sig_pr_kegg <- gsea_pr_kegg$pathway[gsea_pr_kegg$padj < 0.25]
shared_kegg <- intersect(sig_tc_kegg, sig_pr_kegg)

cat(sprintf("  KEGG enriched (padj<0.25): TC=%d, PR=%d, shared=%d\n",
            length(sig_tc_kegg), length(sig_pr_kegg), length(shared_kegg)))

if (length(shared_kegg) > 0) {
  conv_df <- data.frame(
    pathway = shared_kegg,
    TC_NES = gsea_tc_kegg$NES[match(shared_kegg, gsea_tc_kegg$pathway)],
    TC_padj = gsea_tc_kegg$padj[match(shared_kegg, gsea_tc_kegg$pathway)],
    PR_NES = gsea_pr_kegg$NES[match(shared_kegg, gsea_pr_kegg$pathway)],
    PR_padj = gsea_pr_kegg$padj[match(shared_kegg, gsea_pr_kegg$pathway)],
    stringsAsFactors = FALSE
  )
  conv_df$concordant <- sign(conv_df$TC_NES) == sign(conv_df$PR_NES)
  conv_df <- conv_df[order(conv_df$TC_padj), ]
  write.csv(conv_df, file.path(OUT_DIR, "cross_omics_convergent_KEGG.csv"), row.names = FALSE)

  cat(sprintf("  KEGG concordant direction: %d / %d (%.1f%%)\n",
              sum(conv_df$concordant), nrow(conv_df),
              100 * sum(conv_df$concordant) / nrow(conv_df)))

  conv_df$pathway_short <- gsub("^KEGG_MEDICUS_", "", conv_df$pathway)
  conv_df$pathway_short <- gsub("_", " ", conv_df$pathway_short)
  conv_df$pathway_short <- substr(conv_df$pathway_short, 1, 45)

  p_conv <- ggplot(conv_df, aes(TC_NES, PR_NES)) +
    geom_point(aes(color = concordant), size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = pathway_short), size = 2.5, max.overlaps = 20) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey70") +
    scale_color_manual(values = c("TRUE" = "#00A087", "FALSE" = "#E64B35")) +
    labs(title = "Cross-omics KEGG Pathway Convergence",
         x = "Transcriptomics NES", y = "Proteomics NES",
         color = "Concordant") +
    theme_bw(base_size = 12) +
    coord_equal()
  ggsave(file.path(FIG_DIR, "cross_omics_KEGG_convergence.pdf"), p_conv,
         width = 9, height = 8)
}

# Venn for KEGG
if (length(sig_tc_kegg) > 0 || length(sig_pr_kegg) > 0) {
  pdf(file.path(FIG_DIR, "Venn_KEGG_enriched.pdf"), width = 6, height = 5)
  grid.newpage()
  draw.pairwise.venn(
    area1 = length(sig_tc_kegg), area2 = length(sig_pr_kegg),
    cross.area = length(shared_kegg),
    category = c("Transcriptomics", "Proteomics"),
    fill = c("#4DBBD5", "#E64B35"), alpha = 0.5,
    cat.cex = 1.2, cex = 1.5, fontfamily = "sans"
  )
  dev.off()
  cat("  KEGG Venn saved.\n")
}

# ---- 4b. Hallmark convergence ----
sig_tc_hm <- gsea_tc_hallmark$pathway[gsea_tc_hallmark$padj < 0.25]
sig_pr_hm <- gsea_pr_hallmark$pathway[gsea_pr_hallmark$padj < 0.25]
shared_hm <- intersect(sig_tc_hm, sig_pr_hm)

cat(sprintf("  Hallmark enriched (padj<0.25): TC=%d, PR=%d, shared=%d\n",
            length(sig_tc_hm), length(sig_pr_hm), length(shared_hm)))

if (length(shared_hm) > 0) {
  conv_hm <- data.frame(
    pathway = shared_hm,
    TC_NES = gsea_tc_hallmark$NES[match(shared_hm, gsea_tc_hallmark$pathway)],
    TC_padj = gsea_tc_hallmark$padj[match(shared_hm, gsea_tc_hallmark$pathway)],
    PR_NES = gsea_pr_hallmark$NES[match(shared_hm, gsea_pr_hallmark$pathway)],
    PR_padj = gsea_pr_hallmark$padj[match(shared_hm, gsea_pr_hallmark$pathway)]
  )
  conv_hm$concordant <- sign(conv_hm$TC_NES) == sign(conv_hm$PR_NES)
  write.csv(conv_hm, file.path(OUT_DIR, "cross_omics_convergent_Hallmark.csv"), row.names = FALSE)
  cat(sprintf("  Hallmark concordant: %d / %d\n", sum(conv_hm$concordant), nrow(conv_hm)))

  conv_hm$pathway_short <- gsub("^HALLMARK_", "", conv_hm$pathway)
  conv_hm$pathway_short <- gsub("_", " ", conv_hm$pathway_short)

  p_conv_hm <- ggplot(conv_hm, aes(TC_NES, PR_NES)) +
    geom_point(aes(color = concordant), size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = pathway_short), size = 2.5, max.overlaps = 25) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey70") +
    scale_color_manual(values = c("TRUE" = "#00A087", "FALSE" = "#E64B35")) +
    labs(title = "Cross-omics Hallmark Pathway Convergence",
         x = "Transcriptomics NES", y = "Proteomics NES", color = "Concordant") +
    theme_bw(base_size = 12) + coord_equal()
  ggsave(file.path(FIG_DIR, "cross_omics_Hallmark_convergence.pdf"), p_conv_hm,
         width = 9, height = 8)

  pdf(file.path(FIG_DIR, "Venn_Hallmark_enriched.pdf"), width = 6, height = 5)
  grid.newpage()
  draw.pairwise.venn(
    area1 = length(sig_tc_hm), area2 = length(sig_pr_hm),
    cross.area = length(shared_hm),
    category = c("Transcriptomics", "Proteomics"),
    fill = c("#4DBBD5", "#E64B35"), alpha = 0.5,
    cat.cex = 1.2, cex = 1.5, fontfamily = "sans"
  )
  dev.off()
  cat("  Hallmark convergence plots saved.\n")
}

# ---- 4c. GO-BP convergence ----
sig_tc_gobp <- gsea_tc_gobp$pathway[gsea_tc_gobp$padj < 0.25]
sig_pr_gobp <- gsea_pr_gobp$pathway[gsea_pr_gobp$padj < 0.25]
shared_gobp <- intersect(sig_tc_gobp, sig_pr_gobp)
cat(sprintf("  GO-BP enriched (padj<0.25): TC=%d, PR=%d, shared=%d\n",
            length(sig_tc_gobp), length(sig_pr_gobp), length(shared_gobp)))

if (length(shared_gobp) > 0) {
  conv_gobp <- data.frame(
    pathway = shared_gobp,
    TC_NES = gsea_tc_gobp$NES[match(shared_gobp, gsea_tc_gobp$pathway)],
    TC_padj = gsea_tc_gobp$padj[match(shared_gobp, gsea_tc_gobp$pathway)],
    PR_NES = gsea_pr_gobp$NES[match(shared_gobp, gsea_pr_gobp$pathway)],
    PR_padj = gsea_pr_gobp$padj[match(shared_gobp, gsea_pr_gobp$pathway)]
  )
  conv_gobp$concordant <- sign(conv_gobp$TC_NES) == sign(conv_gobp$PR_NES)
  write.csv(conv_gobp, file.path(OUT_DIR, "cross_omics_convergent_GOBP.csv"), row.names = FALSE)
  cat(sprintf("  GO-BP concordant: %d / %d\n", sum(conv_gobp$concordant), nrow(conv_gobp)))
}


# ============================================================================
# 5. COMBINED ENRICHMENT HEATMAP (NES across omics)
# ============================================================================
cat("\n>>> 5. Combined NES Heatmaps\n")

make_nes_heatmap <- function(gsea_tc, gsea_pr, sig_tc, sig_pr, label, prefix_strip, fig_dir) {
  all_pw <- union(sig_tc, sig_pr)
  if (length(all_pw) < 2) { cat(sprintf("  [skip] %s: <2 pathways\n", label)); return(NULL) }

  nes_mat <- data.frame(
    pathway = all_pw,
    TC_NES = gsea_tc$NES[match(all_pw, gsea_tc$pathway)],
    PR_NES = gsea_pr$NES[match(all_pw, gsea_pr$pathway)],
    stringsAsFactors = FALSE
  )
  rn <- gsub(prefix_strip, "", nes_mat$pathway)

  # If GO-BP, use pathway_name if available
  if ("pathway_name" %in% colnames(gsea_tc)) {
    tc_names <- setNames(gsea_tc$pathway_name, gsea_tc$pathway)
    pr_names <- setNames(gsea_pr$pathway_name, gsea_pr$pathway)
    all_names <- c(tc_names, pr_names)
    rn <- all_names[all_pw]
    rn[is.na(rn)] <- all_pw[is.na(rn)]
  }

  rn <- gsub("_", " ", rn)
  rn <- substr(rn, 1, 50)
  rn <- make.unique(rn)
  rownames(nes_mat) <- rn
  mat <- as.matrix(nes_mat[, c("TC_NES", "PR_NES")])
  colnames(mat) <- c("Transcriptomics", "Proteomics")
  mat <- mat[complete.cases(mat), , drop = FALSE]

  if (nrow(mat) > 50) mat <- mat[1:50, ]
  if (nrow(mat) < 2) { cat(sprintf("  [skip] %s: <2 complete rows\n", label)); return(NULL) }

  pdf(file.path(fig_dir, paste0("heatmap_NES_", label, "_cross_omics.pdf")),
      width = 7, height = max(5, nrow(mat) * 0.3 + 2))
  pheatmap(
    mat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_cols = FALSE, cluster_rows = TRUE,
    fontsize_row = 8,
    main = paste0(label, " NES (Transcriptomics vs Proteomics)")
  )
  dev.off()
  cat(sprintf("  %s NES heatmap saved (%d pathways).\n", label, nrow(mat)))
}

make_nes_heatmap(gsea_tc_kegg, gsea_pr_kegg, sig_tc_kegg, sig_pr_kegg,
                 "KEGG", "^KEGG_MEDICUS_", FIG_DIR)
make_nes_heatmap(gsea_tc_hallmark, gsea_pr_hallmark, sig_tc_hm, sig_pr_hm,
                 "Hallmark", "^HALLMARK_", FIG_DIR)
make_nes_heatmap(gsea_tc_gobp, gsea_pr_gobp, sig_tc_gobp, sig_pr_gobp,
                 "GOBP", "^GO:", FIG_DIR)


# ============================================================================
# 6. SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("Phase 2 SUMMARY\n")
cat("========================================\n")
cat("GSEA results (padj < 0.05 / padj < 0.25):\n\n")

for (label in c("Transcriptomics", "Proteomics")) {
  if (label == "Transcriptomics") {
    h <- gsea_tc_hallmark; k <- gsea_tc_kegg; r <- gsea_tc_reactome; g <- gsea_tc_gobp
  } else {
    h <- gsea_pr_hallmark; k <- gsea_pr_kegg; r <- gsea_pr_reactome; g <- gsea_pr_gobp
  }
  cat(sprintf("  %s:\n", label))
  cat(sprintf("    Hallmark: %d / %d (padj<0.05 / padj<0.25)\n",
              sum(h$padj<0.05), sum(h$padj<0.25)))
  cat(sprintf("    KEGG:     %d / %d\n", sum(k$padj<0.05), sum(k$padj<0.25)))
  cat(sprintf("    Reactome: %d / %d\n", sum(r$padj<0.05), sum(r$padj<0.25)))
  cat(sprintf("    GO-BP:    %d / %d\n",
              sum(g$padj<0.05, na.rm=TRUE), sum(g$padj<0.25, na.rm=TRUE)))
}

cat(sprintf("\nCross-omics convergent pathways:\n"))
cat(sprintf("  Hallmark: %d shared\n", length(shared_hm)))
cat(sprintf("  KEGG:     %d shared\n", length(shared_kegg)))
cat(sprintf("  GO-BP:    %d shared\n", length(shared_gobp)))
cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("========================================\n")
cat("Phase 2 COMPLETE.\n")
