#!/usr/bin/env Rscript
# ==============================================================================
# Enhancement M12: TCR/BCR Immune Repertoire Analysis for HAE Multi-omics Study
# ==============================================================================
# This script performs immune receptor repertoire analysis from bulk RNA-seq data
# Two analysis modes:
#   Mode 1: TRUST4 output available - full clonotype analysis with immunarch
#   Mode 2: Expression proxy - IG/TR gene family expression from transcriptome
# ==============================================================================

cat("=== Enhancement M12: TCR/BCR Immune Repertoire Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -------------- Setup and Configuration --------------
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(RColorBrewer)
    library(pheatmap)
    library(ComplexHeatmap)
    library(circlize)
    library(gridExtra)
})

# Set paths
BASE_DIR <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
DATA_DIR <- file.path(BASE_DIR, "analysis/data/processed")
RESULTS_DIR <- file.path(BASE_DIR, "analysis/results/enhancement_tcr_bcr")
IMMUNE_DIR <- file.path(BASE_DIR, "analysis/results/enhancement7_immune_deconvolution")
BAM_DIR <- "/Volumes/阿热-实验数据备份/工作/1.肝包虫/1.泡型包虫/1.2024-11-泡球蚴/1.转录组学测序/02.Bam"
TRUST4_OUTPUT <- file.path(RESULTS_DIR, "trust4_output")

# Create output directories
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TRUST4_OUTPUT, recursive = TRUE, showWarnings = FALSE)

# Color scheme
COLORS <- list(
    group = c("Adjacent" = "#E41A1C", "Normal" = "#377EB8"),
    subtype = c("DS1" = "#FF7F00", "DS2" = "#4DAF4A"),
    chain = c("TRA" = "#E41A1C", "TRB" = "#377EB8", "TRD" = "#4DAF4A", "TRG" = "#984EA3",
              "IGH" = "#FF7F00", "IGK" = "#FFFF33", "IGL" = "#A65628")
)

# -------------- Helper Functions --------------
run_trust4_on_bam <- function(bam_file, output_prefix, ref_dir = NULL) {
    # Check if TRUST4 is available
    trust4_path <- "/Users/rishat/miniforge3/envs/multiomics/bin/trust4"
    
    if (!file.exists(trust4_path)) {
        warning("TRUST4 not found at expected path")
        return(FALSE)
    }
    
    # Reference files location (downloaded via download_trust4_references.sh)
    ref_genome <- file.path(BASE_DIR, "analysis/data/reference/trust4/hg38_bcrtcr.fa")
    ref_imgt <- file.path(BASE_DIR, "analysis/data/reference/trust4/human_IMGT+C.fa")
    
    if (!file.exists(ref_genome) || !file.exists(ref_imgt)) {
        warning("TRUST4 reference files not found. Run download_trust4_references.sh first.")
        return(FALSE)
    }
    
    # TRUST4 command
    cmd <- sprintf("%s -b %s -f %s --ref %s -o %s --od %s -t 4 2>&1",
                   trust4_path, bam_file, ref_genome, ref_imgt,
                   basename(output_prefix), dirname(output_prefix))
    
    cat("Running:", cmd, "\n")
    result <- system(cmd, intern = TRUE)
    return(TRUE)
}

parse_trust4_report <- function(report_file) {
    if (!file.exists(report_file)) return(NULL)
    # TRUST4 v1.1.9 output has 10 columns: count, frequency, CDR3nt, CDR3aa, V, D, J, C, cid, cid_full_length
    df <- read.table(report_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     comment.char = "#", fill = TRUE,
                     col.names = c("count", "frequency", "CDR3_nt", "CDR3_aa", "V", "D", "J", "C", "cid", "full_length"))
    # Filter out any rows with NA count (from malformed lines)
    df <- df[!is.na(df$count) & df$count > 0, ]
    return(df)
}

extract_ig_tr_genes <- function(expr_matrix, gene_annotation = NULL) {
    # Extract immunoglobulin (IG) and T cell receptor (TR) gene families
    # Gene naming patterns: IGHV, IGKV, IGLV, TRAV, TRBV, TRDV, TRGV, etc.
    
    gene_ids <- rownames(expr_matrix)
    
    # If no annotation provided, try to load from DEG file
    if (is.null(gene_annotation)) {
        deg_file <- file.path(BASE_DIR, "analysis/results/phase1_diff/DEGs_Adjacent_vs_Normal.csv")
        if (file.exists(deg_file)) {
            deg_data <- read.csv(deg_file)
            gene_annotation <- data.frame(
                gene_id = deg_data$gene_id,
                gene_name = deg_data$gene_name,
                stringsAsFactors = FALSE
            )
            cat("Loaded gene annotation from DEG file:", nrow(gene_annotation), "genes\n")
        }
    }
    
    # Create mapping from gene_id to gene_name
    if (!is.null(gene_annotation)) {
        id_to_name <- setNames(gene_annotation$gene_name, gene_annotation$gene_id)
        gene_names <- sapply(gene_ids, function(x) {
            if (x %in% names(id_to_name)) id_to_name[x] else x
        })
    } else {
        gene_names <- gene_ids
    }
    
    # Define gene family patterns
    ig_patterns <- c("^IGH[VDJCGAMDE]", "^IGK[VJC]", "^IGL[LVJC]", "^IGLL")
    tr_patterns <- c("^TRA[VJC]", "^TRB[VDJC]", "^TRD[VJC]", "^TRG[VJC]")
    
    ig_idx <- grep(paste(ig_patterns, collapse = "|"), gene_names)
    tr_idx <- grep(paste(tr_patterns, collapse = "|"), gene_names)
    
    cat("  Pattern matching found", length(ig_idx), "IG genes and", length(tr_idx), "TR genes\n")
    
    list(
        IG_genes = gene_ids[ig_idx],
        TR_genes = gene_ids[tr_idx],
        IG_names = gene_names[ig_idx],
        TR_names = gene_names[tr_idx],
        all_gene_names = gene_names
    )
}

calculate_diversity_metrics <- function(clonotype_counts) {
    # Shannon entropy
    freq <- clonotype_counts / sum(clonotype_counts)
    freq <- freq[freq > 0]
    shannon <- -sum(freq * log2(freq))
    
    # Simpson index (1 - D)
    simpson <- 1 - sum(freq^2)
    
    # Inverse Simpson
    inv_simpson <- 1 / sum(freq^2)
    
    # Clonality (normalized entropy)
    n <- length(freq)
    max_entropy <- if(n > 1) log2(n) else 1
    clonality <- 1 - (shannon / max_entropy)
    
    # Chao1 estimate (for richness)
    f1 <- sum(clonotype_counts == 1)
    f2 <- sum(clonotype_counts == 2)
    observed <- length(clonotype_counts)
    chao1 <- observed + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    
    # Top clone frequency
    top1_freq <- max(freq)
    top10_freq <- sum(sort(freq, decreasing = TRUE)[1:min(10, length(freq))])
    
    c(shannon = shannon, simpson = simpson, inv_simpson = inv_simpson,
      clonality = clonality, chao1 = chao1, richness = observed,
      top1_freq = top1_freq, top10_freq = top10_freq)
}

# -------------- Check External HDD and TRUST4 Reference Files --------------
cat("=== Checking data availability ===\n")

bam_available <- dir.exists(BAM_DIR) && length(list.files(BAM_DIR, pattern = "\\.bam$")) > 0

# Check for reference files in project directory
trust4_ref_dir <- file.path(BASE_DIR, "analysis/data/reference/trust4")
trust4_ref_available <- file.exists(file.path(trust4_ref_dir, "hg38_bcrtcr.fa")) && 
                        file.exists(file.path(trust4_ref_dir, "human_IMGT+C.fa"))

cat("External HDD BAM files available:", bam_available, "\n")
cat("TRUST4 reference files available:", trust4_ref_available, "\n")
if (!trust4_ref_available) {
    cat("  To download references, run: bash analysis/scripts/download_trust4_references.sh\n")
}

# Check if TRUST4 output already exists
trust4_results_exist <- length(list.files(TRUST4_OUTPUT, pattern = "_report.tsv$")) > 0
cat("TRUST4 results already exist:", trust4_results_exist, "\n\n")

# Determine analysis mode
ANALYSIS_MODE <- if (trust4_results_exist || (bam_available && trust4_ref_available)) {
    "TRUST4"
} else {
    "EXPRESSION_PROXY"
}
cat("Analysis mode:", ANALYSIS_MODE, "\n\n")

# -------------- MODE 1: TRUST4 Analysis --------------
if (ANALYSIS_MODE == "TRUST4" && !trust4_results_exist && bam_available && trust4_ref_available) {
    cat("=== Running TRUST4 on BAM files ===\n")
    cat("Note: This may take several hours for all samples\n")
    cat("TRUST4 extracts TCR/BCR sequences from aligned reads\n\n")
    
    bam_files <- list.files(BAM_DIR, pattern = "\\.bam$", full.names = TRUE)
    bam_files <- bam_files[!grepl("\\.bai$", bam_files)]
    
    cat("Found", length(bam_files), "BAM files\n")
    
    for (bam_file in bam_files) {
        sample_name <- gsub("\\.bam$", "", basename(bam_file))
        output_prefix <- file.path(TRUST4_OUTPUT, sample_name)
        
        # Check if already processed
        if (file.exists(paste0(output_prefix, "_report.tsv"))) {
            cat("  Skipping", sample_name, "- already processed\n")
            next
        }
        
        cat("  Processing:", sample_name, "\n")
        tryCatch({
            run_trust4_on_bam(bam_file, output_prefix)
        }, error = function(e) {
            cat("    Error:", conditionMessage(e), "\n")
        })
    }
    
    trust4_results_exist <- length(list.files(TRUST4_OUTPUT, pattern = "_report.tsv$")) > 0
}

# -------------- Load Expression Data for Proxy Analysis --------------
cat("=== Loading transcriptomics data ===\n")

expr_file <- file.path(DATA_DIR, "transcriptomics_norm_counts.csv")
expr_matrix <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")

# Sample metadata
sample_info <- data.frame(
    sample = colnames(expr_matrix),
    group = ifelse(grepl("^Normal", colnames(expr_matrix)), "Normal", "Adjacent"),
    patient_id = as.numeric(gsub("Normal|Adjacent", "", colnames(expr_matrix)))
)
rownames(sample_info) <- sample_info$sample

# Try to load subtype information
subtype_file <- file.path(BASE_DIR, "analysis/results/phase6_subtyping/subtype_assignments.csv")
if (file.exists(subtype_file)) {
    subtype_data <- read.csv(subtype_file)
    sample_info$subtype <- NA
    for (i in 1:nrow(subtype_data)) {
        sample_info$subtype[sample_info$sample == subtype_data$sample[i]] <- subtype_data$subtype[i]
    }
    cat("Subtype information loaded\n")
} else {
    sample_info$subtype <- NA
}

# -------------- Extract IG/TR Gene Expression (for both modes) --------------
cat("\n=== Extracting IG/TR gene families ===\n")

ig_tr_genes <- extract_ig_tr_genes(expr_matrix)
cat("Found", length(ig_tr_genes$IG_genes), "IG genes\n")
cat("Found", length(ig_tr_genes$TR_genes), "TR genes\n")

# Check if we have any IG/TR genes
has_ig_tr_genes <- (length(ig_tr_genes$IG_genes) > 0 || length(ig_tr_genes$TR_genes) > 0)

if (!has_ig_tr_genes) {
    cat("\nWARNING: No IG/TR genes found in expression matrix.\n")
    cat("This may indicate:\n")
    cat("  1. Different gene naming convention in the data\n")
    cat("  2. IG/TR genes were filtered out during preprocessing\n")
    cat("  3. Expression levels below detection threshold\n")
    cat("\nProceeding with alternative analysis using all immune-related genes...\n")
}

if (has_ig_tr_genes) {
    # Combine all immune receptor genes
    all_ir_genes <- c(ig_tr_genes$IG_genes, ig_tr_genes$TR_genes)
    all_ir_names <- c(ig_tr_genes$IG_names, ig_tr_genes$TR_names)
    
    ir_expr <- expr_matrix[all_ir_genes, , drop = FALSE]
    # Use gene names as rownames for easier interpretation
    rownames(ir_expr) <- all_ir_names
    
    # Log2 transform for visualization
    ir_expr_log <- log2(ir_expr + 1)
    
    # Classify genes by family - use the gene names directly
    gene_family <- sapply(all_ir_names, function(g) {
        if (grepl("^IGHV|^IGHD|^IGHJ|^IGHC|^IGHG|^IGHA|^IGHM|^IGHE", g)) return("IGH")
        if (grepl("^IGKV|^IGKJ|^IGKC", g)) return("IGK")
        if (grepl("^IGLV|^IGLJ|^IGLC|^IGLL", g)) return("IGL")
        if (grepl("^TRAV|^TRAJ|^TRAC", g)) return("TRA")
        if (grepl("^TRBV|^TRBD|^TRBJ|^TRBC", g)) return("TRB")
        if (grepl("^TRDV|^TRDD|^TRDJ|^TRDC", g)) return("TRD")
        if (grepl("^TRGV|^TRGJ|^TRGC", g)) return("TRG")
        return("Other")
    })
    names(gene_family) <- all_ir_names
    
    # Remove "Other" category
    valid_genes <- names(gene_family)[gene_family != "Other"]
    gene_family <- gene_family[valid_genes]
    
    cat("Gene family classification:\n")
    print(table(gene_family))
    
    # Calculate family-level expression summary
    family_expr <- data.frame(
        sample = colnames(ir_expr),
        group = sample_info[colnames(ir_expr), "group"],
        patient_id = sample_info[colnames(ir_expr), "patient_id"]
    )
    
    for (fam in unique(gene_family)) {
        fam_genes <- names(gene_family)[gene_family == fam]
        if (length(fam_genes) > 0) {
            # ir_expr_log now uses gene names as rownames
            fam_expr_vals <- ir_expr_log[fam_genes, , drop = FALSE]
            family_expr[[fam]] <- colMeans(fam_expr_vals, na.rm = TRUE)
        }
    }
    
    # Save family expression data
    write.csv(family_expr, file.path(RESULTS_DIR, "ig_tr_family_expression.csv"), row.names = FALSE)
    cat("Saved IG/TR family expression summary\n")
    
    # -------------- Statistical Tests for IG/TR Expression --------------
    cat("\n=== Statistical analysis of IG/TR expression ===\n")
    
    # Paired comparison between Adjacent and Normal
    paired_samples <- intersect(
        gsub("Adjacent", "", sample_info$sample[sample_info$group == "Adjacent"]),
        gsub("Normal", "", sample_info$sample[sample_info$group == "Normal"])
    )
    cat("Paired samples:", length(paired_samples), "\n")
    
    stats_results <- data.frame()
    for (fam in setdiff(colnames(family_expr), c("sample", "group", "patient_id"))) {
        adj_vals <- family_expr[[fam]][family_expr$group == "Adjacent"]
        norm_vals <- family_expr[[fam]][family_expr$group == "Normal"]
        
        # Paired test - match by patient_id
        adj_data <- family_expr[family_expr$group == "Adjacent", c("patient_id", fam)]
        norm_data <- family_expr[family_expr$group == "Normal", c("patient_id", fam)]
        
        # Merge to get paired data
        paired_data <- merge(adj_data, norm_data, by = "patient_id", suffixes = c("_adj", "_norm"))
        
        if (nrow(paired_data) >= 3 && 
            sum(!is.na(paired_data[[paste0(fam, "_adj")]])) >= 3 &&
            sum(!is.na(paired_data[[paste0(fam, "_norm")]])) >= 3) {
            
            tryCatch({
                test_result <- wilcox.test(paired_data[[paste0(fam, "_adj")]], 
                                           paired_data[[paste0(fam, "_norm")]], 
                                           paired = TRUE)
                
                stats_results <- rbind(stats_results, data.frame(
                    family = fam,
                    mean_adjacent = mean(adj_vals, na.rm = TRUE),
                    mean_normal = mean(norm_vals, na.rm = TRUE),
                    log2FC = mean(adj_vals - norm_vals, na.rm = TRUE),
                    p_value = test_result$p.value,
                    n_pairs = nrow(paired_data),
                    direction = ifelse(mean(adj_vals, na.rm = TRUE) > mean(norm_vals, na.rm = TRUE), 
                                      "Up in Adjacent", "Down in Adjacent")
                ))
            }, error = function(e) {
                cat("  Warning: Could not test family", fam, "-", conditionMessage(e), "\n")
            })
        }
    }
    
    if (nrow(stats_results) > 0) {
        stats_results$p_adj <- p.adjust(stats_results$p_value, method = "BH")
        stats_results <- stats_results[order(stats_results$p_value), ]
        write.csv(stats_results, file.path(RESULTS_DIR, "ig_tr_differential_stats.csv"), row.names = FALSE)
        cat("Saved differential expression statistics\n")
        print(stats_results)
    }
}

# -------------- MODE 2: Expression Proxy Analysis --------------
if (ANALYSIS_MODE == "EXPRESSION_PROXY" && has_ig_tr_genes) {
    cat("\n=== Expression Proxy Analysis ===\n")
    
    # Calculate pseudo-diversity based on gene expression patterns
    # Use variance and entropy of IG/TR gene expression as diversity proxy
    
    diversity_proxy <- data.frame(
        sample = colnames(ir_expr),
        group = sample_info[colnames(ir_expr), "group"],
        patient_id = sample_info[colnames(ir_expr), "patient_id"]
    )
    
    for (i in 1:ncol(ir_expr)) {
        sample <- colnames(ir_expr)[i]
        expr_vals <- ir_expr[, i]
        expr_vals <- expr_vals[expr_vals > 0]  # Non-zero expression
        
        if (length(expr_vals) > 1) {
            # Treat expression values as pseudo-abundances
            metrics <- calculate_diversity_metrics(expr_vals)
            diversity_proxy$shannon[i] <- metrics["shannon"]
            diversity_proxy$simpson[i] <- metrics["simpson"]
            diversity_proxy$richness[i] <- metrics["richness"]
            diversity_proxy$clonality[i] <- metrics["clonality"]
        }
    }
    
    write.csv(diversity_proxy, file.path(RESULTS_DIR, "expression_diversity_proxy.csv"), row.names = FALSE)
    cat("Saved expression-based diversity proxy metrics\n")
}

# -------------- TRUST4 Results Processing (if available) --------------
if (ANALYSIS_MODE == "TRUST4" && trust4_results_exist) {
    cat("\n=== Processing TRUST4 Results ===\n")
    
    suppressPackageStartupMessages(library(immunarch))
    
    trust4_files <- list.files(TRUST4_OUTPUT, pattern = "_report.tsv$", full.names = TRUE)
    cat("Found", length(trust4_files), "TRUST4 report files\n")
    
    # Parse all TRUST4 reports
    all_clonotypes <- list()
    diversity_metrics <- data.frame()
    
    for (f in trust4_files) {
        sample_name <- gsub("_report.tsv$", "", basename(f))
        df <- parse_trust4_report(f)
        
        if (!is.null(df) && nrow(df) > 0) {
            all_clonotypes[[sample_name]] <- df
            
            # Calculate diversity for each chain type
            for (chain in c("TRA", "TRB", "IGH", "IGK", "IGL")) {
                chain_df <- df[grepl(paste0("^", chain), df$V), ]
                if (nrow(chain_df) > 0) {
                    metrics <- calculate_diversity_metrics(chain_df$count)
                    diversity_metrics <- rbind(diversity_metrics, data.frame(
                        sample = sample_name,
                        chain = chain,
                        t(metrics)
                    ))
                }
            }
        }
    }
    
    if (nrow(diversity_metrics) > 0) {
        diversity_metrics$group <- ifelse(grepl("^Normal", diversity_metrics$sample), "Normal", "Adjacent")
        diversity_metrics$patient_id <- as.numeric(gsub("Normal|Adjacent", "", diversity_metrics$sample))
        
        write.csv(diversity_metrics, file.path(RESULTS_DIR, "trust4_diversity_metrics.csv"), row.names = FALSE)
        cat("Saved TRUST4 diversity metrics\n")
    }
}

# -------------- Load Immune Deconvolution Results for Correlation --------------
cat("\n=== Loading immune deconvolution results ===\n")

immune_file <- file.path(IMMUNE_DIR, "ssGSEA_immune_scores.csv")
if (file.exists(immune_file)) {
    immune_scores <- read.csv(immune_file, row.names = 1)
    cat("Loaded immune scores:", nrow(immune_scores), "cell types\n")
    
    # Prepare for correlation with IG/TR expression
    if (exists("family_expr") && nrow(family_expr) > 0) {
        # Match samples
        common_samples <- intersect(colnames(immune_scores), family_expr$sample)
        cat("Common samples:", length(common_samples), "\n")
        
        if (length(common_samples) > 5) {
            # Correlation analysis
            ir_families <- setdiff(colnames(family_expr), c("sample", "group", "patient_id"))
            immune_cells <- rownames(immune_scores)
            
            cor_matrix <- matrix(NA, nrow = length(ir_families), ncol = length(immune_cells),
                                 dimnames = list(ir_families, immune_cells))
            pval_matrix <- cor_matrix
            
            for (ir in ir_families) {
                for (cell in immune_cells) {
                    ir_vals <- family_expr[[ir]][match(common_samples, family_expr$sample)]
                    cell_vals <- as.numeric(immune_scores[cell, common_samples])
                    
                    if (sum(!is.na(ir_vals) & !is.na(cell_vals)) > 5) {
                        test <- cor.test(ir_vals, cell_vals, method = "spearman")
                        cor_matrix[ir, cell] <- test$estimate
                        pval_matrix[ir, cell] <- test$p.value
                    }
                }
            }
            
            # Save correlation results
            cor_df <- as.data.frame(cor_matrix)
            cor_df$IR_family <- rownames(cor_df)
            write.csv(cor_df, file.path(RESULTS_DIR, "ir_immune_correlation.csv"), row.names = FALSE)
            cat("Saved IR-immune correlation matrix\n")
        }
    }
} else {
    cat("Immune deconvolution file not found\n")
    immune_scores <- NULL
}

# -------------- Visualization --------------
cat("\n=== Generating visualizations ===\n")

# 1. IG/TR Family Expression Heatmap
if (exists("ir_expr_log") && has_ig_tr_genes && nrow(ir_expr_log) > 0) {
    cat("Creating IG/TR expression heatmap...\n")
    
    # Filter for highly variable genes
    gene_vars <- apply(ir_expr_log, 1, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(50, length(gene_vars))]
    
    # Annotation for columns
    col_anno <- HeatmapAnnotation(
        Group = sample_info[colnames(ir_expr_log), "group"],
        col = list(Group = COLORS$group),
        annotation_name_side = "left"
    )
    
    # Row annotation for gene family
    row_family <- gene_family[top_genes]
    row_anno <- rowAnnotation(
        Family = row_family,
        col = list(Family = COLORS$chain),
        annotation_name_side = "top"
    )
    
    pdf(file.path(RESULTS_DIR, "ig_tr_expression_heatmap.pdf"), width = 12, height = 10)
    ht <- Heatmap(
        as.matrix(ir_expr_log[top_genes, ]),
        name = "log2(CPM+1)",
        col = colorRamp2(c(0, 5, 10), c("blue", "white", "red")),
        top_annotation = col_anno,
        left_annotation = row_anno,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        column_title = "IG/TR Gene Expression in HAE Samples"
    )
    draw(ht)
    dev.off()
    
    # PNG version
    png(file.path(RESULTS_DIR, "ig_tr_expression_heatmap.png"), width = 1200, height = 1000, res = 120)
    draw(ht)
    dev.off()
}

# 2. Family Expression Boxplots
if (exists("family_expr") && has_ig_tr_genes && nrow(family_expr) > 0) {
    cat("Creating family expression boxplots...\n")
    
    # Get only numeric columns (families)
    family_cols <- setdiff(colnames(family_expr), c("sample", "group", "patient_id"))
    
    if (length(family_cols) > 0) {
        family_long <- family_expr %>%
            pivot_longer(cols = all_of(family_cols),
                         names_to = "family", values_to = "expression") %>%
            filter(!is.na(expression))
        
        # Separate TCR and BCR
        family_long$receptor_type <- ifelse(grepl("^IG", family_long$family), "BCR (Immunoglobulin)", "TCR")
        
        # Check if we have both types
        if (length(unique(family_long$receptor_type)) > 1) {
            p <- ggplot(family_long, aes(x = family, y = expression, fill = group)) +
                geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
                geom_point(aes(color = group), position = position_jitterdodge(jitter.width = 0.1),
                           size = 1, alpha = 0.6) +
                facet_wrap(~receptor_type, scales = "free_x", ncol = 2) +
                scale_fill_manual(values = COLORS$group) +
                scale_color_manual(values = COLORS$group) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = "top") +
                labs(title = "IG/TR Gene Family Expression: Adjacent vs Normal",
                     subtitle = "HAE Multi-omics Study - Expression Proxy for Immune Repertoire Activity",
                     x = "Gene Family", y = "Mean log2(CPM+1)")
        } else {
            # Single receptor type
            p <- ggplot(family_long, aes(x = family, y = expression, fill = group)) +
                geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
                geom_point(aes(color = group), position = position_jitterdodge(jitter.width = 0.1),
                           size = 1, alpha = 0.6) +
                scale_fill_manual(values = COLORS$group) +
                scale_color_manual(values = COLORS$group) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = "top") +
                labs(title = "IG/TR Gene Family Expression: Adjacent vs Normal",
                     subtitle = "HAE Multi-omics Study - Expression Proxy for Immune Repertoire Activity",
                     x = "Gene Family", y = "Mean log2(CPM+1)")
        }
        
        ggsave(file.path(RESULTS_DIR, "ig_tr_family_boxplot.pdf"), p, width = 10, height = 6)
        ggsave(file.path(RESULTS_DIR, "ig_tr_family_boxplot.png"), p, width = 10, height = 6, dpi = 150)
    }
}

# 3. Paired comparison plot
if (exists("family_expr") && has_ig_tr_genes && length(paired_samples) > 0) {
    cat("Creating paired comparison plot...\n")
    
    family_cols <- setdiff(colnames(family_expr), c("sample", "group", "patient_id"))
    
    if (length(family_cols) > 0) {
        paired_data <- family_expr %>%
            filter(patient_id %in% as.numeric(paired_samples)) %>%
            pivot_longer(cols = all_of(family_cols),
                         names_to = "family", values_to = "expression")
        
        if (nrow(paired_data) > 0) {
            p_paired <- ggplot(paired_data, aes(x = group, y = expression, group = patient_id)) +
                geom_line(alpha = 0.3, color = "gray50") +
                geom_point(aes(color = group), size = 2) +
                facet_wrap(~family, scales = "free_y", ncol = 4) +
                scale_color_manual(values = COLORS$group) +
                theme_bw() +
                theme(legend.position = "top") +
                labs(title = "Paired Comparison of IG/TR Expression",
                     subtitle = "Lines connect matched Adjacent-Normal pairs",
                     x = "", y = "Mean log2(CPM+1)")
            
            ggsave(file.path(RESULTS_DIR, "ig_tr_paired_comparison.pdf"), p_paired, width = 12, height = 8)
            ggsave(file.path(RESULTS_DIR, "ig_tr_paired_comparison.png"), p_paired, width = 12, height = 8, dpi = 150)
        }
    }
}

# 4. Correlation with Immune Infiltration
if (exists("cor_matrix") && !is.null(immune_scores)) {
    cat("Creating correlation heatmap...\n")
    
    # Filter to significant correlations
    sig_mask <- pval_matrix < 0.05
    cor_matrix_filtered <- cor_matrix
    cor_matrix_filtered[!sig_mask] <- 0
    
    # Select immune cells with at least one significant correlation
    sig_cells <- colnames(cor_matrix_filtered)[colSums(sig_mask, na.rm = TRUE) > 0]
    
    if (length(sig_cells) > 0) {
        pdf(file.path(RESULTS_DIR, "ir_immune_correlation_heatmap.pdf"), width = 12, height = 6)
        ht_cor <- Heatmap(
            cor_matrix[, sig_cells, drop = FALSE],
            name = "Spearman r",
            col = colorRamp2(c(-0.8, 0, 0.8), c("blue", "white", "red")),
            cell_fun = function(j, i, x, y, w, h, fill) {
                if (!is.na(pval_matrix[i, sig_cells[j]]) && pval_matrix[i, sig_cells[j]] < 0.05) {
                    grid.text("*", x, y, gp = gpar(fontsize = 14))
                }
            },
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_names_gp = gpar(fontsize = 8),
            row_names_gp = gpar(fontsize = 10),
            column_title = "Correlation: IG/TR Expression vs Immune Cell Infiltration\n(* p < 0.05)"
        )
        draw(ht_cor)
        dev.off()
        
        png(file.path(RESULTS_DIR, "ir_immune_correlation_heatmap.png"), width = 1200, height = 600, res = 120)
        draw(ht_cor)
        dev.off()
    }
}

# 5. Diversity Proxy Visualization
if (exists("diversity_proxy") && nrow(diversity_proxy) > 0) {
    cat("Creating diversity proxy plots...\n")
    
    div_long <- diversity_proxy %>%
        pivot_longer(cols = c(shannon, simpson, richness, clonality),
                     names_to = "metric", values_to = "value") %>%
        filter(!is.na(value))
    
    p_div <- ggplot(div_long, aes(x = group, y = value, fill = group)) +
        geom_boxplot(alpha = 0.8) +
        geom_point(position = position_jitter(width = 0.1), size = 2, alpha = 0.6) +
        facet_wrap(~metric, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = COLORS$group) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(title = "Expression-based Diversity Proxy Metrics",
             subtitle = "Based on IG/TR gene expression patterns",
             x = "", y = "Metric Value")
    
    ggsave(file.path(RESULTS_DIR, "diversity_proxy_boxplot.pdf"), p_div, width = 10, height = 4)
    ggsave(file.path(RESULTS_DIR, "diversity_proxy_boxplot.png"), p_div, width = 10, height = 4, dpi = 150)
}

# 6. Summary Statistics Plot
if (exists("stats_results") && nrow(stats_results) > 0) {
    cat("Creating summary statistics plot...\n")
    
    stats_results$significance <- ifelse(stats_results$p_adj < 0.05, "Significant", "Not Significant")
    
    p_stats <- ggplot(stats_results, aes(x = reorder(family, -log10(p_value)), y = -log10(p_value))) +
        geom_bar(aes(fill = direction), stat = "identity", alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        coord_flip() +
        scale_fill_manual(values = c("Up in Adjacent" = "#E41A1C", "Down in Adjacent" = "#377EB8")) +
        theme_bw() +
        labs(title = "Differential IG/TR Gene Family Expression",
             subtitle = "Adjacent vs Normal (Paired Wilcoxon test)",
             x = "Gene Family", y = "-log10(p-value)",
             fill = "Direction")
    
    ggsave(file.path(RESULTS_DIR, "differential_ig_tr_barplot.pdf"), p_stats, width = 8, height = 5)
    ggsave(file.path(RESULTS_DIR, "differential_ig_tr_barplot.png"), p_stats, width = 8, height = 5, dpi = 150)
}

# -------------- Generate Summary Report --------------
cat("\n=== Generating Summary Report ===\n")

output_files <- list.files(RESULTS_DIR, full.names = FALSE)

summary_text <- paste0(
    "=================================================================\n",
    "Enhancement M12: TCR/BCR Immune Repertoire Analysis - Summary\n",
    "=================================================================\n\n",
    "Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "ANALYSIS MODE: ", ANALYSIS_MODE, "\n",
    if (ANALYSIS_MODE == "EXPRESSION_PROXY") {
        paste0(
            "  - External HDD BAM files: ", ifelse(bam_available, "Available", "Not mounted"), "\n",
            "  - TRUST4 reference files: ", ifelse(trust4_ref_available, "Available", "Not found"), "\n",
            "  - Using expression proxy from transcriptome data\n\n"
        )
    } else {
        "  - Full TRUST4 clonotype analysis\n\n"
    },
    "DATA SUMMARY:\n",
    "  - Expression matrix: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " samples\n",
    "  - IG genes found: ", length(ig_tr_genes$IG_genes), "\n",
    "  - TR genes found: ", length(ig_tr_genes$TR_genes), "\n",
    "  - Paired samples: ", length(paired_samples), "\n\n",
    "KEY FINDINGS:\n"
)

if (exists("stats_results") && nrow(stats_results) > 0) {
    sig_families <- stats_results$family[stats_results$p_adj < 0.05]
    summary_text <- paste0(summary_text,
        "  Differentially expressed IG/TR families (FDR < 0.05): ", length(sig_families), "\n")
    if (length(sig_families) > 0) {
        for (fam in sig_families) {
            row <- stats_results[stats_results$family == fam, ]
            summary_text <- paste0(summary_text,
                "    - ", fam, ": ", row$direction, " (log2FC=", round(row$log2FC, 2), 
                ", p.adj=", format(row$p_adj, digits = 3), ")\n")
        }
    }
}

summary_text <- paste0(summary_text,
    "\nOUTPUT FILES:\n")
for (f in output_files) {
    summary_text <- paste0(summary_text, "  - ", f, "\n")
}

summary_text <- paste0(summary_text,
    "\nINTERPRETATION:\n",
    "  This analysis provides a proxy assessment of immune receptor repertoire\n",
    "  activity in HAE (Hepatic Alveolar Echinococcosis) based on bulk RNA-seq.\n",
    "  While full clonotype-level analysis requires specialized tools like TRUST4,\n",
    "  the expression patterns of IG/TR gene families provide valuable insights\n",
    "  into B cell and T cell receptor activity in peri-lesional vs normal tissue.\n\n",
    "  Key biological implications:\n",
    "  - IGH/IGK/IGL expression reflects B cell receptor (BCR) activity and\n",
    "    antibody-mediated immune responses\n",
    "  - TRA/TRB expression reflects T cell receptor (TCR) diversity and\n",
    "    T cell-mediated immunity\n",
    "  - Changes between Adjacent and Normal tissue indicate local immune\n",
    "    activation or suppression around parasitic lesions\n\n",
    "=================================================================\n"
)

cat(summary_text)
writeLines(summary_text, file.path(RESULTS_DIR, "analysis_summary.txt"))

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", RESULTS_DIR, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
