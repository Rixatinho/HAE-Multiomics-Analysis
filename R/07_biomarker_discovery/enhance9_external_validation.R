#!/usr/bin/env Rscript
# Enhancement 9: External Validation using GSE124362
# Efficient strategy: PDict-based bulk probe-to-gene mapping, then validation
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biostrings)
  library(rentrez)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(ggpubr)
})

PROJECT <- normalizePath("~/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学")
DATA_DIR <- file.path(PROJECT, "analysis/data/processed")
OUT_DIR <- file.path(PROJECT, "analysis/results/enhancement9_external_validation")
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

MAPPING_FILE <- file.path(OUT_DIR, "probe_gene_mapping.csv")

cat("===== Enhancement 9: External Validation (GSE124362) =====\n")

# ============================================================
# PART A: Load GSE124362
# ============================================================
cat("\n--- Part A: Load GSE124362 ---\n")
gse <- getGEO("GSE124362", destdir = "/tmp", GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]
expr_mat <- exprs(eset)
pd <- pData(eset)

peri_idx <- grep("Periparasitic", pd$source_name_ch1)
dist_idx <- grep("distal", pd$source_name_ch1)
patient_ids <- gsub(".*(rep\\d+).*", "\\1", pd$source_name_ch1)
cat(sprintf("Samples: %d periparasitic, %d distal, %d patients\n",
            length(peri_idx), length(dist_idx), length(unique(patient_ids))))

# ============================================================
# PART B: Probe-to-gene mapping via PDict + RefSeq sequences
# ============================================================
if (file.exists(MAPPING_FILE) && file.size(MAPPING_FILE) > 500) {
  cat("\n--- Part B: Loading cached probe-gene mapping ---\n")
  mapping_df <- read.csv(MAPPING_FILE, stringsAsFactors = FALSE)
  cat(sprintf("Loaded: %d probes -> %d genes\n", nrow(mapping_df), length(unique(mapping_df$gene))))
} else {
  cat("\n--- Part B: Build probe-gene mapping (PDict approach) ---\n")

  gpl <- getGEO("GPL16956", destdir = "/tmp")
  gpl_tab <- Table(gpl)
  coding <- gpl_tab[gpl_tab$TRANSCRIPT_TYPE == "protein_coding", ]
  cat(sprintf("Protein-coding probes: %d\n", nrow(coding)))

  # Clean sequences: remove any with non-ACGT characters
  valid <- grepl("^[ACGTacgt]+$", coding$SEQUENCE) & nchar(coding$SEQUENCE) >= 55
  coding <- coding[valid, ]
  cat(sprintf("Valid probes after filtering: %d\n", nrow(coding)))

  # All probes must be same length for PDict; pad/trim to modal length
  seq_lens <- nchar(coding$SEQUENCE)
  modal_len <- as.integer(names(sort(table(seq_lens), decreasing = TRUE))[1])
  cat(sprintf("Modal probe length: %d\n", modal_len))

  # Keep only probes with modal length (vast majority)
  coding <- coding[seq_lens == modal_len, ]
  cat(sprintf("Probes with length %d: %d\n", modal_len, nrow(coding)))

  probe_dna <- DNAStringSet(coding$SEQUENCE)
  names(probe_dna) <- coding$ID

  # Build PDict for forward strand
  pdict_fwd <- PDict(probe_dna)
  # Build PDict for reverse complement
  probe_rc <- reverseComplement(probe_dna)
  pdict_rc <- PDict(probe_rc)

  # Get key genes
  convergent <- read.csv(file.path(PROJECT,
    "analysis/results/optimization_figures/cross_omics_convergent_genes.csv"))
  conv_genes <- convergent$gene[convergent$evidence_count >= 3]

  # Hallmark pathway genes
  gmt_file <- "/tmp/h.all.v2024.1.Hs.symbols.gmt"
  if (!file.exists(gmt_file)) {
    download.file(
      "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
      gmt_file, quiet = TRUE)
  }
  read_gmt <- function(f) {
    lines <- readLines(f)
    gs <- list()
    for (line in lines) { parts <- strsplit(line, "\t")[[1]]; gs[[parts[1]]] <- parts[-(1:2)] }
    gs
  }
  hallmark <- read_gmt(gmt_file)

  # Cell type signature genes
  cell_genes <- c("ALB","APOA1","APOB","SERPINA1","TF","TTR","HP","FGB","FGA","FGG",
    "AHSG","CYP2E1","CYP3A4","PCK1","KRT19","KRT7","EPCAM","SOX9","SPP1","CFTR",
    "CLEC4G","CLEC4M","STAB2","LYVE1","FCN2","FCN3","CD163","MARCO","CD5L","TIMD4",
    "NKG7","KLRD1","KLRB1","GNLY","GZMA","GZMB","PRF1","NCAM1","CD3D","CD3E","CD3G",
    "CD8A","CD8B","CD4","IL7R","TCF7","LEF1","FOXP3","CD79A","CD79B","MS4A1","CD19",
    "JCHAIN","IGHG1","MZB1","COL1A1","COL1A2","COL3A1","ACTA2","DCN","LUM","BGN",
    "TNF","IL1B","IL6","CCL2","CCL3","CXCL10","IDO1","NOS2","PTGS2","CD80","MRC1",
    "MSR1","TGFB1","IL10","CCL18","CLEC7A","S100A8","S100A9","S100A12","FCGR3B",
    "CXCR2","CSF3R","MMP9","PECAM1","CDH5","VWF","ENG","KDR","FLT1","EMCN","PLVAP")

  all_genes <- unique(c(conv_genes, unique(unlist(hallmark)), cell_genes))
  cat(sprintf("Total genes for mapping: %d\n", length(all_genes)))

  # Map genes to RefSeq accessions directly
  eg2refseq <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = all_genes, keytype = "SYMBOL",
           column = "REFSEQ", multiVals = "list"))
  # Flatten: gene -> list of NM_ accessions
  gene_refseq <- list()
  for (g in names(eg2refseq)) {
    accs <- eg2refseq[[g]]
    accs <- accs[!is.na(accs) & grepl("^NM_", accs)]
    if (length(accs) > 0) gene_refseq[[g]] <- accs[1]  # Take first NM_
  }
  refseq_vec <- unlist(gene_refseq)  # named: gene -> NM_accession
  cat(sprintf("Mapped %d genes to RefSeq NM_ accessions\n", length(refseq_vec)))

  # Reverse map: NM_ -> gene symbol
  refseq2gene <- setNames(names(refseq_vec), unname(refseq_vec))

  # Download mRNA sequences using entrez_fetch in small batches (by accession)
  probe_to_gene <- character(0)
  all_accs <- unique(unname(refseq_vec))

  # Process a FASTA text block: parse sequences, match probes via PDict
  process_fasta <- function(fasta_text) {
    entries <- strsplit(fasta_text, "\n>|\n(?=>)", perl = TRUE)[[1]]
    for (entry in entries) {
      lns <- strsplit(entry, "\n")[[1]]
      header <- lns[1]
      seq_str <- paste(gsub("[^ACGTacgt]", "", lns[-1]), collapse = "")
      if (nchar(seq_str) < modal_len) next

      # Extract gene symbol: try "(GENE)," pattern or use RefSeq->gene map
      gene_sym <- sub(".*\\(([A-Za-z0-9_-]+)\\),.*", "\\1", header)
      if (gene_sym == header || nchar(gene_sym) > 25) {
        # Try matching by accession
        acc <- regmatches(header, regexpr("NM_[0-9]+", header))
        if (length(acc) > 0 && acc %in% names(refseq2gene)) {
          gene_sym <- refseq2gene[acc]
        } else { next }
      }
      if (!(gene_sym %in% all_genes)) next

      target <- tryCatch(DNAString(toupper(seq_str)), error = function(e) NULL)
      if (is.null(target)) next

      # Forward PDict match (exact, all probes at once)
      hits <- tryCatch(matchPDict(pdict_fwd, target, max.mismatch = 0),
                       error = function(e) NULL)
      if (!is.null(hits)) {
        for (i in which(elementNROWS(hits) > 0)) {
          pid <- names(probe_dna)[i]
          if (!(pid %in% names(probe_to_gene))) probe_to_gene[pid] <<- gene_sym
        }
      }
      # Reverse complement match
      hits_rc <- tryCatch(matchPDict(pdict_rc, target, max.mismatch = 0),
                          error = function(e) NULL)
      if (!is.null(hits_rc)) {
        for (i in which(elementNROWS(hits_rc) > 0)) {
          pid <- names(probe_dna)[i]
          if (!(pid %in% names(probe_to_gene))) probe_to_gene[pid] <<- gene_sym
        }
      }
    }
  }

  batch_size <- 30  # Small batches to stay under URL limit
  for (b_start in seq(1, length(all_accs), by = batch_size)) {
    b_end <- min(b_start + batch_size - 1, length(all_accs))
    batch_accs <- all_accs[b_start:b_end]
    cat(sprintf("  RefSeq batch %d-%d / %d (mapped %d probes)...\n",
                b_start, b_end, length(all_accs), length(probe_to_gene)))

    tryCatch({
      # Search by accession
      query <- paste(batch_accs, "[ACCN]", collapse = " OR ")
      sr <- entrez_search(db = "nucleotide", term = query, retmax = length(batch_accs),
                          use_history = TRUE)
      if (sr$count > 0) {
        # Fetch using web history (avoids URL length issue)
        fasta <- entrez_fetch(db = "nucleotide", web_history = sr$web_history,
                              rettype = "fasta", retmax = length(batch_accs))
        process_fasta(fasta)
      }
      Sys.sleep(0.4)
    }, error = function(e) {
      cat(sprintf("    Error: %s\n", e$message))
      Sys.sleep(2)
    })
  }

  cat(sprintf("\n==> Mapped %d probes to %d unique genes\n",
              length(probe_to_gene), length(unique(probe_to_gene))))

  mapping_df <- data.frame(probe_id = names(probe_to_gene), gene = unname(probe_to_gene),
                           stringsAsFactors = FALSE)
  write.csv(mapping_df, MAPPING_FILE, row.names = FALSE)
}

# ============================================================
# PART C: Gene-level expression & differential analysis
# ============================================================
cat("\n--- Part C: Gene-level expression from GSE124362 ---\n")

mapped_probes <- mapping_df$probe_id[mapping_df$probe_id %in% rownames(expr_mat)]
cat(sprintf("Mapped probes in expression data: %d\n", length(mapped_probes)))

if (length(mapped_probes) < 20) {
  stop("Too few probes mapped. Cannot proceed with validation.")
}

expr_mapped <- expr_mat[mapped_probes, , drop = FALSE]
gene_labels <- mapping_df$gene[match(rownames(expr_mapped), mapping_df$probe_id)]

# Collapse to gene level
gene_expr <- do.call(rbind, lapply(unique(gene_labels), function(g) {
  idx <- which(gene_labels == g)
  if (length(idx) == 1) return(expr_mapped[idx, , drop = FALSE])
  colMeans(expr_mapped[idx, , drop = FALSE], na.rm = TRUE)
}))
rownames(gene_expr) <- unique(gene_labels)
cat(sprintf("Gene-level matrix: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

# Paired limma analysis
group <- factor(ifelse(seq_len(ncol(gene_expr)) %in% peri_idx, "Peri", "Dist"),
                levels = c("Dist", "Peri"))
patient <- factor(patient_ids)
design <- model.matrix(~ patient + group)
fit <- lmFit(gene_expr, design)
fit <- eBayes(fit)
deg_ext <- topTable(fit, coef = "groupPeri", number = Inf, sort.by = "P")
deg_ext$gene <- rownames(deg_ext)

cat(sprintf("DEGs: %d total, P<0.05: %d, adj.P<0.05: %d\n",
            nrow(deg_ext), sum(deg_ext$P.Value < 0.05), sum(deg_ext$adj.P.Val < 0.05)))
write.csv(deg_ext, file.path(OUT_DIR, "GSE124362_DEG_results.csv"), row.names = FALSE)

cat("Top 20 DEGs:\n")
print(head(deg_ext[, c("gene", "logFC", "P.Value", "adj.P.Val")], 20))

# ============================================================
# PART D: DEG direction concordance with our transcriptomics
# ============================================================
cat("\n--- Part D: DEG direction concordance ---\n")

# Find our DEG results
tc_files <- list.files(file.path(PROJECT, "analysis/results"),
                       pattern = "DEG|deseq2|differential",
                       recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
tc_files <- tc_files[grepl("TC|transcript", tc_files, ignore.case = TRUE)]
cat("TC DEG files found:\n"); print(tc_files)

our_deg <- NULL
for (f in tc_files) {
  tryCatch({
    d <- read.csv(f)
    gene_col <- intersect(c("gene", "Gene", "gene_name", "symbol"), colnames(d))[1]
    lfc_col <- grep("log2|logFC|lfc", colnames(d), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(gene_col) && !is.na(lfc_col) && nrow(d) > 100) {
      our_deg <- d; cat(sprintf("Using: %s\n", basename(f))); break
    }
  }, error = function(e) NULL)
}

if (!is.null(our_deg)) {
  gene_col <- intersect(c("gene", "Gene", "gene_name", "symbol"), colnames(our_deg))[1]
  lfc_col <- grep("log2|logFC|lfc", colnames(our_deg), ignore.case = TRUE, value = TRUE)[1]
  common <- intersect(our_deg[[gene_col]], deg_ext$gene)
  cat(sprintf("Common genes: %d\n", length(common)))

  if (length(common) >= 10) {
    concordance <- data.frame(
      gene = common,
      our_logFC = our_deg[[lfc_col]][match(common, our_deg[[gene_col]])],
      ext_logFC = deg_ext$logFC[match(common, deg_ext$gene)],
      stringsAsFactors = FALSE
    )
    concordance$same_dir <- sign(concordance$our_logFC) == sign(concordance$ext_logFC)
    concordance <- concordance[complete.cases(concordance), ]
    write.csv(concordance, file.path(OUT_DIR, "DEG_direction_concordance.csv"), row.names = FALSE)

    concord_rate <- mean(concordance$same_dir)
    cor_test <- cor.test(concordance$our_logFC, concordance$ext_logFC, method = "spearman")
    cat(sprintf("Concordance: %.1f%%, Spearman rho=%.3f, P=%.2e\n",
                concord_rate * 100, cor_test$estimate, cor_test$p.value))

    # Scatter plot
    p1 <- ggplot(concordance, aes(x = our_logFC, y = ext_logFC)) +
      geom_point(aes(color = same_dir), alpha = 0.6, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                         labels = c("Discordant", "Concordant"), name = "") +
      labs(x = "Our study log2FC (Adjacent vs Normal, n=12)",
           y = "GSE124362 log2FC (Periparasitic vs Distal, n=6)",
           title = sprintf("DEG Direction Concordance (n=%d genes)\nrho=%.3f, P=%.2e, concordance=%.0f%%",
                           nrow(concordance), cor_test$estimate, cor_test$p.value, concord_rate*100)) +
      theme_bw(base_size = 11)
    ggsave(file.path(OUT_DIR, "figures/DEG_concordance_scatter.pdf"), p1, width = 7, height = 6)
    cat("Saved scatter plot.\n")
  }
}

# ============================================================
# PART E: Pathway validation (ssGSEA)
# ============================================================
cat("\n--- Part E: Pathway-level validation (ssGSEA) ---\n")

# Load Hallmark sets
if (!exists("hallmark")) {
  gmt_file <- "/tmp/h.all.v2024.1.Hs.symbols.gmt"
  read_gmt <- function(f) {
    lines <- readLines(f)
    gs <- list()
    for (line in lines) { parts <- strsplit(line, "\t")[[1]]; gs[[parts[1]]] <- parts[-(1:2)] }
    gs
  }
  hallmark <- read_gmt(gmt_file)
}

# Filter to genes present in our matrix
hm_filt <- lapply(hallmark, function(gs) intersect(gs, rownames(gene_expr)))
hm_filt <- hm_filt[sapply(hm_filt, length) >= 5]
cat(sprintf("Hallmark pathways with >=5 genes: %d\n", length(hm_filt)))

# ssGSEA implementation
ssgsea <- function(expr, gsets) {
  scores <- matrix(NA, length(gsets), ncol(expr),
                   dimnames = list(names(gsets), colnames(expr)))
  for (i in seq_along(gsets)) {
    idx <- which(rownames(expr) %in% gsets[[i]])
    if (length(idx) < 3) next
    for (j in seq_len(ncol(expr))) {
      r <- rank(expr[, j])
      scores[i, j] <- mean(r[idx]) - mean(r[-idx])
    }
  }
  scores[complete.cases(scores), , drop = FALSE]
}

pw_scores <- ssgsea(gene_expr, hm_filt)
cat(sprintf("ssGSEA computed for %d pathways\n", nrow(pw_scores)))

# Paired test per pathway
pw_diff <- data.frame()
for (pw in rownames(pw_scores)) {
  tt <- tryCatch(t.test(pw_scores[pw, peri_idx], pw_scores[pw, dist_idx], paired = TRUE),
                 error = function(e) NULL)
  if (!is.null(tt)) {
    pw_diff <- rbind(pw_diff, data.frame(
      pathway = pw,
      mean_peri = mean(pw_scores[pw, peri_idx]),
      mean_dist = mean(pw_scores[pw, dist_idx]),
      diff = mean(pw_scores[pw, peri_idx]) - mean(pw_scores[pw, dist_idx]),
      pvalue = tt$p.value))
  }
}
pw_diff$padj <- p.adjust(pw_diff$pvalue, method = "BH")
pw_diff <- pw_diff[order(pw_diff$pvalue), ]
write.csv(pw_diff, file.path(OUT_DIR, "GSE124362_pathway_ssGSEA.csv"), row.names = FALSE)
cat(sprintf("Pathways: P<0.05=%d, padj<0.05=%d\n",
            sum(pw_diff$pvalue < 0.05), sum(pw_diff$padj < 0.05)))
cat("Top 10:\n"); print(head(pw_diff[, c("pathway","diff","pvalue","padj")], 10))

# Compare with our pathway results
our_pw <- tryCatch(
  read.csv(file.path(PROJECT, "analysis/results/phase2_enrichment/cross_omics_convergent_Hallmark.csv")),
  error = function(e) NULL)

if (!is.null(our_pw)) {
  # Attempt to match pathway names
  pw_diff$pw_short <- gsub("^HALLMARK_", "", pw_diff$pathway)
  if ("pathway" %in% colnames(our_pw)) {
    our_pw$pw_short <- gsub("^HALLMARK_", "", our_pw$pathway)
  }
  common_pw <- intersect(pw_diff$pw_short, our_pw$pw_short)
  cat(sprintf("Common pathways for comparison: %d\n", length(common_pw)))

  if (length(common_pw) >= 5) {
    # Determine direction columns in our data
    our_dir_col <- grep("direction|NES|diff|logFC|TC_dir", colnames(our_pw), ignore.case=TRUE, value=TRUE)
    cat("Our pathway direction columns:", paste(our_dir_col, collapse=", "), "\n")

    pw_compare <- data.frame(
      pathway = common_pw,
      ext_diff = pw_diff$diff[match(common_pw, pw_diff$pw_short)],
      ext_pval = pw_diff$pvalue[match(common_pw, pw_diff$pw_short)]
    )
    if (length(our_dir_col) > 0) {
      pw_compare$our_value <- our_pw[[our_dir_col[1]]][match(common_pw, our_pw$pw_short)]
    }
    write.csv(pw_compare, file.path(OUT_DIR, "pathway_concordance.csv"), row.names = FALSE)

    # Bar plot of pathway scores
    plot_pw <- pw_diff[pw_diff$pw_short %in% common_pw, ]
    plot_pw$pw_label <- gsub("_", " ", plot_pw$pw_short)
    plot_pw$sig <- ifelse(plot_pw$pvalue < 0.05, "*", "")

    p2 <- ggplot(plot_pw, aes(x = reorder(pw_label, diff), y = diff)) +
      geom_bar(stat = "identity", aes(fill = diff > 0), width = 0.7) +
      geom_text(aes(label = sig), hjust = ifelse(plot_pw$diff > 0, -0.3, 1.3), size = 5) +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
      labs(x = "", y = "ssGSEA score difference (Periparasitic - Distal)",
           title = "Hallmark Pathway Activity in GSE124362 (n=6 paired)") +
      theme_bw(base_size = 10)
    ggsave(file.path(OUT_DIR, "figures/pathway_validation_barplot.pdf"), p2,
           width = 10, height = max(5, length(common_pw) * 0.3))
    cat("Saved pathway bar plot.\n")
  }
}

# ============================================================
# PART F: Cell type deconvolution validation
# ============================================================
cat("\n--- Part F: Cell type deconvolution validation ---\n")

cell_sigs <- list(
  Hepatocytes = c("ALB","APOA1","APOB","SERPINA1","TF","TTR","HP","FGB","FGA","FGG","AHSG","CYP2E1","CYP3A4","PCK1"),
  Cholangiocytes = c("KRT19","KRT7","EPCAM","SOX9","SPP1","CFTR","ANXA4"),
  LSECs = c("CLEC4G","CLEC4M","STAB2","LYVE1","FCN2","FCN3","CLEC1B"),
  Kupffer_cells = c("CD163","MARCO","CD5L","TIMD4","VSIG4","SIGLEC1"),
  NK_NKT = c("NKG7","KLRD1","KLRB1","GNLY","GZMA","GZMB","PRF1","NCAM1","KLRK1","CD160"),
  T_cells = c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","IL7R","TCF7","LEF1","FOXP3"),
  B_Plasma = c("CD79A","CD79B","MS4A1","CD19","JCHAIN","IGHG1","MZB1","IGKC","IGLC2"),
  Stellate_Fibro = c("COL1A1","COL1A2","COL3A1","ACTA2","DCN","LUM","BGN","TAGLN"),
  Macrophages_M1 = c("TNF","IL1B","IL6","CCL2","CCL3","CXCL10","IDO1","NOS2","PTGS2","CD80"),
  Macrophages_M2 = c("CD163","MRC1","MSR1","TGFB1","IL10","CCL18","CLEC7A","CD209","STAB1"),
  Neutrophils = c("S100A8","S100A9","S100A12","FCGR3B","CXCR2","CSF3R","MMP9","CEACAM8"),
  Endothelial = c("PECAM1","CDH5","VWF","ENG","KDR","FLT1","EMCN","PLVAP")
)

ct_filt <- lapply(cell_sigs, function(gs) intersect(gs, rownames(gene_expr)))
ct_filt <- ct_filt[sapply(ct_filt, length) >= 3]
cat(sprintf("Cell types with >=3 marker genes: %d/%d\n", length(ct_filt), length(cell_sigs)))

if (length(ct_filt) >= 3) {
  ct_scores <- ssgsea(gene_expr, ct_filt)

  ct_diff <- data.frame()
  for (ct in rownames(ct_scores)) {
    tt <- tryCatch(t.test(ct_scores[ct, peri_idx], ct_scores[ct, dist_idx], paired = TRUE),
                   error = function(e) NULL)
    if (!is.null(tt)) {
      ct_diff <- rbind(ct_diff, data.frame(
        cell_type = ct,
        mean_peri = mean(ct_scores[ct, peri_idx]),
        mean_dist = mean(ct_scores[ct, dist_idx]),
        diff = mean(ct_scores[ct, peri_idx]) - mean(ct_scores[ct, dist_idx]),
        pvalue = tt$p.value))
    }
  }
  ct_diff$padj <- p.adjust(ct_diff$pvalue, method = "BH")
  ct_diff <- ct_diff[order(ct_diff$pvalue), ]
  write.csv(ct_diff, file.path(OUT_DIR, "GSE124362_celltype_deconv.csv"), row.names = FALSE)
  cat("External cell type results:\n"); print(ct_diff[, c("cell_type","diff","pvalue","padj")])

  # Compare with our Enhancement 8 results
  our_ct <- tryCatch(
    read.csv(file.path(PROJECT, "analysis/results/enhancement8_validation/cell_type_differential.csv")),
    error = function(e) NULL)

  if (!is.null(our_ct)) {
    common_ct <- intersect(our_ct$cell_type, ct_diff$cell_type)
    if (length(common_ct) >= 3) {
      ct_compare <- data.frame(
        cell_type = common_ct,
        our_diff = our_ct$diff[match(common_ct, our_ct$cell_type)],
        ext_diff = ct_diff$diff[match(common_ct, ct_diff$cell_type)],
        our_pval = our_ct$pvalue[match(common_ct, our_ct$cell_type)],
        ext_pval = ct_diff$pvalue[match(common_ct, ct_diff$cell_type)])
      ct_compare$concordant <- sign(ct_compare$our_diff) == sign(ct_compare$ext_diff)
      write.csv(ct_compare, file.path(OUT_DIR, "celltype_concordance.csv"), row.names = FALSE)
      cat(sprintf("Cell type concordance: %d/%d (%.0f%%)\n",
                  sum(ct_compare$concordant), nrow(ct_compare), 100*mean(ct_compare$concordant)))

      # Side-by-side bar plot
      pdata <- rbind(
        data.frame(cell_type = ct_compare$cell_type, diff = ct_compare$our_diff,
                   dataset = "Our study (n=12)"),
        data.frame(cell_type = ct_compare$cell_type, diff = ct_compare$ext_diff,
                   dataset = "GSE124362 (n=6)"))
      p3 <- ggplot(pdata, aes(x = reorder(cell_type, diff), y = diff, fill = dataset)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = c("Our study (n=12)" = "#2166AC", "GSE124362 (n=6)" = "#B2182B")) +
        labs(x = "", y = "Score difference (Disease - Control)",
             title = "Cell Type Deconvolution Validation", fill = "") +
        theme_bw(base_size = 11) + theme(legend.position = "bottom")
      ggsave(file.path(OUT_DIR, "figures/celltype_validation_barplot.pdf"), p3, width = 9, height = 6)
      cat("Saved cell type validation plot.\n")
    }
  }
}

# ============================================================
# PART G: Convergent gene validation
# ============================================================
cat("\n--- Part G: Convergent gene validation ---\n")

conv_in_ext <- intersect(conv_genes, rownames(gene_expr))
cat(sprintf("Convergent genes in GSE124362: %d/%d\n", length(conv_in_ext), length(conv_genes)))

if (length(conv_in_ext) >= 3) {
  conv_val <- data.frame()
  for (g in conv_in_ext) {
    tt <- tryCatch(t.test(gene_expr[g, peri_idx], gene_expr[g, dist_idx], paired = TRUE),
                   error = function(e) NULL)
    if (!is.null(tt)) {
      conv_val <- rbind(conv_val, data.frame(
        gene = g, logFC = mean(gene_expr[g, peri_idx]) - mean(gene_expr[g, dist_idx]),
        pvalue = tt$p.value))
    }
  }
  conv_val$padj <- p.adjust(conv_val$pvalue, method = "BH")
  conv_val <- conv_val[order(conv_val$pvalue), ]
  write.csv(conv_val, file.path(OUT_DIR, "convergent_gene_validation.csv"), row.names = FALSE)
  cat("Convergent gene validation:\n"); print(conv_val)

  # Heatmap
  if (length(conv_in_ext) >= 3) {
    ann_col <- data.frame(
      Group = factor(ifelse(seq_len(ncol(gene_expr)) %in% peri_idx, "Periparasitic", "Distal")),
      row.names = colnames(gene_expr))
    ann_colors <- list(Group = c(Periparasitic = "#B2182B", Distal = "#2166AC"))
    mat <- gene_expr[conv_in_ext, , drop = FALSE]
    pdf(file.path(OUT_DIR, "figures/convergent_gene_heatmap_GSE124362.pdf"),
        width = 8, height = max(4, length(conv_in_ext) * 0.4))
    print(pheatmap(mat, scale = "row", annotation_col = ann_col,
                   annotation_colors = ann_colors, cluster_cols = FALSE,
                   fontsize = 10, main = "Convergent Genes in GSE124362"))
    dev.off()
    cat("Saved convergent gene heatmap.\n")
  }
}

# ============================================================
# Summary
# ============================================================
cat("\n========== Enhancement 9 COMPLETE ==========\n")
out_files <- list.files(OUT_DIR, recursive = TRUE)
cat(sprintf("Generated %d files in %s:\n", length(out_files), OUT_DIR))
for (f in out_files) cat(sprintf("  %s\n", f))
