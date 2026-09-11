#!/usr/bin/env Rscript
# enhancement39: 跨方法免疫去卷积共识 (EPIC + MCP-counter) ========================
# 用两种基于标记基因的独立方法 (与 BayesPrism/ssGSEA 原理不同) 对
# 24 个 bulk 样本做免疫细胞定量, 评估 Adjacent vs Normal 位移方向,
# 与 BayesPrism 结果形成三角共识.
# ==============================================================================
suppressMessages({
  library(matrixStats)
})
set.seed(42)

ROOT <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
OUT <- file.path(ROOT, "02_analysis/results/enhancement39_cross_deconv")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 1. bulk 表达矩阵 (线性 scale: 反 log2 CPM -> CPM) -----------------------
bulk <- read.csv(file.path(ROOT, "02_analysis/data/processed/transcriptomics_logcpm_paired.csv"),
                 row.names = 1, check.names = FALSE)
deg <- read.csv(file.path(ROOT, "02_analysis/results/phase1_diff/DEGs_Adjacent_vs_Normal.csv"),
                stringsAsFactors = FALSE)
gmap <- setNames(deg$gene_name, deg$gene_id)
rn <- rownames(bulk)
sym <- ifelse(!is.na(gmap[rn]) & gmap[rn] != "", gmap[rn], sub("_\\d+$", "", rn))
expr <- rowsum(bulk, group = sym)  # aggregate duplicated symbols
cpm <- 2^expr - 1  # 反 log2
samples <- colnames(cpm)
pats <- sort(unique(gsub("^(Normal|Adjacent)", "", samples)))
cat(sprintf("bulk: %d genes x %d samples (%d patients)\n", nrow(cpm), ncol(cpm), length(pats)))

paired_test <- function(mat, label) {
  # mat: samples x cell types (任意量纲), 输出配对差值 Wilcoxon (符号检验)
  # pats 顺序在 Adjacent/Normal 两个子集间一致 -> 按位置配对相减
  out <- list()
  for (ct in colnames(mat)) {
    a <- mat[paste0("Adjacent", pats), ct]
    n <- mat[paste0("Normal", pats), ct]
    if (length(a) != length(n) || length(a) < 3) next
    d <- a - n
    p <- tryCatch(wilcox.test(d)$p.value, error = function(e) NA)
    out[[ct]] <- data.frame(method = label, cell = ct, n = length(a),
                            delta = mean(d), wilcoxon_p = p, row.names = NULL)
  }
  do.call(rbind, out)
}

## ---- 2. MCP-counter (10 谱系, 绝对丰度分数) ----------------------------------
EXT <- file.path(ROOT, "02_analysis/data/external")
mcp_ok <- file.exists(file.path(EXT, "MCPcounter", "MCPcounter.R"))
if (mcp_ok) {
  source(file.path(EXT, "MCPcounter", "MCPcounter.R"))
  mcp_genes <- read.table(file.path(EXT, "MCPcounter", "genes.txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                          colClasses = "character", check.names = FALSE)
  mcp_probes <- read.table(file.path(EXT, "MCPcounter", "probesets.txt"),
                           sep = "\t", stringsAsFactors = FALSE,
                           colClasses = "character")
  mcp <- MCPcounter.estimate(as.matrix(cpm), featuresType = "HUGO_symbols",
                             probesets = mcp_probes, genes = mcp_genes)
  mcp <- t(mcp)
  write.csv(mcp, file.path(OUT, "mcpcounter_scores.csv"))
  mcp_res <- paired_test(mcp, "MCP-counter")
  write.csv(mcp_res, file.path(OUT, "mcpcounter_paired.csv"), row.names = FALSE)
  cat("\n-- MCP-counter paired Adjacent vs Normal --\n"); print(mcp_res)
} else cat("MCPcounter source not found\n")

## ---- 3. EPIC (TIL 参考, 输出细胞比例) ----------------------------------------
epic_ok <- file.exists(file.path(EXT, "EPIC", "EPIC_fun.R"))
if (epic_ok) {
  source(file.path(EXT, "EPIC", "EPIC_fun.R"))
  epic_env <- new.env()
  load(file.path(EXT, "EPIC", "TRef.rda"), envir = epic_env)
  load(file.path(EXT, "EPIC", "mRNA_cell_default.rda"), envir = epic_env)
  e <- EPIC(bulk = as.matrix(cpm),
            reference = epic_env$TRef,
            mRNA_cell = epic_env$mRNA_cell_default)
  ep <- e$cellFractions
  rownames(ep) <- colnames(cpm)
  write.csv(ep, file.path(OUT, "epic_fractions.csv"))
  ep_res <- paired_test(ep, "EPIC")
  write.csv(ep_res, file.path(OUT, "epic_paired.csv"), row.names = FALSE)
  cat("\n-- EPIC paired Adjacent vs Normal --\n"); print(ep_res)
} else cat("EPIC source not found\n")

## ---- 4. 与 BayesPrism 三角共识 ----------------------------------------------
bp <- read.csv(file.path(ROOT, "02_analysis/results/enhancement16_deconvolution/bayesprism_cell_proportions.csv"),
               row.names = 1)
consensus <- function(m2, map2, label2) {
  rows <- list()
  for (ct in names(map2)) {
    if (!ct %in% colnames(m2)) next
    a <- m2[paste0("Adjacent", pats), ct]
    n <- m2[paste0("Normal", pats), ct]
    d2 <- mean(a - n)
    bpc <- map2[[ct]]
    if (!bpc %in% colnames(bp)) next
    bpa <- bp[paste0("Adjacent", pats), bpc]
    bpn <- bp[paste0("Normal", pats), bpc]
    d1 <- mean(bpa - bpn)
    rows[[ct]] <- data.frame(cell = ct, bayesprism_delta = d1, other_delta = d2, method = label2)
  }
  do.call(rbind, rows)
}

if (mcp_ok) {
  cm <- consensus(mcp, c(`T cells` = "T_CD8", `CD8 T cells` = "T_CD8",
                          `B lineage` = "B_cell",
                          `Monocytic lineage` = "Monocyte_derived_Mac",
                          `NK cells` = "NK_cell", Neutrophils = "Neutrophil"), "MCP-counter")
  write.csv(cm, file.path(OUT, "consensus_mcp_vs_bayesprism.csv"), row.names = FALSE)
  cat("\n-- MCP vs BayesPrism direction consensus --\n"); print(cm)
}
if (epic_ok) {
  # EPIC TIL reference cell names detected at runtime
  epic_cells <- colnames(ep)
  cat("\nEPIC cell types:", paste(epic_cells, collapse = ", "), "\n")
  map_epic <- c(Bcells = "B_cell", CD4_Tcells = "T_CD4",
                CD8_Tcells = "T_CD8", NKcells = "NK_cell",
                Macrophages = "Monocyte_derived_Mac",
                CAFs = "Stellate_cell", Endothelial = "Endothelial_LSEC")
  map_epic <- map_epic[names(map_epic) %in% epic_cells]
  ce <- consensus(ep, map_epic, "EPIC")
  write.csv(ce, file.path(OUT, "consensus_epic_vs_bayesprism.csv"), row.names = FALSE)
  cat("\n-- EPIC vs BayesPrism direction consensus --\n"); print(ce)
}

## ---- 5. summary ---------------------------------------------------------------
sm <- c("# enhancement39: cross-method deconvolution consensus (MCP-counter + EPIC vs BayesPrism)",
        "", sprintf("Date: %s", format(Sys.Date())), "")
if (mcp_ok) sm <- c(sm, "## MCP-counter (marker geometric mean, absolute scores)",
  "See mcpcounter_paired.csv - paired signed-rank test Adjacent vs Normal.")
if (epic_ok) sm <- c(sm, "## EPIC (TIL reference, cell fractions)",
  "See epic_paired.csv - paired signed-rank test Adjacent vs Normal.")
sm <- c(sm, "## Triangulation",
  "Direction of peri-lesional immune shifts compared between MCP-counter, EPIC and the in-house BayesPrism estimates; a cell population is called consensus-shifted when all available methods agree in sign.",
  "")
writeLines(sm, file.path(OUT, "summary.md"))
