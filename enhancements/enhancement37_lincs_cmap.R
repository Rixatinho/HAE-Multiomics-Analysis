#!/usr/bin/env Rscript
# enhancement37: LINCS L1000 / CMap 真实数据库连通性查询 ==========================
# 用配对疾病签名 (Adjacent vs Normal, 名义 P<0.05 DEG) 查询:
#   1) CMap build-2 (refdb="cmap", ~7,000 signatures, 1,309 compounds)
#   2) LINCS L1000 (refdb="lincs", >1M signatures, 真实药物扰动谱)
# 输出反向签名的候选药物 (负 ES = 逆转疾病表型), FDR 校正,
# 并核查 pirfenidone / nintedanib / albendazole 是否在库及得分.
# ==============================================================================
suppressMessages({
  library(signatureSearch)
  library(SummarizedExperiment)
})
set.seed(42)

ROOT <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
OUT <- file.path(ROOT, "02_analysis/results/enhancement37_lincs_cmap")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 1. 疾病签名 (名义 P<0.05) ----------------------------------------------
deg <- read.csv(file.path(ROOT, "02_analysis/results/phase1_diff/DEGs_Adjacent_vs_Normal.csv"),
                stringsAsFactors = FALSE)
deg$gene_name <- trimws(deg$gene_name)
deg <- deg[deg$gene_name != "" & !grepl("^ENSG", deg$gene_name) & !is.na(deg$P.Value), ]
deg <- deg[!duplicated(deg$gene_name), ]
up <- deg$gene_name[deg$significance == "Up"]
dn <- deg$gene_name[deg$significance == "Down"]
cat(sprintf("signature: %d up / %d down (nominal P<0.05)\n", length(up), length(dn)))

run_db <- function(refdb, tag) {
  cat(sprintf("\n===== querying %s =====\n", tag))
  qsig <- try(qSig(q = list(up = up, down = dn), gess_method = "LINCS", refdb = refdb))
  if (inherits(qsig, "try-error")) { cat("qSig failed\n"); return(invisible(NULL)) }
  gess <- try(gess_lincs(qsig, tau = FALSE, sortby = "NCS", workers = 4))
  if (inherits(gess, "try-error")) { cat("gess_lincs failed\n"); return(invisible(NULL)) }
  res <- as.data.frame(result(gess))
  cat(sprintf("total perturbation instances: %d\n", nrow(res)))
  res$padj <- p.adjust(res$pval, method = "BH")

  ## 按化合物聚合 (跨细胞系)
  agg <- aggregate(cbind(WTCS, NCS) ~ pert + cell_type, data = res, FUN = mean)
  cmp <- aggregate(NCS ~ pert, data = res, FUN = function(x) mean(x, na.rm = TRUE))
  cmp$n_instances <- as.integer(table(res$pert)[cmp$pert])
  cmp$min_padj <- tapply(res$padj, res$pert, min)[cmp$pert]
  cmp <- cmp[order(cmp$NCS), ]
  write.csv(res, file.path(OUT, sprintf("%s_all_instances.csv", tag)), row.names = FALSE)
  write.csv(cmp, file.path(OUT, sprintf("%s_compound_level.csv", tag)), row.names = FALSE)

  ## 关键药物核查
  keydrugs <- c("pirfenidone", "nintedanib", "albendazole", "imatinib", "ponatinib",
                "rifampicin", "dexamethasone", "sirolimus", "everolimus", "prednisolone")
  hit <- cmp[tolower(cmp$pert) %in% keydrugs, ]
  cat("\n-- key drug check --\n"); print(hit)

  ## 肝系特异 (HEPG2/HEPG2-lids)
  hep <- res[grepl("hepg2|hep_g2|hepatocyte", res$cell_type, ignore.case = TRUE), ]
  if (nrow(hep)) {
    cmp_h <- aggregate(NCS ~ pert, data = hep, FUN = mean)
    cmp_h <- cmp_h[order(cmp_h$NCS), ]
    write.csv(cmp_h, file.path(OUT, sprintf("%s_hepg2_compound_level.csv", tag)), row.names = FALSE)
    cat(sprintf("\nHEPG2 instances: %d; top reversers:\n", nrow(hep)))
    print(head(cmp_h, 15))
  }
  cat(sprintf("\ntop 20 global reversers (%s):\n", tag)); print(head(cmp, 20))
  invisible(res)
}

## ---- 2. CMap build2 (小库先行) ----------------------------------------------
r_cmap <- run_db("cmap", "cmap")

## ---- 3. LINCS L1000 (大库, 自动经 ExperimentHub 下载缓存) ---------------------
r_lincs <- run_db("lincs", "lincs")

## ---- 4. 汇总 -----------------------------------------------------------------
cat("\n[done] outputs in ", OUT, "\n", sep = "")
