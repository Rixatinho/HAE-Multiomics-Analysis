#!/usr/bin/env Rscript
# ============================================================================
# MOFA因子生物学注释分析
# HAE多组学项目 - Task #4
# ============================================================================

# 设置随机种子
set.seed(42)

# 设置工作目录
proj_root <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
setwd(proj_root)

# 创建输出目录
out_dir <- "analysis/results/mofa_annotation"
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("MOFA因子生物学注释分析\n")
cat("============================================================\n\n")

# ============================================================================
# 加载所需的R包
# ============================================================================
cat("[Step 0] 加载R包...\n")

suppressPackageStartupMessages({
  library(MOFA2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# NC主题 - 使用默认sans字体
nc_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 9)
  )

# ============================================================================
# Step 1: 加载MOFA模型和数据
# ============================================================================
cat("\n[Step 1] 加载MOFA模型和数据...\n")

# 加载MOFA模型
mofa_model <- readRDS("analysis/results/phase5_integration/mofa2_model.rds")
cat("  - MOFA模型加载完成\n")
cat("  - 因子数量:", get_dimensions(mofa_model)$K, "\n")
cat("  - 样本数量:", get_dimensions(mofa_model)$N, "\n")
cat("  - 视图:", paste(views_names(mofa_model), collapse = ", "), "\n")

# 提取因子值（样本×因子矩阵）
factor_values <- get_factors(mofa_model)[[1]]
cat("  - 因子值维度:", nrow(factor_values), "x", ncol(factor_values), "\n")

# 提取所有权重
weights_list <- get_weights(mofa_model, as.data.frame = TRUE)
cat("  - 权重数据行数:", nrow(weights_list), "\n")

# 加载临床数据
clinical <- read.csv("analysis/results/phase6_subtyping/clinical_parsed_corrected.csv", 
                     stringsAsFactors = FALSE)
cat("  - 临床数据: ", nrow(clinical), " 样本, ", ncol(clinical), " 变量\n")

# 加载代谢物注释
metab_annot <- read.csv("analysis/data/processed/metabolomics_annotation.csv",
                        stringsAsFactors = FALSE)
cat("  - 代谢物注释: ", nrow(metab_annot), " 条目\n")

# 加载Hallmark基因集（从本地GMT文件）
hallmark_gmt_path <- "analysis/data/gmt/h.all.v2024.1.Hs.symbols.gmt"
if (file.exists(hallmark_gmt_path)) {
  # 读取GMT文件
  gmt_lines <- readLines(hallmark_gmt_path)
  hallmark_list <- list()
  hallmark_t2g <- data.frame(gs_name = character(), gene_symbol = character(), stringsAsFactors = FALSE)
  
  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      pathway_name <- parts[1]
      genes <- parts[3:length(parts)]
      hallmark_list[[pathway_name]] <- genes
      hallmark_t2g <- rbind(hallmark_t2g, data.frame(
        gs_name = pathway_name,
        gene_symbol = genes,
        stringsAsFactors = FALSE
      ))
    }
  }
  cat("  - Hallmark基因集(本地): ", length(hallmark_list), " 个通路\n")
} else {
  # 尝试msigdbr
  tryCatch({
    hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
    hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
    hallmark_t2g <- hallmark_df %>% dplyr::select(gs_name, gene_symbol)
    cat("  - Hallmark基因集(msigdbr): ", length(hallmark_list), " 个通路\n")
  }, error = function(e) {
    hallmark_list <- list()
    hallmark_t2g <- data.frame(gs_name = character(), gene_symbol = character())
    cat("  - 警告: 无法加载Hallmark基因集\n")
  })
}

# ============================================================================
# Step 2: 因子权重基因的通路富集分析
# ============================================================================
cat("\n[Step 2] 因子权重基因通路富集分析...\n")

# 获取所有视图的权重
all_weights <- weights_list %>%
  mutate(
    # 清理feature名称，移除后缀
    gene_symbol = gsub("_transcriptomics|_proteomics", "", feature),
    # 标记代谢物
    is_metabolite = grepl("^Com_", feature)
  )

# 获取因子名称
factor_names <- unique(all_weights$factor)
view_names <- unique(all_weights$view)

cat("  - 因子:", paste(factor_names, collapse = ", "), "\n")
cat("  - 视图:", paste(view_names, collapse = ", "), "\n")

# Gene symbol到ENTREZ ID映射函数
symbol_to_entrez <- function(symbols) {
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

# 富集分析函数
run_enrichment <- function(genes, universe_genes = NULL, direction = "positive") {
  if (length(genes) < 5) {
    return(list(GO = NULL, KEGG = NULL, Hallmark = NULL))
  }
  
  # 转换为ENTREZ ID
  entrez_ids <- symbol_to_entrez(genes)
  
  if (length(entrez_ids) < 3) {
    return(list(GO = NULL, KEGG = NULL, Hallmark = NULL))
  }
  
  # GO BP富集
  go_result <- tryCatch({
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }, error = function(e) NULL)
  
  # KEGG富集
  kegg_result <- tryCatch({
    enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) NULL)
  
  # Hallmark富集
  hallmark_result <- tryCatch({
    enricher(
      gene = genes,
      TERM2GENE = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) NULL)
  
  return(list(GO = go_result, KEGG = kegg_result, Hallmark = hallmark_result))
}

# 对每个因子进行富集分析
enrichment_results <- list()
all_enrichment_df <- data.frame()

for (fac in factor_names) {
  cat("  - 分析", fac, "...\n")
  
  # 获取转录组和蛋白组权重（排除代谢物）
  fac_weights <- all_weights %>%
    filter(factor == fac, !is_metabolite, view %in% c("transcriptomics", "proteomics")) %>%
    arrange(desc(abs(value)))
  
  # 按绝对权重排序，取Top 100
  top_features <- fac_weights %>%
    head(100)
  
  # 分正负权重
  pos_genes <- top_features %>%
    filter(value > 0) %>%
    head(50) %>%
    pull(gene_symbol) %>%
    unique()
  
  neg_genes <- top_features %>%
    filter(value < 0) %>%
    head(50) %>%
    pull(gene_symbol) %>%
    unique()
  
  cat("    - 正权重基因:", length(pos_genes), "\n")
  cat("    - 负权重基因:", length(neg_genes), "\n")
  
  # 正权重富集
  pos_enrich <- run_enrichment(pos_genes, direction = "positive")
  neg_enrich <- run_enrichment(neg_genes, direction = "negative")
  
  enrichment_results[[fac]] <- list(
    positive = pos_enrich,
    negative = neg_enrich,
    pos_genes = pos_genes,
    neg_genes = neg_genes
  )
  
  # 合并富集结果到数据框
  for (dir in c("positive", "negative")) {
    enrich <- if (dir == "positive") pos_enrich else neg_enrich
    
    for (type in c("GO", "KEGG", "Hallmark")) {
      if (!is.null(enrich[[type]]) && nrow(as.data.frame(enrich[[type]])) > 0) {
        df <- as.data.frame(enrich[[type]]) %>%
          mutate(
            factor = fac,
            direction = dir,
            enrichment_type = type
          ) %>%
          head(10)  # 每类取Top 10
        all_enrichment_df <- bind_rows(all_enrichment_df, df)
      }
    }
  }
}

# 保存富集结果
if (nrow(all_enrichment_df) > 0) {
  write.csv(all_enrichment_df, file.path(out_dir, "factor_pathway_enrichment.csv"),
            row.names = FALSE)
  cat("  - 富集结果保存: factor_pathway_enrichment.csv (", nrow(all_enrichment_df), "行)\n")
} else {
  cat("  - 警告: 没有显著的富集结果\n")
}

# ============================================================================
# Step 3: 因子与临床特征的关联分析
# ============================================================================
cat("\n[Step 3] 因子与临床特征关联分析...\n")

# 确保样本ID匹配
# MOFA样本名格式: Normal1, Normal2, ... Adjacent1, Adjacent2, ...
sample_ids <- rownames(factor_values)

# 提取样本ID中的数字部分
# 取唯一的patient ID（Normal和Adjacent同一个患者）
patient_nums <- unique(as.numeric(gsub("Normal|Adjacent", "", sample_ids)))
cat("  - MOFA患者ID: ", paste(patient_nums, collapse = ", "), "\n")

clinical_matched <- clinical %>%
  filter(patient_id %in% patient_nums)

# 如果仍然没有匹配，用complete_3omics过滤
if (nrow(clinical_matched) == 0) {
  clinical_matched <- clinical %>%
    filter(complete_3omics == TRUE)
}

# 准备因子值与临床数据的匹配
# 对于每个患者，取Adjacent样本的因子值（肿瘤组织）
adjacent_samples <- grep("Adjacent", sample_ids, value = TRUE)
factor_for_clinical <- factor_values[adjacent_samples, , drop = FALSE]
# 提取patient_id并排序
factor_patient_ids <- as.numeric(gsub("Adjacent", "", rownames(factor_for_clinical)))
factor_for_clinical <- factor_for_clinical[order(factor_patient_ids), , drop = FALSE]
clinical_matched <- clinical_matched[order(clinical_matched$patient_id), ]

cat("  - 匹配样本数:", nrow(clinical_matched), "\n")

# 选择数值型临床变量
numeric_vars <- clinical_matched %>%
  select(where(is.numeric)) %>%
  select(-patient_id) %>%
  select(-matches("has_|complete_")) %>%
  colnames()

cat("  - 数值型临床变量:", length(numeric_vars), "\n")

# 计算相关性和p值
cor_matrix <- matrix(NA, nrow = ncol(factor_values), ncol = length(numeric_vars),
                     dimnames = list(colnames(factor_values), numeric_vars))
pval_matrix <- cor_matrix

for (i in 1:ncol(factor_for_clinical)) {
  fac_vals <- factor_for_clinical[, i]
  fac_name <- colnames(factor_for_clinical)[i]
  
  for (var in numeric_vars) {
    clin_vals <- clinical_matched[[var]]
    
    # 移除NA
    valid_idx <- !is.na(clin_vals)
    if (sum(valid_idx) >= 5) {
      cor_test <- cor.test(fac_vals[valid_idx], clin_vals[valid_idx], 
                           method = "spearman", exact = FALSE)
      cor_matrix[fac_name, var] <- cor_test$estimate
      pval_matrix[fac_name, var] <- cor_test$p.value
    }
  }
}

# BH校正
pval_adj <- matrix(p.adjust(as.vector(pval_matrix), method = "BH"),
                   nrow = nrow(pval_matrix), ncol = ncol(pval_matrix),
                   dimnames = dimnames(pval_matrix))

# 保存关联结果
write.csv(cor_matrix, file.path(out_dir, "factor_clinical_correlation.csv"))
write.csv(pval_adj, file.path(out_dir, "factor_clinical_pvalue_adj.csv"))
cat("  - 关联结果保存完成\n")

# ============================================================================
# Step 4: 因子跨组学方差贡献
# ============================================================================
cat("\n[Step 4] 因子跨组学方差贡献...\n")

# 获取方差解释
var_explained <- get_variance_explained(mofa_model)

# 整理为数据框 - r2_per_factor[[1]] 是一个矩阵 [factors x views]
var_mat <- var_explained$r2_per_factor[[1]]
var_df <- data.frame()

for (i in 1:nrow(var_mat)) {
  for (j in 1:ncol(var_mat)) {
    var_df <- bind_rows(var_df, data.frame(
      Factor = rownames(var_mat)[i],
      View = colnames(var_mat)[j],
      Variance = var_mat[i, j]
    ))
  }
}

# 保存方差解释
write.csv(var_df, file.path(out_dir, "factor_variance_explained.csv"), row.names = FALSE)
cat("  - 方差解释保存完成\n")
cat("  - 总方差解释:\n")
print(var_explained$r2_total[[1]])

# ============================================================================
# Step 5: 因子生物学注释总结
# ============================================================================
cat("\n[Step 5] 因子生物学注释总结...\n")

# 基于富集结果为每个因子生成生物学命名
annotation_summary <- data.frame(
  Factor = factor_names,
  stringsAsFactors = FALSE
)

annotation_summary$Top_GO_Positive <- ""
annotation_summary$Top_GO_Negative <- ""
annotation_summary$Top_Hallmark <- ""
annotation_summary$Biological_Theme <- ""
annotation_summary$Key_Genes_Positive <- ""
annotation_summary$Key_Genes_Negative <- ""

for (i in seq_along(factor_names)) {
  fac <- factor_names[i]
  
  # 获取Top富集通路
  pos_go <- all_enrichment_df %>%
    filter(factor == fac, direction == "positive", enrichment_type == "GO") %>%
    head(3)
  neg_go <- all_enrichment_df %>%
    filter(factor == fac, direction == "negative", enrichment_type == "GO") %>%
    head(3)
  hallmark <- all_enrichment_df %>%
    filter(factor == fac, enrichment_type == "Hallmark") %>%
    head(3)
  
  annotation_summary$Top_GO_Positive[i] <- paste(pos_go$Description, collapse = "; ")
  annotation_summary$Top_GO_Negative[i] <- paste(neg_go$Description, collapse = "; ")
  annotation_summary$Top_Hallmark[i] <- paste(gsub("HALLMARK_", "", hallmark$ID), collapse = "; ")
  
  # 关键基因
  if (!is.null(enrichment_results[[fac]])) {
    annotation_summary$Key_Genes_Positive[i] <- paste(
      head(enrichment_results[[fac]]$pos_genes, 10), collapse = ", ")
    annotation_summary$Key_Genes_Negative[i] <- paste(
      head(enrichment_results[[fac]]$neg_genes, 10), collapse = ", ")
  }
  
  # 基于通路推断生物学主题
  all_terms <- c(pos_go$Description, neg_go$Description, hallmark$ID)
  
  if (length(all_terms) > 0) {
    # 简单的关键词匹配
    if (any(grepl("metabol|lipid|fatty|oxidation", all_terms, ignore.case = TRUE))) {
      annotation_summary$Biological_Theme[i] <- "Metabolic reprogramming"
    } else if (any(grepl("immune|inflam|cytokine|T cell|B cell", all_terms, ignore.case = TRUE))) {
      annotation_summary$Biological_Theme[i] <- "Immune response"
    } else if (any(grepl("cell cycle|prolif|DNA|mitot", all_terms, ignore.case = TRUE))) {
      annotation_summary$Biological_Theme[i] <- "Cell proliferation"
    } else if (any(grepl("signal|pathway|recept", all_terms, ignore.case = TRUE))) {
      annotation_summary$Biological_Theme[i] <- "Signaling pathways"
    } else if (any(grepl("fibro|ECM|collagen|matrix", all_terms, ignore.case = TRUE))) {
      annotation_summary$Biological_Theme[i] <- "Fibrosis/ECM remodeling"
    } else {
      annotation_summary$Biological_Theme[i] <- "Multi-functional"
    }
  }
}

write.csv(annotation_summary, file.path(out_dir, "factor_biological_annotation.csv"), 
          row.names = FALSE)
cat("  - 因子注释总结保存完成\n")

# ============================================================================
# Step 6: 可视化
# ============================================================================
cat("\n[Step 6] 生成可视化...\n")

# 6.1 因子×通路富集热图
cat("  - 6.1 生成因子通路热图...\n")

if (nrow(all_enrichment_df) > 0) {
  # 准备热图数据：选择Hallmark通路
  hallmark_enrich <- all_enrichment_df %>%
    filter(enrichment_type == "Hallmark") %>%
    mutate(
      pathway = gsub("HALLMARK_", "", ID),
      neg_log_padj = -log10(p.adjust + 1e-10)
    ) %>%
    select(factor, pathway, neg_log_padj, direction) %>%
    mutate(
      neg_log_padj = ifelse(direction == "negative", -neg_log_padj, neg_log_padj)
    )
  
  if (nrow(hallmark_enrich) > 0) {
    # 透视为宽格式
    heatmap_data <- hallmark_enrich %>%
      select(factor, pathway, neg_log_padj) %>%
      pivot_wider(names_from = factor, values_from = neg_log_padj, values_fill = 0)
    
    heatmap_mat <- as.matrix(heatmap_data[, -1])
    rownames(heatmap_mat) <- heatmap_data$pathway
    
    # 生成热图
    pdf(file.path(fig_dir, "factor_pathway_heatmap.pdf"), width = 8, height = 10)
    
    col_fun <- colorRamp2(
      c(-3, 0, 3),
      c("blue", "white", "red")
    )
    
    ht <- Heatmap(
      heatmap_mat,
      name = "-log10(padj)",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 10),
      column_title = "MOFA Factor - Hallmark Pathway Enrichment",
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      heatmap_legend_param = list(
        title = "-log10(padj)",
        title_gp = gpar(fontsize = 9),
        labels_gp = gpar(fontsize = 8)
      )
    )
    draw(ht)
    dev.off()
    cat("    - factor_pathway_heatmap.pdf 生成完成\n")
  }
}

# 6.2 因子跨组学方差贡献堆叠条图
cat("  - 6.2 生成方差贡献图...\n")

p_var <- ggplot(var_df, aes(x = Factor, y = Variance, fill = View)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "MOFA Factor Variance Contribution by View",
    x = "Factor",
    y = "Variance Explained (%)",
    fill = "Omics Layer"
  ) +
  nc_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(fig_dir, "factor_variance_contribution.pdf"), p_var,
       width = 6, height = 4)
cat("    - factor_variance_contribution.pdf 生成完成\n")

# 6.3 因子×临床变量关联热图
cat("  - 6.3 生成临床关联热图...\n")

# 筛选有显著关联的变量
sig_vars <- colnames(pval_adj)[apply(pval_adj, 2, function(x) any(x < 0.1, na.rm = TRUE))]

if (length(sig_vars) > 0) {
  cor_subset <- cor_matrix[, sig_vars, drop = FALSE]
  pval_subset <- pval_adj[, sig_vars, drop = FALSE]
} else {
  # 如果没有显著变量，取相关性最强的10个
  mean_abs_cor <- apply(abs(cor_matrix), 2, mean, na.rm = TRUE)
  top_vars <- names(sort(mean_abs_cor, decreasing = TRUE))[1:min(15, length(mean_abs_cor))]
  cor_subset <- cor_matrix[, top_vars, drop = FALSE]
  pval_subset <- pval_adj[, top_vars, drop = FALSE]
}

# 替换NA为0用于可视化
cor_subset[is.na(cor_subset)] <- 0

pdf(file.path(fig_dir, "factor_clinical_correlation.pdf"), width = 10, height = 5)

col_fun2 <- colorRamp2(c(-0.8, 0, 0.8), c("#2166AC", "white", "#B2182B"))

# 创建星号标记矩阵
sig_marks <- matrix("", nrow = nrow(pval_subset), ncol = ncol(pval_subset))
sig_marks[pval_subset < 0.05] <- "*"
sig_marks[pval_subset < 0.01] <- "**"
sig_marks[pval_subset < 0.001] <- "***"

ht2 <- Heatmap(
  cor_subset,
  name = "Spearman r",
  col = col_fun2,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  column_title = "MOFA Factor - Clinical Variable Correlation",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_marks[i, j], x, y, gp = gpar(fontsize = 8))
  },
  heatmap_legend_param = list(
    title = "Spearman r",
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  )
)
draw(ht2)
dev.off()
cat("    - factor_clinical_correlation.pdf 生成完成\n")

# 6.4 每个因子的Top权重特征条形图
cat("  - 6.4 生成Top权重特征图...\n")

# 只展示前3个因子（方差贡献最高），每个因子top 8特征
factors_to_show <- factor_names[1:min(3, length(factor_names))]

# 增大PDF高度以容纳3个垂直排列的图
pdf(file.path(fig_dir, "factor_top_features.pdf"), width = 8, height = 12)

plot_list <- list()
for (fac in factors_to_show) {
  # 获取该因子的权重 - 只取top 8 (减少以提高可读性)
  fac_weights <- all_weights %>%
    filter(factor == fac) %>%
    arrange(desc(abs(value))) %>%
    head(8) %>%
    mutate(
      feature_label = ifelse(is_metabolite, feature, gene_symbol),
      direction = ifelse(value > 0, "Positive", "Negative")
    )
  
  p <- ggplot(fac_weights, aes(x = reorder(feature_label, value), y = value, fill = view)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = paste(fac, "- Top Features"),
      x = "",
      y = "Weight",
      fill = "View"
    ) +
    nc_theme +
    theme(
      axis.text.y = element_text(size = 9),  # 增大Y轴字体
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom"
    )
  
  plot_list[[fac]] <- p
}

# 组合图 - 3个因子垂直排列，每个有足够空间
combined <- ggarrange(plotlist = plot_list, ncol = 1, nrow = 3, 
                      common.legend = TRUE, legend = "bottom",
                      heights = c(1, 1, 1))
print(combined)
dev.off()
cat("    - factor_top_features.pdf 生成完成\n")

# ============================================================================
# 完成
# ============================================================================
cat("\n============================================================\n")
cat("MOFA因子生物学注释分析完成!\n")
cat("============================================================\n")
cat("\n输出文件:\n")
cat("  CSV文件 (", out_dir, "):\n")
cat("    - factor_pathway_enrichment.csv\n")
cat("    - factor_clinical_correlation.csv\n")
cat("    - factor_clinical_pvalue_adj.csv\n")
cat("    - factor_variance_explained.csv\n")
cat("    - factor_biological_annotation.csv\n")
cat("\n  图片文件 (", fig_dir, "):\n")
cat("    - factor_pathway_heatmap.pdf\n")
cat("    - factor_variance_contribution.pdf\n")
cat("    - factor_clinical_correlation.pdf\n")
cat("    - factor_top_features.pdf\n")

# 打印因子注释摘要
cat("\n因子生物学注释摘要:\n")
print(annotation_summary[, c("Factor", "Biological_Theme", "Top_Hallmark")])
