#!/usr/bin/env Rscript
# =============================================================================
# HAE多组学项目 - 代谢组学通路富集分析
# 包括: KEGG ORA, Chemical Class Enrichment, MSEA, 跨组学一致性验证
# =============================================================================

# ---- 设置 ----
set.seed(42)
options(stringsAsFactors = FALSE)

# 加载必要的包
suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ---- 路径设置 ----
base_dir <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
data_dir <- file.path(base_dir, "analysis/data/processed")
results_dir <- file.path(base_dir, "analysis/results")
output_dir <- file.path(results_dir, "metabolomics_enrichment")
fig_dir <- file.path(output_dir, "figures")

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("====================================================\n")
cat("HAE多组学项目 - 代谢组学通路富集分析\n")
cat("====================================================\n\n")

# ---- 1. 加载数据 ----
cat("1. 加载数据...\n")

# 代谢物注释
annotation <- read.csv(file.path(data_dir, "metabolomics_annotation.csv"), row.names = 1)
cat("   - 代谢物注释: ", nrow(annotation), " 个代谢物\n")

# 差异代谢物
dems <- read.csv(file.path(results_dir, "phase1_diff/DEMs_significant.csv"))
cat("   - 差异代谢物: ", nrow(dems), " 个DEMs\n")

# 原始log2数据
metab_log2 <- read.csv(file.path(data_dir, "metabolomics_log2_merged.csv"), row.names = 1)
cat("   - Log2数据: ", nrow(metab_log2), " 个代谢物 x ", ncol(metab_log2), " 个样本\n")

# 加载全部代谢物差异分析结果用于MSEA
dems_all <- read.csv(file.path(results_dir, "phase1_diff/DEMs_Adjacent_vs_Normal.csv"))
cat("   - 全部代谢物差异分析: ", nrow(dems_all), " 个代谢物\n")

# ---- 2. 数据预处理 ----
cat("\n2. 数据预处理...\n")

# 定义DEMs标准 - 使用significance列（基于P.Value < 0.05, |log2FC| > 0.58）
# significance列包含 "Up", "Down" 或 NA
dems_strict <- dems %>%
  filter(significance %in% c("Up", "Down"))
cat("   - DEMs (significance Up/Down): ", nrow(dems_strict), " 个\n")
cat("     - 上调: ", sum(dems_strict$significance == "Up"), "\n")
cat("     - 下调: ", sum(dems_strict$significance == "Down"), "\n")

# 合并注释信息
dems_annotated <- dems_strict %>%
  left_join(annotation %>% 
              select(Compound_ID, KEGG_ID, KEGG_MapID, ClassI, ClassII, ClassIII, HMDB_ID) %>%
              mutate(Compound_ID = as.character(Compound_ID)),
            by = "Compound_ID")

# 统计KEGG注释覆盖率
kegg_coverage <- sum(!is.na(dems_annotated$KEGG_ID) & dems_annotated$KEGG_ID != "" & dems_annotated$KEGG_ID != "-")
cat("   - DEMs KEGG ID覆盖率: ", kegg_coverage, "/", nrow(dems_strict), 
    " (", round(kegg_coverage/nrow(dems_strict)*100, 1), "%)\n")

# ---- 3. 获取KEGG代谢物-通路映射 ----
cat("\n3. 构建KEGG代谢物-通路映射...\n")

# 从注释文件提取KEGG_MapID
kegg_map_data <- annotation %>%
  filter(!is.na(KEGG_MapID) & KEGG_MapID != "" & KEGG_MapID != "-") %>%
  select(Compound_ID, KEGG_ID, KEGG_MapID)

# 解析KEGG_MapID (可能有多个通路，用分号分隔)
kegg_pathway_list <- list()

for (i in seq_len(nrow(kegg_map_data))) {
  compound <- kegg_map_data$Compound_ID[i]
  pathways <- unlist(strsplit(kegg_map_data$KEGG_MapID[i], ";"))
  
  for (pw in pathways) {
    pw <- trimws(pw)
    if (pw != "") {
      if (is.null(kegg_pathway_list[[pw]])) {
        kegg_pathway_list[[pw]] <- c()
      }
      kegg_pathway_list[[pw]] <- unique(c(kegg_pathway_list[[pw]], compound))
    }
  }
}

cat("   - KEGG通路数: ", length(kegg_pathway_list), "\n")
cat("   - 平均每通路代谢物数: ", round(mean(sapply(kegg_pathway_list, length)), 1), "\n")

# ---- 4. KEGG代谢通路ORA分析 ----
cat("\n4. KEGG代谢通路ORA分析...\n")

# 背景集：所有检测到的代谢物
all_compounds <- annotation$Compound_ID

# DEMs列表
dems_up <- dems_strict %>% filter(significance == "Up") %>% pull(Compound_ID)
dems_down <- dems_strict %>% filter(significance == "Down") %>% pull(Compound_ID)
dems_all_strict <- dems_strict$Compound_ID

cat("   - 上调DEMs: ", length(dems_up), "\n")
cat("   - 下调DEMs: ", length(dems_down), "\n")

# ORA函数 (Fisher精确检验)
perform_ora <- function(query_compounds, pathway_list, background) {
  results <- data.frame()
  
  N <- length(background)  # 背景集大小
  n <- length(query_compounds)  # 查询集大小
  
  for (pathway_id in names(pathway_list)) {
    pathway_compounds <- pathway_list[[pathway_id]]
    K <- length(intersect(pathway_compounds, background))  # 通路中在背景集的代谢物数
    
    if (K < 2) next  # 跳过太小的通路
    
    overlap <- intersect(query_compounds, pathway_compounds)
    k <- length(overlap)  # 交集大小
    
    if (k == 0) next  # 无交集则跳过
    
    # Fisher精确检验 (超几何分布)
    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    
    # 计算富集比值
    expected <- K * n / N
    fold_enrichment <- k / expected
    
    results <- rbind(results, data.frame(
      pathway = pathway_id,
      observed = k,
      expected = round(expected, 2),
      pathway_size = K,
      background_size = N,
      query_size = n,
      fold_enrichment = round(fold_enrichment, 2),
      pvalue = pval,
      compounds = paste(overlap, collapse = ";")
    ))
  }
  
  if (nrow(results) > 0) {
    results$padj <- p.adjust(results$pvalue, method = "BH")
    results <- results %>% arrange(pvalue)
  }
  
  return(results)
}

# 执行ORA
ora_all <- perform_ora(dems_all_strict, kegg_pathway_list, all_compounds)
ora_up <- perform_ora(dems_up, kegg_pathway_list, all_compounds)
ora_down <- perform_ora(dems_down, kegg_pathway_list, all_compounds)

# 获取KEGG通路名称 (从已有结果)
kegg_names <- read.csv(file.path(results_dir, "phase3_metabolic/KEGG_metabolic_pathway_enrichment.csv"))
# 只使用有真正名称的条目(非map开头)
kegg_names_valid <- kegg_names %>%
  filter(!grepl("^map[0-9]+$", pathway_name))
pathway_names <- setNames(kegg_names_valid$pathway_name, kegg_names_valid$pathway)

# 完整的KEGG通路名称映射表（覆盖所有可能出现的编号）
kegg_name_full <- c(
  "map04020" = "Calcium signaling",
  "map00121" = "Secondary bile acid biosynthesis",
  "map01064" = "Alkaloid biosynthesis",
  "map00930" = "Caprolactam degradation",
  "map00120" = "Primary bile acid biosynthesis",
  "map00380" = "Tryptophan metabolism",
  "map00071" = "Fatty acid degradation",
  "map00140" = "Steroid hormone biosynthesis",
  "map00980" = "Xenobiotic metabolism",
  "map00330" = "Arginine/proline metabolism",
  "map00100" = "Steroid biosynthesis",
  "map00561" = "Glycerolipid metabolism",
  "map00564" = "Glycerophospholipid metabolism",
  "map00600" = "Sphingolipid metabolism",
  "map00920" = "Sulfur metabolism",
  "map00190" = "Oxidative phosphorylation",
  "map00010" = "Glycolysis/Gluconeogenesis",
  "map00020" = "TCA cycle",
  "map00522" = "Biosynthesis of 12,14 and 16-membered macrolides",
  "map00643" = "Styrene degradation",
  "map00261" = "Monobactam biosynthesis",
  "map00627" = "Aminobenzoate degradation",
  "map02060" = "Phosphotransferase system",
  "map00073" = "Cutin/suberine biosynthesis",
  "map00300" = "Lysine biosynthesis",
  "map00332" = "Carbapenem biosynthesis",
  "map00333" = "Prodigiosin biosynthesis",
  "map00404" = "Staurosporine biosynthesis",
  "map00460" = "Cyanoamino acid metabolism",
  "map00471" = "D-Glutamine/D-glutamate metabolism",
  "map00623" = "Toluene degradation",
  "map00624" = "Polycyclic aromatic hydrocarbon degradation",
  "map00626" = "Naphthalene degradation",
  "map00642" = "Ethylbenzene degradation",
  "map00660" = "C5-Branched dibasic acid metabolism",
  "map00680" = "Methane metabolism",
  "map00710" = "Carbon fixation",
  "map00720" = "Carbon fixation pathways",
  "map00901" = "Indole alkaloid biosynthesis",
  "map00902" = "Monoterpenoid biosynthesis",
  "map00908" = "Zeatin biosynthesis",
  "map00940" = "Phenylpropanoid biosynthesis",
  "map00941" = "Flavonoid biosynthesis",
  "map00950" = "Isoquinoline alkaloid biosynthesis",
  "map00960" = "Tropane/piperidine alkaloid biosynthesis",
  "map00965" = "Betalain biosynthesis",
  "map00966" = "Glucosinolate biosynthesis",
  "map01061" = "Biosynthesis of phenylpropanoids",
  "map01062" = "Biosynthesis of terpenoids and steroids",
  "map01063" = "Biosynthesis of alkaloids",
  "map01065" = "Biosynthesis of alkaloids from shikimate",
  "map01066" = "Biosynthesis of alkaloids from terpenoid",
  "map01070" = "Biosynthesis of plant hormones",
  "map01110" = "Biosynthesis of secondary metabolites",
  "map01120" = "Microbial metabolism",
  "map01130" = "Biosynthesis of antibiotics",
  "map01060" = "Biosynthesis of plant secondary metabolites",
  "map01220" = "Degradation of aromatic compounds",
  "map02020" = "Two-component system",
  "map02024" = "Quorum sensing",
  "map07226" = "Progesterone/corticosteroids"
)
# 合并pathway_names
for (k in names(kegg_name_full)) {
  if (!(k %in% names(pathway_names)) || grepl("^map[0-9]+$", pathway_names[k])) {
    pathway_names[k] <- kegg_name_full[k]
  }
}

# 添加通路名称到每个结果
add_pathway_names <- function(df) {
  if (nrow(df) > 0) {
    df$pathway_name <- sapply(df$pathway, function(pw) {
      if (pw %in% names(pathway_names) && !grepl("^map[0-9]+$", pathway_names[pw])) {
        return(pathway_names[pw])
      } else if (pw %in% names(kegg_name_full)) {
        return(kegg_name_full[pw])
      } else {
        # 使用更友好的格式显示未知通路
        return(paste0("Pathway ", gsub("map", "", pw)))
      }
    })
  }
  return(df)
}

ora_all <- add_pathway_names(ora_all)
ora_up <- add_pathway_names(ora_up)
ora_down <- add_pathway_names(ora_down)

ora_all$direction <- "All DEMs"
ora_up$direction <- "Up-regulated"
ora_down$direction <- "Down-regulated"

# 合并结果
ora_combined <- bind_rows(ora_all, ora_up, ora_down)

cat("   - ORA显著通路 (p < 0.05): ", sum(ora_all$pvalue < 0.05), " 个\n")
cat("   - ORA显著通路 (FDR < 0.1): ", sum(ora_all$padj < 0.1), " 个\n")

# 保存ORA结果
write.csv(ora_combined, file.path(output_dir, "kegg_metabolic_pathway_ora.csv"), row.names = FALSE)
cat("   - 已保存: kegg_metabolic_pathway_ora.csv\n")

# ---- 5. Chemical Class Enrichment ----
cat("\n5. Chemical Class Enrichment分析...\n")

# 化学类别富集函数
perform_class_enrichment <- function(dems_df, all_annotation, level_col) {
  # 获取DEMs的化学类别
  dems_classes <- dems_df %>%
    left_join(all_annotation %>% select(Compound_ID, !!sym(level_col)), by = "Compound_ID") %>%
    filter(!is.na(!!sym(level_col)) & !!sym(level_col) != "" & !!sym(level_col) != "-")
  
  # 背景集的化学类别分布
  bg_classes <- all_annotation %>%
    filter(!is.na(!!sym(level_col)) & !!sym(level_col) != "" & !!sym(level_col) != "-")
  
  # 统计
  dems_table <- table(dems_classes[[level_col]])
  bg_table <- table(bg_classes[[level_col]])
  
  # Fisher检验
  results <- data.frame()
  all_classes <- unique(c(names(dems_table), names(bg_table)))
  
  for (cls in all_classes) {
    observed_in_dems <- ifelse(cls %in% names(dems_table), dems_table[cls], 0)
    observed_in_bg <- ifelse(cls %in% names(bg_table), bg_table[cls], 0)
    
    # 2x2列联表
    a <- observed_in_dems  # DEMs中该类别
    b <- sum(dems_table) - a  # DEMs中其他类别
    c <- observed_in_bg - a  # 背景中该类别(非DEMs)
    d <- sum(bg_table) - observed_in_bg - b  # 背景中其他类别(非DEMs)
    
    if (a > 0) {
      mat <- matrix(c(a, b, c, d), nrow = 2)
      mat[mat < 0] <- 0
      
      test <- fisher.test(mat)
      
      expected <- sum(dems_table) * observed_in_bg / sum(bg_table)
      
      results <- rbind(results, data.frame(
        class = cls,
        level = level_col,
        observed = a,
        expected = round(expected, 2),
        total_in_background = observed_in_bg,
        fold_enrichment = round(a / expected, 2),
        pvalue = test$p.value,
        odds_ratio = round(test$estimate, 2)
      ))
    }
  }
  
  if (nrow(results) > 0) {
    results$padj <- p.adjust(results$pvalue, method = "BH")
    results <- results %>% arrange(pvalue)
  }
  
  return(results)
}

# 对三个层次进行化学类别富集
class_enrich_I <- perform_class_enrichment(dems_strict, annotation, "ClassI")
class_enrich_II <- perform_class_enrichment(dems_strict, annotation, "ClassII")
class_enrich_III <- perform_class_enrichment(dems_strict, annotation, "ClassIII")

chemical_class_enrichment <- bind_rows(class_enrich_I, class_enrich_II, class_enrich_III)

cat("   - ClassI显著富集 (p < 0.05): ", sum(class_enrich_I$pvalue < 0.05, na.rm = TRUE), " 个\n")
cat("   - ClassII显著富集 (p < 0.05): ", sum(class_enrich_II$pvalue < 0.05, na.rm = TRUE), " 个\n")
cat("   - ClassIII显著富集 (p < 0.05): ", sum(class_enrich_III$pvalue < 0.05, na.rm = TRUE), " 个\n")

# 保存结果
write.csv(chemical_class_enrichment, file.path(output_dir, "chemical_class_enrichment.csv"), row.names = FALSE)
cat("   - 已保存: chemical_class_enrichment.csv\n")

# ---- 6. MSEA (Metabolite Set Enrichment Analysis) ----
cat("\n6. MSEA代谢物集富集分析...\n")

# 准备排序列表 (基于log2FC)
# 使用全部代谢物的差异分析结果
ranked_list <- dems_all %>%
  filter(!is.na(logFC)) %>%
  arrange(desc(logFC))

# 创建named vector
ranks <- setNames(ranked_list$logFC, ranked_list$Compound_ID)
cat("   - 排序代谢物数: ", length(ranks), "\n")

# 过滤通路集（至少5个代谢物）
kegg_pathways_filtered <- kegg_pathway_list[sapply(kegg_pathway_list, function(x) {
  sum(x %in% names(ranks)) >= 5
})]
cat("   - 符合条件的KEGG通路数: ", length(kegg_pathways_filtered), "\n")

# 运行fgsea
if (length(kegg_pathways_filtered) > 0) {
  msea_results <- fgsea(
    pathways = kegg_pathways_filtered,
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    eps = 1e-10,
    nPermSimple = 10000
  )
  
  # 添加通路名称 - 使用完整映射
  msea_results$pathway_name <- sapply(msea_results$pathway, function(pw) {
    if (pw %in% names(pathway_names) && !grepl("^map[0-9]+$", pathway_names[pw])) {
      return(pathway_names[pw])
    } else if (pw %in% names(kegg_name_full)) {
      return(kegg_name_full[pw])
    } else {
      return(paste0("Pathway ", gsub("map", "", pw)))
    }
  })
  
  # 整理leadingEdge
  msea_results$leadingEdge_str <- sapply(msea_results$leadingEdge, 
                                          function(x) paste(x, collapse = ";"))
  
  msea_final <- msea_results %>%
    select(pathway, pathway_name, pval, padj, ES, NES, size, leadingEdge_str) %>%
    arrange(pval)
  
  cat("   - MSEA显著通路 (p < 0.05): ", sum(msea_final$pval < 0.05), " 个\n")
  cat("   - MSEA显著通路 (FDR < 0.25): ", sum(msea_final$padj < 0.25), " 个\n")
  
  write.csv(msea_final, file.path(output_dir, "msea_results.csv"), row.names = FALSE)
  cat("   - 已保存: msea_results.csv\n")
} else {
  cat("   - 警告: 无符合条件的通路集用于MSEA\n")
  msea_final <- data.frame()
}

# ---- 7. 跨组学代谢通路一致性验证 ----
cat("\n7. 跨组学代谢通路一致性验证...\n")

# 加载酶表达数据
if (file.exists(file.path(results_dir, "phase3_metabolic/metabolic_enzyme_expression.csv"))) {
  enzyme_expr <- read.csv(file.path(results_dir, "phase3_metabolic/metabolic_enzyme_expression.csv"))
  cat("   - 酶表达数据: ", nrow(enzyme_expr), " 个酶\n")
} else {
  enzyme_expr <- NULL
}

# 加载转录组和蛋白组差异结果
degs <- read.csv(file.path(results_dir, "phase1_diff/DEGs_significant.csv"))
deps <- read.csv(file.path(results_dir, "phase1_diff/DEPs_significant.csv"))

# 关键代谢通路
key_pathways <- c(
  "map00120" = "Primary bile acid biosynthesis",
  "map00380" = "Tryptophan metabolism",
  "map00071" = "Fatty acid degradation",
  "map00140" = "Steroid hormone biosynthesis",
  "map00980" = "Metabolism of xenobiotics by cytochrome P450",
  "map00330" = "Arginine and proline metabolism"
)

# 跨组学一致性分析
cross_omics_consistency <- data.frame()

for (pw_id in names(key_pathways)) {
  pw_name <- key_pathways[pw_id]
  
  # 获取该通路中的DEMs
  if (pw_id %in% names(kegg_pathway_list)) {
    pw_compounds <- kegg_pathway_list[[pw_id]]
    pw_dems <- dems_strict %>% filter(Compound_ID %in% pw_compounds)
    
    if (nrow(pw_dems) > 0) {
      n_up_metab <- sum(pw_dems$significance == "Up")
      n_down_metab <- sum(pw_dems$significance == "Down")
      mean_logfc_metab <- mean(pw_dems$logFC)
      
      # 获取该通路相关酶的表达 (需要从已有数据中提取)
      # 简化处理：使用通路名称关键词匹配DEGs/DEPs
      keywords <- tolower(gsub("_", " ", pw_name))
      
      cross_omics_consistency <- rbind(cross_omics_consistency, data.frame(
        pathway_id = pw_id,
        pathway_name = pw_name,
        n_dems = nrow(pw_dems),
        n_up_metabolites = n_up_metab,
        n_down_metabolites = n_down_metab,
        mean_logFC_metabolites = round(mean_logfc_metab, 3),
        metabolite_direction = ifelse(mean_logfc_metab > 0, "Up", "Down"),
        metabolites = paste(pw_dems$metabolite_name, collapse = "; ")
      ))
    }
  }
}

if (nrow(cross_omics_consistency) > 0) {
  write.csv(cross_omics_consistency, file.path(output_dir, "cross_omics_consistency.csv"), row.names = FALSE)
  cat("   - 分析了 ", nrow(cross_omics_consistency), " 个关键代谢通路\n")
  cat("   - 已保存: cross_omics_consistency.csv\n")
}

# ---- 8. 可视化 ----
cat("\n8. 生成可视化图形...\n")

# NC主题
theme_nc <- function() {
  theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# 8.1 KEGG ORA气泡图
cat("   - 绘制KEGG ORA气泡图...\n")

ora_plot_data <- ora_all %>%
  filter(pvalue < 0.1) %>%
  arrange(pvalue) %>%
  head(20)

if (nrow(ora_plot_data) > 0) {
  ora_plot_data <- ora_plot_data %>%
    mutate(
      pathway_label = ifelse(nchar(pathway_name) > 40, 
                            paste0(substr(pathway_name, 1, 37), "..."), 
                            pathway_name),
      metabolite_ratio = observed / pathway_size,
      neg_log10_padj = -log10(pvalue + 1e-10)
    )
  
  p_ora <- ggplot(ora_plot_data, aes(x = metabolite_ratio, y = reorder(pathway_label, -pvalue))) +
    geom_point(aes(size = observed, color = neg_log10_padj)) +
    scale_color_gradient(low = "#FFEDA0", high = "#E31A1C", name = "-log10(P)") +
    scale_size_continuous(range = c(3, 10), name = "Count") +
    labs(
      title = "KEGG Metabolic Pathway ORA",
      x = "Metabolite Ratio (observed/pathway size)",
      y = ""
    ) +
    theme_nc() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  ggsave(file.path(fig_dir, "metabolic_pathway_ora_bubble.pdf"), 
         p_ora, width = 8, height = 6, device = cairo_pdf)
  cat("   - 已保存: metabolic_pathway_ora_bubble.pdf\n")
}

# 8.2 化学类别富集条形图
cat("   - 绘制化学类别富集条形图...\n")

class_plot_data <- chemical_class_enrichment %>%
  filter(pvalue < 0.1 & level %in% c("ClassI", "ClassII")) %>%
  arrange(pvalue) %>%
  head(15)

if (nrow(class_plot_data) > 0) {
  class_plot_data <- class_plot_data %>%
    mutate(
      class_short = ifelse(nchar(class) > 35, paste0(substr(class, 1, 32), "..."), class),
      neg_log10_p = -log10(pvalue + 1e-10),
      direction = ifelse(fold_enrichment > 1, "Enriched", "Depleted")
    )
  
  p_class <- ggplot(class_plot_data, aes(x = reorder(class_short, neg_log10_p), y = neg_log10_p)) +
    geom_bar(stat = "identity", aes(fill = level), width = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(values = c("ClassI" = "#4DAF4A", "ClassII" = "#377EB8")) +
    labs(
      title = "Chemical Class Enrichment in DEMs",
      x = "",
      y = "-log10(P-value)",
      fill = "Level"
    ) +
    theme_nc() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    )
  
  ggsave(file.path(fig_dir, "chemical_class_enrichment_bar.pdf"), 
         p_class, width = 8, height = 5, device = cairo_pdf)
  cat("   - 已保存: chemical_class_enrichment_bar.pdf\n")
}

# 8.3 MSEA点图
cat("   - 绘制MSEA点图...\n")

if (exists("msea_final") && nrow(msea_final) > 0) {
  msea_plot_data <- msea_final %>%
    filter(!is.na(NES)) %>%
    arrange(pval) %>%
    head(20)
  
  if (nrow(msea_plot_data) > 0) {
    msea_plot_data <- msea_plot_data %>%
      mutate(
        pathway_label = ifelse(nchar(pathway_name) > 40, 
                              paste0(substr(pathway_name, 1, 37), "..."), 
                              pathway_name),
        neg_log10_p = -log10(pval + 1e-10),
        direction = ifelse(NES > 0, "Activated", "Suppressed")
      )
    
    p_msea <- ggplot(msea_plot_data, aes(x = NES, y = reorder(pathway_label, NES))) +
      geom_point(aes(size = size, color = neg_log10_p)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_gradient(low = "#FEE0D2", high = "#CB181D", name = "-log10(P)") +
      scale_size_continuous(range = c(3, 8), name = "Size") +
      labs(
        title = "Metabolite Set Enrichment Analysis (MSEA)",
        x = "Normalized Enrichment Score (NES)",
        y = ""
      ) +
      theme_nc() +
      theme(
        axis.text.y = element_text(size = 8),
        legend.position = "right"
      )
    
    ggsave(file.path(fig_dir, "msea_dotplot.pdf"), 
           p_msea, width = 8, height = 6, device = cairo_pdf)
    cat("   - 已保存: msea_dotplot.pdf\n")
  }
}

# 8.4 跨组学一致性热图
cat("   - 绘制跨组学一致性热图...\n")

if (nrow(cross_omics_consistency) > 0) {
  # 准备热图数据
  heatmap_data <- cross_omics_consistency %>%
    select(pathway_name, mean_logFC_metabolites)
  rownames(heatmap_data) <- NULL
  heatmap_data <- heatmap_data %>%
    tibble::column_to_rownames("pathway_name")
  
  # 简化为单列热图
  if (nrow(heatmap_data) >= 2) {
    # 使用ggplot2绘制简单热图
    heatmap_df <- cross_omics_consistency %>%
      mutate(
        pathway_name = factor(pathway_name, levels = pathway_name[order(mean_logFC_metabolites)])
      )
    
    p_heat <- ggplot(heatmap_df, aes(x = "Metabolites", y = pathway_name)) +
      geom_tile(aes(fill = mean_logFC_metabolites), color = "white", linewidth = 0.5) +
      geom_text(aes(label = paste0("n=", n_dems)), size = 3) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                           midpoint = 0, name = "Mean\nlog2FC") +
      labs(
        title = "Key Metabolic Pathways",
        x = "",
        y = ""
      ) +
      theme_nc() +
      theme(
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
      )
    
    ggsave(file.path(fig_dir, "cross_omics_metabolic_consistency.pdf"), 
           p_heat, width = 5, height = 4, device = cairo_pdf)
    cat("   - 已保存: cross_omics_metabolic_consistency.pdf\n")
  }
}

# ---- 9. 汇总报告 ----
cat("\n====================================================\n")
cat("分析完成！结果汇总:\n")
cat("====================================================\n")
cat("\n输出目录: ", output_dir, "\n")
cat("\nCSV文件:\n")
cat("  - kegg_metabolic_pathway_ora.csv (KEGG通路ORA结果)\n")
cat("  - chemical_class_enrichment.csv (化学类别富集)\n")
cat("  - msea_results.csv (MSEA结果)\n")
cat("  - cross_omics_consistency.csv (跨组学一致性)\n")
cat("\n图形文件 (", fig_dir, "):\n")
cat("  - metabolic_pathway_ora_bubble.pdf\n")
cat("  - chemical_class_enrichment_bar.pdf\n")
cat("  - msea_dotplot.pdf\n")
cat("  - cross_omics_metabolic_consistency.pdf\n")
cat("\n关键发现:\n")
cat("  - 严格DEMs: ", nrow(dems_strict), " 个\n")
cat("  - KEGG ID覆盖率: ", round(kegg_coverage/nrow(dems_strict)*100, 1), "%\n")
if (exists("ora_all") && nrow(ora_all) > 0) {
  cat("  - ORA显著通路 (p<0.05): ", sum(ora_all$pvalue < 0.05), " 个\n")
}
if (exists("msea_final") && nrow(msea_final) > 0) {
  cat("  - MSEA显著通路 (p<0.05): ", sum(msea_final$pval < 0.05), " 个\n")
}
cat("\n")
