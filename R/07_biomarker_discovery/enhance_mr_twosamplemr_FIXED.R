#!/usr/bin/env Rscript
# =============================================================================
# enhance_mr_twosamplemr_FIXED.R -- 修复版TwoSampleMR分析 (无超时问题)
# =============================================================================
# 问题：available_outcomes()会超时
# 解决：采用增量查询策略，直接查询所需GWAS而不获取完整列表
# =============================================================================

cat("
╔═══════════════════════════════════════════════════════════════════════════════╗
║      HAE Multi-Omics: TwoSampleMR 修复版 (API超时已解决)                       ║
║      使用增量查询策略 + ieugwasr::gwasinfo()                                   ║
╚═══════════════════════════════════════════════════════════════════════════════╝
\n")

# =============================================================================
# 0. 环境配置
# =============================================================================
cat("=== Step 0: 环境配置 ===\n")

project_root <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
setwd(project_root)

output_dir <- file.path(project_root, "analysis/results/enhancement_mr_fixed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("  输出目录: %s\n", output_dir))

# 关键：增加R超时选项（备用）
options(timeout = 180)
cat("  R超时: 180秒\n")

# JWT Token配置
JWT_TOKEN <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJyaXNoYXQucnV6aUB4am11LmVkdS5jbiIsImlhdCI6MTc3NDI1OTI5NCwiZXhwIjoxNzc1NDY4ODk0fQ.GU2kv_jH3jBhYPavXHpImroVYZSdwc8LgF17bpo4VcU7nPGeW9dBfqzQplmg_JM1XEdaV0FuTz4Fj6CtnKXfPVsQUzAjIqUHfwW8cqLU0bvud9xNVaRoxBYQqA5kZ262fWhhbttr5uYUuFKViA9CB00tDpbRjR2hp-ABrFCQ5bFMq6rf74mrvGoJFzfU25IsScLXEeTBuj_AvRBuv0nMJ7aSeaVQC2AkPcbgnwc99s8VN0UyiT_NH0I6JemFA9EP4BnCR5zhoL8fkdCPfhFmTo813P9qJNh02oHx6qmXq61fhIr9MMVLxqES03y3a7gyEvwdhAQMYNl-jmzBQJndgg"
Sys.setenv(OPENGWAS_JWT = JWT_TOKEN)
cat("  JWT Token已配置 (有效期至 2026-04-04)\n")

# 加载基础包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

source("analysis/scripts/nc_theme.R")

# 加载TwoSampleMR及相关包
cat("\n检查MR相关包...\n")
pkg_loaded <- list()

if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(TwoSampleMR))
  pkg_loaded$TwoSampleMR <- packageVersion("TwoSampleMR")
  cat(sprintf("  [OK] TwoSampleMR v%s 已加载\n", pkg_loaded$TwoSampleMR))
} else {
  stop("TwoSampleMR包未安装")
}

if (requireNamespace("ieugwasr", quietly = TRUE)) {
  suppressPackageStartupMessages(suppressWarnings(library(ieugwasr)))
  cat("  [OK] ieugwasr 已加载 (采用增量查询模式)\n")
}

if (requireNamespace("MRPRESSO", quietly = TRUE)) {
  suppressPackageStartupMessages(library(MRPRESSO))
  cat("  [OK] MRPRESSO 已加载\n")
}

if (requireNamespace("RadialMR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(RadialMR))
  cat("  [OK] RadialMR 已加载\n")
}

# =============================================================================
# 1. 定义GWAS数据集（关键创新：不调用available_outcomes()）
# =============================================================================
cat("\n=== Step 1: 定义所需GWAS数据集 ===\n")

# HAE相关的UK Biobank GWAS datasets
gwas_datasets <- list(
  # 免疫细胞比例
  Lymphocyte_pct = "ieu-b-9760",      # 修改：使用正确的IEU ID
  Monocyte_pct = "ieu-b-8189",
  Neutrophil_pct = "ieu-b-8190",
  
  # 脂质代谢
  CRP = "ieu-b-9715",                 # C-reactive protein
  HDL = "ieu-b-109",                  # HDL cholesterol
  LDL = "ieu-b-110",                  # LDL cholesterol
  Triglycerides = "ieu-b-111",        # Triglycerides
  
  # 肝酶
  ALT = "ieu-b-30",                   # Alanine aminotransferase (代理)
  GGT = "ieu-b-39",                   # GGT (代理)
  AST = "ieu-b-31"                    # Aspartate aminotransferase (代理)
)

cat(sprintf("  已定义 %d 个GWAS datasets\n", length(gwas_datasets)))
for (i in seq_along(gwas_datasets)) {
  cat(sprintf("    %2d. %s = %s\n", i, names(gwas_datasets)[i], gwas_datasets[[i]]))
}

# =============================================================================
# 2. 快速验证GWAS可用性（增量查询，不会超时）
# =============================================================================
cat("\n=== Step 2: 验证GWAS数据集可用性 ===\n")

gwas_available <- list()
for (name in names(gwas_datasets)) {
  id <- gwas_datasets[[name]]
  cat(sprintf("  检查 %s (%s)...", name, id))
  
  tryCatch({
    info <- gwasinfo(id)
    if (nrow(info) > 0) {
      gwas_available[[name]] <- list(
        id = id,
        n_snp = info$snp[1],
        sample_size = info$sample_size[1],
        year = info$year[1]
      )
      cat(sprintf(" ✓ OK (%s SNPs, n=%s)\n", 
                  gwas_available[[name]]$n_snp,
                  gwas_available[[name]]$sample_size))
    } else {
      cat(" ✗ 无效\n")
    }
  }, error = function(e) {
    cat(sprintf(" ✗ 错误: %s\n", substr(e$message, 1, 50)))
  })
}

cat(sprintf("\n  可用GWAS数: %d/%d\n", length(gwas_available), length(gwas_datasets)))

# =============================================================================
# 3. 配置参数
# =============================================================================
cat("\n=== Step 3: 配置MR参数 ===\n")

MR_CONFIG <- list(
  pval_threshold = 5e-8,
  clump_r2 = 0.001,
  clump_kb = 10000,
  f_stat_threshold = 10,
  dpi = 300,
  width_mm = 183,
  height_mm = 150
)

cat(sprintf("  p值阈值: %s\n", format(MR_CONFIG$pval_threshold, scientific = TRUE)))
cat(sprintf("  LD clumping: r² < %.3f, window = %d kb\n", MR_CONFIG$clump_r2, MR_CONFIG$clump_kb))
cat(sprintf("  F统计量阈值: > %d\n", MR_CONFIG$f_stat_threshold))

# =============================================================================
# 4. 定义MR分析方案
# =============================================================================
cat("\n=== Step 4: 定义MR分析方案 ===\n")

schemes <- list(
  list(
    name = "Lipids_to_LiverEnzymes",
    exposures = list(LDL = "ieu-b-110", Triglycerides = "ieu-b-111", HDL = "ieu-b-109"),
    outcomes = list(ALT = "ieu-b-30", GGT = "ieu-b-39")
  ),
  list(
    name = "CRP_to_LiverEnzymes", 
    exposures = list(CRP = "ieu-b-9715"),
    outcomes = list(ALT = "ieu-b-30", GGT = "ieu-b-39")
  )
)

for (s in schemes) {
  cat(sprintf("  %s: %d exposures × %d outcomes\n", 
              s$name, length(s$exposures), length(s$outcomes)))
}

# =============================================================================
# 5. MR分析函数
# =============================================================================
cat("\n=== Step 5: 定义MR分析函数 ===\n")

run_mr_analysis <- function(exp_id, out_id, exp_name, out_name) {
  cat(sprintf("\n  >>> %s → %s\n", exp_name, out_name))
  
  result <- list(success = FALSE, data = NULL, error = NULL)
  
  tryCatch({
    # 提取工具变量
    cat(sprintf("      提取工具变量 (extract_instruments %s)...\n", exp_id))
    exp_dat <- extract_instruments(
      outcomes = exp_id,
      p1 = MR_CONFIG$pval_threshold,
      clump = TRUE,
      r2 = MR_CONFIG$clump_r2,
      kb = MR_CONFIG$clump_kb
    )
    
    if (is.null(exp_dat) || nrow(exp_dat) == 0) {
      result$error <- "无显著工具变量"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    cat(sprintf("      获取 %d 个SNP\n", nrow(exp_dat)))
    
    # F统计量过滤
    exp_dat$F_stat <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
    exp_dat <- exp_dat[exp_dat$F_stat >= MR_CONFIG$f_stat_threshold, ]
    
    if (nrow(exp_dat) < 3) {
      result$error <- "F统计量过滤后工具变量不足"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    cat(sprintf("      F>%d后: %d SNPs (平均F=%.1f)\n", 
                MR_CONFIG$f_stat_threshold, nrow(exp_dat), mean(exp_dat$F_stat)))
    
    # 提取结局数据
    cat(sprintf("      提取结局数据 (extract_outcome_data %s)...\n", out_id))
    out_dat <- extract_outcome_data(
      snps = exp_dat$SNP,
      outcomes = out_id
    )
    
    if (is.null(out_dat) || nrow(out_dat) == 0) {
      result$error <- "无匹配的结局数据"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    cat(sprintf("      获取 %d 个SNP\n", nrow(out_dat)))
    
    # 协调
    dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
    cat(sprintf("      协调后: %d 个SNP\n", nrow(dat)))
    
    if (nrow(dat) < 3) {
      result$error <- "协调后工具变量不足"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    
    # MR分析
    mr_results <- mr(dat)
    cat(sprintf("      MR分析: %d 种方法\n", nrow(mr_results)))
    
    result$success <- TRUE
    result$data <- list(
      mr = mr_results,
      dat = dat,
      exposure = exp_name,
      outcome = out_name
    )
    
  }, error = function(e) {
    result$error <<- e$message
    cat(sprintf("      [错误] %s\n", substr(e$message, 1, 80)))
  })
  
  return(result)
}

# =============================================================================
# 6. 执行MR分析
# =============================================================================
cat("\n=== Step 6: 执行MR分析 ===\n")

all_results <- list()
result_count <- 0

for (scheme in schemes) {
  cat(sprintf("\n方案: %s\n", scheme$name))
  
  for (exp_name in names(scheme$exposures)) {
    for (out_name in names(scheme$outcomes)) {
      exp_id <- scheme$exposures[[exp_name]]
      out_id <- scheme$outcomes[[out_name]]
      
      res <- run_mr_analysis(exp_id, out_id, exp_name, out_name)
      
      if (res$success) {
        result_count <- result_count + 1
        all_results[[sprintf("%s_%s", exp_name, out_name)]] <- res
      }
    }
  }
}

cat(sprintf("\n=== 成功完成 %d 个分析 ===\n", result_count))

# =============================================================================
# 7. 汇总结果
# =============================================================================
cat("\n=== Step 7: 汇总结果 ===\n")

if (result_count > 0) {
  # 提取IVW结果
  ivw_results <- data.frame()
  
  for (name in names(all_results)) {
    res <- all_results[[name]]
    ivw <- res$data$mr %>%
      filter(method == "Inverse variance weighted") %>%
      select(exposure, outcome, nsnp, b = b, se, pval, OR = or, OR_lci95 = or_lci95, OR_uci95 = or_uci95)
    
    if (nrow(ivw) > 0) {
      ivw_results <- bind_rows(ivw_results, ivw)
    }
  }
  
  cat("\n主要MR结果 (IVW方法):\n")
  print(ivw_results)
  
  # 保存结果
  if (nrow(ivw_results) > 0) {
    write.csv(ivw_results, 
              file.path(output_dir, "mr_ivw_results.csv"),
              row.names = FALSE)
    cat(sprintf("\n✓ 结果已保存: %s\n", file.path(output_dir, "mr_ivw_results.csv")))
  }
}

# =============================================================================
# 完成
# =============================================================================
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║      分析完成！(修复版 - 无超时问题)                                           ║\n")
cat("║      输出目录: ")
cat(output_dir)
cat("                                              ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════════════╝\n")

cat("\n=== 运行信息 ===\n")
cat(sprintf("R版本: %s\n", R.version$version.string))
cat(sprintf("运行时间: %s\n", Sys.time()))
cat(sprintf("工作目录: %s\n", getwd()))

