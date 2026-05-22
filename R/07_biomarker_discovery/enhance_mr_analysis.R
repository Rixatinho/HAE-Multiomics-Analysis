#!/usr/bin/env Rscript
# =============================================================================
# enhance_mr_analysis.R -- 孟德尔随机化分析脚本
# =============================================================================
# HAE多组学项目 - 免疫细胞/代谢物 → 肝纤维化的因果推断
# =============================================================================

cat("
╔═══════════════════════════════════════════════════════════════════════════════╗
║      HAE Multi-Omics: Mendelian Randomization Analysis                        ║
║      免疫细胞/代谢物 → 肝纤维化 因果推断                                          ║
╚═══════════════════════════════════════════════════════════════════════════════╝
\n")

# =============================================================================
# 0. 环境配置
# =============================================================================
cat("=== Step 0: 环境配置 ===\n")

# 设置工作目录
project_root <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
setwd(project_root)

# 创建输出目录
output_dir <- file.path(project_root, "analysis/results/enhancement_mr")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("  输出目录: %s\n", output_dir))

# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

# 加载NC制图主题
source("analysis/scripts/nc_theme.R")

# 检查并加载MR相关包
mr_packages <- list(
  ieugwasr = FALSE,
  MendelianRandomization = FALSE,
  TwoSampleMR = FALSE,
  RadialMR = FALSE,
  coloc = FALSE
)

for (pkg in names(mr_packages)) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    mr_packages[[pkg]] <- TRUE
    cat(sprintf("  [OK] %s 已加载\n", pkg))
  } else {
    cat(sprintf("  [--] %s 未安装\n", pkg))
  }
}

# 确定使用哪个MR分析包
USE_TWOSAMPLEMR <- mr_packages$TwoSampleMR
USE_MR_PKG <- mr_packages$MendelianRandomization
cat(sprintf("\n  MR分析方案: %s\n", 
            ifelse(USE_TWOSAMPLEMR, "TwoSampleMR", 
                   ifelse(USE_MR_PKG, "MendelianRandomization", "手动实现"))))

# =============================================================================
# 1. 配置参数
# =============================================================================
cat("\n=== Step 1: 配置参数 ===\n")

# MR分析参数
MR_CONFIG <- list(
  # 工具变量选择标准
  pval_threshold = 5e-8,      # genome-wide significance
  clump_r2 = 0.001,           # LD clumping r²阈值
  clump_kb = 10000,           # LD clumping窗口 (kb)
  f_stat_threshold = 10,      # F统计量阈值
  
  # 可视化参数
  dpi = 300,
  width_mm = 183,
  height_mm = 150,
  
  # API配置（如果使用IEU OpenGWAS）
  use_api = TRUE,  # 尝试使用API
  api_timeout = 300
)

cat(sprintf("  p值阈值: %s\n", format(MR_CONFIG$pval_threshold, scientific = TRUE)))
cat(sprintf("  LD clumping: r² < %.3f, window = %d kb\n", MR_CONFIG$clump_r2, MR_CONFIG$clump_kb))
cat(sprintf("  F统计量阈值: > %d\n", MR_CONFIG$f_stat_threshold))

# =============================================================================
# 2. 辅助函数定义
# =============================================================================
cat("\n=== Step 2: 定义辅助函数 ===\n")

# --- 2.1 计算F统计量 ---
calc_f_stat <- function(beta, se, n = NULL, maf = NULL) {
  # F = beta²/se²
  # 或者 F = (n-2) * r² / (1 - r²)，其中 r² = 2*beta²*maf*(1-maf)
  f_stat <- (beta/se)^2
  return(f_stat)
}

# --- 2.2 计算OR和95%CI ---
calc_or_ci <- function(beta, se) {
  or <- exp(beta)
  or_lci <- exp(beta - 1.96 * se)
  or_uci <- exp(beta + 1.96 * se)
  return(data.frame(OR = or, OR_LCI = or_lci, OR_UCI = or_uci))
}

# --- 2.3 IVW方法（固定效应） ---
mr_ivw <- function(beta_exp, se_exp, beta_out, se_out) {
  # 权重
  w <- 1 / se_out^2
  
  # beta估计 (Wald ratio的加权平均)
  beta_iv <- beta_out / beta_exp
  
  # IVW估计
  beta_ivw <- sum(w * beta_iv) / sum(w)
  se_ivw <- sqrt(1 / sum(w * beta_exp^2 / se_out^2))
  
  # p值
  z <- beta_ivw / se_ivw
  pval <- 2 * pnorm(-abs(z))
  
  # 异质性检验 (Cochran's Q)
  q_stat <- sum(w * (beta_iv - beta_ivw)^2)
  q_df <- length(beta_exp) - 1
  q_pval <- pchisq(q_stat, q_df, lower.tail = FALSE)
  
  return(list(
    method = "IVW",
    beta = beta_ivw,
    se = se_ivw,
    pval = pval,
    nsnp = length(beta_exp),
    Q = q_stat,
    Q_df = q_df,
    Q_pval = q_pval
  ))
}

# --- 2.4 MR-Egger回归 ---
mr_egger <- function(beta_exp, se_exp, beta_out, se_out) {
  # 对数据进行方向一致性处理
  sign_exp <- sign(beta_exp)
  beta_exp_abs <- abs(beta_exp)
  beta_out_adj <- beta_out * sign_exp
  
  # 加权回归
  w <- 1 / se_out^2
  
  # 回归: beta_out ~ intercept + slope * beta_exp
  X <- cbind(1, beta_exp_abs)
  W <- diag(w)
  
  # 加权最小二乘
  XtWX <- t(X) %*% W %*% X
  XtWX_inv <- solve(XtWX)
  XtWy <- t(X) %*% W %*% beta_out_adj
  coef <- as.vector(XtWX_inv %*% XtWy)
  
  # 残差
  fitted <- X %*% coef
  resid <- beta_out_adj - fitted
  
  # 标准误（使用残差方差）
  n <- length(beta_exp)
  sigma2 <- sum(w * resid^2) / (n - 2)
  se_coef <- sqrt(diag(XtWX_inv) * sigma2)
  
  # 截距（多效性检验）
  intercept <- coef[1]
  se_intercept <- se_coef[1]
  pval_intercept <- 2 * pnorm(-abs(intercept / se_intercept))
  
  # 斜率（因果效应）
  slope <- coef[2]
  se_slope <- se_coef[2]
  pval_slope <- 2 * pnorm(-abs(slope / se_slope))
  
  return(list(
    method = "MR-Egger",
    beta = slope,
    se = se_slope,
    pval = pval_slope,
    nsnp = length(beta_exp),
    intercept = intercept,
    se_intercept = se_intercept,
    pval_intercept = pval_intercept
  ))
}

# --- 2.5 加权中位数方法 ---
mr_weighted_median <- function(beta_exp, se_exp, beta_out, se_out, nboot = 1000) {
  # Wald ratio估计
  beta_iv <- beta_out / beta_exp
  se_iv <- sqrt((se_out^2 / beta_exp^2) + (beta_out^2 * se_exp^2 / beta_exp^4))
  
  # 权重
  w <- 1 / se_iv^2
  w_norm <- w / sum(w)
  
  # 加权中位数
  weighted_median <- function(b, weights) {
    ord <- order(b)
    b_sorted <- b[ord]
    w_sorted <- weights[ord]
    cum_w <- cumsum(w_sorted)
    idx <- min(which(cum_w >= 0.5))
    return(b_sorted[idx])
  }
  
  beta_wm <- weighted_median(beta_iv, w_norm)
  
  # Bootstrap估计标准误
  beta_boot <- numeric(nboot)
  for (i in 1:nboot) {
    beta_iv_boot <- rnorm(length(beta_iv), beta_iv, se_iv)
    beta_boot[i] <- weighted_median(beta_iv_boot, w_norm)
  }
  se_wm <- sd(beta_boot)
  
  # p值
  pval <- 2 * pnorm(-abs(beta_wm / se_wm))
  
  return(list(
    method = "Weighted median",
    beta = beta_wm,
    se = se_wm,
    pval = pval,
    nsnp = length(beta_exp)
  ))
}

# --- 2.6 加权众数方法 (简化版) ---
mr_weighted_mode <- function(beta_exp, se_exp, beta_out, se_out, bandwidth = 0.5) {
  # Wald ratio估计
  beta_iv <- beta_out / beta_exp
  se_iv <- sqrt((se_out^2 / beta_exp^2) + (beta_out^2 * se_exp^2 / beta_exp^4))
  
  # 权重
  w <- 1 / se_iv^2
  
  # 核密度估计找众数
  dens <- density(beta_iv, weights = w / sum(w), bw = bandwidth, n = 1024)
  mode_idx <- which.max(dens$y)
  beta_mode <- dens$x[mode_idx]
  
  # 使用bootstrap估计标准误
  nboot <- 1000
  beta_boot <- numeric(nboot)
  for (i in 1:nboot) {
    beta_iv_boot <- rnorm(length(beta_iv), beta_iv, se_iv)
    dens_boot <- density(beta_iv_boot, weights = w / sum(w), bw = bandwidth, n = 512)
    beta_boot[i] <- dens_boot$x[which.max(dens_boot$y)]
  }
  se_mode <- sd(beta_boot)
  
  # p值
  pval <- 2 * pnorm(-abs(beta_mode / se_mode))
  
  return(list(
    method = "Weighted mode",
    beta = beta_mode,
    se = se_mode,
    pval = pval,
    nsnp = length(beta_exp)
  ))
}

# --- 2.7 Leave-one-out分析 ---
mr_loo <- function(beta_exp, se_exp, beta_out, se_out, snp_names = NULL) {
  n <- length(beta_exp)
  if (is.null(snp_names)) snp_names <- paste0("SNP", 1:n)
  
  loo_results <- data.frame(
    SNP = snp_names,
    beta = numeric(n),
    se = numeric(n),
    pval = numeric(n)
  )
  
  for (i in 1:n) {
    # 排除第i个SNP
    idx <- setdiff(1:n, i)
    res <- mr_ivw(beta_exp[idx], se_exp[idx], beta_out[idx], se_out[idx])
    loo_results$beta[i] <- res$beta
    loo_results$se[i] <- res$se
    loo_results$pval[i] <- res$pval
  }
  
  # 添加全部SNP的结果
  res_all <- mr_ivw(beta_exp, se_exp, beta_out, se_out)
  loo_results <- rbind(
    loo_results,
    data.frame(
      SNP = "All SNPs",
      beta = res_all$beta,
      se = res_all$se,
      pval = res_all$pval
    )
  )
  
  return(loo_results)
}

# --- 2.8 综合MR分析 ---
run_mr_analysis <- function(exposure_data, outcome_data, exposure_name, outcome_name) {
  cat(sprintf("\n  >>> 分析: %s → %s\n", exposure_name, outcome_name))
  
  # 合并数据（通过SNP匹配）
  dat <- merge(exposure_data, outcome_data, by = "SNP", suffixes = c("_exp", "_out"))
  
  if (nrow(dat) < 3) {
    cat(sprintf("    [跳过] SNP数量不足: %d < 3\n", nrow(dat)))
    return(NULL)
  }
  
  # 统一效应方向（确保暴露效应为正）
  flip <- dat$beta_exp < 0
  dat$beta_exp[flip] <- -dat$beta_exp[flip]
  dat$beta_out[flip] <- -dat$beta_out[flip]
  
  # 计算F统计量并过滤弱工具变量
  dat$f_stat <- calc_f_stat(dat$beta_exp, dat$se_exp)
  dat <- dat[dat$f_stat >= MR_CONFIG$f_stat_threshold, ]
  
  if (nrow(dat) < 3) {
    cat(sprintf("    [跳过] F统计量过滤后SNP不足: %d < 3\n", nrow(dat)))
    return(NULL)
  }
  
  cat(sprintf("    有效SNP数: %d (平均F = %.1f)\n", nrow(dat), mean(dat$f_stat)))
  
  # 运行各种MR方法
  results <- list()
  
  # IVW
  results$ivw <- tryCatch(
    mr_ivw(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
    error = function(e) NULL
  )
  
  # MR-Egger
  results$egger <- tryCatch(
    mr_egger(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
    error = function(e) NULL
  )
  
  # 加权中位数
  results$wmedian <- tryCatch(
    mr_weighted_median(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
    error = function(e) NULL
  )
  
  # 加权众数
  results$wmode <- tryCatch(
    mr_weighted_mode(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
    error = function(e) NULL
  )
  
  # Leave-one-out
  results$loo <- tryCatch(
    mr_loo(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out, dat$SNP),
    error = function(e) NULL
  )
  
  # 整理结果
  mr_results <- do.call(rbind, lapply(names(results)[1:4], function(m) {
    res <- results[[m]]
    if (is.null(res)) return(NULL)
    
    or_ci <- calc_or_ci(res$beta, res$se)
    
    data.frame(
      exposure = exposure_name,
      outcome = outcome_name,
      method = res$method,
      nsnp = res$nsnp,
      beta = res$beta,
      se = res$se,
      pval = res$pval,
      OR = or_ci$OR,
      OR_LCI = or_ci$OR_LCI,
      OR_UCI = or_ci$OR_UCI,
      Q = ifelse(!is.null(res$Q), res$Q, NA),
      Q_pval = ifelse(!is.null(res$Q_pval), res$Q_pval, NA),
      egger_intercept = ifelse(m == "egger", res$intercept, NA),
      egger_pval = ifelse(m == "egger", res$pval_intercept, NA),
      stringsAsFactors = FALSE
    )
  }))
  
  if (!is.null(results$ivw)) {
    cat(sprintf("    IVW: beta = %.3f (95%%CI: %.3f-%.3f), p = %.2e\n",
                results$ivw$beta, 
                results$ivw$beta - 1.96*results$ivw$se,
                results$ivw$beta + 1.96*results$ivw$se,
                results$ivw$pval))
  }
  
  return(list(
    results = mr_results,
    loo = results$loo,
    harmonized_data = dat
  ))
}

cat("  辅助函数定义完成\n")

# =============================================================================
# 3. GWAS数据获取函数
# =============================================================================
cat("\n=== Step 3: 定义GWAS数据获取函数 ===\n")

# --- 3.1 从IEU OpenGWAS获取数据 ---
get_gwas_data_api <- function(gwas_id, pval_threshold = 5e-8) {
  if (!mr_packages$ieugwasr) {
    cat("    [错误] ieugwasr包未安装\n")
    return(NULL)
  }
  
  cat(sprintf("    从IEU API获取: %s\n", gwas_id))
  
  tryCatch({
    # 获取显著SNP
    tophits <- ieugwasr::tophits(gwas_id, pval = pval_threshold)
    
    if (is.null(tophits) || nrow(tophits) == 0) {
      cat(sprintf("    [警告] %s 未获取到显著SNP\n", gwas_id))
      return(NULL)
    }
    
    # 格式化
    dat <- data.frame(
      SNP = tophits$rsid,
      chr = tophits$chr,
      pos = tophits$position,
      effect_allele = tophits$ea,
      other_allele = tophits$nea,
      beta = tophits$beta,
      se = tophits$se,
      pval = tophits$p,
      eaf = tophits$eaf,
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("    获取到 %d 个显著SNP\n", nrow(dat)))
    return(dat)
    
  }, error = function(e) {
    cat(sprintf("    [错误] API获取失败: %s\n", e$message))
    return(NULL)
  })
}

# --- 3.2 LD clumping ---
perform_ld_clumping <- function(dat, r2 = 0.001, kb = 10000) {
  if (!mr_packages$ieugwasr) {
    cat("    [跳过] LD clumping (ieugwasr不可用)\n")
    return(dat)
  }
  
  cat(sprintf("    LD clumping: r² < %.3f, window = %d kb\n", r2, kb))
  
  tryCatch({
    clumped <- ieugwasr::ld_clump(
      dat = data.frame(rsid = dat$SNP, pval = dat$pval),
      clump_r2 = r2,
      clump_kb = kb,
      plink_bin = NULL  # 使用在线服务
    )
    
    dat_clumped <- dat[dat$SNP %in% clumped$rsid, ]
    cat(sprintf("    Clumping后SNP数: %d → %d\n", nrow(dat), nrow(dat_clumped)))
    return(dat_clumped)
    
  }, error = function(e) {
    cat(sprintf("    [警告] LD clumping失败: %s\n", e$message))
    return(dat)
  })
}

# --- 3.3 从本地文件读取GWAS数据 ---
read_gwas_local <- function(file_path, snp_col = "SNP", beta_col = "beta", 
                            se_col = "se", pval_col = "pval",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele") {
  
  cat(sprintf("    从本地文件读取: %s\n", basename(file_path)))
  
  if (!file.exists(file_path)) {
    cat("    [错误] 文件不存在\n")
    return(NULL)
  }
  
  # 读取文件
  if (grepl("\\.gz$", file_path)) {
    dat <- data.table::fread(cmd = paste("gzip -dc", file_path))
  } else if (grepl("\\.csv$", file_path)) {
    dat <- data.table::fread(file_path)
  } else {
    dat <- data.table::fread(file_path)
  }
  
  # 标准化列名
  dat <- as.data.frame(dat)
  
  # 检查必需列
  required_cols <- c(snp_col, beta_col, se_col, pval_col)
  missing <- setdiff(required_cols, names(dat))
  
  if (length(missing) > 0) {
    cat(sprintf("    [错误] 缺少必需列: %s\n", paste(missing, collapse = ", ")))
    return(NULL)
  }
  
  # 格式化
  result <- data.frame(
    SNP = dat[[snp_col]],
    beta = as.numeric(dat[[beta_col]]),
    se = as.numeric(dat[[se_col]]),
    pval = as.numeric(dat[[pval_col]]),
    stringsAsFactors = FALSE
  )
  
  # 添加可选列
  if (effect_allele_col %in% names(dat)) {
    result$effect_allele <- dat[[effect_allele_col]]
  }
  if (other_allele_col %in% names(dat)) {
    result$other_allele <- dat[[other_allele_col]]
  }
  
  cat(sprintf("    读取到 %d 个SNP\n", nrow(result)))
  return(result)
}

cat("  GWAS数据获取函数定义完成\n")

# =============================================================================
# 4. 可视化函数
# =============================================================================
cat("\n=== Step 4: 定义可视化函数 ===\n")

# --- 4.1 森林图 ---
plot_mr_forest <- function(results, title = "Mendelian Randomization Forest Plot") {
  
  # 创建标签
  results$label <- paste0(results$exposure, " → ", results$outcome)
  
  # 只绘制IVW结果
  plot_data <- results %>%
    filter(method == "IVW") %>%
    mutate(
      significant = pval < 0.05,
      color = ifelse(significant, ifelse(beta > 0, COL_UP, COL_DOWN), COL_NS)
    )
  
  if (nrow(plot_data) == 0) {
    cat("    [警告] 无IVW结果可绘制\n")
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = OR, y = reorder(label, OR))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = OR_LCI, xmax = OR_UCI), height = 0.2, linewidth = 0.5) +
    geom_point(aes(color = color), size = 3) +
    scale_color_identity() +
    scale_x_log10() +
    labs(
      x = "Odds Ratio (95% CI)",
      y = NULL,
      title = title
    ) +
    theme_pub +
    theme(
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 10, face = "bold")
    )
  
  return(p)
}

# --- 4.2 散点图 ---
plot_mr_scatter <- function(harmonized_data, results, exposure_name, outcome_name) {
  
  dat <- harmonized_data
  
  p <- ggplot(dat, aes(x = beta_exp, y = beta_out)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_errorbar(aes(ymin = beta_out - 1.96*se_out, ymax = beta_out + 1.96*se_out),
                  width = 0, alpha = 0.3) +
    geom_errorbarh(aes(xmin = beta_exp - 1.96*se_exp, xmax = beta_exp + 1.96*se_exp),
                   height = 0, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
  
  # 添加各方法的回归线
  colors <- c("IVW" = COL_UP, "MR-Egger" = COL_DOWN, 
              "Weighted median" = PAL_CAT[3], "Weighted mode" = PAL_CAT[4])
  
  for (i in 1:nrow(results)) {
    method <- results$method[i]
    slope <- results$beta[i]
    if (method %in% names(colors)) {
      p <- p + geom_abline(intercept = 0, slope = slope, 
                          color = colors[method], linewidth = 0.8, alpha = 0.7)
    }
  }
  
  p <- p +
    labs(
      x = paste0("SNP effect on ", exposure_name, " (β)"),
      y = paste0("SNP effect on ", outcome_name, " (β)"),
      title = paste0(exposure_name, " → ", outcome_name)
    ) +
    theme_pub
  
  return(p)
}

# --- 4.3 漏斗图 ---
plot_mr_funnel <- function(harmonized_data, results_ivw, exposure_name, outcome_name) {
  
  dat <- harmonized_data
  
  # 计算Wald ratio
  dat$beta_iv <- dat$beta_out / dat$beta_exp
  dat$se_iv <- sqrt((dat$se_out^2 / dat$beta_exp^2) + 
                     (dat$beta_out^2 * dat$se_exp^2 / dat$beta_exp^4))
  dat$precision <- 1 / dat$se_iv
  
  ivw_beta <- results_ivw$beta
  
  p <- ggplot(dat, aes(x = beta_iv, y = precision)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_vline(xintercept = ivw_beta, linetype = "dashed", color = COL_UP, linewidth = 0.8) +
    labs(
      x = "Causal estimate (β)",
      y = "Precision (1/SE)",
      title = paste0(exposure_name, " → ", outcome_name, " (Funnel Plot)")
    ) +
    theme_pub
  
  return(p)
}

# --- 4.4 Leave-one-out图 ---
plot_mr_loo <- function(loo_data, exposure_name, outcome_name) {
  
  loo_data$lci <- loo_data$beta - 1.96 * loo_data$se
  loo_data$uci <- loo_data$beta + 1.96 * loo_data$se
  loo_data$is_all <- loo_data$SNP == "All SNPs"
  
  p <- ggplot(loo_data, aes(x = beta, y = factor(SNP, levels = rev(SNP)))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.2, linewidth = 0.5) +
    geom_point(aes(color = is_all), size = 2) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = COL_UP), guide = "none") +
    labs(
      x = "IVW causal estimate (β)",
      y = NULL,
      title = paste0(exposure_name, " → ", outcome_name, " (Leave-one-out)")
    ) +
    theme_pub +
    theme(axis.text.y = element_text(size = 6))
  
  return(p)
}

# --- 4.5 因果网络图 ---
plot_causal_network <- function(results, pval_threshold = 0.05) {
  
  sig_results <- results %>%
    filter(method == "IVW", pval < pval_threshold)
  
  if (nrow(sig_results) == 0) {
    cat("    [警告] 无显著因果关系可绘制\n")
    return(NULL)
  }
  
  # 简单网络图
  nodes <- unique(c(sig_results$exposure, sig_results$outcome))
  
  # 创建节点坐标
  n_nodes <- length(nodes)
  angles <- seq(0, 2*pi, length.out = n_nodes + 1)[1:n_nodes]
  node_df <- data.frame(
    name = nodes,
    x = cos(angles),
    y = sin(angles)
  )
  
  # 创建边
  edge_df <- sig_results %>%
    left_join(node_df, by = c("exposure" = "name")) %>%
    rename(x_start = x, y_start = y) %>%
    left_join(node_df, by = c("outcome" = "name")) %>%
    rename(x_end = x, y_end = y)
  
  p <- ggplot() +
    geom_segment(data = edge_df, 
                aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                arrow = arrow(length = unit(0.3, "cm")),
                linewidth = abs(edge_df$beta) * 2,
                color = ifelse(edge_df$beta > 0, COL_UP, COL_DOWN),
                alpha = 0.7) +
    geom_point(data = node_df, aes(x = x, y = y), size = 8, color = PAL_CAT[4]) +
    geom_text(data = node_df, aes(x = x, y = y, label = name), size = 3) +
    coord_fixed() +
    theme_void() +
    labs(title = "Causal Network")
  
  return(p)
}

cat("  可视化函数定义完成\n")

# =============================================================================
# 5. 主分析流程
# =============================================================================
cat("\n=== Step 5: 开始MR分析 ===\n")

# 存储所有结果
all_results <- list()
all_harmonized <- list()
all_loo <- list()

# --- 5.1 检测API可用性 ---
api_available <- FALSE
if (mr_packages$ieugwasr && MR_CONFIG$use_api) {
  cat("\n检测IEU OpenGWAS API可用性...\n")
  tryCatch({
    # 尝试简单查询
    test_info <- ieugwasr::gwasinfo("ieu-b-30")
    if (!is.null(test_info) && nrow(test_info) > 0) {
      api_available <- TRUE
      cat("  [OK] API可用\n")
    }
  }, error = function(e) {
    cat(sprintf("  [警告] API不可用: %s\n", e$message))
  })
}

# --- 5.2 定义分析方案 ---
cat("\n定义分析方案...\n")

# 暴露GWAS ID（血细胞计数 - 与免疫细胞相关）
exposure_ids <- list(
  # 血细胞计数 (UK Biobank)
  "Lymphocyte" = "ieu-b-30",
  "Monocyte" = "ieu-b-31",
  "Neutrophil" = "ieu-b-32",
  "Eosinophil" = "ieu-b-33",
  "Basophil" = "ieu-b-34"
)

# 结局GWAS ID
outcome_ids <- list(
  "Liver_cirrhosis" = "finn-b-CIRRHOSIS",
  "ALT" = "ieu-b-107",
  "AST" = "ieu-b-108",
  "GGT" = "ieu-b-109"
)

cat(sprintf("  暴露数量: %d\n", length(exposure_ids)))
cat(sprintf("  结局数量: %d\n", length(outcome_ids)))

# --- 5.3 执行分析 ---
if (api_available) {
  cat("\n开始API模式分析...\n")
  
  for (exp_name in names(exposure_ids)) {
    exp_id <- exposure_ids[[exp_name]]
    
    cat(sprintf("\n--- 处理暴露: %s (%s) ---\n", exp_name, exp_id))
    
    # 获取暴露数据
    exp_data <- get_gwas_data_api(exp_id, MR_CONFIG$pval_threshold)
    
    if (is.null(exp_data) || nrow(exp_data) == 0) {
      cat("  [跳过] 无法获取暴露数据\n")
      next
    }
    
    # LD clumping
    exp_data <- perform_ld_clumping(exp_data, MR_CONFIG$clump_r2, MR_CONFIG$clump_kb)
    
    for (out_name in names(outcome_ids)) {
      out_id <- outcome_ids[[out_name]]
      
      cat(sprintf("\n  结局: %s (%s)\n", out_name, out_id))
      
      # 获取结局数据（只获取暴露SNP的效应）
      tryCatch({
        out_data <- ieugwasr::associations(exp_data$SNP, out_id)
        
        if (is.null(out_data) || nrow(out_data) == 0) {
          cat("    [跳过] 无法获取结局数据\n")
          next
        }
        
        # 格式化结局数据
        out_data_fmt <- data.frame(
          SNP = out_data$rsid,
          beta = out_data$beta,
          se = out_data$se,
          pval = out_data$p,
          effect_allele = out_data$ea,
          other_allele = out_data$nea,
          stringsAsFactors = FALSE
        )
        
        # 运行MR分析
        analysis_result <- run_mr_analysis(exp_data, out_data_fmt, exp_name, out_name)
        
        if (!is.null(analysis_result)) {
          key <- paste0(exp_name, "_", out_name)
          all_results[[key]] <- analysis_result$results
          all_harmonized[[key]] <- analysis_result$harmonized_data
          all_loo[[key]] <- analysis_result$loo
        }
        
      }, error = function(e) {
        cat(sprintf("    [错误] %s\n", e$message))
      })
    }
  }
  
} else {
  cat("\n=== API不可用，使用模拟数据演示 ===\n")
  cat("
提示：要使用真实GWAS数据，您需要：
1. 设置IEU OpenGWAS API token
2. 或手动下载GWAS summary statistics

手动下载指南：
- IEU OpenGWAS: https://gwas.mrcieu.ac.uk/
- FinnGen: https://finngen.fi/en/access_results
- UK Biobank GWAS: https://pan.ukbb.broadinstitute.org/

下载后，使用 read_gwas_local() 函数读取本地文件。
\n")
  
  # 生成模拟数据进行演示
  cat("生成模拟数据进行方法演示...\n")
  
  set.seed(42)
  n_snps <- 50
  
  # 模拟暴露数据
  sim_exposure <- data.frame(
    SNP = paste0("rs", sample(1e6:9e6, n_snps)),
    beta = rnorm(n_snps, 0.1, 0.05),
    se = runif(n_snps, 0.01, 0.03),
    pval = runif(n_snps, 1e-20, 1e-8),
    effect_allele = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE),
    other_allele = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE)
  )
  
  # 模拟结局数据（假设真实因果效应 = 0.3）
  true_effect <- 0.3
  sim_outcome <- data.frame(
    SNP = sim_exposure$SNP,
    beta = sim_exposure$beta * true_effect + rnorm(n_snps, 0, 0.02),
    se = runif(n_snps, 0.02, 0.04),
    pval = NA,
    effect_allele = sim_exposure$effect_allele,
    other_allele = sim_exposure$other_allele
  )
  sim_outcome$pval <- 2 * pnorm(-abs(sim_outcome$beta / sim_outcome$se))
  
  # 运行分析
  cat("\n运行模拟数据MR分析...\n")
  analysis_result <- run_mr_analysis(sim_exposure, sim_outcome, 
                                     "Simulated_Immune_Cell", "Simulated_Liver_Trait")
  
  if (!is.null(analysis_result)) {
    all_results[["Simulated_Analysis"]] <- analysis_result$results
    all_harmonized[["Simulated_Analysis"]] <- analysis_result$harmonized_data
    all_loo[["Simulated_Analysis"]] <- analysis_result$loo
  }
}

# =============================================================================
# 6. 汇总结果
# =============================================================================
cat("\n=== Step 6: 汇总结果 ===\n")

if (length(all_results) > 0) {
  # 合并所有结果
  combined_results <- do.call(rbind, all_results)
  rownames(combined_results) <- NULL
  
  cat(sprintf("  总分析数: %d\n", length(all_results)))
  cat(sprintf("  总结果行数: %d\n", nrow(combined_results)))
  
  # 显示主要结果
  cat("\n=== 主要MR结果 (IVW方法) ===\n")
  ivw_results <- combined_results %>%
    filter(method == "IVW") %>%
    arrange(pval)
  
  print(ivw_results[, c("exposure", "outcome", "nsnp", "beta", "se", "pval", "OR", "OR_LCI", "OR_UCI")])
  
  # 保存结果
  cat("\n=== Step 7: 保存结果 ===\n")
  
  # 主结果表
  write.csv(combined_results, 
            file.path(output_dir, "mr_results_all_methods.csv"),
            row.names = FALSE)
  cat("  已保存: mr_results_all_methods.csv\n")
  
  # IVW结果
  write.csv(ivw_results, 
            file.path(output_dir, "mr_results_ivw.csv"),
            row.names = FALSE)
  cat("  已保存: mr_results_ivw.csv\n")
  
  # 敏感性分析汇总
  sensitivity <- combined_results %>%
    select(exposure, outcome, method, nsnp, beta, se, pval, Q, Q_pval, egger_intercept, egger_pval) %>%
    pivot_wider(
      id_cols = c(exposure, outcome, nsnp),
      names_from = method,
      values_from = c(beta, se, pval),
      names_glue = "{method}_{.value}"
    )
  write.csv(sensitivity, 
            file.path(output_dir, "mr_sensitivity_analysis.csv"),
            row.names = FALSE)
  cat("  已保存: mr_sensitivity_analysis.csv\n")
  
  # =============================================================================
  # 7. 生成可视化
  # =============================================================================
  cat("\n=== Step 8: 生成可视化 ===\n")
  
  # 森林图
  p_forest <- plot_mr_forest(combined_results, "MR Forest Plot: Immune Cells/Metabolites → Liver Traits")
  if (!is.null(p_forest)) {
    ggsave(file.path(output_dir, "mr_forest_plot.pdf"), p_forest,
           width = mm2in(183), height = mm2in(120), dpi = 300)
    cat("  已保存: mr_forest_plot.pdf\n")
  }
  
  # 散点图集
  scatter_plots <- list()
  for (key in names(all_harmonized)) {
    if (!is.null(all_harmonized[[key]]) && !is.null(all_results[[key]])) {
      parts <- strsplit(key, "_")[[1]]
      exp_name <- parts[1]
      out_name <- paste(parts[-1], collapse = "_")
      
      p <- plot_mr_scatter(all_harmonized[[key]], all_results[[key]], exp_name, out_name)
      scatter_plots[[key]] <- p
    }
  }
  
  if (length(scatter_plots) > 0) {
    pdf(file.path(output_dir, "mr_scatter_plots.pdf"), 
        width = mm2in(183), height = mm2in(150))
    for (p in scatter_plots) {
      print(p)
    }
    dev.off()
    cat("  已保存: mr_scatter_plots.pdf\n")
  }
  
  # 漏斗图集
  funnel_plots <- list()
  for (key in names(all_harmonized)) {
    if (!is.null(all_harmonized[[key]]) && !is.null(all_results[[key]])) {
      parts <- strsplit(key, "_")[[1]]
      exp_name <- parts[1]
      out_name <- paste(parts[-1], collapse = "_")
      
      ivw_res <- all_results[[key]] %>% filter(method == "IVW")
      if (nrow(ivw_res) > 0) {
        p <- plot_mr_funnel(all_harmonized[[key]], 
                           list(beta = ivw_res$beta[1]),
                           exp_name, out_name)
        funnel_plots[[key]] <- p
      }
    }
  }
  
  if (length(funnel_plots) > 0) {
    pdf(file.path(output_dir, "mr_funnel_plots.pdf"),
        width = mm2in(183), height = mm2in(150))
    for (p in funnel_plots) {
      print(p)
    }
    dev.off()
    cat("  已保存: mr_funnel_plots.pdf\n")
  }
  
  # Leave-one-out图集
  loo_plots <- list()
  for (key in names(all_loo)) {
    if (!is.null(all_loo[[key]])) {
      parts <- strsplit(key, "_")[[1]]
      exp_name <- parts[1]
      out_name <- paste(parts[-1], collapse = "_")
      
      p <- plot_mr_loo(all_loo[[key]], exp_name, out_name)
      loo_plots[[key]] <- p
    }
  }
  
  if (length(loo_plots) > 0) {
    pdf(file.path(output_dir, "mr_loo_plots.pdf"),
        width = mm2in(183), height = mm2in(200))
    for (p in loo_plots) {
      print(p)
    }
    dev.off()
    cat("  已保存: mr_loo_plots.pdf\n")
  }
  
  # 因果网络图
  p_network <- plot_causal_network(combined_results, 0.05)
  if (!is.null(p_network)) {
    ggsave(file.path(output_dir, "mr_causal_network.pdf"), p_network,
           width = mm2in(150), height = mm2in(150), dpi = 300)
    cat("  已保存: mr_causal_network.pdf\n")
  }
  
  # =============================================================================
  # 8. 生成摘要报告
  # =============================================================================
  cat("\n=== Step 9: 生成摘要报告 ===\n")
  
  report_lines <- c(
    "================================================================================",
    "HAE Multi-Omics: Mendelian Randomization Analysis Summary Report",
    paste0("Generated: ", Sys.time()),
    "================================================================================",
    "",
    "1. ANALYSIS OVERVIEW",
    "--------------------",
    sprintf("Total exposure-outcome pairs analyzed: %d", length(all_results)),
    sprintf("MR methods applied: IVW, MR-Egger, Weighted Median, Weighted Mode"),
    "",
    "2. SIGNIFICANT CAUSAL RELATIONSHIPS (IVW p < 0.05)",
    "---------------------------------------------------"
  )
  
  sig_ivw <- ivw_results %>% filter(pval < 0.05)
  if (nrow(sig_ivw) > 0) {
    for (i in 1:nrow(sig_ivw)) {
      report_lines <- c(report_lines,
                       sprintf("  %s → %s:", sig_ivw$exposure[i], sig_ivw$outcome[i]),
                       sprintf("    OR = %.2f (95%%CI: %.2f-%.2f), p = %.2e",
                              sig_ivw$OR[i], sig_ivw$OR_LCI[i], sig_ivw$OR_UCI[i], sig_ivw$pval[i]))
    }
  } else {
    report_lines <- c(report_lines, "  No significant causal relationships detected")
  }
  
  report_lines <- c(report_lines,
                   "",
                   "3. SENSITIVITY ANALYSIS SUMMARY",
                   "-------------------------------")
  
  for (key in names(all_results)) {
    res <- all_results[[key]]
    egger_res <- res %>% filter(method == "MR-Egger")
    ivw_res <- res %>% filter(method == "IVW")
    
    if (nrow(egger_res) > 0 && nrow(ivw_res) > 0) {
      report_lines <- c(report_lines,
                       sprintf("  %s:", key),
                       sprintf("    Heterogeneity (Cochran's Q): p = %.3f", ivw_res$Q_pval[1]),
                       sprintf("    Pleiotropy (Egger intercept): p = %.3f", egger_res$egger_pval[1]))
    }
  }
  
  report_lines <- c(report_lines,
                   "",
                   "4. OUTPUT FILES",
                   "---------------",
                   "  - mr_results_all_methods.csv: All MR results",
                   "  - mr_results_ivw.csv: IVW results only",
                   "  - mr_sensitivity_analysis.csv: Sensitivity analysis summary",
                   "  - mr_forest_plot.pdf: Forest plot",
                   "  - mr_scatter_plots.pdf: Scatter plots",
                   "  - mr_funnel_plots.pdf: Funnel plots",
                   "  - mr_loo_plots.pdf: Leave-one-out plots",
                   "  - mr_causal_network.pdf: Causal network visualization",
                   "",
                   "================================================================================",
                   "END OF REPORT",
                   "================================================================================")
  
  writeLines(report_lines, file.path(output_dir, "mr_summary_report.txt"))
  cat("  已保存: mr_summary_report.txt\n")
  
} else {
  cat("\n[警告] 无MR分析结果\n")
}

# =============================================================================
# 完成
# =============================================================================
cat("\n
╔═══════════════════════════════════════════════════════════════════════════════╗
║      MR分析完成！                                                             ║
╚═══════════════════════════════════════════════════════════════════════════════╝
\n")

cat(sprintf("输出目录: %s\n", output_dir))
cat("\n建议后续步骤:\n")
cat("1. 检查 mr_summary_report.txt 了解分析结果概要\n")
cat("2. 查看森林图识别显著因果关系\n")
cat("3. 检查敏感性分析结果验证结果稳健性\n")
cat("4. 如需使用真实GWAS数据，请配置IEU API token或下载本地数据\n")

# 打印会话信息
cat("\n=== Session Info ===\n")
print(sessionInfo())
