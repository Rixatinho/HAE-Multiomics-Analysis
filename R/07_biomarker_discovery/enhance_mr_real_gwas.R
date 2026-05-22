#!/usr/bin/env Rscript
# =============================================================================
# enhance_mr_real_gwas.R -- 真实GWAS数据孟德尔随机化分析
# =============================================================================
# HAE多组学项目 - 使用IEU OpenGWAS真实数据
# =============================================================================

cat("
╔═══════════════════════════════════════════════════════════════════════════════╗
║      HAE Multi-Omics: Real GWAS Mendelian Randomization Analysis              ║
║      使用IEU OpenGWAS API真实数据                                              ║
╚═══════════════════════════════════════════════════════════════════════════════╝
\n")

# =============================================================================
# 0. 环境配置与JWT Token设置
# =============================================================================
cat("=== Step 0: 环境配置 ===\n")

# 设置工作目录
project_root <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
setwd(project_root)

# 创建输出目录
output_dir <- file.path(project_root, "analysis/results/enhancement_mr")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("  输出目录: %s\n", output_dir))

# *** 配置IEU OpenGWAS JWT Token ***
JWT_TOKEN <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJyaXNoYXQucnV6aUB4am11LmVkdS5jbiIsImlhdCI6MTc3NDI1OTI5NCwiZXhwIjoxNzc1NDY4ODk0fQ.GU2kv_jH3jBhYPavXHpImroVYZSdwc8LgF17bpo4VcU7nPGeW9dBfqzQplmg_JM1XEdaV0FuTz4Fj6CtnKXfPVsQUzAjIqUHfwW8cqLU0bvud9xNVaRoxBYQqA5kZ262fWhhbttr5uYUuFKViA9CB00tDpbRjR2hp-ABrFCQ5bFMq6rf74mrvGoJFzfU25IsScLXEeTBuj_AvRBuv0nMJ7aSeaVQC2AkPcbgnwc99s8VN0UyiT_NH0I6JemFA9EP4BnCR5zhoL8fkdCPfhFmTo813P9qJNh02oHx6qmXq61fhIr9MMVLxqES03y3a7gyEvwdhAQMYNl-jmzBQJndgg"

# 设置环境变量（ieugwasr会自动读取）
Sys.setenv(OPENGWAS_JWT = JWT_TOKEN)
# 设置更长的超时时间
options(ieugwasr_api = "https://api.opengwas.io/api/")
options(ieugwasr_timeout = 600)
cat("  JWT Token已配置\n")

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
cat("\n检查MR相关包...\n")
pkg_status <- list()

# ieugwasr
if (requireNamespace("ieugwasr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ieugwasr))
  pkg_status$ieugwasr <- TRUE
  cat("  [OK] ieugwasr 已加载\n")
  # 设置JWT token
  tryCatch({
    ieugwasr::set_jwt(JWT_TOKEN)
    cat("  [OK] JWT Token已设置到ieugwasr\n")
  }, error = function(e) {
    cat(sprintf("  [警告] 设置JWT失败: %s\n", e$message))
  })
} else {
  pkg_status$ieugwasr <- FALSE
  cat("  [--] ieugwasr 未安装\n")
}

# MendelianRandomization
if (requireNamespace("MendelianRandomization", quietly = TRUE)) {
  suppressPackageStartupMessages(library(MendelianRandomization))
  pkg_status$MendelianRandomization <- TRUE
  cat("  [OK] MendelianRandomization 已加载\n")
} else {
  pkg_status$MendelianRandomization <- FALSE
  cat("  [--] MendelianRandomization 未安装\n")
}

# TwoSampleMR - 尝试安装
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  cat("\n尝试安装TwoSampleMR...\n")
  tryCatch({
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never", quiet = TRUE)
    cat("  [OK] TwoSampleMR 安装成功\n")
  }, error = function(e) {
    cat(sprintf("  [警告] TwoSampleMR 安装失败: %s\n", e$message))
  })
}

if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(TwoSampleMR))
  pkg_status$TwoSampleMR <- TRUE
  cat("  [OK] TwoSampleMR 已加载\n")
} else {
  pkg_status$TwoSampleMR <- FALSE
  cat("  [--] TwoSampleMR 不可用，将使用手动实现\n")
}

# coloc
if (requireNamespace("coloc", quietly = TRUE)) {
  suppressPackageStartupMessages(library(coloc))
  pkg_status$coloc <- TRUE
  cat("  [OK] coloc 已加载\n")
} else {
  pkg_status$coloc <- FALSE
  cat("  [--] coloc 未安装\n")
}

# =============================================================================
# 1. 测试API连通性
# =============================================================================
cat("\n=== Step 1: 测试API连通性 ===\n")

api_available <- FALSE
tryCatch({
  status <- ieugwasr::api_status()
  cat(sprintf("  API状态: %s\n", status$status %||% "connected"))
  api_available <- TRUE
  cat("  [OK] API连接成功!\n")
}, error = function(e) {
  cat(sprintf("  [错误] API连接失败: %s\n", e$message))
})

if (!api_available) {
  cat("\n尝试备用方法验证API...\n")
  tryCatch({
    test_gwas <- ieugwasr::gwasinfo("ieu-b-30")
    if (!is.null(test_gwas) && nrow(test_gwas) > 0) {
      api_available <- TRUE
      cat(sprintf("  [OK] API验证成功! 测试GWAS: %s\n", test_gwas$trait[1]))
    }
  }, error = function(e) {
    cat(sprintf("  [错误] 备用验证失败: %s\n", e$message))
  })
}

if (!api_available) {
  stop("无法连接IEU OpenGWAS API，请检查网络和JWT Token")
}

# =============================================================================
# 2. 配置参数
# =============================================================================
cat("\n=== Step 2: 配置参数 ===\n")

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

# =============================================================================
# 3. 辅助函数定义
# =============================================================================
cat("\n=== Step 3: 定义辅助函数 ===\n")

# F统计量计算
calc_f_stat <- function(beta, se) (beta/se)^2

# OR和CI计算
calc_or_ci <- function(beta, se) {
  data.frame(OR = exp(beta), OR_LCI = exp(beta - 1.96*se), OR_UCI = exp(beta + 1.96*se))
}

# IVW方法
mr_ivw <- function(beta_exp, se_exp, beta_out, se_out) {
  w <- 1 / se_out^2
  beta_iv <- beta_out / beta_exp
  beta_ivw <- sum(w * beta_iv) / sum(w)
  se_ivw <- sqrt(1 / sum(w * beta_exp^2 / se_out^2))
  z <- beta_ivw / se_ivw
  pval <- 2 * pnorm(-abs(z))
  q_stat <- sum(w * (beta_iv - beta_ivw)^2)
  q_df <- length(beta_exp) - 1
  q_pval <- pchisq(q_stat, q_df, lower.tail = FALSE)
  list(method = "IVW", beta = beta_ivw, se = se_ivw, pval = pval,
       nsnp = length(beta_exp), Q = q_stat, Q_df = q_df, Q_pval = q_pval)
}

# MR-Egger回归
mr_egger <- function(beta_exp, se_exp, beta_out, se_out) {
  sign_exp <- sign(beta_exp)
  beta_exp_abs <- abs(beta_exp)
  beta_out_adj <- beta_out * sign_exp
  w <- 1 / se_out^2
  X <- cbind(1, beta_exp_abs)
  W <- diag(w)
  XtWX <- t(X) %*% W %*% X
  XtWX_inv <- solve(XtWX)
  XtWy <- t(X) %*% W %*% beta_out_adj
  coef <- as.vector(XtWX_inv %*% XtWy)
  fitted <- X %*% coef
  resid <- beta_out_adj - fitted
  n <- length(beta_exp)
  sigma2 <- sum(w * resid^2) / (n - 2)
  se_coef <- sqrt(diag(XtWX_inv) * sigma2)
  list(method = "MR-Egger", beta = coef[2], se = se_coef[2],
       pval = 2 * pnorm(-abs(coef[2] / se_coef[2])), nsnp = n,
       intercept = coef[1], se_intercept = se_coef[1],
       pval_intercept = 2 * pnorm(-abs(coef[1] / se_coef[1])))
}

# 加权中位数方法
mr_weighted_median <- function(beta_exp, se_exp, beta_out, se_out, nboot = 1000) {
  beta_iv <- beta_out / beta_exp
  se_iv <- sqrt((se_out^2 / beta_exp^2) + (beta_out^2 * se_exp^2 / beta_exp^4))
  w <- 1 / se_iv^2
  w_norm <- w / sum(w)
  weighted_median <- function(b, weights) {
    ord <- order(b)
    b_sorted <- b[ord]
    w_sorted <- weights[ord]
    cum_w <- cumsum(w_sorted)
    idx <- min(which(cum_w >= 0.5))
    b_sorted[idx]
  }
  beta_wm <- weighted_median(beta_iv, w_norm)
  beta_boot <- sapply(1:nboot, function(i) {
    beta_iv_boot <- rnorm(length(beta_iv), beta_iv, se_iv)
    weighted_median(beta_iv_boot, w_norm)
  })
  se_wm <- sd(beta_boot)
  list(method = "Weighted median", beta = beta_wm, se = se_wm,
       pval = 2 * pnorm(-abs(beta_wm / se_wm)), nsnp = length(beta_exp))
}

# 加权众数方法
mr_weighted_mode <- function(beta_exp, se_exp, beta_out, se_out, bandwidth = 0.5) {
  beta_iv <- beta_out / beta_exp
  se_iv <- sqrt((se_out^2 / beta_exp^2) + (beta_out^2 * se_exp^2 / beta_exp^4))
  w <- 1 / se_iv^2
  dens <- density(beta_iv, weights = w / sum(w), bw = bandwidth, n = 1024)
  beta_mode <- dens$x[which.max(dens$y)]
  beta_boot <- sapply(1:1000, function(i) {
    beta_iv_boot <- rnorm(length(beta_iv), beta_iv, se_iv)
    dens_boot <- density(beta_iv_boot, weights = w / sum(w), bw = bandwidth, n = 512)
    dens_boot$x[which.max(dens_boot$y)]
  })
  se_mode <- sd(beta_boot)
  list(method = "Weighted mode", beta = beta_mode, se = se_mode,
       pval = 2 * pnorm(-abs(beta_mode / se_mode)), nsnp = length(beta_exp))
}

# Leave-one-out分析
mr_loo <- function(beta_exp, se_exp, beta_out, se_out, snp_names = NULL) {
  n <- length(beta_exp)
  if (is.null(snp_names)) snp_names <- paste0("SNP", 1:n)
  loo_results <- do.call(rbind, lapply(1:n, function(i) {
    idx <- setdiff(1:n, i)
    res <- mr_ivw(beta_exp[idx], se_exp[idx], beta_out[idx], se_out[idx])
    data.frame(SNP = snp_names[i], beta = res$beta, se = res$se, pval = res$pval)
  }))
  res_all <- mr_ivw(beta_exp, se_exp, beta_out, se_out)
  rbind(loo_results, data.frame(SNP = "All SNPs", beta = res_all$beta,
                                 se = res_all$se, pval = res_all$pval))
}

cat("  辅助函数定义完成\n")

# =============================================================================
# 4. 搜索可用GWAS数据 (跳过大规模搜索，直接使用已知ID)
# =============================================================================
cat("\n=== Step 4: 使用预定义的可靠GWAS ID ===\n")
cat("  (跳过gwasinfo大规模搜索以避免超时)\n")

# =============================================================================
# 5. 定义分析方案 (精简版，使用UK Biobank可靠数据)
# =============================================================================
cat("\n=== Step 5: 定义分析方案 ===\n")

# 使用UK Biobank和经过验证的GWAS ID
# 方案1: 免疫细胞 → 肝酶 (UK Biobank)
scheme1 <- list(
  name = "Scheme1_BloodCells_to_LiverEnzymes",
  exposures = list(
    "Lymphocyte_pct" = "ukb-d-30180_irnt",  # Lymphocyte percentage
    "Monocyte_pct" = "ukb-d-30190_irnt"     # Monocyte percentage
  ),
  outcomes = list(
    "ALT" = "ukb-d-30620_irnt",  # ALT (UK Biobank)
    "GGT" = "ukb-d-30730_irnt"   # GGT (UK Biobank)
  )
)

# 方案2: 代谢物/血脂 → 肝酶
scheme2 <- list(
  name = "Scheme2_Lipids_to_LiverEnzymes",
  exposures = list(
    "LDL" = "ukb-d-30780_irnt",        # LDL cholesterol
    "Triglycerides" = "ukb-d-30870_irnt"  # Triglycerides
  ),
  outcomes = list(
    "ALT" = "ukb-d-30620_irnt",
    "AST" = "ukb-d-30650_irnt"   # AST (UK Biobank)
  )
)

# 方案3: 炎症标志物 → 肝酶
scheme3 <- list(
  name = "Scheme3_Inflammation_to_LiverEnzymes",
  exposures = list(
    "CRP" = "ukb-d-30710_irnt"  # C-reactive protein
  ),
  outcomes = list(
    "ALT" = "ukb-d-30620_irnt",
    "GGT" = "ukb-d-30730_irnt"
  )
)

# 方案4: 反向MR (肝酶 → 血细胞)
scheme4 <- list(
  name = "Scheme4_Reverse_LiverEnzymes_to_BloodCells",
  exposures = list(
    "ALT" = "ukb-d-30620_irnt"
  ),
  outcomes = list(
    "Lymphocyte_pct" = "ukb-d-30180_irnt",
    "Monocyte_pct" = "ukb-d-30190_irnt"
  )
)

all_schemes <- list(scheme1, scheme2, scheme3, scheme4)

cat(sprintf("  定义了 %d 个分析方案\n", length(all_schemes)))
for (s in all_schemes) {
  cat(sprintf("    - %s: %d exposures × %d outcomes\n", 
              s$name, length(s$exposures), length(s$outcomes)))
}

# =============================================================================
# 6. 执行MR分析
# =============================================================================
cat("\n=== Step 6: 执行MR分析 ===\n")

# 存储结果
all_results <- list()
all_harmonized <- list()
all_loo <- list()
analysis_log <- list()

# 执行单个MR分析的函数
run_single_mr <- function(exp_id, out_id, exp_name, out_name, max_snps = 100) {
  cat(sprintf("\n  >>> %s → %s\n", exp_name, out_name))
  cat(sprintf("      Exposure: %s, Outcome: %s\n", exp_id, out_id))
  
  result <- list(success = FALSE, error = NULL, data = NULL)
  
  tryCatch({
    # 1. 获取暴露数据
    cat("      获取暴露工具变量...\n")
    exp_data <- ieugwasr::tophits(exp_id, pval = MR_CONFIG$pval_threshold)
    
    if (is.null(exp_data) || nrow(exp_data) == 0) {
      result$error <- "无显著暴露SNP"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    cat(sprintf("      获取到 %d 个显著SNP\n", nrow(exp_data)))
    
    # 2. LD clumping
    cat("      执行LD clumping...\n")
    clumped <- tryCatch({
      ieugwasr::ld_clump(
        dat = data.frame(rsid = exp_data$rsid, pval = exp_data$p),
        clump_r2 = MR_CONFIG$clump_r2,
        clump_kb = MR_CONFIG$clump_kb
      )
    }, error = function(e) {
      cat(sprintf("      [警告] Clumping失败，使用原始数据: %s\n", e$message))
      data.frame(rsid = exp_data$rsid)
    })
    
    exp_data <- exp_data[exp_data$rsid %in% clumped$rsid, ]
    cat(sprintf("      Clumping后: %d SNPs\n", nrow(exp_data)))
    
    if (nrow(exp_data) < 3) {
      result$error <- "Clumping后SNP不足"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    
    # 限制SNP数量以避免API超时
    if (nrow(exp_data) > max_snps) {
      cat(sprintf("      限制SNP数量: %d → %d (按p值排序)\n", nrow(exp_data), max_snps))
      exp_data <- exp_data[order(exp_data$p), ][1:max_snps, ]
    }
    
    # 3. 获取结局数据 (分批查询)
    cat("      获取结局数据...\n")
    snp_batches <- split(exp_data$rsid, ceiling(seq_along(exp_data$rsid) / 50))
    out_data_list <- list()
    
    for (i in seq_along(snp_batches)) {
      cat(sprintf("        批次 %d/%d...\n", i, length(snp_batches)))
      Sys.sleep(2)  # API限速
      batch_result <- tryCatch({
        ieugwasr::associations(snp_batches[[i]], out_id)
      }, error = function(e) {
        cat(sprintf("        [警告] 批次%d失败: %s\n", i, e$message))
        NULL
      })
      if (!is.null(batch_result) && nrow(batch_result) > 0) {
        out_data_list[[i]] <- batch_result
      }
    }
    
    if (length(out_data_list) == 0) {
      result$error <- "无法获取结局数据"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    
    out_data <- do.call(rbind, out_data_list)
    cat(sprintf("      获取到 %d 个结局SNP效应\n", nrow(out_data)))
    
    # 4. 数据协调
    exp_fmt <- data.frame(
      SNP = exp_data$rsid, beta = exp_data$beta, se = exp_data$se,
      pval = exp_data$p, effect_allele = exp_data$ea, other_allele = exp_data$nea
    )
    out_fmt <- data.frame(
      SNP = out_data$rsid, beta = out_data$beta, se = out_data$se,
      pval = out_data$p, effect_allele = out_data$ea, other_allele = out_data$nea
    )
    
    dat <- merge(exp_fmt, out_fmt, by = "SNP", suffixes = c("_exp", "_out"))
    
    if (nrow(dat) < 3) {
      result$error <- "匹配后SNP不足"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    
    # 效应方向统一
    flip <- dat$beta_exp < 0
    dat$beta_exp[flip] <- -dat$beta_exp[flip]
    dat$beta_out[flip] <- -dat$beta_out[flip]
    
    # F统计量过滤
    dat$f_stat <- calc_f_stat(dat$beta_exp, dat$se_exp)
    dat <- dat[dat$f_stat >= MR_CONFIG$f_stat_threshold, ]
    
    if (nrow(dat) < 3) {
      result$error <- "F统计量过滤后SNP不足"
      cat(sprintf("      [跳过] %s\n", result$error))
      return(result)
    }
    
    cat(sprintf("      有效SNP: %d (平均F = %.1f)\n", nrow(dat), mean(dat$f_stat)))
    
    # 5. 执行MR分析
    mr_results <- list()
    
    # IVW
    mr_results$ivw <- tryCatch(mr_ivw(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
                                error = function(e) NULL)
    # MR-Egger
    mr_results$egger <- tryCatch(mr_egger(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
                                  error = function(e) NULL)
    # 加权中位数
    mr_results$wmedian <- tryCatch(mr_weighted_median(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
                                    error = function(e) NULL)
    # 加权众数
    mr_results$wmode <- tryCatch(mr_weighted_mode(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out),
                                  error = function(e) NULL)
    # LOO
    loo <- tryCatch(mr_loo(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out, dat$SNP),
                    error = function(e) NULL)
    
    # 整理结果
    results_df <- do.call(rbind, lapply(names(mr_results), function(m) {
      res <- mr_results[[m]]
      if (is.null(res)) return(NULL)
      or_ci <- calc_or_ci(res$beta, res$se)
      data.frame(
        exposure = exp_name, outcome = out_name, method = res$method,
        nsnp = res$nsnp, beta = res$beta, se = res$se, pval = res$pval,
        OR = or_ci$OR, OR_LCI = or_ci$OR_LCI, OR_UCI = or_ci$OR_UCI,
        Q = ifelse(!is.null(res$Q), res$Q, NA),
        Q_pval = ifelse(!is.null(res$Q_pval), res$Q_pval, NA),
        egger_intercept = ifelse(m == "egger", res$intercept, NA),
        egger_pval = ifelse(m == "egger", res$pval_intercept, NA),
        stringsAsFactors = FALSE
      )
    }))
    
    if (!is.null(mr_results$ivw)) {
      cat(sprintf("      IVW: beta = %.3f, OR = %.2f (%.2f-%.2f), p = %.2e\n",
                  mr_results$ivw$beta, exp(mr_results$ivw$beta),
                  exp(mr_results$ivw$beta - 1.96*mr_results$ivw$se),
                  exp(mr_results$ivw$beta + 1.96*mr_results$ivw$se),
                  mr_results$ivw$pval))
    }
    
    result$success <- TRUE
    result$data <- list(results = results_df, harmonized = dat, loo = loo)
    
  }, error = function(e) {
    result$error <- e$message
    cat(sprintf("      [错误] %s\n", e$message))
  })
  
  return(result)
}

# 执行所有方案
for (scheme in all_schemes) {
  cat(sprintf("\n\n======== %s ========\n", scheme$name))
  
  for (exp_name in names(scheme$exposures)) {
    exp_id <- scheme$exposures[[exp_name]]
    
    for (out_name in names(scheme$outcomes)) {
      out_id <- scheme$outcomes[[out_name]]
      key <- paste0(scheme$name, "__", exp_name, "__", out_name)
      
      result <- run_single_mr(exp_id, out_id, exp_name, out_name)
      
      analysis_log[[key]] <- list(
        scheme = scheme$name, exposure = exp_name, outcome = out_name,
        exp_id = exp_id, out_id = out_id,
        success = result$success, error = result$error
      )
      
      if (result$success && !is.null(result$data)) {
        all_results[[key]] <- result$data$results
        all_harmonized[[key]] <- result$data$harmonized
        all_loo[[key]] <- result$data$loo
      }
      
      Sys.sleep(3)  # API限速 - 增加到3秒
    }
  }
  Sys.sleep(5)  # 方案之间的额外延迟
}

# =============================================================================
# 7. 汇总结果
# =============================================================================
cat("\n\n=== Step 7: 汇总结果 ===\n")

# 分析日志
log_df <- do.call(rbind, lapply(analysis_log, as.data.frame))
cat(sprintf("  总分析数: %d\n", nrow(log_df)))
cat(sprintf("  成功: %d, 失败: %d\n", sum(log_df$success), sum(!log_df$success)))

if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)
  rownames(combined_results) <- NULL
  
  cat(sprintf("  总结果行数: %d\n", nrow(combined_results)))
  
  # 显示IVW结果
  ivw_results <- combined_results %>%
    filter(method == "IVW") %>%
    arrange(pval)
  
  cat("\n=== 主要MR结果 (IVW方法，按p值排序) ===\n")
  print(ivw_results[, c("exposure", "outcome", "nsnp", "beta", "OR", "OR_LCI", "OR_UCI", "pval")])
  
  # 显著结果
  sig_results <- ivw_results %>% filter(pval < 0.05)
  cat(sprintf("\n显著因果关系 (p < 0.05): %d 个\n", nrow(sig_results)))
  if (nrow(sig_results) > 0) {
    print(sig_results)
  }
  
  # =============================================================================
  # 8. 保存结果
  # =============================================================================
  cat("\n=== Step 8: 保存结果 ===\n")
  
  write.csv(combined_results, file.path(output_dir, "mr_real_gwas_all_methods.csv"), row.names = FALSE)
  cat("  已保存: mr_real_gwas_all_methods.csv\n")
  
  write.csv(ivw_results, file.path(output_dir, "mr_real_gwas_ivw.csv"), row.names = FALSE)
  cat("  已保存: mr_real_gwas_ivw.csv\n")
  
  write.csv(log_df, file.path(output_dir, "mr_analysis_log.csv"), row.names = FALSE)
  cat("  已保存: mr_analysis_log.csv\n")
  
  # =============================================================================
  # 9. 生成可视化
  # =============================================================================
  cat("\n=== Step 9: 生成可视化 ===\n")
  
  # 森林图
  if (nrow(ivw_results) > 0) {
    ivw_results$label <- paste0(ivw_results$exposure, " → ", ivw_results$outcome)
    ivw_results$significant <- ivw_results$pval < 0.05
    
    p_forest <- ggplot(ivw_results, aes(x = OR, y = reorder(label, OR))) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_errorbarh(aes(xmin = OR_LCI, xmax = OR_UCI), height = 0.2, linewidth = 0.5) +
      geom_point(aes(color = significant), size = 3) +
      scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = COL_UP), guide = "none") +
      scale_x_log10() +
      labs(x = "Odds Ratio (95% CI)", y = NULL,
           title = "MR Forest Plot: Real GWAS Data",
           subtitle = sprintf("n = %d exposure-outcome pairs", nrow(ivw_results))) +
      theme_pub +
      theme(axis.text.y = element_text(size = 7))
    
    ggsave(file.path(output_dir, "mr_real_gwas_forest.pdf"), p_forest,
           width = mm2in(183), height = mm2in(max(80, nrow(ivw_results) * 12)), dpi = 300)
    cat("  已保存: mr_real_gwas_forest.pdf\n")
  }
  
  # 散点图和漏斗图
  if (length(all_harmonized) > 0) {
    pdf(file.path(output_dir, "mr_real_gwas_scatter.pdf"), width = mm2in(183), height = mm2in(150))
    for (key in names(all_harmonized)) {
      dat <- all_harmonized[[key]]
      res <- all_results[[key]]
      parts <- strsplit(key, "__")[[1]]
      exp_name <- parts[2]
      out_name <- parts[3]
      
      p <- ggplot(dat, aes(x = beta_exp, y = beta_out)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_errorbar(aes(ymin = beta_out - 1.96*se_out, ymax = beta_out + 1.96*se_out),
                      width = 0, alpha = 0.3) +
        geom_errorbarh(aes(xmin = beta_exp - 1.96*se_exp, xmax = beta_exp + 1.96*se_exp),
                       height = 0, alpha = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
      
      # 添加回归线
      ivw_beta <- res$beta[res$method == "IVW"][1]
      if (!is.na(ivw_beta)) {
        p <- p + geom_abline(intercept = 0, slope = ivw_beta, color = COL_UP, linewidth = 0.8)
      }
      
      p <- p + labs(x = paste0("SNP effect on ", exp_name),
                    y = paste0("SNP effect on ", out_name),
                    title = paste0(exp_name, " → ", out_name)) +
        theme_pub
      
      print(p)
    }
    dev.off()
    cat("  已保存: mr_real_gwas_scatter.pdf\n")
  }
  
  # =============================================================================
  # 10. 生成摘要报告
  # =============================================================================
  cat("\n=== Step 10: 生成摘要报告 ===\n")
  
  report <- c(
    "================================================================================",
    "HAE Multi-Omics: Real GWAS Mendelian Randomization Analysis",
    paste0("Generated: ", Sys.time()),
    "================================================================================",
    "",
    "1. ANALYSIS OVERVIEW",
    "--------------------",
    sprintf("Total exposure-outcome pairs attempted: %d", nrow(log_df)),
    sprintf("Successful analyses: %d", sum(log_df$success)),
    sprintf("Failed analyses: %d", sum(!log_df$success)),
    "",
    "2. ANALYSIS SCHEMES",
    "-------------------"
  )
  
  for (scheme in all_schemes) {
    report <- c(report, sprintf("  - %s", scheme$name))
  }
  
  report <- c(report, "", "3. SIGNIFICANT CAUSAL RELATIONSHIPS (p < 0.05)", "----------------------------------------------")
  
  if (nrow(sig_results) > 0) {
    for (i in 1:nrow(sig_results)) {
      report <- c(report,
                  sprintf("  %s → %s:", sig_results$exposure[i], sig_results$outcome[i]),
                  sprintf("    nSNP = %d, OR = %.2f (95%%CI: %.2f-%.2f), p = %.2e",
                          sig_results$nsnp[i], sig_results$OR[i],
                          sig_results$OR_LCI[i], sig_results$OR_UCI[i], sig_results$pval[i]))
    }
  } else {
    report <- c(report, "  No significant causal relationships detected at p < 0.05")
  }
  
  report <- c(report, "", "4. SENSITIVITY ANALYSIS", "------------------------")
  
  for (key in names(all_results)) {
    res <- all_results[[key]]
    egger <- res[res$method == "MR-Egger", ]
    ivw <- res[res$method == "IVW", ]
    parts <- strsplit(key, "__")[[1]]
    
    if (nrow(egger) > 0 && nrow(ivw) > 0) {
      report <- c(report,
                  sprintf("  %s → %s:", parts[2], parts[3]),
                  sprintf("    Heterogeneity Q p-value: %.3f", ivw$Q_pval[1]),
                  sprintf("    Egger intercept p-value: %.3f", egger$egger_pval[1]))
    }
  }
  
  report <- c(report, "", "5. OUTPUT FILES", "---------------",
              "  - mr_real_gwas_all_methods.csv",
              "  - mr_real_gwas_ivw.csv",
              "  - mr_analysis_log.csv",
              "  - mr_real_gwas_forest.pdf",
              "  - mr_real_gwas_scatter.pdf",
              "", "================================================================================")
  
  writeLines(report, file.path(output_dir, "mr_real_gwas_summary.txt"))
  cat("  已保存: mr_real_gwas_summary.txt\n")
  
} else {
  cat("\n[警告] 无成功的MR分析结果\n")
  cat("\n失败原因汇总:\n")
  print(table(log_df$error))
}

# =============================================================================
# 完成
# =============================================================================
cat("\n
╔═══════════════════════════════════════════════════════════════════════════════╗
║      真实GWAS MR分析完成！                                                    ║
╚═══════════════════════════════════════════════════════════════════════════════╝
\n")

cat(sprintf("输出目录: %s\n", output_dir))
cat("\n=== Session Info ===\n")
print(sessionInfo())
