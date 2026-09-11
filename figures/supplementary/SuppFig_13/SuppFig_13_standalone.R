#!/usr/bin/env Rscript
# =============================================================================
# SuppFig_13_standalone.R
# Supplementary Figure 13: Machine Learning Classification & Diagnostic Nomogram
# HAE Multi-omics Study | Target: EBioMedicine (Lancet family)
# =============================================================================
# Panels (8):
#   a = Classification metrics (Sens/Spec/PPV/NPV) per ML model (CS1 vs CS2)
#   b = Recursive feature elimination (RFE) accuracy curve with optimal point
#   c = Random Forest top 20 feature importance (gene symbols; UniProt fallback)
#   d = Feature Jaccard similarity matrix between ML methods
#   e = Multi-omics nomogram ROC curves (Full + LOO-CV) + bootstrap-corrected AUC
#   f = Calibration curve (loess smoothed) + Brier score
#   g = Decision curve analysis (clinically relevant 0-0.5 range)
#   h = LASSO-selected nomogram coefficients (mapped metabolite names)
# =============================================================================
# Usage: conda run -n multiomics Rscript SuppFig_13_standalone.R
# =============================================================================

cat("=== Supplementary Figure 13: ML Classification & Nomogram ===\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(ggsci)
  library(jsonlite)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

pdfFonts(Arial = pdfFonts()$Helvetica)
postscriptFonts(Arial = postscriptFonts()$Helvetica)

# =============================================================================
# Paths
# =============================================================================
BASE <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/\u4e2a\u4eba\u6587\u6863/Word/\u6587\u732e\u5199\u4f5c/2.\u809d\u5305\u866b/20260317-\u809d\u5305\u866b\u591a\u7ec4\u5b66"
RES  <- file.path(BASE, "02_analysis/results")
OUT  <- file.path(BASE, "04_figures/supplementary/SuppFig_13")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Constants
# =============================================================================
FONT_FAMILY <- "Arial"
FONT_GRID   <- "Arial"
MM <- 25.4
DPI <- 600

FS_GEOM <- 2.55  # ggplot geom_text size (~7.5pt)
FS_TAG  <- 12    # panel label

COL_UP   <- "#CD534CFF"; COL_DOWN <- "#0073C2FF"; COL_NS <- "#868686FF"
COL_TC   <- "#0073C2FF"; COL_PR   <- "#CD534CFF"; COL_MT <- "#EFC000FF"
COL_NORMAL <- "#7AA6DCFF"; COL_ADJACENT <- "#CD534CFF"
PAL_CAT  <- pal_jco("default")(10)

# =============================================================================
# Theme
# =============================================================================
theme_nc <- theme_bw(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text = element_text(family = FONT_FAMILY, size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.text = element_text(size = 8, color = "black", family = FONT_FAMILY),
    axis.title = element_text(size = 9, face = "bold", family = FONT_FAMILY),
    plot.title = element_text(size = 10, face = "bold", hjust = 0,
                              family = FONT_FAMILY,
                              margin = margin(t = 0, b = 2, l = 6, unit = "mm")),
    plot.title.position = "plot",
    legend.text = element_text(size = 8, family = FONT_FAMILY),
    legend.title = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    legend.key.size = unit(3, "mm"),
    legend.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold", family = FONT_FAMILY),
    plot.margin = margin(3, 4, 3, 4, "mm")
  )
theme_set(theme_nc)

gp_rn <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID)
gp_cn <- function(sz = 8) gpar(fontsize = sz, fontfamily = FONT_GRID, fontface = "bold")
ht_opt$message <- FALSE

save_panel_pdf <- function(filename, w_mm, h_mm, expr) {
  fp <- file.path(OUT, filename)
  cairo_pdf(fp, width = w_mm / MM, height = h_mm / MM, family = FONT_FAMILY)
  tryCatch(force(expr), error = function(e) message("  ERROR: ", e$message))
  dev.off()
  cat("  ->", filename, sprintf("(%.0fx%.0fmm)\n", w_mm, h_mm))
}

# Metabolite display name handler (per project rules)
abbreviate_metabolite <- function(x) {
  if (is.na(x) || x == "") return(x)
  x <- gsub("glycerophosphoethanolamine", "GPE", x, ignore.case = TRUE)
  x <- gsub("glycerophosphocholine",      "GPC", x, ignore.case = TRUE)
  x <- gsub("lysophosphatidylcholine",    "LPC", x, ignore.case = TRUE)
  x <- gsub("lysophosphatidylethanolamine","LPE", x, ignore.case = TRUE)
  x <- gsub("phosphatidylcholine",        "PC",  x, ignore.case = TRUE)
  x <- gsub("phosphatidylethanolamine",   "PE",  x, ignore.case = TRUE)
  x <- gsub("sphingomyelin",              "SM",  x, ignore.case = TRUE)
  x <- gsub("Hydroxy",                    "OH",  x)
  x <- gsub("alpha", "\u03b1", x); x <- gsub("beta", "\u03b2", x)
  x <- gsub("gamma", "\u03b3", x)
  x <- gsub("\\|A", "\u03b1", x); x <- gsub("\\|B", "\u03b2", x)
  x
}

# =============================================================================
# DATA LOADING
# =============================================================================
cat("\n--- Loading data ---\n")
cm   <- read.csv(file.path(RES, "enhancement17_ml_ensemble/confusion_matrix_metrics.csv"))
rfe  <- read.csv(file.path(RES, "enhancement17_ml_ensemble/rfe_accuracy_curve.csv"))
fi   <- read.csv(file.path(RES, "enhancement17_ml_ensemble/rf_top20_features.csv"))
# Project-internal ENSP -> gene-symbol master table (full DEPs, all 8789 proteins)
deps_full <- read.csv(file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"),
                      stringsAsFactors = FALSE)
jac  <- read.csv(file.path(RES, "enhancement17_ml_ensemble/feature_jaccard_overlap.csv"),
                 row.names = 1)
pred <- read.csv(file.path(RES, "clinical_nomogram/nomogram_predictions.csv"))
nomo <- fromJSON(file.path(RES, "clinical_nomogram/nomogram_results.json"))
cal  <- read.csv(file.path(RES, "clinical_nomogram/calibration_data.csv"))
dca  <- read.csv(file.path(RES, "clinical_nomogram/dca_data.csv"))
feat <- read.csv(file.path(RES, "clinical_nomogram/nomogram_features.csv"))
met_annot <- read.csv(file.path(RES, "fix06_metabolite_grading/metabolite_MSI_grading.csv"))
cat("  All data loaded.\n")

# =============================================================================
# PANEL A: Classification metrics per model (grouped bars + value labels)
# =============================================================================
cat("\n--- Panel a: Classification metrics ---\n")
p_A <- NULL
tryCatch({
  cm_long <- cm %>%
    pivot_longer(-Model, names_to = "Metric", values_to = "Value")
  cm_long$Model  <- factor(cm_long$Model,
                           levels = c("XGBoost", "LASSO", "SVM", "RF", "Stacking"))
  cm_long$Metric <- factor(cm_long$Metric,
                           levels = c("Sensitivity", "Specificity", "PPV", "NPV"))

  metric_pal <- c("Sensitivity" = "#0073C2FF",
                  "Specificity" = "#EFC000FF",
                  "PPV"         = "#868686FF",
                  "NPV"         = "#CD534CFF")

  p_A <- ggplot(cm_long, aes(x = Model, y = Value, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72,
             color = "white", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", Value)),
              position = position_dodge(width = 0.78),
              vjust = -0.4, size = FS_GEOM - 0.4,
              family = FONT_FAMILY, color = "black") +
    scale_fill_manual(values = metric_pal, name = NULL) +
    scale_y_continuous(limits = c(0, 1.18),
                       breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    labs(title = "Classification metrics across ML models",
         x = NULL, y = "Score") +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, -2, 0),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  save_panel_pdf("Supp13a_confusion_metrics.pdf", 88, 58, print(p_A))
}, error = function(e) cat("  ERROR Panel a:", e$message, "\n"))

# =============================================================================
# PANEL B: RFE accuracy curve (highlight optimal point + LOO error band)
# =============================================================================
cat("\n--- Panel b: RFE curve ---\n")
p_B <- NULL
tryCatch({
  # SE for LOOCV at n=12 samples per accuracy
  n_samp <- 12
  rfe$SE <- sqrt(rfe$Accuracy * (1 - rfe$Accuracy) / n_samp)
  rfe$lo <- pmax(0, rfe$Accuracy - rfe$SE)
  rfe$hi <- pmin(1, rfe$Accuracy + rfe$SE)
  optimal <- rfe[which.max(rfe$Accuracy), ]

  p_B <- ggplot(rfe, aes(x = NumFeatures, y = Accuracy)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = COL_TC, alpha = 0.18) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.35) +
    geom_line(color = COL_TC, linewidth = 0.9) +
    geom_point(color = COL_TC, fill = "white", shape = 21,
               size = 2.4, stroke = 1.0) +
    geom_point(data = optimal, color = COL_UP, fill = COL_UP,
               shape = 21, size = 3.4, stroke = 1.0) +
    annotate("text",
             x = optimal$NumFeatures, y = optimal$Accuracy + 0.06,
             label = sprintf("Optimal: %d features\nAccuracy = %.3f",
                             optimal$NumFeatures, optimal$Accuracy),
             size = FS_GEOM, family = FONT_FAMILY, hjust = 0.5,
             color = COL_UP, fontface = "bold") +
    annotate("text", x = max(rfe$NumFeatures), y = 0.46,
             label = "Chance (0.5)", hjust = 1, size = FS_GEOM - 0.2,
             family = FONT_FAMILY, color = "grey40") +
    scale_x_continuous(breaks = c(5, 10, 15, 20, 30, 50, 75, 100)) +
    scale_y_continuous(limits = c(0.3, 0.95),
                       breaks = seq(0.3, 0.9, 0.1)) +
    labs(title = "Recursive feature elimination",
         x = "Number of features (log-spaced)",
         y = "LOO-CV accuracy") +
    coord_cartesian(clip = "off")
  save_panel_pdf("Supp13b_RFE_curve.pdf", 95, 58, print(p_B))
}, error = function(e) cat("  ERROR Panel b:", e$message, "\n"))

# =============================================================================
# PANEL C: RF feature importance (top 20; ENSP -> gene symbol, UniProt fallback)
# =============================================================================
cat("\n--- Panel c: RF feature importance ---\n")
p_C <- NULL
tryCatch({
  fi <- fi[order(fi$Importance, decreasing = TRUE), ]
  if (nrow(fi) > 20) fi <- head(fi, 20)
  # 1. Clean ID: strip "_PROT" suffix; keep version-stripped ENSP for matching
  fi$ensp_clean <- gsub("\\.[0-9]+_PROT$", "", fi$Feature)

  # 2. Primary mapping: project DEPs master table (Protein col has ENSP.<ver>)
  deps_full$ensp_clean <- gsub("\\.[0-9]+$", "", deps_full$Protein)
  m <- match(fi$ensp_clean, deps_full$ensp_clean)
  fi$Symbol <- deps_full$gene_name[m]
  # Treat blank / underscore / NA as unresolved
  fi$Symbol[is.na(fi$Symbol) | fi$Symbol %in% c("", "_", "-")] <- NA_character_

  # 3. Fallback: org.Hs.eg.db ENSEMBLPROT -> UNIPROT (then -> SYMBOL)
  miss <- which(is.na(fi$Symbol))
  if (length(miss) > 0L) {
    ok <- requireNamespace("AnnotationDbi", quietly = TRUE) &&
          requireNamespace("org.Hs.eg.db", quietly = TRUE)
    if (ok) {
      valid_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,
                                        keytype = "ENSEMBLPROT")
      q_keys <- intersect(fi$ensp_clean[miss], valid_keys)
      if (length(q_keys) > 0L) {
        suppressMessages({
          ann <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                       keys    = q_keys,
                                       keytype = "ENSEMBLPROT",
                                       columns = c("SYMBOL", "UNIPROT"))
        })
        sym_lookup <- tapply(ann$SYMBOL,  ann$ENSEMBLPROT, function(x) x[!is.na(x)][1])
        uni_lookup <- tapply(ann$UNIPROT, ann$ENSEMBLPROT, function(x) x[!is.na(x)][1])
        sym_hit <- sym_lookup[fi$ensp_clean[miss]]
        uni_hit <- uni_lookup[fi$ensp_clean[miss]]
        fi$Symbol[miss] <- ifelse(!is.na(sym_hit), sym_hit,
                                  ifelse(!is.na(uni_hit), uni_hit, NA_character_))
      }
    }
  }

  # 4. Last resort: keep cleaned ENSP for any still-unresolved entry
  unresolved <- is.na(fi$Symbol) | fi$Symbol == ""
  fi$Symbol[unresolved] <- fi$ensp_clean[unresolved]

  cat("  Mapping summary (ENSP -> display label):\n")
  print(fi[, c("ensp_clean", "Symbol")])
  if (any(unresolved)) {
    cat("  [Legend note] Unresolved IDs kept as ENSP:",
        paste(fi$ensp_clean[unresolved], collapse = ", "), "\n")
  }

  # Disambiguate any duplicate symbols by appending ENSP tail
  dup <- duplicated(fi$Symbol) | duplicated(fi$Symbol, fromLast = TRUE)
  if (any(dup)) {
    fi$Symbol[dup] <- paste0(fi$Symbol[dup], " (",
                             substr(fi$ensp_clean[dup], 12, 15), ")")
  }
  fi$display_id <- factor(fi$Symbol, levels = rev(fi$Symbol))

  p_C <- ggplot(fi, aes(x = Importance, y = display_id)) +
    geom_segment(aes(x = 0, xend = Importance, yend = display_id),
                 color = "grey80", linewidth = 0.5) +
    geom_point(color = COL_PR, fill = COL_PR, shape = 21,
               size = 2.6, stroke = 0.6) +
    geom_text(aes(label = sprintf("%.2f", Importance)),
              hjust = -0.3, size = FS_GEOM - 0.3,
              family = FONT_FAMILY, color = "black") +
    scale_x_continuous(limits = c(0, max(fi$Importance) * 1.18),
                       expand = c(0, 0)) +
    labs(title = "RF top 20 feature importance",
         x = "Mean decrease in accuracy",
         y = NULL) +
    theme(axis.text.y = element_text(size = 7, family = FONT_FAMILY,
                                     face = "italic"))
  save_panel_pdf("Supp13c_feature_importance.pdf", 88, 64, print(p_C))
}, error = function(e) cat("  ERROR Panel c:", e$message, "\n"))

# =============================================================================
# PANEL D: Jaccard similarity matrix between ML methods
# =============================================================================
cat("\n--- Panel d: Jaccard matrix ---\n")
ht_D <- NULL
tryCatch({
  mat <- as.matrix(jac); mode(mat) <- "numeric"
  col_fun <- colorRamp2(c(0, 0.05, 0.5, 1),
                        c("#FAFAFA", "#FFE9C5", "#EFC000FF", COL_UP))

  pdf_path <- file.path(OUT, "Supp13d_jaccard_matrix.pdf")
  cairo_pdf(pdf_path, width = 95/MM, height = 64/MM, family = FONT_FAMILY)
  ht <- Heatmap(
    mat, name = "Jaccard", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = 1.2),
    row_names_gp = gp_rn(8),
    column_names_gp = gp_cn(8),
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = "Cross-method feature overlap (Jaccard)",
    column_title_gp = gpar(fontsize = 10, fontfamily = FONT_GRID,
                           fontface = "bold"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      v <- mat[i, j]
      grid.text(sprintf("%.3f", v), x, y,
                gp = gpar(fontsize = 7.5, fontfamily = FONT_GRID,
                          col = ifelse(v > 0.6, "white", "black"),
                          fontface = ifelse(i == j, "bold", "plain")))
    },
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontfamily = FONT_GRID, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontfamily = FONT_GRID),
      legend_height = unit(28, "mm")
    ))
  draw(ht, padding = unit(c(2, 4, 2, 2), "mm"))
  ht_D <<- ht
  dev.off()
  cat("  -> Supp13d_jaccard_matrix.pdf (95x64mm)\n")
}, error = function(e) cat("  ERROR Panel d:", e$message, "\n"))

# =============================================================================
# PANEL E: Nomogram ROC curves with bootstrap-corrected AUC annotation
# =============================================================================
cat("\n--- Panel e: ROC curves ---\n")
p_E <- NULL
tryCatch({
  compute_roc <- function(probs, labels) {
    th <- sort(unique(c(-Inf, probs, Inf)), decreasing = TRUE)
    tpr <- sapply(th, function(t) sum(probs >= t & labels == 1) / sum(labels == 1))
    fpr <- sapply(th, function(t) sum(probs >= t & labels == 0) / sum(labels == 0))
    data.frame(FPR = fpr, TPR = tpr)
  }
  roc_full <- compute_roc(pred$Predicted_prob, pred$True_label)
  roc_full$Model <- sprintf("Apparent (AUC = %.3f)", nomo$auc_full)
  roc_loo  <- compute_roc(pred$LOO_prob,       pred$True_label)
  roc_loo$Model  <- sprintf("LOO-CV (AUC = %.3f)", nomo$auc_loo)
  roc_all <- rbind(roc_full, roc_loo)

  bs_label <- sprintf("Bootstrap-corrected AUC = %.3f\n(95%% CI %.3f-%.3f, n = 1000)",
                      nomo$bootstrap_auc_mean,
                      nomo$bootstrap_auc_ci_lower,
                      nomo$bootstrap_auc_ci_upper)

  p_E <- ggplot(roc_all, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.35) +
    geom_step(linewidth = 1.0, direction = "vh") +
    annotate("rect", xmin = 0.32, xmax = 0.99, ymin = 0.05, ymax = 0.30,
             fill = "white", color = "grey70", linewidth = 0.3, alpha = 0.95) +
    annotate("text", x = 0.36, y = 0.18, label = bs_label,
             hjust = 0, size = FS_GEOM - 0.25, family = FONT_FAMILY,
             color = "black") +
    scale_color_manual(values = c(COL_PR, COL_TC), name = NULL) +
    scale_linetype_manual(values = c("solid", "longdash"), name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    coord_equal() +
    labs(title = "Nomogram ROC (n = 14 pairs)",
         x = "1 - Specificity", y = "Sensitivity") +
    theme(legend.position = c(0.65, 0.62),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white"))
  save_panel_pdf("Supp13e_nomogram_ROC.pdf", 88, 64, print(p_E))
}, error = function(e) cat("  ERROR Panel e:", e$message, "\n"))

# =============================================================================
# PANEL F: Calibration curve (loess smoothed) + Brier score
# =============================================================================
cat("\n--- Panel f: Calibration ---\n")
p_F <- NULL
tryCatch({
  cal2 <- cal
  brier <- mean((cal2$predicted_prob - cal2$true_label)^2)

  # Quantile bins (deciles trimmed to 5 because n is small)
  set.seed(7); cal2$bin <- ntile(cal2$predicted_prob, 5)
  cal_summary <- cal2 %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(predicted_prob),
              mean_true = mean(true_label),
              n = n(), .groups = "drop")

  p_F <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.35) +
    geom_smooth(data = cal2,
                aes(x = predicted_prob, y = true_label),
                method = "loess", span = 1.0, se = TRUE,
                color = COL_PR, fill = COL_PR, alpha = 0.18,
                linewidth = 0.9) +
    geom_point(data = cal_summary,
               aes(x = mean_pred, y = mean_true, size = n),
               color = COL_PR, fill = "white", shape = 21, stroke = 1.1) +
    scale_size_continuous(range = c(2.4, 5.5), name = "n") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    annotate("text", x = 0.05, y = 0.93,
             label = sprintf("Brier score = %.3f", brier),
             hjust = 0, size = FS_GEOM, family = FONT_FAMILY,
             fontface = "bold", color = "black") +
    annotate("text", x = 0.95, y = 0.05, label = "Ideal",
             hjust = 1, size = FS_GEOM - 0.3, family = FONT_FAMILY,
             color = "grey40", angle = 45) +
    coord_equal() +
    labs(title = "Calibration (n = 28; loess smoothing)",
         x = "Predicted probability",
         y = "Observed frequency") +
    theme(legend.position = c(0.86, 0.30),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key.size = unit(2.5, "mm"))
  save_panel_pdf("Supp13f_calibration.pdf", 95, 64, print(p_F))
}, error = function(e) cat("  ERROR Panel f:", e$message, "\n"))

# =============================================================================
# PANEL G: Decision curve analysis (clinically relevant 0-0.5 range)
# =============================================================================
cat("\n--- Panel g: DCA ---\n")
p_G <- NULL
tryCatch({
  dca_long <- dca %>%
    pivot_longer(-Threshold, names_to = "Strategy", values_to = "Net_benefit") %>%
    mutate(Strategy = gsub("Net_benefit_", "", Strategy))
  dca_long$Strategy <- factor(dca_long$Strategy,
                              levels = c("model", "treat_all", "treat_none"),
                              labels = c("Nomogram", "Treat all", "Treat none"))

  # Restrict to clinically relevant range
  dca_long <- dca_long %>% filter(Threshold <= 0.5)

  p_G <- ggplot(dca_long,
                aes(x = Threshold, y = Net_benefit,
                    color = Strategy, linetype = Strategy)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c("Nomogram"   = COL_PR,
                                  "Treat all"  = COL_TC,
                                  "Treat none" = "grey50"), name = NULL) +
    scale_linetype_manual(values = c("Nomogram"   = "solid",
                                     "Treat all"  = "longdash",
                                     "Treat none" = "dotted"), name = NULL) +
    scale_x_continuous(limits = c(0, 0.5),
                       breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.05, 0.55),
                       breaks = seq(0, 0.5, 0.1)) +
    labs(title = "Decision curve analysis",
         x = "Threshold probability", y = "Net benefit") +
    guides(color = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1)) +
    theme(legend.position = "top",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, -2, 0))
  save_panel_pdf("Supp13g_DCA.pdf", 88, 59, print(p_G))
}, error = function(e) cat("  ERROR Panel g:", e$message, "\n"))

# =============================================================================
# PANEL H: LASSO-selected nomogram features (mapped metabolite names)
# =============================================================================
cat("\n--- Panel h: LASSO coefficients ---\n")
p_H <- NULL
tryCatch({
  feat$Compound_ID <- gsub("^M:", "", feat$Feature)
  met_map <- setNames(met_annot$Name, met_annot$Compound_ID)
  feat$met_name <- met_map[feat$Compound_ID]
  feat$display_name <- ifelse(!is.na(feat$met_name) & feat$met_name != "",
                              vapply(feat$met_name, abbreviate_metabolite,
                                     character(1)),
                              feat$Compound_ID)
  # Replace ultra-long names with compound ID for legend annotation
  long_idx <- nchar(feat$display_name) > 32
  if (any(long_idx)) {
    cat("[Legend note] Metabolite ID replacements:\n")
    for (i in which(long_idx))
      cat(sprintf("  %s = %s\n", feat$Compound_ID[i], feat$display_name[i]))
    feat$display_name[long_idx] <- feat$Compound_ID[long_idx]
  }
  # Wrap medium names at natural break
  feat$display_name <- ifelse(nchar(feat$display_name) > 22,
                              str_wrap(feat$display_name, width = 22),
                              feat$display_name)

  feat$Direction <- ifelse(feat$Coefficient > 0,
                           "Higher in HAE", "Lower in HAE")
  feat$display_name <- factor(feat$display_name,
                              levels = feat$display_name[order(feat$Coefficient)])

  dir_pal <- c("Higher in HAE" = COL_UP, "Lower in HAE" = COL_DOWN)

  p_H <- ggplot(feat, aes(x = Coefficient, y = display_name, fill = Direction)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
    geom_text(aes(label = sprintf("%+0.2f", Coefficient),
                  hjust = ifelse(Coefficient > 0, -0.15, 1.15)),
              size = FS_GEOM - 0.3, family = FONT_FAMILY, color = "black") +
    scale_fill_manual(values = dir_pal, name = NULL) +
    scale_x_continuous(limits = c(min(feat$Coefficient) * 1.25,
                                   max(feat$Coefficient) * 1.25),
                       expand = c(0, 0)) +
    labs(title = sprintf("LASSO coefficients (%d metabolites)", nrow(feat)),
         x = "Standardized coefficient", y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, -2, 0),
          axis.text.y = element_text(size = 7))
  save_panel_pdf("Supp13h_LASSO_coef.pdf", 95, 59, print(p_H))
}, error = function(e) cat("  ERROR Panel h:", e$message, "\n"))

# =============================================================================
# COMPOSITE: Pure-vector cairo_pdf assembly (183x245mm)
# - Uses grid::viewport() to place each ggplot / ComplexHeatmap object directly
#   on a single cairo_pdf canvas; NO magick raster intermediate.
# - Verification: pdfimages -list SuppFig_13.pdf should output only the header.
# =============================================================================
cat("\n--- Assembling composite SuppFig_13 (pure vector, 183x245mm) ---\n")

tryCatch({
  W_TOTAL <- 183; H_TOTAL <- 245
  H1 <- 58;  H2 <- 64;  H3 <- 64
  H4 <- H_TOTAL - H1 - H2 - H3                 # 59
  W_A <- 88; W_B <- W_TOTAL - W_A              # 95
  W_C <- 88; W_D <- W_TOTAL - W_C              # 95
  W_E <- 88; W_F <- W_TOTAL - W_E              # 95
  W_G <- 88; W_H <- W_TOTAL - W_G              # 95

  cat(sprintf("  Layout: %dx%d mm; rows %d/%d/%d/%d; cols 88/95\n",
              W_TOTAL, H_TOTAL, H1, H2, H3, H4))

  # Panel viewport factory: x_mm/y_mm = top-left corner in mm from page top-left
  panel_vp <- function(x_mm, y_mm, w_mm, h_mm) {
    grid::viewport(x = grid::unit(x_mm, "mm"),
                   y = grid::unit(H_TOTAL - y_mm, "mm"),
                   width  = grid::unit(w_mm, "mm"),
                   height = grid::unit(h_mm, "mm"),
                   just = c("left", "top"))
  }

  vp_A <- panel_vp(0,    0,                 W_A, H1)
  vp_B <- panel_vp(W_A,  0,                 W_B, H1)
  vp_C <- panel_vp(0,    H1,                W_C, H2)
  vp_D <- panel_vp(W_C,  H1,                W_D, H2)
  vp_E <- panel_vp(0,    H1 + H2,           W_E, H3)
  vp_F <- panel_vp(W_E,  H1 + H2,           W_F, H3)
  vp_G <- panel_vp(0,    H1 + H2 + H3,      W_G, H4)
  vp_H <- panel_vp(W_G,  H1 + H2 + H3,      W_H, H4)

  render_composite <- function() {
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(col = NA, fill = "white"))

    # ggplot panels — print() preserves vector graphics
    if (!is.null(p_A)) print(p_A, vp = vp_A)
    if (!is.null(p_B)) print(p_B, vp = vp_B)
    if (!is.null(p_C)) print(p_C, vp = vp_C)
    if (!is.null(p_E)) print(p_E, vp = vp_E)
    if (!is.null(p_F)) print(p_F, vp = vp_F)
    if (!is.null(p_G)) print(p_G, vp = vp_G)
    if (!is.null(p_H)) print(p_H, vp = vp_H)

    # ComplexHeatmap panel — draw inside its own viewport (also vector)
    if (!is.null(ht_D)) {
      grid::pushViewport(vp_D)
      ComplexHeatmap::draw(ht_D, newpage = FALSE,
                           padding = grid::unit(c(2, 4, 2, 2), "mm"))
      grid::popViewport()
    }

    # Panel labels (a-h)
    label_data <- data.frame(
      text = c("a", "b", "c", "d", "e", "f", "g", "h"),
      x_mm = c(2, W_A + 2, 2, W_C + 2, 2, W_E + 2, 2, W_G + 2),
      y_mm = c(2, 2,
               H1 + 2, H1 + 2,
               H1 + H2 + 2, H1 + H2 + 2,
               H1 + H2 + H3 + 2, H1 + H2 + H3 + 2),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(label_data))) {
      grid::grid.text(
        label = label_data$text[i],
        x = grid::unit(label_data$x_mm[i], "mm"),
        y = grid::unit(H_TOTAL - label_data$y_mm[i], "mm"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = FS_TAG, fontface = "bold",
                        fontfamily = FONT_FAMILY)
      )
    }
  }

  # Pure-vector PDF
  cairo_pdf(file.path(OUT, "SuppFig_13.pdf"),
            width = W_TOTAL/MM, height = H_TOTAL/MM, family = FONT_FAMILY)
  render_composite(); dev.off()
  cat("  -> SuppFig_13.pdf  (vector)\n")

  # Cairo PNG (rasterised from same viewport tree)
  px_per_mm <- DPI / 25.4
  px_W <- round(W_TOTAL * px_per_mm)
  px_H <- round(H_TOTAL * px_per_mm)
  grDevices::png(file.path(OUT, "SuppFig_13.png"),
                 width = px_W, height = px_H, res = DPI, type = "cairo")
  render_composite(); dev.off()
  cat("  -> SuppFig_13.png  (raster, 600 DPI)\n")

  # Cairo TIFF for journal submission
  grDevices::tiff(file.path(OUT, "SuppFig_13.tiff"),
                  width = px_W, height = px_H, res = DPI,
                  compression = "lzw", type = "cairo")
  render_composite(); dev.off()
  cat("  -> SuppFig_13.tiff (raster, 600 DPI, LZW)\n")

  cat("  SuppFig_13 DONE (183x245mm; PDF is pure vector)\n")
}, error = function(e) cat("  ERROR assembly:", e$message, "\n"))

cat("\n=== Supplementary Figure 13 rendering complete ===\n")
