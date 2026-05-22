#!/usr/bin/env Rscript
# =============================================================================
# Fig2c_volcano_MT.R -- Metabolomics Volcano Plot
# Figure 2c | 90 x 85 mm | plot_nc_volcano
# =============================================================================

slot <- get_panel_slot("Fig2c_volcano_MT.pdf")
W <- slot$w; H <- slot$h

df <- read.csv(file.path(RES, "phase1_diff/DEMs_Adjacent_vs_Normal.csv"),
               stringsAsFactors = FALSE)

sym_col <- find_col(df, c("metabolite_name", "name", "symbol"))
if (!is.null(sym_col)) df$symbol <- df[[sym_col]] else df$symbol <- rownames(df)

p <- plot_nc_volcano(df, title = "Metabolomics",
                     fc_col = "logFC", p_col = "P.Value",
                     symbol_col = "symbol", fc_cutoff = 1.0,
                     p_cutoff = 0.05, top_n = 8)

save_cairo("Fig2c_volcano_MT.pdf", W, H, print(p))
