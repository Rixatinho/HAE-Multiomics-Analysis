#!/usr/bin/env Rscript
# =============================================================================
# Fig2b_volcano_PR.R -- Proteomics Volcano Plot
# Figure 2b | 90 x 85 mm | plot_nc_volcano
# =============================================================================

slot <- get_panel_slot("Fig2b_volcano_PR.pdf")
W <- slot$w; H <- slot$h

df <- read.csv(file.path(RES, "phase1_diff/DEPs_Adjacent_vs_Normal.csv"),
               stringsAsFactors = FALSE)

sym_col <- find_col(df, c("gene_name", "symbol", "name"))
if (!is.null(sym_col)) df$symbol <- df[[sym_col]] else df$symbol <- rownames(df)
df$symbol <- ensg_to_symbol(df$symbol)
df$symbol <- ensp_to_symbol(df$symbol)

p <- plot_nc_volcano(df, title = "Proteomics",
                     fc_col = "logFC", p_col = "P.Value",
                     symbol_col = "symbol", fc_cutoff = 1.0,
                     p_cutoff = 0.05, top_n = 15)

save_cairo("Fig2b_volcano_PR.pdf", W, H, print(p))
