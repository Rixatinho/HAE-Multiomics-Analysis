#!/usr/bin/env Rscript
# =============================================================================
# Fig2a_volcano_TC.R -- Transcriptomics Volcano Plot
# Figure 2a | 90 x 85 mm | factory_volcano + plot_nc_volcano
# =============================================================================

slot <- get_panel_slot("Fig2a_volcano_TC.pdf")
W <- slot$w; H <- slot$h

df <- read.csv(file.path(RES, "phase1_diff/DEGs_Adjacent_vs_Normal.csv"),
               stringsAsFactors = FALSE)

# Prepare symbol column
sym_col <- find_col(df, c("gene_name", "symbol", "name"))
if (!is.null(sym_col)) df$symbol <- df[[sym_col]] else df$symbol <- rownames(df)
df$symbol <- ensg_to_symbol(df$symbol)

p <- plot_nc_volcano(df, title = "Transcriptomics",
                     fc_col = "logFC", p_col = "P.Value",
                     symbol_col = "symbol", fc_cutoff = 1.0,
                     p_cutoff = 0.05, top_n = 15)

save_cairo("Fig2a_volcano_TC.pdf", W, H, print(p))
