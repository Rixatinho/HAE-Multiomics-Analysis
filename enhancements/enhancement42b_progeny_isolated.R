#!/usr/bin/env Rscript
# enhancement42b: PROGENy pathway footprint activity (isolated process to
#                 avoid OmnipathR-cache memory issues seen in 42)
suppressMessages(library(decoupleR))
ROOT <- "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学/02_analysis"
setwd(ROOT)
OUT <- file.path(ROOT, "results", "enhancement42_dorothea_progeny")

emat <- read.csv("data/processed/transcriptomics_logcpm_paired.csv",
                 row.names = 1, check.names = FALSE)
deg_map <- read.csv("results/phase1_diff/DEGs_Adjacent_vs_Normal.csv")
id2sym <- setNames(deg_map$gene_name, deg_map$gene_id)
rn <- rownames(emat)
sym <- ifelse(!is.na(id2sym[rn]) & id2sym[rn] != "" & !grepl("^ENSG", id2sym[rn]),
             id2sym[rn], sub("_[0-9]+$", "", rn))
keep <- !duplicated(sym)
emat <- emat[keep, ]; sym <- sym[keep]; rownames(emat) <- sym

samples <- colnames(emat)
grp <- ifelse(grepl("^Normal", samples), "Normal", "Adjacent")
pid <- sub("^(Normal|Adjacent)", "", samples)
pairs12 <- intersect(pid[grp == "Normal"], pid[grp == "Adjacent"])
adj_cols <- paste0("Adjacent", pairs12); nor_cols <- paste0("Normal", pairs12)

prog <- decoupleR::get_resource("PROGENy")
prog <- prog[prog$genesymbol %in% rownames(emat) & prog$p_value < 0.05, ]
prog <- prog[, c("pathway", "genesymbol", "weight")]
prog <- prog[!duplicated(prog[, c("pathway", "genesymbol")]), ]
names(prog) <- c("source", "target", "mor")
cat("PROGENy core links:", nrow(prog), "| pathways:", length(unique(prog$source)), "\n")

res <- run_ulm(emat, prog, .source = "source", .target = "target", .mor = "mor")
# robust manual spread (avoid reshape quirks)
pw_names <- unique(res$source)
mat_pw <- sapply(pw_names, function(pw) {
  d <- res[res$source == pw, ]
  v <- d$score; names(v) <- d$condition
  v[colnames(emat)]
})
rownames(mat_pw) <- colnames(emat)

pw_diff <- t(sapply(colnames(mat_pw), function(pw) {
  a <- mat_pw[adj_cols, pw]; n <- mat_pw[nor_cols, pw]
  d <- a - n
  c(mean_d = mean(d),
    wilcox_p = tryCatch(wilcox.test(a, n, paired = TRUE)$p.value,
                        error = function(e) NA))
}))
pw_diff <- as.data.frame(pw_diff)
pw_diff$BH <- p.adjust(pw_diff$wilcox_p, "BH")
pw_diff$pathway <- rownames(pw_diff)
pw_diff <- pw_diff[order(pw_diff$wilcox_p), ]
print(pw_diff)
write.csv(pw_diff, file.path(OUT, "progeny_pathway_paired_diff.csv"), row.names = FALSE)
write.csv(mat_pw, file.path(OUT, "progeny_activity_matrix.csv"))

pdf(file.path(OUT, "progeny_heatmap.pdf"), width = 10, height = 5)
ord <- order(grp, decreasing = TRUE)
sub_pw <- t(mat_pw[ord, ])
sub_pw <- apply(sub_pw, 1, function(x) (x - mean(x)) / (sd(x) + 1e-9)) # z by pathway
stats::heatmap(sub_pw, Colv = NA, Rowv = NA, scale = "none",
               col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
               main = "PROGENy pathway footprints (Adjacent | Normal)")
dev.off()
cat("Done.\n")
