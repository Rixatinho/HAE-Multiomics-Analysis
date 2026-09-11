#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
enhancement46: 完整 GES (GSEA 模式) 独立验证 enhancement37 的 LINCS 逆转打分
=============================================================================
问题: enhancement37 用 Enrichr 简化重叠打分 (combined score) 对 5,403 个化合物
     排序得到 nintedanib rank 12。此处用经典加权 GSEA (weighted Kolmogorov-
     Smirnov enrichment score, 即 signatureSearch gess_lincs 的核心统计量)
     在全基因组排序表 (17,527 基因, t 统计量加权) 上对同一 LINCS 化合物库
     (~10,850 化合物 up/down 基因集) 重新计算逆转分数。

GES 逆转语义:
  查询排序表: 按疾病 t 统计量降序 — 顶部=疾病上调, 底部=疾病下调
  ES_up   = 化合物 UP 基因集在排序表中的加权富集分 (正=位于疾病上调侧)
  ES_down = 化合物 DOWN 基因集的加权富集分
  逆转化合物: 其 UP 基因应落在疾病下调侧 (ES_up<0), 其 DOWN 基因落在疾病
  上调侧 (ES_down>0)  =>  GES_reversal = NES_down - NES_up (越大逆转越强)

输出: 02_analysis/results/enhancement46_ges_validation/
"""
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
DEG = os.path.join(ROOT, "02_analysis/results/phase1_diff/DEGs_Adjacent_vs_Normal.csv")
ENR37 = os.path.join(ROOT, "02_analysis/results/enhancement37_lincs_cmap/lincs_cmap2020_reversal_scores.csv")
LIB = os.path.join(ROOT, "02_analysis/data/external/lincs/LINCS_L1000_Chem_Pert_Consensus_Sigs.gmt")
OUT = os.path.join(ROOT, "02_analysis/results/enhancement46_ges_validation")
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(42)

# ---------------------------------------------------------------- 1. 查询排序表
print("[1] building genome-wide ranked list ...", flush=True)
deg = pd.read_csv(DEG)
deg = deg.dropna(subset=["t", "gene_name"])
deg["gene_name"] = deg["gene_name"].astype(str).str.upper().str.strip()
# 聚合同名基因 (取 |t| 最大者)
deg = deg.loc[deg.groupby("gene_name")["t"].apply(lambda s: s.abs().idxmax())]
deg = deg.sort_values("t", ascending=False).reset_index(drop=True)
genes = deg["gene_name"].values
tstat = deg["t"].values.astype(float)
n = len(genes)
pos = {g: i for i, g in enumerate(genes)}
print(f"    ranked genes: {n}", flush=True)

# ---------------------------------------------------------------- 2. 解析 LINCS 共识库
def parse_lib(path):
    sets = {}
    with open(path) as f:
        for line in f:
            parts = [p for p in line.rstrip("\n").split("\t") if p]
            if len(parts) < 5:
                continue
            sets[parts[0]] = [g.upper() for g in parts[1:]]
    return sets

print("[2] parsing LINCS consensus gene set library ...", flush=True)
all_terms = parse_lib(LIB)
up_map, dn_map = {}, {}
for term, genes in all_terms.items():
    m = re.match(r"^(.*)\s+(Up|Down)$", term)
    if not m:
        continue
    if m.group(2) == "Up":
        up_map[m.group(1)] = genes
    else:
        dn_map[m.group(1)] = genes
print(f"    terms total: {len(all_terms)}, up: {len(up_map)}, down: {len(dn_map)}", flush=True)

compounds = sorted(set(up_map) & set(dn_map))
print(f"    compounds with both up+down: {len(compounds)}", flush=True)

# ---------------------------------------------------------------- 3. 加权 GSEA ES
def weighted_es(rank_idx, weights_all, n_total):
    """经典加权 ES: 在全部 N 个位置上构建 running sum
    (hit 位置加 w/W, miss 位置减 1/Nm), 取绝对值最大的偏离."""
    m = len(rank_idx)
    if m == 0:
        return 0.0
    w = np.abs(weights_all[rank_idx])
    W = w.sum()
    Nm = n_total - m
    inc = np.full(n_total, -1.0 / Nm)
    inc[rank_idx] = w / W
    run = np.cumsum(inc)
    i = int(np.argmax(np.abs(run)))
    return float(run[i])

print("[3] computing weighted ES for all compounds ...", flush=True)
rows = []
for c in compounds:
    up_idx = np.array(sorted(pos[g] for g in up_map[c] if g in pos), dtype=int)
    dn_idx = np.array(sorted(pos[g] for g in dn_map[c] if g in pos), dtype=int)
    es_up = weighted_es(up_idx, tstat, n)
    es_dn = weighted_es(dn_idx, tstat, n)
    rows.append({"compound": c, "n_up_mapped": len(up_idx), "n_dn_mapped": len(dn_idx),
                 "ES_up": es_up, "ES_down": es_dn, "raw_reversal": es_dn - es_up})
ges = pd.DataFrame(rows)

# ---------------------------------------------------------------- 4. 大小分层置换零模型 -> NES
print("[4] size-stratified permutation null (NES normalization) ...", flush=True)
sizes = pd.concat([ges["n_up_mapped"], ges["n_dn_mapped"]])
sizes = sizes[(sizes >= 10)]
uniq_sizes = np.array(sorted(sizes.unique()))
bins = np.array_split(uniq_sizes, min(10, len(uniq_sizes)))
null_table = {}  # bin center -> (mean, sd) of ES under random gene sets
N_PERM = 500
for b in bins:
    if len(b) == 0:
        continue
    s0 = int(np.median(b))
    null_es = np.empty(N_PERM)
    for i in range(N_PERM):
        idx = RNG.choice(n, size=s0, replace=False)
        null_es[i] = weighted_es(np.sort(idx), tstat, n)
    null_table[s0] = (null_es.mean(), null_es.std(ddof=1) + 1e-12)

def to_nes(es, size):
    best, bestd = None, 1e18
    for k, (mu, sd) in null_table.items():
        d = abs(size - k)
        if d < bestd:
            best, bestd = (mu, sd), d
    mu, sd = best
    return (es - mu) / sd

ges["NES_up"] = ges.apply(lambda r: to_nes(r.ES_up, r.n_up_mapped), axis=1)
ges["NES_down"] = ges.apply(lambda r: to_nes(r.ES_down, r.n_dn_mapped), axis=1)
ges["GES_reversal"] = ges["NES_down"] - ges["NES_up"]
ges = ges.sort_values("GES_reversal", ascending=False).reset_index(drop=True)
ges.insert(0, "rank", np.arange(1, len(ges) + 1))
ges.to_csv(os.path.join(OUT, "ges_reversal_scores.csv"), index=False)

# ---------------------------------------------------------------- 5. 与 enhancement37 交叉
print("[5] cross-validation vs enhancement37 Enrichr scoring ...", flush=True)
enr = pd.read_csv(ENR37)
enr["compound"] = enr["compound"].astype(str)
mg = ges.merge(enr[["compound", "reversal_score"]], on="compound", how="inner")
mg["rank_enr37"] = mg["reversal_score"].rank(ascending=False)
mg["rank_ges"] = mg["GES_reversal"].rank(ascending=False)
rho, p = spearmanr(mg["GES_reversal"], mg["reversal_score"])
mg = mg.sort_values("rank_ges")
mg.to_csv(os.path.join(OUT, "ges_vs_enrichr_merged.csv"), index=False)

# top30 重叠
top_n = 100
top_ges = set(mg.head(top_n)["compound"])
top_enr = set(mg.nsmallest(top_n, "rank_enr37")["compound"])
ov = len(top_ges & top_enr)
# 超几何检验
from scipy.stats import hypergeom
p_hyper = hypergeom.sf(ov - 1, len(mg), top_n, top_n)

# nintedanib 位置
nint = mg[mg["compound"].str.contains("Nintedanib", case=False)]
nint_info = nint[["compound", "rank_ges", "rank_enr37", "GES_reversal", "NES_up", "NES_down"]].to_dict("records")

print(f"    compounds merged: {len(mg)}", flush=True)
print(f"    Spearman rho (GES vs Enrichr reversal): {rho:.3f} (p={p:.2e})", flush=True)
print(f"    top-{top_n} overlap: {ov}/{top_n} (hypergeometric p={p_hyper:.2e})", flush=True)
print(f"    nintedanib: {nint_info}", flush=True)

# ---------------------------------------------------------------- 6. summary
with open(os.path.join(OUT, "SUMMARY.md"), "w") as f:
    f.write("# enhancement46: 完整 GES (加权 GSEA) 独立验证 LINCS 逆转打分\n\n")
    f.write("## 方法\n")
    f.write("- 查询签名: 全基因组 17,527 基因按 limma t 统计量排序 (非仅 top150), 加权 KS (w=|t|)\n")
    f.write("- 参考库: LINCS L1000 Chem_Pert_Consensus_Sigs 完整基因集 (Enrichr 下载, 与 enhancement37 同库)\n")
    f.write("- 统计量: 加权 ES, 经大小分层置换零模型 (500 perms/层) 归一化为 NES\n")
    f.write("- GES_reversal = NES_down − NES_up (越大逆转疾病签名越强)\n\n")
    f.write("## 与 enhancement37 (Enrichr 简化打分) 一致性\n\n")
    f.write(f"- 化合物数: {len(mg)}\n")
    f.write(f"- Spearman ρ = {rho:.3f} (p = {p:.2e})\n")
    f.write(f"- top-100 重叠: {ov}/100 (超几何 p = {p_hyper:.2e})\n\n")
    f.write("## nintedanib (核心候选)\n\n")
    for r in nint_info:
        f.write(f"- {r['compound']}: GES rank {int(r['rank_ges'])}/{len(mg)} "
                f"(Enrichr rank {int(r['rank_enr37'])}), GES_reversal={r['GES_reversal']:.2f}, "
                f"NES_up={r['NES_up']:.2f}, NES_down={r['NES_down']:.2f}\n")
    f.write("\n## GES top-20 化合物\n\n")
    f.write(mg.head(20)[["compound", "GES_reversal", "NES_up", "NES_down", "rank_enr37"]].to_markdown(index=False))
    f.write("\n")
print("[done]", flush=True)
