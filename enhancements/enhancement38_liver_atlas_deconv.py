#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
enhancement38: 真实人肝单细胞图谱 (MacParland et al. 2018, GSE115469) 参考去卷积
=============================================================================
目标: 用真实人肝单细胞参考(8,444 细胞, 10 聚类)替换先前合成参考,
     以 CIBERSORT 式核 (ν-SVR 不可用时用 NNLS + 秩归一化) 对 24 个 bulk
     样本 (12 Normal / 12 Adjacent) 做去卷积, 与 BayesPrism (合成参考)
     结果交叉验证细胞类型位移方向.

参考: MacParland JC, et al. Nat Commun 2018;9:4383. Single cell RNA sequencing
     of human liver reveals distinct intrahepatic macrophage populations.
数据: GSE115469_Data.csv.gz (UMI counts), GSE115469_CellClusterType.txt.gz
输出: 02_analysis/results/enhancement38_liver_atlas_deconv/
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import wilcoxon, spearmanr

ROOT = "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
EXT = os.path.join(ROOT, "02_analysis/data/external/GSE115469")
PROC = os.path.join(ROOT, "02_analysis/data/processed")
OUT = os.path.join(ROOT, "02_analysis/results/enhancement38_liver_atlas_deconv")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- 1. 读单细胞参考
print("[1] loading GSE115469 atlas ...", flush=True)
sc_data = pd.read_csv(os.path.join(EXT, "GSE115469_Data.csv.gz"), index_col=0, compression="gzip")
# 行=基因, 列=细胞 (MacParland 释放格式为基因x细胞)
print("    raw matrix:", sc_data.shape, flush=True)
annot = pd.read_csv(os.path.join(EXT, "GSE115469_CellClusterType.txt.gz"), sep="\t", compression="gzip")
# 列: CellName / Sample / Cell# / Cluster# / CellType
annot.columns = [c.strip() for c in annot.columns]
cell_labels = annot.set_index("CellName")["CellType"].astype(str).str.strip()
print("    cells annotated:", len(cell_labels), "clusters:", cell_labels.nunique(), flush=True)
print("    cluster sizes:\n", cell_labels.value_counts().to_string(), flush=True)

# 对齐细胞
common = [c for c in sc_data.columns if c in cell_labels.index]
print("    matched cells:", len(common), "/", sc_data.shape[1], flush=True)
X = sc_data[common]
labels = cell_labels.loc[common]

# 基因符号清理: MacParland 用 Symbol, 可能含 "HBG2" 等; 去重取均值
X.index = X.index.astype(str).str.strip()
X = X[~(X.index == "") & ~X.index.str.startswith("MT-")]
X = X.groupby(level=0).mean()

# 聚类级表达矩阵 (CPM) — signature matrix
groups = sorted(labels.unique())
sig_rows = []
for g in groups:
    cells = labels.index[labels == g].intersection(X.columns)
    sub = X[list(cells)]
    tot = sub.sum(axis=0)
    tot = tot.replace(0, np.nan)
    sub = sub.div(tot, axis=1) * 1e6
    # 聚类 signature = 组内平均 CPM
    sig_rows.append(sub.mean(axis=1))
S = pd.concat(sig_rows, axis=1)
S.columns = groups
S = np.log2(S + 1.0)
S = S.replace([np.inf, -np.inf], np.nan).fillna(0.0)
print("[2] signature matrix:", S.shape, flush=True)

# ---------------------------------------------------------------- 2. 读 bulk
bulk = pd.read_csv(os.path.join(PROC, "transcriptomics_logcpm_paired.csv"), index_col=0)
deg = pd.read_csv(os.path.join(ROOT, "02_analysis/results/phase1_diff/DEGs_Adjacent_vs_Normal.csv"))
gmap = dict(zip(deg["gene_id"].astype(str), deg["gene_name"].astype(str)))
bulk.index = [gmap.get(g, g.split("_")[0]) for g in bulk.index.astype(str)]
bulk = bulk.groupby(level=0).mean()
bulk = np.log2(bulk + 1.0)
bulk = bulk.replace([np.inf, -np.inf], np.nan).dropna()

genes = S.index.intersection(bulk.index)
print("[3] shared genes:", len(genes), flush=True)
Sg = S.loc[genes]
Bg = bulk.loc[genes]

# 变异筛选 (signature matrix 选信息基因: 各聚类间方差 top 3000)
v = Sg.var(axis=1).sort_values(ascending=False)
sel = v.head(3000).index
Sg = Sg.loc[sel]; Bg = Bg.loc[sel]

# ---------------------------------------------------------------- 3. NNLS 去卷积
props = {}
for s in Bg.columns:
    w, _ = nnls(Sg.values, Bg[s].values)
    props[s] = w / w.sum() if w.sum() > 0 else w
P = pd.DataFrame(props).T
P.columns = groups
P.to_csv(os.path.join(OUT, "atlas_nnls_proportions.csv"))
print("[4] deconvolution done:", P.shape, flush=True)

# ---------------------------------------------------------------- 4. 配对比较 (Adjacent vs Normal)
samples = list(Bg.columns)
pats = sorted(set(s.replace("Normal", "").replace("Adjacent", "") for s in samples))
rows = []
for g in groups:
    adj = [f"Adjacent{p}" for p in pats if f"Adjacent{p}" in P.index]
    nor = [f"Normal{p}" for p in pats if f"Normal{p}" in P.index]
    pairs = [(a, n) for a, n in zip(adj, nor) if a in P.index and n in P.index]
    if not pairs:
        continue
    d = np.array([P.loc[a, g] - P.loc[n, g] for a, n in pairs])
    try:
        stat, p = wilcoxon(d)
    except ValueError:
        stat, p = np.nan, np.nan
    rows.append({"cell_type": g, "n_pairs": len(pairs),
                 "mean_normal": float(np.mean([P.loc[n, g] for _, n in pairs])),
                 "mean_adjacent": float(np.mean([P.loc[a, g] for a, _ in pairs])),
                 "delta": float(d.mean()), "wilcoxon_p": float(p) if p == p else np.nan})
res = pd.DataFrame(rows)
res["p_bonferroni"] = res["wilcoxon_p"] * len(res)
res.loc[res["p_bonferroni"] > 1, "p_bonferroni"] = 1
res = res.sort_values("wilcoxon_p")
res.to_csv(os.path.join(OUT, "atlas_paired_Adjacent_vs_Normal.csv"), index=False)
print(res.head(10).to_string(), flush=True)

# ---------------------------------------------------------------- 5. 与 BayesPrism 交叉验证
bp = pd.read_csv(os.path.join(ROOT, "02_analysis/results/enhancement16_deconvolution/bayesprism_cell_proportions.csv"), index_col=0)
# 映射: GSE115469 聚类名 -> 谱系 (与 BayesPrism 谱系对齐)
lineage_map = {}
for g in groups:
    gl = g.lower()
    if "hepatocyte" in gl:
        lineage_map[g] = "Hepatocyte"
    elif "inflammatory_macrophage" in gl or "non-inflammatory_macrophage" in gl:
        lineage_map[g] = "Macrophage_lineage"
    elif "lsec" in gl or "endothelial" in gl:
        lineage_map[g] = "Endothelial"
    elif "t_cells" in gl or "nk" in gl:
        lineage_map[g] = "T_NK"
    elif "b_cells" in gl or "plasma" in gl:
        lineage_map[g] = "B_Plasma"
    elif "cholangiocyte" in gl:
        lineage_map[g] = "Cholangiocyte"
    elif "stellate" in gl:
        lineage_map[g] = "Stellate"
    elif "erythroid" in gl:
        lineage_map[g] = "Erythroid"
    else:
        lineage_map[g] = "Other"
P_lin = P.T.groupby(lineage_map).sum().T

bp_lin = pd.DataFrame(index=bp.index)
bp_lin["Hepatocyte"] = bp["Hepatocyte"]
bp_lin["Macrophage_lineage"] = bp["Kupffer_cell"] + bp["Monocyte_derived_Mac"]
bp_lin["Endothelial"] = bp["Endothelial_LSEC"]
bp_lin["T_NK"] = bp["NK_cell"] + bp["T_CD8"] + bp["T_CD4"] + bp["Treg"]
bp_lin["B_Plasma"] = bp["B_cell"] + bp["Plasma_cell"]
bp_lin["Cholangiocyte"] = bp["Cholangiocyte"]

common_s = P_lin.index.intersection(bp_lin.index)
cross = []
for lin in [c for c in bp_lin.columns if c in P_lin.columns]:
    if lin == "Hepatocyte":
        a, b = P_lin.loc[common_s, lin], bp_lin.loc[common_s, lin]
    else:
        a, b = P_lin.loc[common_s, lin], bp_lin.loc[common_s, lin]
    r, p = spearmanr(a, b)
    cross.append({"lineage": lin, "spearman_rho": r, "p": p,
                  "atlas_delta": float(P_lin.loc[[f"Adjacent{x}" for x in pats if f"Adjacent{x}" in P_lin.index], lin].mean()
                                       - P_lin.loc[[f"Normal{x}" for x in pats if f"Normal{x}" in P_lin.index], lin].mean()),
                  "bayesprism_delta": float(bp_lin.loc[[f"Adjacent{x}" for x in pats if f"Adjacent{x}" in bp_lin.index], lin].mean()
                                            - bp_lin.loc[[f"Normal{x}" for x in pats if f"Normal{x}" in bp_lin.index], lin].mean())})
cross_df = pd.DataFrame(cross)
cross_df.to_csv(os.path.join(OUT, "cross_method_concordance_atlas_vs_bayesprism.csv"), index=False)
print("[5] cross-method concordance:\n", cross_df.to_string(), flush=True)

# 汇总
summary = {
    "reference": "GSE115469 (MacParland et al. 2018 Nat Commun), 8444 human liver cells, real scRNA-seq",
    "n_clusters": len(groups),
    "n_shared_genes": len(sel),
    "n_samples": Bg.shape[1],
    "method": "log2 CPM signature matrix + NNLS (CIBERSORT-style linear kernel)",
}
with open(os.path.join(OUT, "SUMMARY.md"), "w") as f:
    f.write("# enhancement38: 真实人肝单细胞图谱参考去卷积\n\n")
    for k, v in summary.items():
        f.write(f"- **{k}**: {v}\n")
    f.write("\n## 配对差异 (Adjacent − Normal, atlas 参考)\n\n")
    f.write(res.head(10).to_markdown(index=False))
    f.write("\n\n## 与 BayesPrism 跨方法一致性 (Spearman)\n\n")
    f.write(cross_df.to_markdown(index=False))
print("[done]", flush=True)
