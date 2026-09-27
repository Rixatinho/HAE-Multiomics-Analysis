#!/usr/bin/env python3
"""
Enhancement 32: Healthy-liver reference validation + hepatic zonation collapse analysis
=====================================================================================
目的（回应 CommsBio 审稿 R1.1/R2.2 对"远端正常肝"对照质量的质疑）：
  A. 健康肝参考验证（替代 GTEx：门户被网络代理拦截，改用 GSE135251 Govaere et al. 10 例
     健康肝活检 bulk RNA-seq）
     A1. 跨数据集 PCA：我们的 Normal(远端≥1cm)/Adjacent(病灶旁<1cm) 与健康肝对照的
         全局转录组定位
     A2. 全基因谱相关性：Normal/Adjacent 与健康肝均值的 Spearman 相关
     A3. 肝功能/药物代谢基因面板：CYP/FMO/ALB 等在 Normal vs 健康肝 vs Adjacent 的比较
     A4. Field-effect 筛查：Normal 相对健康肝的差异基因（field effect 候选）与
         Adjacent vs Normal 的 DEG 取交集，量化"远端正常肝"的受污染程度
  B. 肝小叶分区崩塌分析（GSE83990 Saito et al. LCM 分区转录组，3 例正常供体 Z1/Z2/Z3+WL）
     B1. 从 GSE83990 提取分区标志基因（Z1 门周 vs Z3 中央）
     B2. 计算每样本 zonation 评分（Z1 score / Z3 score / zonation balance）
     B3. 检验 perilesional(Adjacent) 是否出现分区特异性崩塌（如中央区 CYP2E1/GLUL
         特征丧失）——为代谢酶下降提供空间生物学解释

输入:
  data/processed/transcriptomics_logcpm_paired.csv        # 我们的 24 样本 logCPM
  external_data/GSE135251_controls/*.counts.txt.gz        # 10 例健康肝 counts
  external_data/GSE83990_Saito_LiverZonation_NormalizedReads.txt.gz

输出:
  results/enhancement32_healthy_reference/
    healthy_reference_pca.csv / pca_plot.pdf
    global_correlation.csv
    liver_function_panel.csv / panel_plot.pdf
    field_effect_genes.csv
    zonation_markers.csv / zonation_scores.csv / zonation_plot.pdf
    enhancement32_report.md

运行: /Users/rishat/miniforge3/envs/multiomics/bin/python scripts/enhancement32_healthy_reference.py
"""

import gzip
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ---------------------------------------------------------------- 路径
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(BASE, "data", "processed")
EXT = os.path.join(BASE, "external_data")
CTRL_DIR = os.path.join(EXT, "GSE135251_controls")
ZON_F = os.path.join(EXT, "GSE83990_Saito_LiverZonation_NormalizedReads.txt.gz")
OUT = os.path.join(BASE, "results", "enhancement32_healthy_reference")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- 工具
def read_gz_tsv(path):
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", index_col=0)


def cpm_from_counts(counts_df):
    lib = counts_df.sum(axis=0)
    return np.log2(counts_df.div(lib, axis=1) * 1e6 + 1)


def strip_version(x):
    return x.split(".")[0]


def clean_ensg(x):
    """ENSG00000198804_2 -> ENSG00000198804"""
    base = x.split(".")[0]
    if "_" in base:
        base = base.rsplit("_", 1)[0]
    return base


def load_gene_map():
    """Ensembl -> Symbol 映射，从 DESeq2 结果文件的 gene_name 列构建。"""
    mapping = {}
    import glob
    for f in glob.glob(os.path.join(BASE, "results", "*", "DESeq2_paired_results.csv")):
        try:
            df = pd.read_csv(f)
        except Exception:
            continue
        if "gene_id" in df.columns and "gene_name" in df.columns:
            for gid, gname in zip(df["gene_id"], df["gene_name"]):
                if isinstance(gname, str) and gname and gname != "NA":
                    mapping[clean_ensg(str(gid))] = gname
    return mapping


# ---------------------------------------------------------------- 1. 载入我们的数据
ours = pd.read_csv(os.path.join(PROC, "transcriptomics_logcpm_paired.csv"),
                   index_col=0)
ours.index = [clean_ensg(i) for i in ours.index]
# 合并同一 Ensembl 的重复行（取均值）
ours = ours.groupby(level=0).mean()
normal_cols = [c for c in ours.columns if c.startswith("Normal")]
adjacent_cols = [c for c in ours.columns if c.startswith("Adjacent")]
print(f"[1] our data: {ours.shape[0]} genes, "
      f"{len(normal_cols)} Normal + {len(adjacent_cols)} Adjacent")

gene_map = load_gene_map()
print(f"[1b] gene map entries: {len(gene_map)}")

# ---------------------------------------------------------------- 2. 载入健康对照
ctrl_files = sorted(os.listdir(CTRL_DIR))
ctrl_list, ctrl_names = [], []
for f in ctrl_files:
    df = read_gz_tsv(os.path.join(CTRL_DIR, f))
    df = df[~df.index.duplicated(keep="first")]
    ctrl_list.append(df.iloc[:, 0])
    ctrl_names.append(f.replace(".counts.txt.gz", "").split("_")[0])
ctrl_counts = pd.concat(ctrl_list, axis=1)
ctrl_counts.columns = ctrl_names
ctrl_counts.index = [strip_version(str(i)) for i in ctrl_counts.index]
ctrl_logcpm = cpm_from_counts(ctrl_counts)
ctrl_logcpm = ctrl_logcpm.groupby(level=0).mean()
print(f"[2] GSE135251 controls: {ctrl_logcpm.shape[1]} samples, "
      f"{ctrl_logcpm.shape[0]} genes")

# ---------------------------------------------------------------- 3. A1 跨数据集 PCA
shared = sorted(set(ours.index) & set(ctrl_logcpm.index))
mat = pd.concat([ours.loc[shared], ctrl_logcpm.loc[shared]], axis=1)
mat = mat.dropna(axis=0)
X = mat.T.values
Xc = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
pcs = (Xc @ Vt.T)[:, :4]
var_exp = (S ** 2 / (S ** 2).sum())[:4]

pca_df = pd.DataFrame(pcs, index=mat.columns,
                      columns=[f"PC{i+1}" for i in range(4)])
pca_df["group"] = ["Normal" if c.startswith("Normal") else
                   "Adjacent" if c.startswith("Adjacent") else "HealthyCtrl"
                   for c in mat.columns]
pca_df.to_csv(os.path.join(OUT, "healthy_reference_pca.csv"))

fig, ax = plt.subplots(figsize=(5.2, 4.4))
colors = {"Normal": "#2166AC", "Adjacent": "#B2182B", "HealthyCtrl": "#7F7F7F"}
for g in ["HealthyCtrl", "Normal", "Adjacent"]:
    sub = pca_df[pca_df.group == g]
    ax.scatter(sub.PC1, sub.PC2, c=colors[g], s=42 if g != "HealthyCtrl" else 30,
               alpha=0.85, label=f"{g} (n={len(sub)})", edgecolor="white",
               linewidth=0.5, zorder=3 if g != "HealthyCtrl" else 2)
ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)")
ax.set_title("Our cohort vs healthy liver reference (GSE135251)")
ax.legend(frameon=False, fontsize=8)
sns_style = dict(spines={"top": False, "right": False})
for s_ in ["top", "right"]:
    ax.spines[s_].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "pca_plot.pdf"))
plt.close()

# 质心距离: 我们的样本到健康对照质心的 PC1-PC2 欧氏距离（以对照内部距离的 MAD 标准化）
ctrl_pts = pca_df.loc[pca_df.group == "HealthyCtrl", ["PC1", "PC2"]].values
ctrl_cen = ctrl_pts.mean(axis=0)
ctrl_spread = np.median(np.linalg.norm(ctrl_pts - ctrl_cen, axis=1))
dist = {}
for g in ["Normal", "Adjacent"]:
    pts = pca_df.loc[pca_df.group == g, ["PC1", "PC2"]].values
    d = np.linalg.norm(pts - ctrl_cen, axis=1) / max(ctrl_spread, 1e-9)
    dist[g] = d
    print(f"[3] {g}: median centroid distance = {np.median(d):.2f} "
          f"(ctrl MAD units)")
dist_df = pd.DataFrame({g: pd.Series(v) for g, v in dist.items()})
dist_df.to_csv(os.path.join(OUT, "centroid_distance.csv"), index=False)

# ---------------------------------------------------------------- 4. A2 全局相关
ctrl_mean = ctrl_logcpm.loc[shared].mean(axis=1)
corr_rows = []
for c in ours.columns:
    rho = stats.spearmanr(ours.loc[shared, c], ctrl_mean)[0]
    grp = "Normal" if c.startswith("Normal") else "Adjacent"
    corr_rows.append({"sample": c, "group": grp, "spearman_vs_healthy": rho})
corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(os.path.join(OUT, "global_correlation.csv"), index=False)
med_n = corr_df[corr_df.group == "Normal"].spearman_vs_healthy.median()
med_a = corr_df[corr_df.group == "Adjacent"].spearman_vs_healthy.median()
u, p = stats.mannwhitneyu(corr_df[corr_df.group == "Normal"].spearman_vs_healthy,
                          corr_df[corr_df.group == "Adjacent"].spearman_vs_healthy)
print(f"[4] Spearman vs healthy mean: Normal {med_n:.3f}, "
      f"Adjacent {med_a:.3f}, MW p={p:.3f}")

# ---------------------------------------------------------------- 5. A3 肝功能/药代基因面板
panel_genes = {
    "CYP3A4": "drug metabolism", "CYP3A5": "drug metabolism",
    "CYP2E1": "drug metabolism", "CYP2C9": "drug metabolism",
    "CYP2C19": "drug metabolism", "CYP1A2": "drug metabolism",
    "CYP2D6": "drug metabolism", "FMO3": "drug metabolism",
    "ALB": "liver function", "TTR": "liver function",
    "AFP": "liver function", "ASGR1": "liver function",
    "PCK1": "gluconeogenesis", "G6PC": "gluconeogenesis",
    "CPS1": "urea cycle", "OTC": "urea cycle",
    "GLUL": "ammonia detox", "AADAT": "amino acid metabolism",
}
# ENSG 查找: 优先用我们的 DESeq2 gene_name 映射
inv_map = {}
for ens, sym in gene_map.items():
    inv_map.setdefault(sym, ens)
panel_rows = []
ctrl_idx = ctrl_logcpm.loc[shared].index
for sym, func in panel_genes.items():
    ens = inv_map.get(sym)
    if ens is None or ens not in ours.index or ens not in ctrl_idx:
        # 尝试直接用 symbol 在 GSE83990/无法定位则跳过
        panel_rows.append({"gene": sym, "function": func,
                           "Normal_mean": np.nan, "Adjacent_mean": np.nan,
                           "Healthy_mean": np.nan,
                           "log2FC_Adj_vs_Nor": np.nan,
                           "log2FC_Nor_vs_Healthy": np.nan, "p_Nor_vs_Healthy": np.nan,
                           "note": "not mapped"})
        continue
    nor = ours.loc[ens, normal_cols].mean()
    adj = ours.loc[ens, adjacent_cols].mean()
    hea = ctrl_logcpm.loc[ens].mean()
    # Normal vs Healthy: 每样本秩检验
    x = ours.loc[ens, normal_cols].values
    y = ctrl_logcpm.loc[ens].values
    try:
        pv = stats.mannwhitneyu(x, y)[1]
    except ValueError:
        pv = np.nan
    panel_rows.append({
        "gene": sym, "function": func,
        "Normal_mean": round(nor, 3), "Adjacent_mean": round(adj, 3),
        "Healthy_mean": round(hea, 3),
        "log2FC_Adj_vs_Nor": round(adj - nor, 3),
        "log2FC_Nor_vs_Healthy": round(nor - hea, 3),
        "p_Nor_vs_Healthy": pv,
    })
panel_df = pd.DataFrame(panel_rows)
panel_df.to_csv(os.path.join(OUT, "liver_function_panel.csv"), index=False)

# 面板热图式 dot plot
sub = panel_df.dropna(subset=["Normal_mean"]).copy()
if len(sub):
    fig, ax = plt.subplots(figsize=(6.4, 0.34 * len(sub) + 1.2))
    y = np.arange(len(sub))[::-1]
    w = 0.27
    ax.barh(y + w, sub.Healthy_mean, height=w, color="#7F7F7F", label="Healthy (GSE135251)")
    ax.barh(y, sub.Normal_mean, height=w, color="#2166AC", label="Normal (ours)")
    ax.barh(y - w, sub.Adjacent_mean, height=w, color="#B2182B", label="Adjacent (ours)")
    ax.set_yticks(y)
    ax.set_yticklabels(sub.gene, fontsize=8)
    ax.set_xlabel("mean log2(CPM+1)")
    ax.set_title("Hepatic function & drug-metabolism genes", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "panel_plot.pdf"))
    plt.close()

# ---------------------------------------------------------------- 6. A4 Field effect 筛查
# Normal vs Healthy (Mann-Whitney per gene, FDR) 与 Adjacent vs Normal 已知 DEG 求交集
from statsmodels.stats.multitest import multipletests  # noqa: E402

fe_rows = []
for ens in shared:
    if ens not in ours.index or ens not in ctrl_logcpm.index:
        continue
    x = ours.loc[ens, normal_cols].values
    y = ctrl_logcpm.loc[ens].values
    try:
        u2, pv = stats.mannwhitneyu(x, y)
    except ValueError:
        continue
    fe_rows.append({"ensembl": ens, "gene": gene_map.get(ens, ens),
                    "log2FC_Nor_vs_Healthy": float(np.mean(x) - np.mean(y)),
                    "p": pv})
fe_df = pd.DataFrame(fe_rows)
fe_df["fdr"] = multipletests(fe_df["p"], method="fdr_bh")[1]
fe_sig = fe_df[(fe_df.fdr < 0.05) & (fe_df.log2FC_Nor_vs_Healthy.abs() > 1)].copy()
fe_sig.sort_values("fdr", inplace=True)
fe_sig.to_csv(os.path.join(OUT, "field_effect_genes.csv"), index=False)
print(f"[6] field-effect genes (Normal vs healthy, FDR<0.05, |log2FC|>1): "
      f"{len(fe_sig)} / {len(fe_df)} tested")

# 与 Adjacent vs Normal DEG 交集（主分析: phase1_diff, limma）
deg_path = os.path.join(BASE, "results", "phase1_diff",
                        "DEGs_Adjacent_vs_Normal.csv")
if os.path.exists(deg_path):
    adj_deg = pd.read_csv(deg_path)
    adj_deg["ensembl"] = adj_deg["gene_id"].astype(str).map(clean_ensg)
    adj_sig = set(adj_deg.loc[(adj_deg["P.Value"] < 0.05) &
                              (adj_deg["logFC"].abs() > 1), "ensembl"])
    tested_ens = set(fe_df.ensembl) & set(adj_deg.ensembl)
    adj_sig = adj_sig & tested_ens
    fe_sig_ens = set(fe_sig.ensembl) & tested_ens
    overlap = fe_sig_ens & adj_sig
    ov_df = fe_sig[fe_sig.ensembl.isin(overlap)]
    ov_df.to_csv(os.path.join(OUT, "field_effect_DEG_overlap.csv"), index=False)
    # 超几何检验
    M = len(tested_ens)
    n = len(adj_sig)
    N = len(fe_sig_ens)
    k = len(overlap)
    from scipy.stats import hypergeom  # noqa: E402
    pv_hyper = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
    print(f"[6b] internal DEGs: {len(adj_sig)}; overlap with field-effect "
          f"genes: {k}; hypergeometric p={pv_hyper:.2e}")

# ---------------------------------------------------------------- 7. B 分区崩塌分析
with gzip.open(ZON_F, "rt") as fh:
    zon = pd.read_csv(fh, sep="\t", index_col=1)  # index = gene_id (ENSG)
# 列 LM-A..LM-L 对应 3 个供体 × (WL,Z1,Z2,Z3)
# soft 元数据: LM1-WL/GSM2224921 ... 顺序 = LM-A..LM-L
zone_order = []
for d in ["LM1", "LM2", "LM3"]:
    for z in ["WL", "Z1", "Z2", "Z3"]:
        zone_order.append((d, z))
zon_lm_cols = [c for c in zon.columns if c.startswith("LM-")]
zone_map = {c: f"{d}-{z}" for c, (d, z) in zip(zon_lm_cols, zone_order)}
assert len(zon_lm_cols) == 12, f"expected 12 LM columns, got {zon_lm_cols}"
zon = zon.rename(columns=zone_map)
z1_cols = [c for c in zon.columns if c.endswith("-Z1")]
z3_cols = [c for c in zon.columns if c.endswith("-Z3")]
wl_cols = [c for c in zon.columns if c.endswith("-WL")]
zon = zon[~zon.index.duplicated(keep="first")]
zon_g = zon.groupby(level=0).mean(numeric_only=True)
zon_g.index = [strip_version(str(i)) for i in zon_g.index]
zon_g = zon_g.groupby(level=0).mean(numeric_only=True)

# 分区标志基因: Z1 vs Z3 差异 (每组 n=3, 用 fold change 排序)
z1_mean = zon_g[z1_cols].mean(axis=1)
z3_mean = zon_g[z3_cols].mean(axis=1)
zon_lfc = z1_mean - z3_mean
zon_markers = pd.DataFrame({
    "ensembl": zon_g.index,
    "gene": [gene_map.get(e, "") for e in zon_g.index],
    "Z1_mean": z1_mean.values, "Z3_mean": z3_mean.values,
    "log2FC_Z1_vs_Z3": zon_lfc.values,
})
zon_markers = zon_markers[(zon_markers.Z1_mean > 2) | (zon_markers.Z3_mean > 2)]
z1_top = zon_markers.sort_values("log2FC_Z1_vs_Z3", ascending=False).head(60)
z3_top = zon_markers.sort_values("log2FC_Z1_vs_Z3").head(60)
z1_set = set(z1_top.ensembl)
z3_set = set(z3_top.ensembl)
print(f"[7] zonation markers: {len(z1_set)} Z1(periportal), {len(z3_set)} Z3(pericentral)")
z1_top.to_csv(os.path.join(OUT, "zonation_markers_Z1.csv"), index=False)
z3_top.to_csv(os.path.join(OUT, "zonation_markers_Z3.csv"), index=False)

# 每样本 Z1/Z3 评分（相对本数据集所有样本的 mean-centered）
def zone_score(expr_df, gene_set):
    gs = [g for g in gene_set if g in expr_df.index]
    sub = expr_df.loc[gs]
    return sub.mean(axis=0) - expr_df.mean(axis=0)

# 在我们的数据 + GSE83990 WL 上统一计算（用我们的矩阵尺度）
z1_score = zone_score(ours, z1_set)
z3_score = zone_score(ours, z3_set)
zs = pd.DataFrame({
    "sample": ours.columns,
    "group": ["Normal" if c.startswith("Normal") else "Adjacent" for c in ours.columns],
    "Z1_score": z1_score.values,
    "Z3_score": z3_score.values,
    "zonation_balance": (z1_score - z3_score).values,
})
zs.to_csv(os.path.join(OUT, "zonation_scores.csv"), index=False)

u1, p1 = stats.mannwhitneyu(zs[zs.group == "Normal"].Z1_score,
                            zs[zs.group == "Adjacent"].Z1_score)
u3, p3 = stats.mannwhitneyu(zs[zs.group == "Normal"].Z3_score,
                            zs[zs.group == "Adjacent"].Z3_score)
print(f"[7b] Z1 score Normal vs Adjacent p={p1:.4f}; Z3 score p={p3:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
for ax, (col, ttl, pv) in zip(axes, [
        ("Z1_score", "Periportal (Z1) score", p1),
        ("Z3_score", "Pericentral (Z3) score", p3)]):
    data = [zs[zs.group == "Normal"][col], zs[zs.group == "Adjacent"][col]]
    bp = ax.boxplot(data, tick_labels=["Normal", "Adjacent"], widths=0.55,
                    patch_artist=True)
    for patch, c in zip(bp["boxes"], ["#2166AC", "#B2182B"]):
        patch.set_facecolor(c); patch.set_alpha(0.55)
    for i, d_ in enumerate(data):
        ax.scatter(np.random.normal(i + 1, 0.05, len(d_)), d_, s=14,
                   color="black", alpha=0.6, zorder=3)
    ax.set_title(f"{ttl}\nMW p={pv:.3f}", fontsize=9)
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "zonation_plot.pdf"))
plt.close()

# ---------------------------------------------------------------- 8. 报告
report = []
report.append("# Enhancement 32: Healthy-liver reference & zonation collapse\n")
report.append(f"Date: 2026-08-22\n")
report.append("## A. Healthy reference (GSE135251, n=10 controls; GTEx blocked by network)\n")
report.append(f"- Shared genes for PCA: {len(shared)}\n")
report.append(f"- PC variance explained: "
              f"{', '.join(f'{v*100:.1f}%' for v in var_exp)}\n")
report.append(f"- Centroid distance (ctrl MAD units): Normal median "
              f"{np.median(dist['Normal']):.2f}, Adjacent median "
              f"{np.median(dist['Adjacent']):.2f}\n")
report.append(f"- Global Spearman vs healthy mean: Normal {med_n:.3f}, "
              f"Adjacent {med_a:.3f} (MW p={p:.3f})\n")
report.append(f"- Field-effect genes (Normal vs healthy, FDR<0.05 & |lfc|>1): "
              f"{len(fe_sig)} of {len(fe_df)}\n")
report.append(f"- Liver-function panel genes with Normal-vs-Healthy shift "
              f"(|lfc|>0.5): "
              f"{(panel_df.dropna(subset=['Normal_mean']).log2FC_Nor_vs_Healthy.abs() > 0.5).sum()}\n")
report.append("## B. Zonation collapse (GSE83990, LCM zones)\n")
report.append(f"- Z1 (periportal) markers: {len(z1_set)}; Z3 (pericentral) markers: {len(z3_set)}\n")
report.append(f"- Z1 score Normal vs Adjacent: MW p={p1:.4f}\n")
report.append(f"- Z3 score Normal vs Adjacent: MW p={p3:.4f}\n")
with open(os.path.join(OUT, "enhancement32_report.md"), "w") as fh:
    fh.write("\n".join(report))

print("\n[DONE] outputs ->", OUT)
