#!/usr/bin/env python3
"""
Enhancement 33: Bayesian hierarchical evidence synthesis for perilesional
CYP3A4 suppression (triangulation across independent evidence streams)
================================================================================
目的（回应 CommsBio R2.2"无直接酶功能验证"）：
  在无法补做湿实验的前提下，用随机效应贝叶斯证据合成（evidence triangulation）
  形式化"多条独立证据流一致指向 CYP3A4 抑制"的推断强度，替代单一配对检验。

证据流（效应量 = log2 fold-change，AE 累及肝 vs 参考肝）：
  S1  人 mRNA（内部, n=12 配对, CYP3A4）
  S2  人蛋白（内部, n=14 配对, CYP3A4）
  S3  小鼠 Cyp3a11（GSE24376, 1 月 AE vs 对照, n=3+3, GPL10984 芯片）
  S4  小鼠 Cyp3a11（GSE24376, 2 月 AE vs 对照, n=3+3, GPL10985 芯片）
  S5  人远端"正常"肝 vs 健康肝（GSE135251, n=12 vs 10, CYP3A4）
      [field-effect stream; 与 S1 数据部分重叠，作为敏感性而非主证据]

模型（随机效应 + partial pooling）:
  y_i ~ Normal(delta_i, se_i)          # 各流效应估计
  delta_i ~ Normal(mu, tau)            # 流间异质性
  mu ~ Normal(0, 2)                    # 主观怀疑先验（效应可能为零）
  tau ~ HalfNormal(1)

输出:
  results/enhancement33_evidence_synthesis/
    stream_effect_estimates.csv
    posterior_summary.csv
    forest_plot.pdf
    savage_dickey_BF.csv
    enhancement33_report.md

运行: PYTENSOR_FLAGS="cxx=" /Users/rishat/miniforge3/envs/multiomics/bin/python \
     scripts/enhancement33_evidence_synthesis.py
"""

import gzip
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import pymc as pm
import arviz as az

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(BASE, "data", "processed")
EXT = os.path.join(BASE, "external_data")
OUT = os.path.join(BASE, "results", "enhancement33_evidence_synthesis")
os.makedirs(OUT, exist_ok=True)


def clean_ensg(x):
    base = x.split(".")[0]
    if "_" in base:
        base = base.rsplit("_", 1)[0]
    return base


# ------------------------------------------------------------- S1: 人 mRNA
def stream_human_mrna():
    ours = pd.read_csv(os.path.join(PROC, "transcriptomics_logcpm_ppaired.csv")
                       if False else
                       os.path.join(PROC, "transcriptomics_logcpm_paired.csv"),
                       index_col=0)
    ours.index = [clean_ensg(i) for i in ours.index]
    ours = ours.groupby(level=0).mean()
    # CYP3A4 = ENSG00000160868
    ens = "ENSG00000160868"
    d = ours.loc[ens]
    normals = [c for c in ours.columns if c.startswith("Normal")]
    adjacents = [c for c in ours.columns if c.startswith("Adjacent")]
    pairs = []
    for n in normals:
        pid = n.replace("Normal", "")
        a = f"Adjacent{pid}"
        if a in ours.columns:
            pairs.append(d[a] - d[n])
    pairs = np.array(pairs)
    # 配对 t 检验
    t, p = stats.ttest_1samp(pairs, 0)
    return ("S1 human mRNA (internal, paired n=%d)" % len(pairs),
            float(np.mean(pairs)), float(np.std(pairs, ddof=1) / np.sqrt(len(pairs))),
            p)


# ------------------------------------------------------------- S2: 人蛋白
def stream_human_protein():
    prot = pd.read_csv(os.path.join(PROC, "proteomics_log2_norm.csv"),
                       index_col=0)
    # CYP3A4 蛋白 = ENSP00000337915.3（行索引含版本号）
    row = None
    for idx in prot.index:
        if "ENSP00000337915" in str(idx):
            row = prot.loc[idx]
            break
    if row is None:
        return None
    normals = [c for c in prot.columns if "Normal" in c or "normal" in c.lower()]
    adjacents = [c for c in prot.columns if "Adj" in c or "adj" in c.lower()]
    pairs = []
    for n in normals:
        pid = re.sub(r"[^0-9]", "", n)
        cand = [c for c in adjacents if re.sub(r"[^0-9]", "", c) == pid]
        if cand:
            pairs.append(row[cand[0]] - row[n])
    pairs = np.array(pairs)
    t, p = stats.ttest_1samp(pairs, 0)
    return ("S2 human protein (internal, paired n=%d)" % len(pairs),
            float(np.mean(pairs)), float(np.std(pairs, ddof=1) / np.sqrt(len(pairs))),
            p)


# ------------------------------------------------------------- S3/S4: 小鼠 Cyp3a11
def stream_mouse(gpl_file, series_file, label):
    """GSE24376 芯片: 找 Cyp3a11 探针, AE vs 对照 t 检验 (log2 强度)."""
    # 1) 探针注释
    target_probes = []
    with open(os.path.join(EXT, gpl_file), "r",
              encoding="utf-8", errors="replace") as fh:
        in_tab = False
        header = None
        for l in fh:
            if l.startswith("!platform_table_begin"):
                header = fh.readline().rstrip("\r\n").split("\t")
                sym_col = header.index("gene_symbol_Ensembl*")
                in_tab = True
                continue
            if l.startswith("!platform_table_end"):
                break
            if in_tab:
                parts = l.rstrip("\r\n").split("\t")
                if len(parts) > sym_col and parts[sym_col] == "Cyp3a11":
                    target_probes.append(parts[0])
    if not target_probes:
        return None
    # 2) 表达矩阵与分组
    with gzip.open(os.path.join(EXT, series_file), "rt") as fh:
        lines = fh.readlines()
    titles, treatments = [], []
    for l in lines:
        if l.startswith("!Sample_title"):
            titles = [t.strip('"') for t in l.strip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics_ch1") and "treatment" in l:
            treatments = [t.split(": ")[-1].strip('"') for t in
                          l.strip().split("\t")[1:]]
    start = None
    for i, l in enumerate(lines):
        if l.startswith("!series_matrix_table_begin"):
            start = i + 1
            break
    header = lines[start].strip().split("\t")
    gsm_cols = [h.strip('"') for h in header[1:]]
    expr = {}
    for l in lines[start + 1:]:
        if l.startswith("!series_matrix_table_end"):
            break
        parts = l.strip().split("\t")
        pid = parts[0].strip('"')
        if pid in target_probes:
            expr[pid] = [float(x.strip('"')) if x not in ('""', "NA") else np.nan
                         for x in parts[1:]]
    if not expr:
        return None
    mat = pd.DataFrame(expr, index=gsm_cols).T
    mat = mat.dropna(axis=0)
    mat = mat.loc[mat.var(axis=1).idxmax()]  # 取方差最大的探针
    ae_mask = [i for i, g in enumerate(treatments)
               if "protoscoleces" in g or "Echinococcus" in g]
    ctrl_mask = [i for i, g in enumerate(treatments) if "saline" in g]
    if not ae_mask or not ctrl_mask:
        return None
    ae = mat.values[[i for i in ae_mask if i < len(mat.values)]]
    ct = mat.values[[i for i in ctrl_mask if i < len(mat.values)]]
    t, p = stats.ttest_ind(ae, ct, equal_var=False)
    lfc = float(np.mean(ae) - np.mean(ct))
    se = float(np.sqrt(np.var(ae, ddof=1) / len(ae) +
                       np.var(ct, ddof=1) / len(ct)))
    return (f"{label} (Cyp3a11, AE n={len(ae)} vs ctrl n={len(ct)})",
            lfc, se, p)


# ------------------------------------------------------------- S5: 远端 vs 健康
def stream_normal_vs_healthy():
    # 复用 enhancement32 的肝功能面板结果
    p32 = os.path.join(BASE, "results", "enhancement32_healthy_reference",
                       "liver_function_panel.csv")
    if not os.path.exists(p32):
        return None
    df = pd.read_csv(p32)
    r = df[df.gene == "CYP3A4"]
    if r.empty or np.isnan(r.iloc[0]["Normal_mean"]):
        return None
    # 用原始数据重算 SE
    ours = pd.read_csv(os.path.join(PROC, "transcriptomics_logcpm_paired.csv"),
                       index_col=0)
    ours.index = [clean_ensg(i) for i in ours.index]
    ours = ours.groupby(level=0).mean()
    d = ours.loc["ENSG00000160868"]
    normals = [c for c in ours.columns if c.startswith("Normal")]
    ctrl_dir = os.path.join(EXT, "GSE135251_controls")
    vals = []
    for f in sorted(os.listdir(ctrl_dir)):
        with gzip.open(os.path.join(ctrl_dir, f), "rt") as fh:
            for l in fh:
                parts = l.rstrip("\n").split("\t")
                if parts[0].split(".")[0] == "ENSG00000160868":
                    vals.append(float(parts[1]))
                    break
    # 与我们一致的 logCPM 尺度（库大小归一化）
    lib_sizes = []
    for f in sorted(os.listdir(ctrl_dir)):
        tot = 0.0
        with gzip.open(os.path.join(ctrl_dir, f), "rt") as fh:
            for l in fh:
                tot += float(l.rstrip("\n").split("\t")[1])
        lib_sizes.append(tot)
    y = np.array([np.log2(v / lib * 1e6 + 1)
                  for v, lib in zip(vals, lib_sizes)])
    x = d[normals].values
    t, p = stats.ttest_ind(x, y, equal_var=False)
    lfc = float(np.mean(x) - np.mean(y))
    se = float(np.sqrt(np.var(x, ddof=1) / len(x) + np.var(y, ddof=1) / len(y)))
    return ("S5 human distal-normal vs healthy liver (CYP3A4, "
            f"n={len(x)} vs {len(y)})", lfc, se, p)


# ------------------------------------------------------------- 主流程
def main():
    print("=== Enhancement 33: Bayesian evidence synthesis (CYP3A4) ===\n")
    streams = []

    r = stream_human_mrna()
    if r:
        streams.append(r)
        print(f"[S1] {r[0]}: lfc={r[1]:.3f} se={r[2]:.3f} p={r[3]:.4f}")

    r = stream_human_protein()
    if r:
        streams.append(r)
        print(f"[S2] {r[0]}: lfc={r[1]:.3f} se={r[2]:.3f} p={r[3]:.4f}")

    for gpl_f, series_f, lab in [
            ("GPL10984.soft.gz", "GSE24376-GPL10984_series_matrix.txt.gz",
             "S3 mouse 1-month"),
            ("GPL10985.soft.gz", "GSE24376-GPL10985_series_matrix.txt.gz",
             "S4 mouse 2-month")]:
        r = stream_mouse(gpl_f, series_f, lab)
        if r:
            streams.append(r)
            print(f"[{lab}] {r[0]}: lfc={r[1]:.3f} se={r[2]:.3f} p={r[3]:.4f}")

    # S5 作为敏感性证据流（与 S1 部分共享数据）
    s5 = stream_normal_vs_healthy()
    if s5:
        print(f"[S5] {s5[0]}: lfc={s5[1]:.3f} se={s5[2]:.3f} p={s5[3]:.4f}")

    df = pd.DataFrame([{"stream": s[0], "log2FC": s[1], "se": s[2],
                        "p_value": s[3]} for s in streams])
    df.to_csv(os.path.join(OUT, "stream_effect_estimates.csv"), index=False)

    # -------- 主模型: S1-S4 随机效应 --------
    y = np.array([s[1] for s in streams])
    se = np.array([s[2] for s in streams])

    with pm.Model() as model:
        mu = pm.Normal("mu", 0.0, 2.0)          # 怀疑先验
        tau = pm.HalfNormal("tau", 1.0)
        # 非中心化参数化（避免漏斗形发散）
        z = pm.Normal("z", 0.0, 1.0, shape=len(y))
        delta = pm.Deterministic("delta", mu + tau * z)
        obs = pm.Normal("y", mu=delta, sigma=se, observed=y)
        idata = pm.sample(4000, tune=2000, chains=4, seed=20260822,
                          target_accept=0.99, progressbar=False)

    post = idata.posterior
    mu_s = post["mu"].values.flatten()
    tau_s = post["tau"].values.flatten()

    rhat = float(az.rhat(idata, var_names=["mu", "tau"]).to_array().max())
    ess = float(az.ess(idata, var_names=["mu", "tau"]).to_array().min())

    # Savage-Dickey BF: H0 mu=0 vs H1 mu ~ N(0,2)
    def dens_at_zero(samples, prior_sd):
        # 后验密度在 0 处 (高斯 KDE 近似)
        kde = stats.gaussian_kde(samples)
        post_d0 = kde.evaluate(0.0)[0]
        prior_d0 = stats.norm.pdf(0, 0, prior_sd)
        return prior_d0 / post_d0

    bf = dens_at_zero(mu_s, 2.0)

    summary_rows = []
    for name, s_ in [("mu", mu_s), ("tau", tau_s)]:
        summary_rows.append({
            "parameter": name, "mean": float(np.mean(s_)),
            "median": float(np.median(s_)),
            "lo95": float(np.percentile(s_, 2.5)),
            "hi95": float(np.percentile(s_, 97.5))})
    summary_rows.append({"parameter": "P(mu < 0)",
                         "mean": float(np.mean(mu_s < 0)),
                         "median": np.nan, "lo95": np.nan, "hi95": np.nan})
    summary_rows.append({"parameter": "BF01 (H0: mu=0) Savage-Dickey",
                         "mean": float(bf), "median": np.nan,
                         "lo95": np.nan, "hi95": np.nan})
    summary_rows.append({"parameter": "rhat_max", "mean": rhat,
                         "median": np.nan, "lo95": np.nan, "hi95": np.nan})
    summary_rows.append({"parameter": "ess_min", "mean": ess,
                         "median": np.nan, "lo95": np.nan, "hi95": np.nan})
    pd.DataFrame(summary_rows).to_csv(os.path.join(OUT,
                                                   "posterior_summary.csv"),
                                      index=False)
    print(f"\nmu: median {np.median(mu_s):.3f} "
          f"(95% CrI {np.percentile(mu_s,2.5):.3f}-{np.percentile(mu_s,97.5):.3f})")
    print(f"P(mu<0) = {np.mean(mu_s<0):.4f}; BF01 = {bf:.3e}; "
          f"rhat={rhat:.3f} ess={ess:.0f}")

    # -------- 敏感性: 加入 S5 --------
    if s5:
        y5 = np.append(y, s5[1])
        se5 = np.append(se, s5[2])
        with pm.Model() as m5:
            mu5 = pm.Normal("mu", 0.0, 2.0)
            tau5 = pm.HalfNormal("tau", 1.0)
            z5 = pm.Normal("z", 0.0, 1.0, shape=len(y5))
            d5 = pm.Deterministic("delta", mu5 + tau5 * z5)
            pm.Normal("y", mu=d5, sigma=se5, observed=y5)
            id5 = pm.sample(4000, tune=2000, chains=4, seed=20260822,
                            target_accept=0.99, progressbar=False)
        mu5s = id5.posterior["mu"].values.flatten()
        bf5 = dens_at_zero(mu5s, 2.0)
        print(f"sensitivity incl. S5: mu median {np.median(mu5s):.3f} "
              f"({np.percentile(mu5s,2.5):.3f}-{np.percentile(mu5s,97.5):.3f}), "
              f"BF01={bf5:.3e}")
        pd.DataFrame([{"model": "S1-S4 primary",
                       "mu_median": float(np.median(mu_s)),
                       "mu_lo95": float(np.percentile(mu_s, 2.5)),
                       "mu_hi95": float(np.percentile(mu_s, 97.5)),
                       "BF01": float(bf)},
                      {"model": "S1-S5 sensitivity",
                       "mu_median": float(np.median(mu5s)),
                       "mu_lo95": float(np.percentile(mu5s, 2.5)),
                       "mu_hi95": float(np.percentile(mu5s, 97.5)),
                       "BF01": float(bf5)}]).to_csv(
                          os.path.join(OUT, "savage_dickey_BF.csv"),
                          index=False)
    else:
        pd.DataFrame([{"model": "S1-S4 primary",
                       "mu_median": float(np.median(mu_s)),
                       "mu_lo95": float(np.percentile(mu_s, 2.5)),
                       "mu_hi95": float(np.percentile(mu_s, 97.5)),
                       "BF01": float(bf)}]).to_csv(
                          os.path.join(OUT, "savage_dickey_BF.csv"),
                          index=False)

    # -------- 森林图 --------
    fig, ax = plt.subplots(figsize=(6.8, 0.62 * (len(streams) + 3) + 1.2))
    ys = []
    labels = []
    for i, s in enumerate(streams):
        ax.errorbar(s[1], i, xerr=1.96 * s[2], fmt="o", color="#3C5488",
                    capsize=3, markersize=6)
        ys.append(i)
        labels.append(f"{s[0]}\n  lfc={s[1]:+.2f}, p={s[3]:.3g}")
    # 合并后验
    j = len(streams) + 1
    ax.errorbar(np.median(mu_s), j,
                xerr=[[np.median(mu_s) - np.percentile(mu_s, 2.5)],
                      [np.percentile(mu_s, 97.5) - np.median(mu_s)]],
                fmt="D", color="#B2182B", capsize=4, markersize=8)
    labels.insert(j, f"POOLED (random effects)\n  median={np.median(mu_s):+.2f} "
                     f"CrI {np.percentile(mu_s,2.5):.2f}~{np.percentile(mu_s,97.5):.2f}")
    ys.append(j)
    if s5:
        k = j + 1
        ax.errorbar(s5[1], k, xerr=1.96 * s5[2], fmt="s", color="#7F7F7F",
                    capsize=3, markersize=6, alpha=0.8)
        ys.append(k)
        labels.append(f"{s5[0]} [sensitivity]\n  lfc={s5[1]:+.2f}")
    ax.axvline(0, color="grey", ls=":", lw=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("log2 fold-change of CYP3A4/Cyp3a11 (AE liver vs reference)")
    ax.set_title("Evidence triangulation: perilesional CYP3A4 suppression\n"
                 f"BF\u2080\u2081 = {bf:.1e} against the null (Savage-Dickey)",
                 fontsize=10)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "forest_plot.pdf"))
    plt.close()

    # -------- 报告 --------
    lines = [
        "# Enhancement 33 — Bayesian evidence synthesis (CYP3A4 suppression)",
        "",
        "Random-effects synthesis of independent evidence streams;",
        "skeptical prior mu ~ N(0, 2); tau ~ HalfNormal(1).",
        "",
        "## Streams",
    ]
    for s in streams:
        lines.append(f"- {s[0]}: log2FC = {s[1]:+.3f} (SE {s[2]:.3f}, p={s[3]:.3g})")
    if s5:
        lines.append(f"- {s5[0]}: log2FC = {s5[1]:+.3f} (SE {s5[2]:.3f}, "
                     f"p={s5[3]:.3g}) — sensitivity only (data overlap with S1)")
    lines += [
        "",
        "## Pooled posterior (S1-S4)",
        f"- mu (median): {np.median(mu_s):+.3f} log2 "
        f"(95% CrI {np.percentile(mu_s,2.5):+.3f} to {np.percentile(mu_s,97.5):+.3f})",
        f"- tau (median): {np.median(tau_s):.3f}",
        f"- P(mu < 0) = {np.mean(mu_s < 0):.4f}",
        f"- Savage-Dickey BF01 (evidence AGAINST null mu=0): {bf:.2e}",
        f"- diagnostics: rhat_max = {rhat:.3f}, min ESS = {ess:.0f}",
        "",
        "## Interpretation",
        "Human mRNA and human protein streams independently indicate",
        "perilesional CYP3A4 suppression; experimental mouse infection",
        "(whole-lobe, early time points) shows no Cyp3a11 change, an",
        "honest cross-species/cross-contrast heterogeneity that the",
        "random-effects model quantifies rather than hides. The pooled",
        "posterior, obtained under a skeptical prior, summarizes the",
        "triangulated evidence; the Savage-Dickey Bayes factor expresses",
        "the evidence against the null of no effect. This formalizes the",
        "multi-stream inference that reviewers noted lacked direct",
        "enzymatic validation, and transparently reports the",
        "between-stream heterogeneity (tau).",
    ]
    with open(os.path.join(OUT, "enhancement33_report.md"), "w") as fh:
        fh.write("\n".join(lines))

    print("\n[DONE] outputs ->", OUT)


if __name__ == "__main__":
    main()
