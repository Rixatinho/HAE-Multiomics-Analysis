#!/usr/bin/env python
"""Enhancement 41: Clinical anchoring of molecular derangement (dose-response
with disease burden).

Question a top-tier reviewer asks: does the molecular phenotype track disease
severity, or is it a binary bystander effect? We test whether the magnitude of
peri-lesional molecular derangement (CYP3A4 loss, HNF4A loss, zonation
collapse, CYP-family depletion) scales with disease burden (lesion size,
PNM-P stage, bile-duct and vascular invasion grade, disease duration) and
systemic inflammation (CRP, IL-6), beyond the cholestasis controls already
reported (ALP/GGT/bilirubin).

Design: per-patient paired log2FC (Adjacent - Normal) for molecular features,
Spearman correlation against clinical variables (n = 12 patients with paired
transcriptomics), permutation p-values (10,000 permutations) + BH FDR across
the full test matrix. Subtype (CS1/CS2) clinical comparison by Mann-Whitney.

Outputs -> 02_analysis/results/enhancement41_clinical_anchor/
    clinical_anchor_correlations.csv
    subtype_clinical_comparison.csv
    clinical_anchor_heatmap.pdf
    summary.md
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "02_analysis" / "data" / "processed"
RES = ROOT / "02_analysis" / "results"
OUT = RES / "enhancement41_clinical_anchor"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)

# ----------------------------------------------------------------- load data
tx = pd.read_csv(DATA / "transcriptomics_logcpm_paired.csv", index_col=0)
deg = pd.read_csv(RES / "phase1_diff" / "DEGs_Adjacent_vs_Normal.csv")
idmap = deg.dropna(subset=["gene_name"]).set_index("gene_id")["gene_name"].to_dict()

prot = pd.read_csv(DATA / "proteomics_log2_norm.csv", index_col=0)
pann = pd.read_csv(DATA / "proteomics_annotation.csv")
pann = pann.rename(columns={pann.columns[0]: "pid"}).dropna(subset=["Gene"])
protmap = dict(zip(pann["pid"].astype(str), pann["Gene"].astype(str)))

clin = pd.read_csv(ROOT / "05_tables" / "SuppTable1_complete_clinical_data.csv")
zon = pd.read_csv(RES / "zonation_collapse" / "zonation_clinical_merged.csv")

# patients with paired transcriptomics (Adjacent & Normal)
def paired_cols(mat, base):
    normals = {int(c.replace(base, "")) for c in mat.columns if c.startswith("Normal")}
    adj = {int(c.replace("Adjacent", "")) for c in mat.columns if c.startswith("Adjacent")}
    return sorted(normals & adj)

tx_patients = paired_cols(tx, "Normal")
prot_patients = paired_cols(prot, "Normal")
print("paired patients: tx", tx_patients, "\n              prot", prot_patients, flush=True)


def tx_log2fc(gene_name):
    gid = [g for g, n in idmap.items() if n == gene_name and g in tx.index]
    if not gid:
        return None
    row = tx.loc[gid[0]]
    return pd.Series(
        {p: row[f"Adjacent{p}"] - row[f"Normal{p}"] for p in tx_patients}, name=gene_name
    )


def prot_log2fc(gene_name):
    pids = [p for p, g in protmap.items() if g == gene_name and p in prot.index]
    if not pids:
        return None
    row = prot.loc[pids[0]]
    return pd.Series(
        {p: row[f"Adjacent{p}"] - row[f"Normal{p}"] for p in prot_patients}, name=gene_name
    )


# CYP family composite (transcriptome)
cyp_rows = [g for g, n in idmap.items() if n.startswith("CYP") and g in tx.index]
cyp_fc = (tx.loc[cyp_rows, [f"Adjacent{p}" for p in tx_patients]].values
          - tx.loc[cyp_rows, [f"Normal{p}" for p in tx_patients]].values)
cyp_fc = pd.Series(cyp_fc.mean(axis=0), index=tx_patients, name="CYP_family_mean")
print(f"CYP genes detected: {len(cyp_rows)}", flush=True)

feats = {}
for name, gene, layer in [
    ("CYP3A4_mRNA", "CYP3A4", "tx"),
    ("CYP3A4_protein", "CYP3A4", "prot"),
    ("HNF4A_mRNA", "HNF4A", "tx"),
    ("HNF1A_mRNA", "HNF1A", "tx"),
    ("FMO3_mRNA", "FMO3", "tx"),
]:
    s = tx_log2fc(gene) if layer == "tx" else prot_log2fc(gene)
    if s is not None:
        feats[name] = s
feats["CYP_family_mean"] = cyp_fc
# zonation delta (periportal collapse) from zonation module
zps = zon.set_index("patient_int")["delta_ZPS"]
feats["Zonation_collapse_deltaZPS"] = zps.reindex(tx_patients)

feat_df = pd.DataFrame(feats)
print("feature matrix:\n", feat_df.round(3).to_string(), flush=True)

# clinical variables
clin_vars = {
    "Lesion size (mm)": "Lesion size (mm)",
    "PNM-P stage": "PNM-P stage",
    "Bile duct invasion": "Bile duct invasion grade",
    "Vascular invasion": "Vascular invasion grade",
    "Disease duration (y)": "Disease duration (years)",
    "CRP (mg/L)": "CRP (mg/L)",
    "IL-6 (pg/mL)": "IL-6 (pg/mL)",
    "ALP (U/L)": "ALP (U/L)",
    "GGT (U/L)": "GGT (U/L)",
    "Total bilirubin": "Total bilirubin (µmol/L)",
}
clin_ix = clin.set_index("Patient ID")

# ------------------------------------------------- correlation with permutation
N_PERM = 10000


def perm_spearman(x, y, n_perm=N_PERM):
    mask = x.notna() & y.notna()
    if mask.sum() < 5:
        return np.nan, np.nan, int(mask.sum())
    rho, p_asym = stats.spearmanr(x[mask], y[mask])
    null = np.array(
        [stats.spearmanr(x[mask].sample(frac=1, random_state=RNG.integers(1e9)).values,
                         y[mask].values)[0]
         for _ in range(n_perm)]
    )
    p_perm = (np.abs(null) >= abs(rho)).mean()
    return rho, min(p_perm, 1.0), int(mask.sum())


rows = []
for fname, fser in feat_df.items():
    for cname, ccol in clin_vars.items():
        if ccol not in clin_ix.columns:
            continue
        y = clin_ix[ccol]
        rho, p, n = perm_spearman(fser, y.reindex(fser.index))
        rows.append({"feature": fname, "clinical": cname, "rho": rho, "p_perm": p, "n": n})
cor = pd.DataFrame(rows)

# BH FDR across the full matrix
m = cor["p_perm"].notna().sum()
order = cor["p_perm"].rank(method="min")
cor["p_BH"] = (cor["p_perm"] * m / order).clip(upper=1.0)
cor = cor.sort_values("p_perm")
cor.to_csv(OUT / "clinical_anchor_correlations.csv", index=False)
print("\ncorrelations (top 15):\n", cor.head(15).round(4).to_string(index=False), flush=True)

# ------------------------------------------------------ subtype clinical comparison
sub = clin.dropna(subset=["Molecular subtype (K=2)"])
sub_rows = []
for ccol in [clin_vars[k] for k in clin_vars if clin_vars[k] in clin.columns] + ["Age (years)"]:
    g1 = sub.loc[sub["Molecular subtype (K=2)"] == "CS1", ccol].dropna()
    g2 = sub.loc[sub["Molecular subtype (K=2)"] == "CS2", ccol].dropna()
    if len(g1) >= 3 and len(g2) >= 3:
        u, p = stats.mannwhitneyu(g1, g2, alternative="two-sided")
        sub_rows.append(
            {"variable": ccol, "CS1_median": g1.median(), "CS1_n": len(g1),
             "CS2_median": g2.median(), "CS2_n": len(g2), "p_MWU": p}
        )
subtab = pd.DataFrame(sub_rows).sort_values("p_MWU")
subtab.to_csv(OUT / "subtype_clinical_comparison.csv", index=False)
print("\nsubtype comparison:\n", subtab.round(4).to_string(index=False), flush=True)

# ------------------------------------------------------------------ heatmap
hm = cor.pivot(index="feature", columns="clinical", values="rho")
pm = cor.pivot(index="feature", columns="clinical", values="p_perm")
fig, ax = plt.subplots(figsize=(10, 4.2))
vmax = max(0.6, np.nanmax(np.abs(hm.values)) * 1.05)
im = ax.imshow(hm.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(hm.shape[1]), hm.columns, rotation=40, ha="right", fontsize=8)
ax.set_yticks(range(hm.shape[0]), hm.index, fontsize=8)
for i in range(hm.shape[0]):
    for j in range(hm.shape[1]):
        r = hm.values[i, j]
        if np.isnan(r):
            continue
        star = ""
        if pm.values[i, j] < 0.05:
            star = "*"
        if pm.values[i, j] < 0.01:
            star = "**"
        ax.text(j, i, f"{r:.2f}{star}", ha="center", va="center",
                fontsize=6.5, color="white" if abs(r) > vmax * 0.6 else "black")
plt.colorbar(im, ax=ax, label="Spearman rho", shrink=0.8)
ax.set_title("Peri-lesional molecular derangement vs disease burden\n"
             "(paired log2FC; permutation p, * <0.05, ** <0.01; n=12)", fontsize=9)
plt.tight_layout()
plt.savefig(OUT / "clinical_anchor_heatmap.pdf")
plt.savefig(OUT / "clinical_anchor_heatmap.png", dpi=200)
print("saved heatmap", flush=True)

# ------------------------------------------------------------------ summary
sig = cor[cor["p_perm"] < 0.05]
lines = [
    "# Enhancement 41: clinical anchoring (dose-response with disease burden)",
    "",
    f"- Features: {list(feat_df.columns)}",
    f"- Clinical variables: {list(clin_vars.keys())}",
    f"- Tests: {len(cor)} Spearman correlations, permutation p (10,000) + BH FDR; n=12 pairs",
    "",
    "## Significant associations (permutation p < 0.05)",
]
if len(sig):
    for _, r in sig.iterrows():
        lines.append(
            f"- {r['feature']} vs {r['clinical']}: rho={r['rho']:.2f}, "
            f"p_perm={r['p_perm']:.4f}, p_BH={r['p_BH']:.3f}"
        )
else:
    lines.append("- None at permutation p < 0.05")
lines += ["", "## Subtype (CS1/CS2) clinical comparison", ""]
for _, r in subtab.iterrows():
    if r["p_MWU"] < 0.1:
        lines.append(
            f"- {r['variable']}: CS1 median {r['CS1_median']:.1f} (n={r['CS1_n']}) vs "
            f"CS2 {r['CS2_median']:.1f} (n={r['CS2_n']}), p={r['p_MWU']:.3f}"
        )
(OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
print("saved summary.md", flush=True)
