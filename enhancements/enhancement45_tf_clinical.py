#!/usr/bin/env python
"""Enhancement 45: DoRothEA TF-activity x clinical-severity anchoring.

Connects the externally-validated TF activities (enhancement42, DoRothEA
curated regulons) to the clinical-severity gradient (enhancement41 framework):
does the activity of the collapsing hepatic TF circuitry track cholestasis /
disease chronicity in the peri-lesional samples?

Permutation Spearman (10,000 perms) + BH; Adjacent samples, n = 12.
Outputs -> 02_analysis/results/enhancement45_tf_clinical/
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(42)
ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "02_analysis" / "results" / "enhancement45_tf_clinical"
RES.mkdir(parents=True, exist_ok=True)
N_PERM = 10_000

tf = pd.read_csv(ROOT / "02_analysis/results/enhancement42_dorothea_progeny/dorothea_tf_activity_matrix.csv",
                 index_col=0)
tf.index = tf.index.astype(str)

clin = pd.read_csv(ROOT / "05_tables/SuppTable1_complete_clinical_data.csv")
clin.columns = [str(c).strip() for c in clin.columns]


def find_col(kws):
    for c in clin.columns:
        if all(k in c.lower() for k in kws):
            return c
    return None


VAR_MAP = {
    "bilirubin": find_col(["bilirubin"]),
    "ALP": find_col(["alp"]) or find_col(["alkaline"]),
    "GGT": find_col(["ggt"]),
    "CRP": find_col(["crp"]),
    "IL6": find_col(["il-6"]) or find_col(["il6"]),
    "lesion_size": find_col(["lesion", "size"]),
    "disease_duration": find_col(["duration"]) or find_col(["course"]),
}
VAR_MAP = {k: v for k, v in VAR_MAP.items() if v is not None}

pid_col = find_col(["patient"]) or clin.columns[0]
clin[pid_col] = clin[pid_col].astype(str).str.strip()
clin["pid"] = clin[pid_col].str.extract(r"(\d+)")[0]
clin_small = clin[["pid"] + list(VAR_MAP.values())].drop_duplicates("pid")

KEY_TFS = ["HNF4A", "HNF1A", "CEBPA", "RXRA", "PPARA", "NR5A2", "HNF4G",
           "FOXA2", "ONECUT1", "SMAD4"]
key_tfs = [t for t in KEY_TFS if t in tf.columns]
print("TFs available:", key_tfs)

adj = tf.loc[[i for i in tf.index if i.startswith("Adjacent")]].copy()
adj["pid"] = [i.replace("Adjacent", "") for i in adj.index]
merged = adj.merge(clin_small, on="pid", how="inner")
print("patients matched:", len(merged))


def perm_spearman(x, y):
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if len(x) < 6:
        return np.nan, np.nan
    rho, _ = stats.spearmanr(x, y)
    if np.isnan(rho):
        return np.nan, np.nan
    null = np.array([stats.spearmanr(RNG.permutation(x), y)[0] for _ in range(N_PERM)])
    p = (np.abs(null) >= abs(rho)).mean()
    return rho, p


rows = []
for t in key_tfs:
    for vname, vcol in VAR_MAP.items():
        rho, p = perm_spearman(merged[t].values.astype(float),
                               pd.to_numeric(merged[vcol], errors="coerce").values)
        rows.append({"TF": t, "clinical_var": vname, "rho": rho, "p_perm": p})

res = pd.DataFrame(rows).dropna().sort_values("p_perm").reset_index(drop=True)
# BH
p = res["p_perm"].values
order = np.argsort(p)
ranked = p[order]
m = len(p)
q = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
res.loc[order, "BH_q"] = np.clip(q, 0, 1)
res.to_csv(RES / "tf_clinical_permutation.csv", index=False)

sig = res[res["BH_q"] < 0.05]
print("\nSignificant (BH<0.05):")
print(sig.to_string(index=False) if len(sig) else "  none")

lines = ["# Enhancement 45: DoRothEA TF activity x clinical severity (peri-lesional, n=12)",
         "", "Date: 2026-08-22", "",
         f"Permutation Spearman ({N_PERM:,}) + BH; TF activities from enhancement42 (DoRothEA/decoupleR consensus).",
         f"Tests: {len(res)}; significant (BH<0.05): {len(sig)}", "",
         "## Key TFs x clinical variables (all tests)"]
for r in res.itertuples():
    lines.append(f"- {r.TF} x {r.clinical_var}: rho={r.rho:+.2f}, p={r.p_perm:.4f}, q={r.BH_q:.3f}")
lines += ["", "## Interpretation"]
hep = res[res["TF"].isin(["HNF4A", "HNF1A", "CEBPA", "RXRA", "PPARA"])]
bil = hep[hep["clinical_var"] == "bilirubin"]
if len(bil):
    neg = (bil["rho"] < 0).sum()
    lines.append(f"- Hepatic-circuitry TFs (HNF4A/HNF1A/CEBPA/RXRA/PPARA) x bilirubin: {neg}/{len(bil)} negative "
                 + "; ".join(f"{r.TF} {r.rho:+.2f} (p={r.p_perm:.3f})" for r in bil.itertuples())
                 + (" - activity falls as cholestasis deepens, extending the molecular dose-response (enhancement41) "
                    "to the externally-validated TF-activity layer." if neg == len(bil) else ""))
(RES / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("saved:", RES)
