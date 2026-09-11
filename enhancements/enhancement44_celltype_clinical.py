#!/usr/bin/env python
"""Enhancement 44: cell-fraction x clinical-severity anchoring across three
deconvolution methods (BayesPrism / MCP-counter / EPIC).

Extends enhancement41 (molecular-feature anchoring) to the cellular level:
do the immune/stromal shifts - especially the triangulated fibrotic/stromal
expansion - scale with clinical severity (lesion size, bilirubin, ALP/GGT,
disease duration, CRP/IL-6, PNM stage)?

Permutation Spearman (10,000 perms) + BH across tests.
Outputs -> 02_analysis/results/enhancement44_celltype_clinical/
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(42)
ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "02_analysis" / "results" / "enhancement44_celltype_clinical"
RES.mkdir(parents=True, exist_ok=True)
N_PERM = 10_000

# clinical variables (same set as enhancement41)
CLIN = ROOT / "05_tables" / "SuppTable1_complete_clinical_data.csv"
clin = pd.read_csv(CLIN)
clin.columns = [str(c).strip() for c in clin.columns]


def find_col(keywords):
    for c in clin.columns:
        lc = c.lower()
        if all(k in lc for k in keywords):
            return c
    return None


VAR_MAP = {
    "lesion_size_mm": find_col(["lesion", "size"]),
    "PNM_P": find_col(["pnm", "p"]) or find_col(["pnm"]),
    "bilirubin": find_col(["bilirubin"]),
    "ALP": find_col(["alp"]) or find_col(["alkaline"]),
    "GGT": find_col(["ggt"]),
    "CRP": find_col(["crp"]),
    "IL6": find_col(["il-6"]) or find_col(["il6"]),
    "disease_duration": find_col(["duration"]) or find_col(["course"]),
}
VAR_MAP = {k: v for k, v in VAR_MAP.items() if v is not None}
print("clinical variables:", VAR_MAP)

# patient id join key
pid_col = find_col(["patient"]) or find_col(["id"]) or clin.columns[0]
print("patient id column:", pid_col)
clin[pid_col] = clin[pid_col].astype(str).str.strip()
# patient ids like "1".."12" matching sample suffixes
clin["pid"] = clin[pid_col].str.extract(r"(\d+)")[0]

# ---- cell-fraction matrices (samples as rows) ------------------------------
E16 = ROOT / "02_analysis" / "results" / "enhancement16_deconvolution"
E39 = ROOT / "02_analysis" / "results" / "enhancement39_cross_deconv"

mats = {}
bp = pd.read_csv(E16 / "bayesprism_cell_proportions.csv", index_col=0)
mats["BayesPrism"] = bp
mcp = pd.read_csv(E39 / "mcpcounter_scores.csv", index_col=0)
mats["MCPcounter"] = mcp
ep = pd.read_csv(E39 / "epic_fractions.csv", index_col=0)
mats["EPIC"] = ep

for k, m in mats.items():
    m.index = m.index.astype(str)
    mats[k] = m

# subset to Adjacent (peri-lesional) samples for clinical anchoring
def adj_subset(m):
    idx = [i for i in m.index if i.startswith("Adjacent")]
    out = m.loc[idx].copy()
    out["pid"] = [i.replace("Adjacent", "") for i in out.index]
    return out


def perm_spearman(x, y):
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if len(x) < 6:
        return np.nan, np.nan
    rho, _ = stats.spearmanr(x, y)
    null = np.array([
        stats.spearmanr(RNG.permutation(x), y)[0] for _ in range(N_PERM)
    ])
    p = (np.abs(null) >= abs(rho)).mean() if not np.isnan(rho) else np.nan
    return rho, p


rows = []
for method, m in mats.items():
    A = adj_subset(m)
    merged = A.merge(clin[["pid"] + list(VAR_MAP.values())].drop_duplicates("pid"),
                     on="pid", how="inner")
    print(f"{method}: {len(merged)} patients matched")
    for cell in m.columns:
        for vname, vcol in VAR_MAP.items():
            if vcol not in merged.columns:
                continue
            rho, p = perm_spearman(merged[cell].values.astype(float),
                                   pd.to_numeric(merged[vcol], errors="coerce").values)
            rows.append({"method": method, "cell": cell, "clinical_var": vname,
                         "rho": rho, "p_perm": p, "n": int((~merged[cell].isna()).sum())})

res = pd.DataFrame(rows)
res["p_perm"] = res["p_perm"].fillna(1.0)
res["BH_q"] = stats.false_discovery_control if False else None
# manual BH
p = res["p_perm"].values
order = np.argsort(p)
ranked = p[order]
m = len(p)
q = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
q = np.clip(q, 0, 1)
res_sorted = res.loc[order].copy()
res_sorted["BH_q"] = q
res = res_sorted.sort_values("p_perm").reset_index(drop=True)
res.to_csv(RES / "celltype_clinical_permutation.csv", index=False)

sig = res[res["BH_q"] < 0.05]
print("\nSignificant (BH<0.05) cell x clinical associations:")
print(sig.to_string(index=False) if len(sig) else "  none")

# fibrosis focus: stromal/fibrotic features across methods
fibro_keys = {
    "BayesPrism": ["Stellate_cell"],
    "MCPcounter": ["Fibroblasts"],
    "EPIC": ["CAFs"],
}
lines = ["# Enhancement 44: cell-fraction x clinical-severity anchoring (3 deconvolution methods)",
         "", "Date: 2026-08-22", "",
         f"Permutation Spearman ({N_PERM:,} perms) + BH across all method x cell x clinical tests; peri-lesional (Adjacent) samples only.", "",
         f"Total tests: {len(res)}; significant (BH<0.05): {len(sig)}", "",
         "## Fibrotic/stromal features vs clinical severity"]
for method, keys in fibro_keys.items():
    for cell in keys:
        sub = res[(res["method"] == method) & (res["cell"] == cell)]
        for r in sub.itertuples():
            if not np.isnan(r.rho):
                lines.append(f"- {method} {cell} x {r.clinical_var}: rho={r.rho:+.2f}, p={r.p_perm:.4f}, q={r.BH_q:.3f}")
lines += ["", "## Top 15 associations overall"]
for r in res.head(15).itertuples():
    lines.append(f"- {r.method} {r.cell} x {r.clinical_var}: rho={r.rho:+.2f}, p={r.p_perm:.4f}, q={r.BH_q:.3f}")
(RES / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\nsaved:", RES)
