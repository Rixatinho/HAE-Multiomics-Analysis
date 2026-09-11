#!/usr/bin/env python3
"""
Enhancement 35: Bayesian simulation of rifampicin-mediated enzymatic rescue.

Purpose
-------
The manuscript (Results, rescue framework paragraph) currently states a
deterministic claim: "rifampicin-mediated induction (5-fold, upper bound of
the reported 3-5-fold range) restores all 12 modelled patients above
threshold (Fig. 2K)".  Fig. 2K was built in revision_v2_CYP_FMO_PK.R with
    activity_rescue = min(activity * 5, 1.0)          # capped at normal
    exposure = baseline * (act*f_cyp + fmo3*(1-f_cyp))
using SIMULATED patient activities ~ N(0.356, 0.12).

This script replaces the deterministic claim with a posterior probability
under the same structural model as Enhancement 31 (Bayesian MCMC), with:
  - the induction multiplier as an uncertain quantity drawn from the
    reported 3-5-fold range (ind ~ TruncNormal(4, 0.6, [3, 5]));
  - observed per-patient CYP3A4 mRNA log2FC (hierarchical model, n = 12);
  - induced activity capped at physiological normal (1.0), consistent with
    the published Fig. 2K logic (local activation cannot exceed the
    systemic ceiling in this two-compartment approximation).

Sensitivity scenarios
---------------------
  S1  primary: no change in systemic ABZ-SO (baseline prior as model A)
  S2  conservative: systemic ABZ-SO reduced 25% (baseline * 0.75),
      reflecting rifampicin-accelerated sulfoxide inactivation
  S3  supra-normal induction allowed (cap at 1.5 instead of 1.0)

Outputs
-------
  02_analysis/results/enhancement35_rifampicin_rescue/
    rifampicin_rescue_summary.csv
    rifampicin_rescue_patient.csv
    rifampicin_rescue_report.md
"""

import os
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az

BASE = "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
OUT = os.path.join(BASE, "02_analysis/results/enhancement35_rifampicin_rescue")
os.makedirs(OUT, exist_ok=True)

THRESHOLD = 250.0  # ng/mL

# Observed per-patient CYP3A4 mRNA log2FC (identical to enhancement31 model A)
MRNA_LOG2FC = np.array([
    -1.840323, -2.883263, -8.711582, -1.344386,  0.685132, -0.658449,
    -2.085768, -0.811944, -1.194944,  0.272610,  0.132242,  0.368645,
])
MRNA_PATIENTS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 14]


def fit_rescue(label, sys_factor=1.0, act_cap=1.0, draws=4000, tune=2000,
               chains=4, seed=20260823):
    """Bayesian rifampicin-rescue model (structure mirrors enhancement31)."""
    with pm.Model() as model:
        # Literature-grounded priors (identical to enhancement31 model A)
        baseline = pm.TruncatedNormal("baseline", mu=600, sigma=100,
                                      lower=350, upper=1000)
        f_cyp = pm.TruncatedNormal("f_cyp", mu=0.85, sigma=0.08,
                                   lower=0.50, upper=0.98)
        fmo3_resid = pm.TruncatedNormal("fmo3_resid", mu=0.87, sigma=0.08,
                                        lower=0.55, upper=1.10)

        # Rifampicin induction multiplier: reported 3-5-fold range
        # (Niemi et al. 2003 Clin Pharmacol Ther; manuscript Fig. 2K anchor)
        ind = pm.TruncatedNormal("induction_fold", mu=4.0, sigma=0.6,
                                 lower=3.0, upper=5.0)

        # Hierarchical per-patient activity (fitted to observed log2FC)
        mu_act = pm.Normal("mu_log2FC", mu=np.mean(MRNA_LOG2FC), sigma=1.5)
        sigma_act = pm.HalfNormal("sigma_log2FC", sigma=1.0)
        log2fc = pm.Normal("log2FC_obs", mu=mu_act, sigma=sigma_act,
                           observed=MRNA_LOG2FC)
        log2fc_new = pm.Normal("log2FC_new", mu=mu_act, sigma=sigma_act)

        act_mean = pm.Deterministic("activity_mean", 2.0 ** mu_act)
        act_new = pm.Deterministic("activity_new", 2.0 ** log2fc_new)

        # Induced (rescue) activity, capped at physiological normal
        act_res_typ = pm.math.minimum(act_mean * ind, act_cap)
        act_res_new = pm.math.minimum(act_new * ind, act_cap)
        act_obs = 2.0 ** MRNA_LOG2FC
        act_res_obs = pm.math.minimum(act_obs * ind, act_cap)

        base_sys = baseline * sys_factor

        # (a) typical patient
        pm.Deterministic(
            "abzso_rescue_typ",
            base_sys * (act_res_typ * f_cyp + fmo3_resid * (1 - f_cyp)))
        # (b) posterior predictive new patient
        abzso_res_new = pm.Deterministic(
            "abzso_rescue_new",
            base_sys * (act_res_new * f_cyp + fmo3_resid * (1 - f_cyp)))
        pm.Deterministic("p_above_new", pm.math.gt(abzso_res_new, THRESHOLD))
        # (c) observed patients
        abzso_res_obs = pm.Deterministic(
            "abzso_rescue_obs",
            base_sys * (act_res_obs * f_cyp + fmo3_resid * (1 - f_cyp)))

        idata = pm.sample(draws, tune=tune, chains=chains, seed=seed,
                          target_accept=0.95, progressbar=False,
                          return_inferencedata=True)

    post = idata.posterior
    typ = post["abzso_rescue_typ"].values.flatten()
    new = post["abzso_rescue_new"].values.flatten()
    obs = post["abzso_rescue_obs"].values.reshape(-1, len(MRNA_PATIENTS))

    ess = float(az.ess(idata, var_names=["baseline", "f_cyp", "mu_log2FC",
                                         "induction_fold"]).to_array().min())
    rhat = float(az.rhat(idata, var_names=["baseline", "f_cyp", "mu_log2FC",
                                           "induction_fold"]).to_array().max())

    pat_rows = []
    for j, p in enumerate(MRNA_PATIENTS):
        pat_rows.append({
            "patient": p,
            "exposure_rescue_median": float(np.median(obs[:, j])),
            "exposure_rescue_lo": float(np.percentile(obs[:, j], 2.5)),
            "exposure_rescue_hi": float(np.percentile(obs[:, j], 97.5)),
            "p_above_threshold": float(np.mean(obs[:, j] > THRESHOLD)),
        })

    return {
        "label": label,
        "typ_median": float(np.median(typ)),
        "typ_lo": float(np.percentile(typ, 2.5)),
        "typ_hi": float(np.percentile(typ, 97.5)),
        "p_typ_above": float(np.mean(typ > THRESHOLD)),
        "new_median": float(np.median(new)),
        "new_lo": float(np.percentile(new, 2.5)),
        "new_hi": float(np.percentile(new, 97.5)),
        "p_new_above": float(np.mean(new > THRESHOLD)),
        "ind_median": float(np.median(post["induction_fold"].values.flatten())),
        "n_patients_above_0.5": int(sum(r["p_above_threshold"] > 0.5
                                        for r in pat_rows)),
        "patients": pd.DataFrame(pat_rows),
        "ess_bulk_min": ess,
        "rhat_max": rhat,
    }


def main():
    print("=== Enhancement 35: rifampicin rescue Bayesian simulation ===\n")

    scenarios = [
        ("primary (cap 1.0, systemic unchanged)", dict(sys_factor=1.0, act_cap=1.0)),
        ("conservative (cap 1.0, systemic -25%)", dict(sys_factor=0.75, act_cap=1.0)),
        ("supra-normal induction (cap 1.5)",     dict(sys_factor=1.0, act_cap=1.5)),
    ]

    results = []
    for label, kw in scenarios:
        print(f"[fit] {label} ...")
        results.append(fit_rescue(label, **kw))

    # ---- summary csv ----
    rows = []
    for r in results:
        rows.append({
            "scenario": r["label"],
            "typical_patient_median_ng_mL": round(r["typ_median"], 1),
            "typical_95CrI": f'{r["typ_lo"]:.0f}-{r["typ_hi"]:.0f}',
            "P_typical_above_threshold": round(r["p_typ_above"], 3),
            "new_patient_median_ng_mL": round(r["new_median"], 1),
            "new_patient_95CrI": f'{r["new_lo"]:.0f}-{r["new_hi"]:.0f}',
            "P_new_patient_above_threshold": round(r["p_new_above"], 3),
            "n_observed_patients_P_above_0.5": r["n_patients_above_0.5"],
            "induction_fold_posterior_median": round(r["ind_median"], 2),
            "rhat_max": round(r["rhat_max"], 3),
            "ess_bulk_min": round(r["ess_bulk_min"]),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "rifampicin_rescue_summary.csv"), index=False)
    print(df.to_string(index=False))

    pat = results[0]["patients"].copy()
    pat.insert(0, "scenario", "primary")
    for r in results[1:]:
        p2 = r["patients"].copy()
        p2.insert(0, "scenario", r["label"])
        pat = pd.concat([pat, p2])
    pat.to_csv(os.path.join(OUT, "rifampicin_rescue_patient.csv"), index=False)

    # ---- report ----
    r0 = results[0]
    lines = [
        "# Enhancement 35 — Bayesian rifampicin rescue simulation",
        "",
        f"Threshold: {THRESHOLD:.0f} ng/mL. Structure identical to Enhancement 31 "
        "model A; induction multiplier ind ~ TruncNormal(4.0, 0.6, [3, 5]) "
        "(reported 3-5-fold range); induced activity capped at physiological "
        "normal, consistent with the published Fig. 2K logic.",
        "",
        "## Primary scenario (cap 1.0, systemic unchanged)",
        f"- Typical-patient rescued exposure: **{r0['typ_median']:.1f} ng/mL** "
        f"(95% CrI {r0['typ_lo']:.0f}-{r0['typ_hi']:.0f}), "
        f"P(above threshold) = **{r0['p_typ_above']:.3f}**",
        f"- New-patient predictive: {r0['new_median']:.1f} ng/mL "
        f"({r0['new_lo']:.0f}-{r0['new_hi']:.0f}), "
        f"P(above) = {r0['p_new_above']:.3f}",
        f"- Observed patients with P(above) > 0.5: "
        f"{r0['n_patients_above_0.5']}/12",
        f"- Posterior median induction fold: {r0['ind_median']:.2f}",
        "",
        "## Sensitivity scenarios",
    ]
    for r in results[1:]:
        lines.append(
            f"- {r['label']}: typical {r['typ_median']:.1f} "
            f"({r['typ_lo']:.0f}-{r['typ_hi']:.0f}), "
            f"P(above) = {r['p_typ_above']:.3f}; new patient "
            f"P(above) = {r['p_new_above']:.3f}; "
            f"{r['n_patients_above_0.5']}/12 observed patients P>0.5")
    lines += [
        "",
        "## Diagnostics",
        f"- Max Rhat: {r0['rhat_max']:.3f}; min bulk ESS: {r0['ess_bulk_min']:.0f}",
        "",
        "## Interpretation",
        "The deterministic claim 'restores all 12 modelled patients above",
        "threshold' is replaced by a posterior probability. Under the primary",
        "scenario the typical patient is rescued with high probability, while",
        "the most severely depleted patient (P3, log2FC = -8.71) remains",
        "sub-threshold even at the upper induction bound, quantifying",
        "patient-to-patient heterogeneity that the deterministic figure",
        "could not express.",
    ]
    with open(os.path.join(OUT, "rifampicin_rescue_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n[done] outputs in", OUT)


if __name__ == "__main__":
    main()
