#!/usr/bin/env python3
"""
Enhancement 31: Bayesian MCMC uncertainty propagation for the
two-compartment peri-lesional ABZ-SO pharmacokinetic model.

Purpose
-------
Replaces the point estimates (181.6 / 259.9 ng/mL) with full posterior
distributions. All model inputs are assigned informative priors grounded
in their primary literature; per-patient CYP3A4 activity enters through
a hierarchical model fitted to the observed paired log2FC values.

Responds to: CommsBio R3.1 residual concern (literature-parameter-based
model) and R1.2 (no direct functional validation) by quantifying the
probability of sub-therapeutic peri-lesional exposure.

Model
-----
  ABZ-SO_peri,i = baseline * ( activity_i * f_cyp + fmo3_resid * (1 - f_cyp) )
  ABZ-SO_cyp,i  = baseline * ( activity_i * f_cyp )

Priors
------
  baseline        ~ TruncNormal(600, 100, [350, 1000])   systemic ABZ-SO Cmax
                     (Duthaler 2019 BJC; Mingjie 2014; sensitivity range 400-800)
  f_cyp           ~ TruncNormal(0.85, 0.08, [0.50, 0.98]) CYP3A4 share of
                     ABZ sulfoxidation (Rawden 2000 EPS; sensitivity 50-85%)
  fmo3_resid      ~ TruncNormal(0.87, 0.08, [0.55, 1.10]) FMO3 residual protein
                     (observed log2FC = -0.225)
  log2FC_i        ~ Normal(mu_act, sigma_act)             hierarchical patient
                     activity: activity_i = 2 ** log2FC_i
                     (fitted to 12 paired mRNA values; protein-based sensitivity)

Outputs
-------
  02_analysis/results/enhancement31_pk_bayesian/
    pk_mcmc_posterior_summary.csv   parameter + derived quantity posteriors
    pk_mcmc_patient_risk.csv        per-patient P(exposure < threshold)
    pk_mcmc_model_comparison.csv    mRNA vs protein activity calibre
    FigS_pk_mcmc_posterior.pdf      posterior density figure
    pk_mcmc_report.md               human-readable summary
"""

import os
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
OUT = os.path.join(BASE, "02_analysis/results/enhancement31_pk_bayesian")
os.makedirs(OUT, exist_ok=True)

THRESHOLD = 250.0  # ng/mL, therapeutic target (Kern 2017)

# ----------------------------------------------------------------------------
# Observed data (extracted from data/processed matrices; see header)
# ----------------------------------------------------------------------------
# Per-patient paired log2FC, CYP3A4 mRNA (logCPM difference, n=12)
MRNA_LOG2FC = np.array([
    -1.840323, -2.883263, -8.711582, -1.344386,  0.685132, -0.658449,
    -2.085768, -0.811944, -1.194944,  0.272610,  0.132242,  0.368645,
])
MRNA_PATIENTS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 14]

# Per-patient paired log2FC, CYP3A4 protein (log2 intensity difference, n=14)
PROT_LOG2FC = np.array([
     0.220154, -0.396241, -1.227346, -0.727717, -0.913271, -0.480279,
    -0.534814,  0.032358, -0.181044, -4.141127, -1.062312,  0.688140,
    -0.622174, -0.047374,
])
PROT_PATIENTS = list(range(1, 15))

# FMO3 protein residual activity observed: 2^-0.225 = 0.854 (geometric).
# Manuscript uses 0.87 (arithmetic of 2^x would differ); we centre the prior
# at 0.87 with SD 0.08 and run a prior-shift sensitivity at 0.854.
FMO3_PRIOR_MU = 0.87


def fit_model(log2fc_obs, label, draws=4000, tune=2000, chains=4,
              fmo3_mu=FMO3_PRIOR_MU, seed=20260822):
    """Hierarchical MCMC for peri-lesional ABZ-SO exposure."""
    with pm.Model() as model:
        # Literature-grounded priors
        baseline = pm.TruncatedNormal("baseline", mu=600, sigma=100,
                                      lower=350, upper=1000)
        f_cyp = pm.TruncatedNormal("f_cyp", mu=0.85, sigma=0.08,
                                   lower=0.50, upper=0.98)
        fmo3_resid = pm.TruncatedNormal("fmo3_resid", mu=fmo3_mu, sigma=0.08,
                                        lower=0.55, upper=1.10)

        # Hierarchical per-patient activity
        mu_act = pm.Normal("mu_log2FC", mu=np.mean(log2fc_obs),
                           sigma=1.5)
        sigma_act = pm.HalfNormal("sigma_log2FC", sigma=1.0)
        log2fc = pm.Normal("log2FC_obs", mu=mu_act, sigma=sigma_act,
                           observed=log2fc_obs)

        # Population-level representative activity (posterior predictive)
        log2fc_new = pm.Normal("log2FC_new", mu=mu_act, sigma=sigma_act)
        activity_new = pm.Deterministic("activity_new", 2.0 ** log2fc_new)
        activity_mean = pm.Deterministic("activity_mean",
                                         2.0 ** mu_act)

        # Derived exposure quantities
        # (a) typical patient: exposure at the population-typical activity 2^mu
        pm.Deterministic("abzso_cyp_only_typ",
                         baseline * activity_mean * f_cyp)
        pm.Deterministic("abzso_cyp_fmo3_typ",
                         baseline * (activity_mean * f_cyp
                                     + fmo3_resid * (1 - f_cyp)))
        # (b) posterior predictive for a new patient (includes sigma heterogeneity)
        abzso_cyp = pm.Deterministic(
            "abzso_cyp_only", baseline * activity_new * f_cyp)
        abzso_fmo = pm.Deterministic(
            "abzso_cyp_fmo3",
            baseline * (activity_new * f_cyp + fmo3_resid * (1 - f_cyp)))

        # Observed-patient exposures (using each patient's measured log2FC)
        act_obs = 2.0 ** log2fc_obs
        abzso_cyp_pat = pm.Deterministic(
            "abzso_cyp_only_pat",
            baseline * (act_obs * f_cyp)[:, None] * np.ones((1,)))
        abzso_fmo_pat = pm.Deterministic(
            "abzso_cyp_fmo3_pat",
            baseline * (act_obs * f_cyp + fmo3_resid * (1 - f_cyp)))

        # Probabilities of sub-therapeutic exposure
        p_below_cyp = pm.Deterministic("p_below_cyp_only",
                                       pm.math.lt(abzso_cyp, THRESHOLD))
        p_below_fmo = pm.Deterministic("p_below_cyp_fmo3",
                                       pm.math.lt(abzso_fmo, THRESHOLD))

        idata = pm.sample(draws, tune=tune, chains=chains, seed=seed,
                          target_accept=0.95, progressbar=False,
                          return_inferencedata=True)
    return idata


def summarise(idata, label, patients):
    post = idata.posterior
    res = {}

    def ci(name):
        x = post[name].values.flatten()
        return np.median(x), np.percentile(x, 2.5), np.percentile(x, 97.5)

    for p in ["baseline", "f_cyp", "fmo3_resid", "activity_mean",
              "abzso_cyp_only_typ", "abzso_cyp_fmo3_typ",
              "abzso_cyp_only", "abzso_cyp_fmo3"]:
        m, lo, hi = ci(p)
        res[p] = (m, lo, hi)

    p_below_cyp = float(np.mean(post["p_below_cyp_only"].values))
    p_below_fmo = float(np.mean(post["p_below_cyp_fmo3"].values))

    # Per-patient posterior P(below threshold), CYP3A4+FMO3 calibre
    pat_fmo = post["abzso_cyp_fmo3_pat"].values.reshape(-1, len(patients))
    pat_cyp = post["abzso_cyp_only_pat"].values.reshape(-1, len(patients))
    pat_rows = []
    for j, p in enumerate(patients):
        pat_rows.append({
            "patient": p,
            "exposure_cyp_fmo3_median": float(np.median(pat_fmo[:, j])),
            "exposure_cyp_fmo3_lo": float(np.percentile(pat_fmo[:, j], 2.5)),
            "exposure_cyp_fmo3_hi": float(np.percentile(pat_fmo[:, j], 97.5)),
            "p_below_threshold": float(np.mean(pat_fmo[:, j] < THRESHOLD)),
            "p_below_threshold_cyp_only": float(np.mean(pat_cyp[:, j] < THRESHOLD)),
        })

    # Diagnostics
    ess = float(az.ess(idata, var_names=["baseline", "f_cyp",
                                  "mu_log2FC"]).to_array().min())
    rhat = float(az.rhat(idata, var_names=["baseline", "f_cyp",
                       "mu_log2FC"]).to_array().max())

    return {
        "label": label,
        "params": res,
        "p_below_cyp_only": p_below_cyp,
        "p_below_cyp_fmo3": p_below_fmo,
        "patients": pd.DataFrame(pat_rows),
        "ess_bulk_min": ess,
        "rhat_max": rhat,
        "idata": idata,
    }


def main():
    print("=== Enhancement 31: PK model Bayesian MCMC ===\n")

    # ---- Model A: mRNA-derived activity (primary, manuscript calibre) ----
    print("[1/3] Fitting model A (mRNA activity calibre, n=12)...")
    idata_a = fit_model(MRNA_LOG2FC, "mRNA")
    summ_a = summarise(idata_a, "mRNA_calibre", MRNA_PATIENTS)

    # ---- Model B: protein-derived activity (sensitivity) ----
    print("[2/3] Fitting model B (protein activity calibre, n=14)...")
    idata_b = fit_model(PROT_LOG2FC, "protein")
    summ_b = summarise(idata_b, "protein_calibre", PROT_PATIENTS)

    # ---- Model C: FMO3 prior shifted to observed geometric residual ----
    print("[3/3] Fitting model C (FMO3 prior 0.854 sensitivity)...")
    idata_c = fit_model(MRNA_LOG2FC, "fmo3_alt", fmo3_mu=0.854)
    summ_c = summarise(idata_c, "fmo3_0.854_sensitivity", MRNA_PATIENTS)

    # ---------------- Reporting ----------------
    rows = []
    for s in (summ_a, summ_b, summ_c):
        for p, (m, lo, hi) in s["params"].items():
            rows.append({"model": s["label"], "quantity": p,
                         "median": m, "lo95": lo, "hi95": hi})
        rows.append({"model": s["label"], "quantity": "P(new patient below threshold | CYP3A4 only)",
                     "median": s["p_below_cyp_only"], "lo95": np.nan, "hi95": np.nan})
        rows.append({"model": s["label"], "quantity": "P(new patient below threshold | CYP3A4+FMO3)",
                     "median": s["p_below_cyp_fmo3"], "lo95": np.nan, "hi95": np.nan})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "pk_mcmc_posterior_summary.csv"),
                              index=False)

    pat_all = []
    for s, cal in [(summ_a, "mRNA"), (summ_b, "protein")]:
        df = s["patients"].copy()
        df.insert(0, "activity_calibre", cal)
        pat_all.append(df)
    pd.concat(pat_all).to_csv(os.path.join(OUT, "pk_mcmc_patient_risk.csv"),
                              index=False)

    comp = pd.DataFrame([{
        "model": s["label"],
        "typical_patient_exposure_cyp_fmo3": s["params"]["abzso_cyp_fmo3_typ"][0],
        "typical_patient_95CrI": f'{s["params"]["abzso_cyp_fmo3_typ"][1]:.0f}-{s["params"]["abzso_cyp_fmo3_typ"][2]:.0f}',
        "typical_patient_exposure_cyp_only": s["params"]["abzso_cyp_only_typ"][0],
        "new_patient_predictive_median": s["params"]["abzso_cyp_fmo3"][0],
        "new_patient_predictive_95CrI": f'{s["params"]["abzso_cyp_fmo3"][1]:.0f}-{s["params"]["abzso_cyp_fmo3"][2]:.0f}',
        "P_below_threshold_cyp_fmo3": s["p_below_cyp_fmo3"],
        "P_below_threshold_cyp_only": s["p_below_cyp_only"],
        "rhat_max": s["rhat_max"], "ess_bulk_min": s["ess_bulk_min"],
    } for s in (summ_a, summ_b, summ_c)])
    comp.to_csv(os.path.join(OUT, "pk_mcmc_model_comparison.csv"), index=False)
    print(comp.to_string(index=False))

    # ---------------- Figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    post_a = summ_a["idata"].posterior

    ax = axes[0]
    for name, colour, lab in [
            ("abzso_cyp_fmo3_typ", "#3C5488", "typical patient (CYP3A4+FMO3)"),
            ("abzso_cyp_only_typ", "#E64B35", "typical patient (CYP3A4 only)"),
            ("abzso_cyp_fmo3", "#9FB2C8", "new-patient predictive")]:
        x = post_a[name].values.flatten()
        ax.hist(x, bins=80, density=True, alpha=0.55, color=colour,
                label=lab)
    ax.axvline(THRESHOLD, color="red", ls="--", lw=1)
    ax.text(THRESHOLD + 8, ax.get_ylim()[1] * 0.9, "threshold\n250 ng/mL",
            color="red", fontsize=8)
    ax.set_xlabel("Peri-lesional ABZ-SO (ng/mL)")
    ax.set_ylabel("Posterior density")
    ax.set_title("A  Posterior exposure (mRNA calibre)")
    ax.legend(fontsize=8)

    ax = axes[1]
    x = post_a["activity_new"].values.flatten()
    ax.hist(x, bins=80, density=True, color="#00A087", alpha=0.7)
    ax.axvline(1.0, color="grey", ls=":", lw=1)
    ax.set_xlabel("Peri-lesional CYP3A4 activity (fraction of normal)")
    ax.set_title("B  Posterior CYP3A4 activity")

    ax = axes[2]
    dfp = summ_a["patients"]
    ys = np.arange(len(dfp))
    meds = dfp["exposure_cyp_fmo3_median"].values
    los = dfp["exposure_cyp_fmo3_lo"].values
    his = dfp["exposure_cyp_fmo3_hi"].values
    risks = dfp["p_below_threshold"].values
    colours = plt.cm.RdYlGn(1 - risks)
    ax.hlines(ys, los, his, color=colours, lw=3)
    ax.scatter(meds, ys, color=colours, s=25, zorder=3)
    ax.axvline(THRESHOLD, color="red", ls="--", lw=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"P{int(p)} ({r:.2f})" for p, r in
                        zip(dfp['patient'], risks)], fontsize=7)
    ax.set_xlabel("ABZ-SO ng/mL (CYP3A4+FMO3)")
    ax.set_title("C  Per-patient posterior (P(below thr))")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "FigS_pk_mcmc_posterior.pdf"))
    plt.close(fig)

    # ---------------- Markdown report ----------------
    a = summ_a["params"]
    lines = [
        "# Enhancement 31 — PK model Bayesian MCMC uncertainty propagation",
        "",
        f"Threshold: {THRESHOLD:.0f} ng/mL. 4 chains x 4000 draws (target_accept 0.95).",
        "All point estimates below are posterior **medians**.",
        "",
        "## Primary model (mRNA activity calibre, n=12)",
        "### Typical patient (population-typical activity 2^mu)",
        f"- Peri-lesional ABZ-SO (CYP3A4 only): **{a['abzso_cyp_only_typ'][0]:.0f} ng/mL** "
        f"(95% CrI {a['abzso_cyp_only_typ'][1]:.0f}-{a['abzso_cyp_only_typ'][2]:.0f})",
        f"- Peri-lesional ABZ-SO (CYP3A4+FMO3): **{a['abzso_cyp_fmo3_typ'][0]:.0f} ng/mL** "
        f"(95% CrI {a['abzso_cyp_fmo3_typ'][1]:.0f}-{a['abzso_cyp_fmo3_typ'][2]:.0f})",
        "### Posterior predictive for a new patient (includes inter-patient heterogeneity sigma)",
        f"- ABZ-SO (CYP3A4 only): median {a['abzso_cyp_only'][0]:.0f} ng/mL "
        f"(95% CrI {a['abzso_cyp_only'][1]:.0f}-{a['abzso_cyp_only'][2]:.0f}), "
        f"P(< threshold) = {summ_a['p_below_cyp_only']:.3f}",
        f"- ABZ-SO (CYP3A4+FMO3): median {a['abzso_cyp_fmo3'][0]:.0f} ng/mL "
        f"(95% CrI {a['abzso_cyp_fmo3'][1]:.0f}-{a['abzso_cyp_fmo3'][2]:.0f}), "
        f"P(< threshold) = {summ_a['p_below_cyp_fmo3']:.3f}",
        f"- Posterior median CYP3A4 activity (typical patient): {a['activity_mean'][0]:.3f} "
        f"({a['activity_mean'][1]:.3f}-{a['activity_mean'][2]:.3f})",
        f"- baseline: {a['baseline'][0]:.0f} ({a['baseline'][1]:.0f}-{a['baseline'][2]:.0f}) "
        f"ng/mL; f_cyp: {a['f_cyp'][0]:.2f}; fmo3_resid: {a['fmo3_resid'][0]:.2f}",
        "",
        "## Sensitivity models",
        f"- Protein calibre (n=14): typical-patient exposure (CYP3A4+FMO3) "
        f"{summ_b['params']['abzso_cyp_fmo3_typ'][0]:.0f} ng/mL "
        f"({summ_b['params']['abzso_cyp_fmo3_typ'][1]:.0f}-{summ_b['params']['abzso_cyp_fmo3_typ'][2]:.0f}), "
        f"new-patient P(< threshold) = {summ_b['p_below_cyp_fmo3']:.3f}",
        f"- FMO3 prior 0.854: typical-patient exposure (CYP3A4+FMO3) "
        f"{summ_c['params']['abzso_cyp_fmo3_typ'][0]:.0f} ng/mL, "
        f"new-patient P(< threshold) = {summ_c['p_below_cyp_fmo3']:.3f}",
        "",
        "## Diagnostics",
        f"- Max Rhat (model A): {summ_a['rhat_max']:.3f}; "
        f"min bulk ESS: {summ_a['ess_bulk_min']:.0f}",
        "",
        "## Interpretation",
        "The point-estimate narrative (181.6 ng/mL CYP3A4-only; 259.9 ng/mL",
        "with FMO3) is replaced by full posterior distributions. The",
        "typical-patient posterior replaces the plug-in point estimate and",
        "quantifies joint uncertainty in systemic Cmax, CYP3A4 contribution",
        "fraction, FMO3 residual activity, and the population-level CYP3A4",
        "depletion. The new-patient posterior predictive additionally",
        "propagates the substantial inter-patient heterogeneity (posterior",
        "sigma of log2FC ~ 2), yielding the probability that a newly",
        "treated AE patient fails to reach the therapeutic threshold.",
    ]
    with open(os.path.join(OUT, "pk_mcmc_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print("\n[done] outputs in", OUT)


if __name__ == "__main__":
    main()
