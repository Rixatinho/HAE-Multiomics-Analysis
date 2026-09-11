#!/usr/bin/env python
"""Enhancement 43: Cross-library CMap validation.

Two INDEPENDENT LINCS libraries queried with the same disease signature
(enhancement37: paired peri-lesional vs distal, t-ranked top/bottom 150 genes):

1. LINCS_L1000_CRISPR_KO_Consensus_Sigs (10,424 gene-knockout consensus
   signatures): genetic-perturbation mirror of the pharmacologic analysis.
   Top reversal hits = genes whose loss reverses the disease signature
   (candidate disease-maintaining genes); direction for HNF4A/HNF1A KOs
   tests the host-metabolic-collapse narrative.
2. LINCS_L1000_Chem_Pert_up / _down (33,132 raw non-consensus chemical
   perturbation signatures): robustness of the chemical reversal ranking
   (incl. nintedanib) outside the consensus-processing pipeline.

Reversal logic identical to enhancement37.

Outputs -> 02_analysis/results/enhancement43_lincs_crosslib/
"""
import json
import math
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "02_analysis" / "results" / "enhancement43_lincs_crosslib"
RES.mkdir(parents=True, exist_ok=True)
SRC = ROOT / "02_analysis" / "results" / "enhancement37_lincs_cmap"
BASE = "https://maayanlab.cloud/Enrichr"
UA = {"User-Agent": "Mozilla/5.0 (research; HAE multi-omics)"}


def add_list(genes, description):
    boundary = "----codebuddy" + str(int(time.time()))
    gene_str = "\n".join(genes)
    body = []
    body.append(f'--{boundary}\r\nContent-Disposition: form-data; name="description"\r\n\r\n{description}\r\n'.encode())
    body.append(f'--{boundary}\r\nContent-Disposition: form-data; name="list"\r\n\r\n{gene_str}\r\n'.encode())
    body.append(f"--{boundary}--\r\n".encode())
    req = urllib.request.Request(
        f"{BASE}/addList", data=b"".join(body),
        headers={**UA, "Content-Type": f"multipart/form-data; boundary={boundary}"},
    )
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)["userListId"]


def fetch_enrichment(user_list_id, library, retries=4):
    url = f"{BASE}/enrich?userListId={user_list_id}&backgroundType={library}"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as r:
                d = json.load(r)
            if library in d and d[library]:
                return d[library]
            print(f"    attempt {attempt+1}: empty, retrying...", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"    attempt {attempt+1} error: {e}", flush=True)
        time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"Enrichr failed for list {user_list_id} lib {library}")


def to_df(rows):
    df = pd.DataFrame(rows, columns=["rank", "term", "p_value", "z_score",
                                     "combined_score", "genes", "adj_p",
                                     "old_p", "old_adj"])
    df["nlp"] = df["p_value"].apply(lambda x: 0.0 if (x is None or x <= 0) else min(50, -math.log10(x)))
    return df.drop(columns=["genes"])


def nlp_map(df):
    return df.set_index("term")["nlp"].to_dict()


def main():
    up = pd.read_csv(SRC / "signature_up150.csv")["gene_name"].dropna().astype(str).tolist()
    dn = pd.read_csv(SRC / "signature_down150.csv")["gene_name"].dropna().astype(str).tolist()
    print(f"signature: {len(up)} up / {len(dn)} down genes", flush=True)

    uid_up = add_list(up, "HAE_UP150_crosslib")
    uid_dn = add_list(dn, "HAE_DOWN150_crosslib")
    print("userListIds:", uid_up, uid_dn, flush=True)

    summary_lines = [
        "# Enhancement 43: cross-library LINCS CMap validation",
        "", "Date: 2026-08-22", "",
        f"Signature: {len(up)} up / {len(dn)} down genes (identical to enhancement37).",
        "",
    ]

    # ---------- (1) CRISPR KO consensus library -----------------------------
    KO = "LINCS_L1000_CRISPR_KO_Consensus_Sigs"
    print("\n===", KO, "===", flush=True)
    up_df = to_df(fetch_enrichment(uid_up, KO)); up_df.to_csv(RES / "ko_up_query_results.csv", index=False)
    print("  UP rows:", len(up_df), flush=True)
    dn_df = to_df(fetch_enrichment(uid_dn, KO)); dn_df.to_csv(RES / "ko_down_query_results.csv", index=False)
    print("  DOWN rows:", len(dn_df), flush=True)
    up_m, dn_m = nlp_map(up_df), nlp_map(dn_df)

    genes = {}
    for term in set(up_m) | set(dn_m):
        m = re.match(r"^(.*)\s+(Up|Down)$", term)
        if m:
            genes.setdefault(m.group(1), {})[m.group(2)] = term
    rows = []
    for g, dirs in genes.items():
        rev = (up_m.get(dirs.get("Down", ""), 0.0) + dn_m.get(dirs.get("Up", ""), 0.0)
               - up_m.get(dirs.get("Up", ""), 0.0) - dn_m.get(dirs.get("Down", ""), 0.0))
        rows.append({"ko_gene": g, "reversal_score": rev})
    ko_df = pd.DataFrame(rows).sort_values("reversal_score", ascending=False).reset_index(drop=True)
    ko_df.to_csv(RES / "ko_reversal_scores.csv", index=False)
    print("  KOs scored:", len(ko_df))
    print("\nTop 20 disease-reversing knockouts:")
    print(ko_df.head(20).round(2).to_string(index=False))

    for q in ["HNF4A", "HNF1A", "CEBPA", "PPARA"]:
        hit = ko_df[ko_df["ko_gene"].str.upper() == q]
        if len(hit):
            rank = hit.index[0] + 1
            summary_lines.append(
                f"- KO {q}: reversal rank {rank}/{len(ko_df)} (score {hit.iloc[0]['reversal_score']:.2f})")
        else:
            summary_lines.append(f"- KO {q}: not scored in this library")
    summary_lines.append(f"- Top 5 disease-reversing KOs: "
                         + ", ".join(f"{r.ko_gene} ({r.reversal_score:.1f})"
                                     for r in ko_df.head(5).itertuples()))
    summary_lines.append("")
    summary_lines.append(
        "Note on interpretation: a high reversal rank for KO of gene G means the transcriptional state "
        "after losing G is opposite to the disease state - i.e. G's activity is required to maintain the "
        "disease-like signature. Conversely, KO signatures that MIRROR the disease (negative reversal) "
        "identify genes whose loss phenocopies the peri-lesional state.")

    # ---------- (2) Raw chemical perturbation library (non-consensus) -------
    print("\n=== LINCS_L1000_Chem_Pert raw up/down libraries ===", flush=True)
    lib_dn = "LINCS_L1000_Chem_Pert_down"   # compound-DOWN gene sets
    lib_up = "LINCS_L1000_Chem_Pert_up"     # compound-UP gene sets
    # disease-UP vs compound-DOWN (reversal evidence)
    q1 = to_df(fetch_enrichment(uid_up, lib_dn)); q1.to_csv(RES / "chempert_diseaseUP_vs_compoundDOWN.csv", index=False)
    print("  diseaseUP x compoundDOWN rows:", len(q1), flush=True)
    # disease-DOWN vs compound-UP (reversal evidence)
    q2 = to_df(fetch_enrichment(uid_dn, lib_up)); q2.to_csv(RES / "chempert_diseaseDOWN_vs_compoundUP.csv", index=False)
    print("  diseaseDOWN x compoundUP rows:", len(q2), flush=True)
    # disease-UP vs compound-UP (similarity penalty)
    q3 = to_df(fetch_enrichment(uid_up, lib_up)); q3.to_csv(RES / "chempert_diseaseUP_vs_compoundUP.csv", index=False)
    print("  diseaseUP x compoundUP rows:", len(q3), flush=True)
    # disease-DOWN vs compound-DOWN (similarity penalty)
    q4 = to_df(fetch_enrichment(uid_dn, lib_dn)); q4.to_csv(RES / "chempert_diseaseDOWN_vs_compoundDOWN.csv", index=False)
    print("  diseaseDOWN x compoundDOWN rows:", len(q4), flush=True)

    m1, m2, m3, m4 = (nlp_map(x) for x in (q1, q2, q3, q4))
    # terms may repeat (same compound, multiple cells/doses) -> take max nlp per compound
    def max_map(df):
        return df.groupby("term")["nlp"].max().to_dict()
    m1, m2, m3, m4 = (max_map(x) for x in (q1, q2, q3, q4))

    all_terms = set(m1) | set(m2) | set(m3) | set(m4)
    rows = []
    for t in all_terms:
        rev = m1.get(t, 0.0) + m2.get(t, 0.0) - m3.get(t, 0.0) - m4.get(t, 0.0)
        rows.append({"compound": t, "reversal_score": rev,
                     "rev_via_compoundDOWN_hits_diseaseUP": m1.get(t, 0.0),
                     "rev_via_compoundUP_hits_diseaseDOWN": m2.get(t, 0.0)})
    chem_df = pd.DataFrame(rows).sort_values("reversal_score", ascending=False).reset_index(drop=True)
    chem_df.to_csv(RES / "chempert_raw_reversal_scores.csv", index=False)
    print("  compounds scored (raw lib):", len(chem_df))

    key = ["nintedanib", "pirfenidone", "losartan", "rifampicin"]
    summary_lines.append("")
    summary_lines.append(f"Raw chemical-perturbation library ({len(chem_df)} signature terms):")
    for k in key:
        hit = chem_df[chem_df["compound"].str.lower().str.contains(k)]
        if len(hit):
            best = hit.iloc[0]
            rank = chem_df.index[chem_df["compound"] == best["compound"]][0] + 1
            summary_lines.append(f"- {best['compound'][:60]}: rank {rank}/{len(chem_df)} (score {best['reversal_score']:.2f})")

    (RES / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\nsaved:", RES)


if __name__ == "__main__":
    main()
