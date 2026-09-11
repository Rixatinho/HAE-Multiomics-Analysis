#!/usr/bin/env python
"""Enhancement 37: Real LINCS L1000 / CMap connectivity query via Enrichr REST API.

Disease signature: paired peri-lesional vs distal normal liver (HAE, n=12 pairs),
t-statistic ranked top/bottom 150 genes (signature_up150/down150.csv).

Database: Enrichr library `LINCS_L1000_Chem_Pert_Consensus_Sigs` — consensus
differential-expression signatures (up/down gene sets) for chemical perturbations
from the LINCS L1000 CMap release.

Reversal logic (CMap paradigm): perturbation X reverses the disease if the genes
it DOWN-regulates overlap disease-UP genes and the genes it UP-regulates overlap
disease-DOWN genes:
    reversal(X) = nlp(UP_query -> "X Down") + nlp(DOWN_query -> "X Up")
                - nlp(UP_query -> "X Up")  - nlp(DOWN_query -> "X Down")

Enrichr result row format (9 fields):
    [rank, term, p_value, z_score, combined_score, overlapping_genes,
     adjusted_p, old_p, old_adjusted]

Outputs -> 02_analysis/results/enhancement37_lincs_cmap/
"""
import json
import math
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "02_analysis" / "results" / "enhancement37_lincs_cmap"
RES.mkdir(parents=True, exist_ok=True)
LIB = "LINCS_L1000_Chem_Pert_Consensus_Sigs"
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


def fetch_enrichment(user_list_id, library=LIB, retries=4):
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
    raise RuntimeError(f"Enrichr enrichment failed for list {user_list_id}")


def to_df(rows):
    df = pd.DataFrame(rows, columns=["rank", "term", "p_value", "z_score",
                                     "combined_score", "genes", "adj_p",
                                     "old_p", "old_adj"])
    df["n_overlap"] = df["genes"].apply(len)
    df["nlp"] = df["p_value"].apply(lambda x: 0.0 if (x is None or x <= 0) else min(50, -math.log10(x)))
    return df.drop(columns=["genes"])


def nlp_map(df):
    return df.set_index("term")["nlp"].to_dict()


def main():
    up = pd.read_csv(RES / "signature_up150.csv")["gene_name"].dropna().astype(str).tolist()
    dn = pd.read_csv(RES / "signature_down150.csv")["gene_name"].dropna().astype(str).tolist()
    print(f"signature: {len(up)} up / {len(dn)} down genes", flush=True)

    uid_up = add_list(up, "HAE_peri_lesional_UP150")
    uid_dn = add_list(dn, "HAE_peri_lesional_DOWN150")
    print("userListIds:", uid_up, uid_dn, flush=True)

    print("querying", LIB, "UP list ...", flush=True)
    up_df = to_df(fetch_enrichment(uid_up))
    up_df.to_csv(RES / "lincs_up_query_results.csv", index=False)
    print("  rows:", len(up_df), flush=True)
    print("querying", LIB, "DOWN list ...", flush=True)
    dn_df = to_df(fetch_enrichment(uid_dn))
    dn_df.to_csv(RES / "lincs_down_query_results.csv", index=False)
    print("  rows:", len(dn_df), flush=True)

    up_m = nlp_map(up_df)
    dn_m = nlp_map(dn_df)

    # parse terms "COMPOUND Up" / "COMPOUND Down"
    compounds = {}
    for term in set(up_m) | set(dn_m):
        m = re.match(r"^(.*)\s+(Up|Down)$", term)
        if m:
            compounds.setdefault(m.group(1), {})[m.group(2)] = term

    rows = []
    for comp, dirs in compounds.items():
        up_in_disease_up = up_m.get(dirs.get("Up", ""), 0.0)     # similarity (bad)
        dn_in_disease_up = up_m.get(dirs.get("Down", ""), 0.0)   # reversal (good)
        up_in_disease_dn = dn_m.get(dirs.get("Up", ""), 0.0)     # reversal (good)
        dn_in_disease_dn = dn_m.get(dirs.get("Down", ""), 0.0)   # similarity (bad)
        rev = dn_in_disease_up + up_in_disease_dn - up_in_disease_up - dn_in_disease_dn
        rows.append({
            "compound": comp, "reversal_score": rev,
            "rev_via_downsig_in_diseaseUP": dn_in_disease_up,
            "rev_via_upsig_in_diseaseDOWN": up_in_disease_dn,
            "sim_up_up": up_in_disease_up, "sim_dn_dn": dn_in_disease_dn,
        })
    rev_df = pd.DataFrame(rows).sort_values("reversal_score", ascending=False).reset_index(drop=True)
    rev_df.to_csv(RES / "lincs_cmap2020_reversal_scores.csv", index=False)

    print("\nTop 30 reversing compounds:")
    print(rev_df.head(30).round(2).to_string(index=False))

    cands = ["pirfenidone", "nintedanib", "rifampicin", "dexamethasone", "sirolimus",
             "imatinib", "metformin", "valproic", "trichostatin", "vorinostat",
             "all-trans", "tretinoin", "fasudil", "losartan", "spironolactone"]
    low = rev_df["compound"].str.lower()
    print("\nNominated candidates in LINCS:")
    for c in cands:
        hit = rev_df[low.str.contains(c)]
        if len(hit):
            best = hit.iloc[0]
            rank = rev_df.index[rev_df["compound"] == best["compound"]][0] + 1
            print(f"  {c}: rank {rank}/{len(rev_df)} reversal={best['reversal_score']:.2f}")

    (RES / "summary.md").write_text(
        "# Enhancement 37: LINCS L1000 CMap real-database query (Enrichr)\n\n"
        f"- Signature: {len(up)} up / {len(dn)} down genes (paired peri-lesional vs distal, t-ranked)\n"
        f"- Library: {LIB} (UP query rows: {len(up_df)}; DOWN query rows: {len(dn_df)})\n"
        f"- Compounds scored: {len(rev_df)}\n"
        f"- Top reversing compound: {rev_df.iloc[0]['compound']} "
        f"(score {rev_df.iloc[0]['reversal_score']:.2f})\n\n"
        "Reversal = nlp(perturbation-Down | disease-UP) + nlp(perturbation-Up | disease-DOWN)\n"
        "         - nlp(perturbation-Up | disease-UP) - nlp(perturbation-Down | disease-DOWN)\n",
        encoding="utf-8",
    )
    print("\nsaved:", RES / "lincs_cmap2020_reversal_scores.csv", flush=True)


if __name__ == "__main__":
    main()
