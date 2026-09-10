"""
Download PubChem 3D structures for peptides, replicating the methodology:
"Molecular docking analyses were carried out when the three-dimensional
structures ... were available in the PubChem Database. Molecular docking
analysis was performed on all ... peptides with reported PubChem or SatPDB
structures."

This script covers the PubChem half of that pipeline:
  1. For each peptide's SMILES, look up its PubChem CID
  2. If a CID exists, try to download its 3D SDF structure
  3. Peptides with NO CID, or a CID but NO 3D conformer, are logged and
     SKIPPED (not backfilled with a generated conformer) - matching the
     paper's approach of only docking peptides with a reported structure.

SatPDB has no public API, so it can't be automated here. Peptides that fail
the PubChem lookup are flagged in the summary so you can check them on
SatPDB manually if you want to extend coverage.

Requirements:
    pip install requests pandas tqdm
"""

import os
import time
import requests
import pandas as pd
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_FILE = "/home/maksym/lupinus/ADMET_inputs/optimised_peptides_for_admet_metadata.csv"
OUTPUT_DIR = "/home/maksym/lupinus/Docking_inputs/pubchem_3d_sdf"
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "pubchem_lookup_summary.csv")

METADATA_COLS = ["Sequence", "Protein name", "Activity", "Enzyme ID", "Peptide ID", "smiles"]
SMILES_COL = "smiles"
ID_COL = "Peptide ID"

# ---------------------------------------------------------------------------
# PubChem PUG REST settings
# ---------------------------------------------------------------------------
PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
SLEEP_BETWEEN_CALLS = 0.25  # PubChem asks for <=5 requests/sec; stay well under that
TIMEOUT = 30
RETRY_ATTEMPTS = 3
RETRY_BACKOFF = 2


def _request_with_retries(method, url, **kwargs):
    """GET/POST with retries on transient errors (429/503/connection issues)."""
    backoff = RETRY_BACKOFF
    last_exc = None
    for attempt in range(1, RETRY_ATTEMPTS + 1):
        try:
            resp = requests.request(method, url, timeout=TIMEOUT, **kwargs)
            if resp.status_code in (429, 503):
                raise requests.exceptions.HTTPError(f"{resp.status_code} rate-limited/unavailable")
            return resp
        except requests.exceptions.RequestException as e:
            last_exc = e
            if attempt < RETRY_ATTEMPTS:
                time.sleep(backoff)
                backoff *= 2
    raise last_exc


def get_cid_from_smiles(smiles):
    """
    Look up a PubChem CID for a given SMILES.
    Returns the CID (int) or None if PubChem has no matching compound.
    Uses POST so unusual SMILES characters don't break URL encoding.
    """
    url = f"{PUG_BASE}/compound/smiles/cids/JSON"
    resp = _request_with_retries("POST", url, data={"smiles": smiles})

    if resp.status_code == 404:
        return None  # no matching compound in PubChem
    resp.raise_for_status()

    data = resp.json()
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def download_3d_sdf(cid):
    """
    Try to download the 3D SDF record for a CID.
    Returns the SDF text, or None if PubChem has no 3D conformer for this CID.
    """
    url = f"{PUG_BASE}/compound/cid/{cid}/SDF"
    resp = _request_with_retries("GET", url, params={"record_type": "3d"})

    if resp.status_code == 404:
        return None  # compound exists but has no 3D conformer
    resp.raise_for_status()
    return resp.text


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_FILE)
    missing_meta = [c for c in METADATA_COLS if c not in df.columns]
    if missing_meta:
        raise ValueError(f"Missing columns {missing_meta}. Available: {df.columns.tolist()}")

    summary_rows = []

    for _, r in tqdm(df.iterrows(), total=len(df), desc="Querying PubChem"):
        pid = r[ID_COL]
        seq = r["Sequence"]
        smi = r[SMILES_COL]

        status = None
        cid = None
        sdf_path = None

        try:
            cid = get_cid_from_smiles(smi)
            time.sleep(SLEEP_BETWEEN_CALLS)

            if cid is None:
                status = "no_pubchem_entry"
            else:
                sdf_text = download_3d_sdf(cid)
                time.sleep(SLEEP_BETWEEN_CALLS)

                if sdf_text is None or not sdf_text.strip():
                    status = "cid_found_no_3d_structure"
                else:
                    fname = f"{pid}_{seq}_CID{cid}_3d.sdf".replace(" ", "_")
                    sdf_path = os.path.join(OUTPUT_DIR, fname)
                    with open(sdf_path, "w") as f:
                        f.write(sdf_text)
                    status = "downloaded"

        except Exception as e:
            status = f"error: {e}"

        summary_rows.append({
            "Peptide ID": pid,
            "Sequence": seq,
            "smiles": smi,
            "pubchem_cid": cid,
            "status": status,
            "sdf_file": sdf_path,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SUMMARY_FILE, index=False)

    n_downloaded = (summary_df["status"] == "downloaded").sum()
    n_no_entry = (summary_df["status"] == "no_pubchem_entry").sum()
    n_no_3d = (summary_df["status"] == "cid_found_no_3d_structure").sum()
    n_errors = summary_df["status"].str.startswith("error").sum()

    print(f"\nDone. Summary saved to {SUMMARY_FILE}")
    print(f"  Downloaded 3D structures : {n_downloaded} / {len(summary_df)}")
    print(f"  No PubChem entry at all  : {n_no_entry}")
    print(f"  CID exists, no 3D record : {n_no_3d}")
    print(f"  Errors                   : {n_errors}")
    print(f"\nPeptides with no PubChem 3D structure ({n_no_entry + n_no_3d}) were "
          f"skipped, matching the paper's approach - check SatPDB manually for "
          f"those if you want to extend coverage (no public API available for it).")
