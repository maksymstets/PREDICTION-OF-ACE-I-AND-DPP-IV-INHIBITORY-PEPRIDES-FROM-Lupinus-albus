"""
ADMETLab 3.0 peptide ADMET pipeline (single-molecule endpoint)
================================================================

CONFIRMED via diagnostic testing (see admetlab3_diagnose.py):
  - Endpoint: POST https://admetlab3.scbdd.com/api/single/admet
  - Payload:  {"SMILES": "<one SMILES string>"}   <- capital SMILES, ONE molecule per call
  - The server does NOT reliably error on malformed payloads - some malformed
    shapes silently return a cached/default example compound's data instead
    of an error. This script guards against that by verifying the SMILES
    echoed back in the response actually corresponds to what was sent.

Input CSV columns expected: Sequence, Protein name, Activity, Enzyme ID,
Peptide ID, smiles.

Output columns: all of the above metadata columns first (so you can always
trace a result back to its sequence, source protein, activity, enzyme, and
database ID), followed by every ADMET endpoint returned by the API.

Workflow:
1. Load peptides + metadata from CSV
2. Call the API once per peptide (this endpoint doesn't batch)
3. Validate each response against the input SMILES
4. Save incrementally so a crash/interruption doesn't lose completed work

This script only fetches raw ADMET data - it does not apply any filtering.
Run admetlab3_apply_filter.py separately on the raw output to flag
"optimal ADMET" peptides (f20 < 0.3 and ROA < 0.3).

Requirements:
    pip install requests pandas tqdm rdkit
"""

import os
import time
import requests
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_FILE = "/home/maksym/lupinus/ADMET_inputs/optimised_peptides_for_admet_metadata.csv"
OUTPUT_DIR = "/home/maksym/lupinus/ADMETLab_3.0_results"

RAW_OUTPUT_FILE = os.path.join(OUTPUT_DIR, "admetlab3_raw_results.csv")
FAILED_LOG_FILE = os.path.join(OUTPUT_DIR, "admetlab3_failed.csv")

# ---------------------------------------------------------------------------
# Input columns
# ---------------------------------------------------------------------------
SMILES_COL = "smiles"
# Metadata columns to carry through to the output, in this exact order.
METADATA_COLS = ["Sequence", "Protein name", "Activity", "Enzyme ID", "Peptide ID", "smiles"]

# ---------------------------------------------------------------------------
# API settings
# ---------------------------------------------------------------------------
BASE_URL = "https://admetlab3.scbdd.com"
ENDPOINT = "/api/single/admet"
SLEEP_BETWEEN_CALLS = 1.5   # seconds - be polite to the free server, one call at a time
TIMEOUT = 60
RETRY_ATTEMPTS = 3
RETRY_BACKOFF = 3           # seconds, doubles each retry

# Known "default fallback" compound the server returns when a payload is
# malformed but doesn't raise an error - used as a tripwire.
KNOWN_DEFAULT_SMILES = {"CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N",
                         "CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"}


def canonical(smi):
    """Canonicalize a SMILES for robust comparison. Returns None if invalid."""
    try:
        mol = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(mol) if mol else None
    except Exception:
        return None


def query_one(smiles):
    """
    Query ADMETLab3 for a single SMILES. Returns the parsed 'data' dict
    (category -> list of endpoint dicts), or raises on failure/validation error.
    """
    url = BASE_URL + ENDPOINT
    payload = {"SMILES": smiles}

    last_exc = None
    backoff = RETRY_BACKOFF
    for attempt in range(1, RETRY_ATTEMPTS + 1):
        try:
            resp = requests.post(url, json=payload, timeout=TIMEOUT)
            resp.raise_for_status()
            body = resp.json()

            returned_smiles = body.get("data", {}).get("smiles")
            returned_raw = body.get("data", {}).get("raw_smiles")

            # Tripwire: catch the known "silently returned default compound" failure mode
            if returned_smiles in KNOWN_DEFAULT_SMILES or returned_raw in KNOWN_DEFAULT_SMILES:
                raise ValueError(
                    f"Server returned the known default/fallback compound instead of "
                    f"a real prediction for input SMILES: {smiles}"
                )

            # Sanity check: returned molecule should match what we sent
            input_canon = canonical(smiles)
            returned_canon = canonical(returned_smiles) if returned_smiles else None
            if input_canon is not None and returned_canon is not None:
                if input_canon != returned_canon:
                    raise ValueError(
                        f"Response SMILES mismatch. Sent: {smiles} "
                        f"(canonical {input_canon}) but got back "
                        f"{returned_smiles} (canonical {returned_canon})"
                    )

            return body["data"]["data"]  # the actual category->endpoints dict

        except (requests.exceptions.RequestException, ValueError, KeyError) as e:
            last_exc = e
            if attempt < RETRY_ATTEMPTS:
                time.sleep(backoff)
                backoff *= 2

    raise RuntimeError(f"Failed after {RETRY_ATTEMPTS} attempts: {last_exc}")


def flatten_one(mol_data):
    """Flatten one molecule's {category: [ {name, value}, ... ]} into a flat dict."""
    row = {}
    for _, endpoints in mol_data.items():
        for ep in endpoints:
            row[ep["name"]] = ep["value"]
    return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_FILE)

    missing_meta = [c for c in METADATA_COLS if c not in df.columns]
    if missing_meta:
        raise ValueError(
            f"Expected metadata columns not found in {INPUT_FILE}: {missing_meta}. "
            f"Available columns: {df.columns.tolist()}"
        )

    valid_df = df.dropna(subset=[SMILES_COL]).copy()
    missing = len(df) - len(valid_df)
    if missing:
        print(f"WARNING: {missing} rows have no SMILES and will be skipped.")

    rows = []
    failed = []

    for _, r in tqdm(valid_df.iterrows(), total=len(valid_df), desc="Querying ADMETLab3"):
        smi = r[SMILES_COL]
        meta = {col: r[col] for col in METADATA_COLS}
        try:
            mol_data = query_one(smi)
            row = dict(meta)  # Sequence, Protein name, Activity, Enzyme ID, Peptide ID, smiles
            row.update(flatten_one(mol_data))
            rows.append(row)
        except Exception as e:
            print(f"  FAILED [{r['Peptide ID']}]: {e}")
            failed.append({**meta, "error": str(e)})

        time.sleep(SLEEP_BETWEEN_CALLS)

        # Save progress incrementally every 20 molecules in case of interruption
        if len(rows) % 20 == 0 and rows:
            pd.DataFrame(rows).to_csv(RAW_OUTPUT_FILE, index=False)

    results_df = pd.DataFrame(rows)
    results_df.to_csv(RAW_OUTPUT_FILE, index=False)
    print(f"\nSaved raw results to {RAW_OUTPUT_FILE}  ({len(results_df)} / {len(valid_df)} succeeded)")

    if failed:
        pd.DataFrame(failed).to_csv(FAILED_LOG_FILE, index=False)
        print(f"Saved {len(failed)} failures to {FAILED_LOG_FILE}")

    admet_cols = [c for c in results_df.columns if c not in METADATA_COLS]
    print("ADMET endpoint columns returned by the API:")
    print(admet_cols)
    