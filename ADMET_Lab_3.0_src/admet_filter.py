"""
Apply the optimal-ADMET filter to an already-fetched ADMETLab 3.0 raw results
CSV. Doesn't touch the API - just re-processes admetlab3_raw_results.csv.

Expects the raw results file to already contain the metadata columns
(Sequence, Protein name, Activity, Enzyme ID, Peptide ID, smiles) followed by
ADMET endpoint columns, as produced by admetlab3_peptide_pipeline.py.
"""

import os
import pandas as pd

OUTPUT_DIR = "/home/maksym/lupinus/ADMETLab_3.0_results"
RAW_FILE = os.path.join(OUTPUT_DIR, "admetlab3_raw_results.csv")
FINAL_FILE = os.path.join(OUTPUT_DIR, "admetlab3_final_results.csv")

METADATA_COLS = ["Sequence", "Protein name", "Activity", "Enzyme ID", "Peptide ID", "smiles"]

# Confirmed column names from the live API response:
#   F(20%) bioavailability   -> f20
#   rat oral acute toxicity  -> ROA
F20_COL = "f20"
ROAT_COL = "ROA"

df = pd.read_csv(RAW_FILE)

missing_meta = [c for c in METADATA_COLS if c not in df.columns]
if missing_meta:
    raise ValueError(
        f"Missing expected metadata columns {missing_meta} in {RAW_FILE}. "
        f"Available: {df.columns.tolist()}\n"
        f"(If this file was produced by an older version of the pipeline, "
        f"it may only contain 'peptide_sequence'/'peptide_id' + 'smiles' - "
        f"rerun the pipeline to regenerate it with full metadata.)"
    )

missing_admet = [c for c in (F20_COL, ROAT_COL) if c not in df.columns]
if missing_admet:
    raise ValueError(f"Missing expected ADMET columns {missing_admet}. Available: {df.columns.tolist()}")

df["optimal_ADMET"] = (df[F20_COL] < 0.3) & (df[ROAT_COL] < 0.3)

# Keep metadata columns first, in order, followed by all ADMET result columns
ordered_cols = METADATA_COLS + [c for c in df.columns if c not in METADATA_COLS]
df = df[ordered_cols]

n_total = len(df)
df = df[df["optimal_ADMET"]].reset_index(drop=True)

df.to_csv(FINAL_FILE, index=False)
print(f"Saved final results to {FINAL_FILE}")
print(f"{len(df)} / {n_total} peptides met the optimal ADMET criteria "
      f"(f20 < 0.3 and ROA < 0.3) and were kept in the final file.")