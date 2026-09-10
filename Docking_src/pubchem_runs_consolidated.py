"""
Consolidate multiple pubchem_structure_download.py run outputs into one.

Since transient network errors caused different peptides to fail in
different runs, this script:
  1. Auto-discovers every subfolder under BASE_DIR containing a
     pubchem_lookup_summary.csv (i.e. each prior run's output)
  2. For each peptide, picks the BEST result found across all runs
     (downloaded > cid_found_no_3d_structure > no_pubchem_entry > error)
  3. Copies the winning SDF file (if any) into one clean merged folder
  4. Writes a consolidated summary CSV
  5. Flags any peptide that failed with "error" in every single run -
     these are the only ones that genuinely still need a retry/investigation

Requirements:
    pip install pandas
"""

import shutil
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path("/home/maksym/lupinus/Docking_inputs")
METADATA_FILE = Path("/home/maksym/lupinus/ADMET_inputs/optimised_peptides_for_admet_metadata.csv")

MERGED_SDF_DIR = BASE_DIR / "pubchem_3d_sdf_merged"
MERGED_SDF_DIR.mkdir(parents=True, exist_ok=True)

CONSOLIDATED_SUMMARY = BASE_DIR / "pubchem_lookup_summary_consolidated.csv"
STILL_FAILING_FILE = BASE_DIR / "peptides_still_failing.csv"

# Priority: lower number = better result, kept when a peptide appears in multiple runs
STATUS_PRIORITY = {
    "downloaded": 0,
    "cid_found_no_3d_structure": 1,   # definitive PubChem answer, just no 3D record
    "no_pubchem_entry": 1,            # definitive PubChem answer, no entry at all
}
DEFAULT_PRIORITY = 2  # anything else (errors) - inconclusive, worth retrying


def status_rank(status):
    if isinstance(status, str) and status.startswith("error"):
        return DEFAULT_PRIORITY
    return STATUS_PRIORITY.get(status, DEFAULT_PRIORITY)


def normalize_id(x):
    try:
        f = float(x)
        return str(int(f)) if f.is_integer() else str(x).strip()
    except (TypeError, ValueError):
        return str(x).strip()


def resolve_sdf_path(raw_path, run_dir):
    """Try the recorded path as-is, then fall back to same filename inside run_dir."""
    if pd.isna(raw_path) or not str(raw_path).strip():
        return None
    candidate = Path(str(raw_path))
    if candidate.exists():
        return candidate
    fallback = run_dir / candidate.name
    return fallback if fallback.exists() else None


# ---------------------------------------------------------------------------
# 1. Discover all prior runs
# ---------------------------------------------------------------------------
run_summaries = sorted(BASE_DIR.rglob("pubchem_lookup_summary.csv"))
if not run_summaries:
    raise SystemExit(f"No pubchem_lookup_summary.csv files found under {BASE_DIR}")

print(f"Found {len(run_summaries)} prior run(s):")
for p in run_summaries:
    print(f"  {p}")

all_rows = []
for summary_path in run_summaries:
    run_dir = summary_path.parent
    df = pd.read_csv(summary_path)
    df["Peptide ID"] = df["Peptide ID"].apply(normalize_id)
    df["_run_dir"] = str(run_dir)
    df["_resolved_sdf"] = df.apply(
        lambda r: resolve_sdf_path(r.get("sdf_file"), run_dir), axis=1
    )
    df["_status_rank"] = df["status"].apply(status_rank)
    all_rows.append(df)

combined = pd.concat(all_rows, ignore_index=True)

# ---------------------------------------------------------------------------
# 2. Pick the best row per peptide across all runs
# ---------------------------------------------------------------------------
# Among rows with the best (lowest) status_rank for a peptide, prefer one
# whose SDF file actually resolves on disk (covers edge case of a
# "downloaded" status whose file went missing in one particular run folder).
combined["_has_file"] = combined["_resolved_sdf"].notna()
combined_sorted = combined.sort_values(
    by=["_status_rank", "_has_file"], ascending=[True, False]
)
best = combined_sorted.drop_duplicates(subset=["Peptide ID"], keep="first").copy()

# ---------------------------------------------------------------------------
# 3. Copy winning SDF files into one merged folder
# ---------------------------------------------------------------------------
merged_paths = []
for _, row in best.iterrows():
    if row["status"] == "downloaded" and row["_resolved_sdf"] is not None:
        dest = MERGED_SDF_DIR / row["_resolved_sdf"].name
        if not dest.exists():
            shutil.copy(row["_resolved_sdf"], dest)
        merged_paths.append(str(dest))
    else:
        merged_paths.append(None)

best["sdf_file"] = merged_paths
best = best.drop(columns=["_run_dir", "_resolved_sdf", "_status_rank", "_has_file"], errors="ignore")

# ---------------------------------------------------------------------------
# 4. Cross-check against the full peptide list in the metadata file
# ---------------------------------------------------------------------------
meta_df = pd.read_csv(METADATA_FILE)
meta_df["Peptide ID"] = meta_df["Peptide ID"].apply(normalize_id)

missing_entirely = meta_df[~meta_df["Peptide ID"].isin(best["Peptide ID"])]

best.to_csv(CONSOLIDATED_SUMMARY, index=False)

still_failing = best[best["status"].astype(str).str.startswith("error")]
if len(missing_entirely):
    # Peptides that appeared in NO run at all (shouldn't normally happen, but check anyway)
    missing_entirely.to_csv(STILL_FAILING_FILE, index=False)
    still_failing = pd.concat([still_failing, missing_entirely], ignore_index=True)
if len(still_failing):
    still_failing.to_csv(STILL_FAILING_FILE, index=False)

# ---------------------------------------------------------------------------
# 5. Report
# ---------------------------------------------------------------------------
n_total = len(meta_df)
n_downloaded = (best["status"] == "downloaded").sum()
n_no_3d = (best["status"] == "cid_found_no_3d_structure").sum()
n_no_entry = (best["status"] == "no_pubchem_entry").sum()
n_still_error = len(still_failing)

print(f"\n=== Consolidated across {len(run_summaries)} runs ===")
print(f"Total peptides in metadata file : {n_total}")
print(f"Best-of-all-runs downloaded     : {n_downloaded}")
print(f"CID found, no 3D (all runs)     : {n_no_3d}")
print(f"No PubChem entry (all runs)     : {n_no_entry}")
print(f"Still failing in every run      : {n_still_error}")
print(f"\nConsolidated summary saved to: {CONSOLIDATED_SUMMARY}")
print(f"Merged SDF files copied to:    {MERGED_SDF_DIR}")
if n_still_error:
    print(f"Peptides needing a retry saved to: {STILL_FAILING_FILE}")
else:
    print("No peptides are stuck on errors across all runs - full coverage achieved.")