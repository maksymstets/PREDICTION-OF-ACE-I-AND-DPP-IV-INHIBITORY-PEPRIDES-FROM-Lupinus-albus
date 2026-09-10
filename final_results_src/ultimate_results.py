"""
Build the ultimate combined results CSV from:
  - ADMET-AI results (file2, primary peptide list)
  - ADMETLab 3.0 results (file1, supplementary columns + extra peptides)
  - Docking results at exhaustiveness=8 and exhaustiveness=32

Requirements:
    pip install pandas
"""

import re
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ADMETLAB3_FILE = Path("/home/maksym/lupinus/ADMETLab_3.0_results/admetlab3_final_results_gid_annotated.csv")
ADMETAI_FILE = Path("/home/maksym/lupinus/ADMET_results/optimal_admet_activity_all_lengths_gid_annotated.csv")

DOCK_EX32_DIR = Path("/home/maksym/lupinus/dock_res_2run")
DOCK_EX8_DIR = Path("/home/maksym/lupinus/docking_results")

OUTPUT_DIR = Path("/home/maksym/lupinus/final_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "ultimate_peptide_results.csv"

# Columns pulled from ADMETLab3 into the merged table (source_name -> output_name)
ADMETLAB3_COLS = {
    "Enzyme ID": "Enzyme ID",
    "hia": "hia",
    "f20": "f20",
    "Lipinski": "Lipinski_ADMETLab3",          # renamed - collides with ADMET-AI's own Lipinski
    "VDss": "VDss",
    "t0.5": "t0.5",
    "ROA": "ROA",
    "GID resistance pH 1.3": "GID resistance pH 1.3_ADMETLab3",   # renamed - collides
    "GID resistance ph>2": "GID resistance ph>2_ADMETLab3",       # renamed - collides
}

FINAL_COLUMN_ORDER = [
    "Sequence", "Length", "Activity_Class", "Activity",
    "Lipinski", "HIA_Hou", "VDss_Lombardo", "Half_Life_Obach", "Bioavailability_Ma", "ClinTox",
    "GID resistance pH 1.3", "GID resistance ph>2",
    "Enzyme ID", "hia", "f20", "Lipinski_ADMETLab3", "VDss", "t0.5", "ROA",
    "GID resistance pH 1.3_ADMETLab3", "GID resistance ph>2_ADMETLab3",
    "exhaustiveness_8_Rank", "exhaustiveness_8_Peptide_ID",
    "exhaustiveness_8_ACE Affinity (kcal/mol)", "exhaustiveness_8_DPP-IV Affinity (kcal/mol)",
    "exhaustiveness_8_Average Affinity",
    "exhaustiveness_32_Rank", "exhaustiveness_32_Peptide_ID",
    "exhaustiveness_32_ACE Affinity (kcal/mol)", "exhaustiveness_32_DPP-IV Affinity (kcal/mol)",
    "exhaustiveness_32_Average Affinity",
]


def classify_activity(act_string):
    """Same classification logic used in the original ADMET-AI script."""
    if not isinstance(act_string, str):
        return "Other"
    act_lower = act_string.lower()
    has_ace = "ace" in act_lower
    has_dpp = "dipeptidyl peptidase" in act_lower or "dpp" in act_lower
    if has_ace and has_dpp:
        return "Dual"
    elif has_ace:
        return "ACE"
    elif has_dpp:
        return "DPP"
    else:
        return "Other"


# ---------------------------------------------------------------------------
# STEP 1: Merge the two ADMET files
# ---------------------------------------------------------------------------
print("=== Step 1: merging ADMET-AI and ADMETLab3 results ===")

admetai_df = pd.read_csv(ADMETAI_FILE)
admetlab3_df = pd.read_csv(ADMETLAB3_FILE)

admetai_df = admetai_df.drop(columns=["smiles"], errors="ignore")

missing_cols = [c for c in ADMETLAB3_COLS if c not in admetlab3_df.columns]
if missing_cols:
    raise ValueError(f"Expected columns not found in ADMETLab3 file: {missing_cols}")

lab3_subset = admetlab3_df[["Sequence"] + list(ADMETLAB3_COLS.keys())].rename(columns=ADMETLAB3_COLS)

# Left join: keep every ADMET-AI peptide, attach ADMETLab3 columns where a Sequence match exists
merged = admetai_df.merge(lab3_subset, on="Sequence", how="left")

# Peptides present in ADMETLab3 but NOT in ADMET-AI - append as extra rows
unmatched_seqs = set(admetlab3_df["Sequence"]) - set(admetai_df["Sequence"])
print(f"ADMET-AI peptides: {len(admetai_df)} | ADMETLab3-only extra peptides: {len(unmatched_seqs)}")

extra_rows = []
for _, row in admetlab3_df[admetlab3_df["Sequence"].isin(unmatched_seqs)].iterrows():
    extra_rows.append({
        "Sequence": row["Sequence"],
        "Length": len(str(row["Sequence"])),
        "Activity_Class": classify_activity(row.get("Activity")),
        "Activity": row.get("Activity"),
        # ADMET-AI-native columns unavailable for these peptides -> left blank
        "Lipinski": None, "HIA_Hou": None, "VDss_Lombardo": None,
        "Half_Life_Obach": None, "Bioavailability_Ma": None, "ClinTox": None,
        "GID resistance pH 1.3": None, "GID resistance ph>2": None,
        # ADMETLab3 columns filled in from the source row
        "Enzyme ID": row.get("Enzyme ID"),
        "hia": row.get("hia"),
        "f20": row.get("f20"),
        "Lipinski_ADMETLab3": row.get("Lipinski"),
        "VDss": row.get("VDss"),
        "t0.5": row.get("t0.5"),
        "ROA": row.get("ROA"),
        "GID resistance pH 1.3_ADMETLab3": row.get("GID resistance pH 1.3"),
        "GID resistance ph>2_ADMETLab3": row.get("GID resistance ph>2"),
    })

if extra_rows:
    merged = pd.concat([merged, pd.DataFrame(extra_rows)], ignore_index=True)

print(f"Combined ADMET table: {len(merged)} peptides total")

# ---------------------------------------------------------------------------
# STEP 2: Load and combine docking results per exhaustiveness level
# ---------------------------------------------------------------------------
def extract_sequence(peptide_id):
    """'10692_VRF' -> 'VRF' (strip the numeric ID prefix before the first underscore)."""
    if not isinstance(peptide_id, str):
        return None
    parts = peptide_id.split("_", 1)
    return parts[1] if len(parts) == 2 else peptide_id


def load_docking_set(dock_dir, ace_file, dpp_file, dual_file, label):
    """
    Load and combine the ACE-only, DPP-only, and Dual ranked-results CSVs for
    one exhaustiveness level into a single dataframe keyed by extracted Sequence.
    """
    ace_df = pd.read_csv(dock_dir / ace_file)
    dpp_df = pd.read_csv(dock_dir / dpp_file)
    dual_df = pd.read_csv(dock_dir / dual_file)

    combined = pd.concat([ace_df, dpp_df, dual_df], ignore_index=True, sort=False)
    combined["Sequence"] = combined["Peptide_ID"].apply(extract_sequence)

    # Handle any Sequence collisions (same sequence from two different numeric IDs):
    # keep the row with the most favorable (lowest/most negative) affinity available.
    def best_affinity(row):
        vals = [v for v in [row.get("Average Affinity"),
                             row.get("ACE Affinity (kcal/mol)"),
                             row.get("DPP-IV Affinity (kcal/mol)")] if pd.notna(v)]
        return min(vals) if vals else float("inf")

    combined["_sort_key"] = combined.apply(best_affinity, axis=1)
    dupes = combined[combined.duplicated(subset=["Sequence"], keep=False)]
    if not dupes.empty:
        print(f"  WARNING [{label}]: {dupes['Sequence'].nunique()} sequence(s) had multiple "
              f"docking entries; keeping the most favorable affinity for each. "
              f"Affected sequences: {sorted(dupes['Sequence'].unique().tolist())}")

    combined = combined.sort_values("_sort_key").drop_duplicates(subset=["Sequence"], keep="first")
    combined = combined.drop(columns=["_sort_key"])
    return combined


print("\n=== Step 2: loading docking results ===")
ex8 = load_docking_set(
    DOCK_EX8_DIR, "ACE_ranked_results.csv", "DPP_ranked_results.csv",
    "Dual_Activity_ranked_results.csv", "exhaustiveness_8"
)
ex32 = load_docking_set(
    DOCK_EX32_DIR, "ACE_ranked_results_2run.csv", "DPP_ranked_results_2run.csv",
    "Dual_Activity_ranked_results_2run.csv", "exhaustiveness_32"
)
print(f"exhaustiveness_8  docking entries: {len(ex8)}")
print(f"exhaustiveness_32 docking entries: {len(ex32)}")

# ---------------------------------------------------------------------------
# STEP 3: Merge docking results onto the combined ADMET table
# ---------------------------------------------------------------------------
def prefix_and_merge(base_df, dock_df, prefix):
    rename_map = {
        "Rank": f"{prefix}_Rank",
        "Peptide_ID": f"{prefix}_Peptide_ID",
        "ACE Affinity (kcal/mol)": f"{prefix}_ACE Affinity (kcal/mol)",
        "DPP-IV Affinity (kcal/mol)": f"{prefix}_DPP-IV Affinity (kcal/mol)",
        "Average Affinity": f"{prefix}_Average Affinity",
    }
    dock_renamed = dock_df.rename(columns=rename_map)
    keep_cols = ["Sequence"] + list(rename_map.values())
    keep_cols = [c for c in keep_cols if c in dock_renamed.columns]
    return base_df.merge(dock_renamed[keep_cols], on="Sequence", how="left")


print("\n=== Step 3: merging docking results into combined table ===")
merged = prefix_and_merge(merged, ex8, "exhaustiveness_8")
merged = prefix_and_merge(merged, ex32, "exhaustiveness_32")

# Flag docking peptides that never matched anything in the ADMET table (unexpected, worth knowing about)
admet_seqs = set(merged["Sequence"])
for dock_df, label in [(ex8, "exhaustiveness_8"), (ex32, "exhaustiveness_32")]:
    orphan = set(dock_df["Sequence"]) - admet_seqs
    if orphan:
        print(f"  WARNING: {len(orphan)} {label} docking peptide(s) not found in the ADMET "
              f"table (not added as new rows): {sorted(orphan)}")

# ---------------------------------------------------------------------------
# Final column order + save
# ---------------------------------------------------------------------------
for col in FINAL_COLUMN_ORDER:
    if col not in merged.columns:
        merged[col] = None

merged = merged[FINAL_COLUMN_ORDER]
merged.to_csv(OUTPUT_FILE, index=False)

print(f"\nDone. {len(merged)} peptides written to {OUTPUT_FILE}")