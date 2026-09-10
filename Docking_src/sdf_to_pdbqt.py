"""
Convert PubChem-derived peptide SDF structures to PDBQT and sort them into
activity-based folders (ACE_only / DPP_only / Dual), all in one pass.

Uses pubchem_lookup_summary.csv (produced by pubchem_structure_download.py)
as the authoritative source for which peptide maps to which SDF file, rather
than reconstructing filenames - this avoids mismatches if files were later
moved/renamed, and cleanly skips peptides that never got a downloaded SDF.

Requirements:
    pip install rdkit meeko pandas
"""

import os
import pandas as pd
from pathlib import Path
from rdkit import Chem
from meeko import MoleculePreparation

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SDF_DIR = Path("/home/maksym/lupinus/Docking_inputs/pubchem_3d_sdf_merged")
SUMMARY_FILE = Path("/home/maksym/lupinus/Docking_inputs/pubchem_lookup_summary_consolidated.csv")
METADATA_FILE = Path("/home/maksym/lupinus/ADMET_inputs/optimised_peptides_for_admet_metadata.csv")

BASE_OUTPUT_DIR = Path("/home/maksym/lupinus/Docking_inputs")
ace_dir = BASE_OUTPUT_DIR / "ACE_only_pdbqt"
dpp_dir = BASE_OUTPUT_DIR / "DPP_only_pdbqt"
dual_dir = BASE_OUTPUT_DIR / "Dual"
for d in (ace_dir, dpp_dir, dual_dir):
    d.mkdir(parents=True, exist_ok=True)

FINAL_CSV = BASE_OUTPUT_DIR / "pdbqt_structures.csv"

TARGET_DIR_BY_CATEGORY = {
    "ACE_only_pdbqt": ace_dir,
    "DPP_only_pdbqt": dpp_dir,
    "Dual": dual_dir,
}


def normalize_id(x):
    """Normalize a Peptide ID for merging, handling int/float('7583.0')/str mismatches."""
    try:
        f = float(x)
        if f.is_integer():
            return str(int(f))
        return str(x).strip()
    except (TypeError, ValueError):
        return str(x).strip()


def classify_activity(activity):
    """Same substring logic as the original sorting script, just remapped to new folder names."""
    if not isinstance(activity, str):
        return None
    has_ace = "ACE inhibitor" in activity
    has_dpp = "dipeptidyl peptidase IV inhibitor" in activity
    if has_ace and has_dpp:
        return "Dual"
    elif has_ace:
        return "ACE_only_pdbqt"
    elif has_dpp:
        return "DPP_only_pdbqt"
    else:
        return None


def resolve_sdf_path(raw_path):
    """
    Given the sdf_file path recorded in the summary CSV, return a Path that
    actually exists, falling back to looking up the same filename inside the
    current SDF_DIR (handles the folder having been moved/renamed since the
    summary CSV was written).
    """
    if pd.isna(raw_path) or not str(raw_path).strip():
        return None
    candidate = Path(str(raw_path))
    if candidate.exists():
        return candidate
    fallback = SDF_DIR / candidate.name
    return fallback if fallback.exists() else None


# ---------------------------------------------------------------------------
# Load & merge
# ---------------------------------------------------------------------------
summary_df = pd.read_csv(SUMMARY_FILE)
meta_df = pd.read_csv(METADATA_FILE)

summary_df["Peptide ID"] = summary_df["Peptide ID"].apply(normalize_id)
meta_df["Peptide ID"] = meta_df["Peptide ID"].apply(normalize_id)

merged = summary_df.merge(
    meta_df[["Peptide ID", "Activity"]],
    on="Peptide ID",
    how="left",
)

# ---------------------------------------------------------------------------
# Convert + sort
# ---------------------------------------------------------------------------
preparer = MoleculePreparation()
locations = []

n_converted = 0
n_no_sdf = 0
n_conv_error = 0
n_unrecognized = 0

print("Beginning batch preparation of peptides...")

for _, row in merged.iterrows():
    pid = row["Peptide ID"]
    seq = row["Sequence"]
    activity = row.get("Activity")

    sdf_path = resolve_sdf_path(row.get("sdf_file"))
    if sdf_path is None:
        n_no_sdf += 1
        locations.append("no_sdf_available")
        continue

    category = classify_activity(activity)
    if category is None:
        print(f"⚠️ Unrecognized/missing activity pattern for {pid} ({seq}): {activity!r}")
        n_unrecognized += 1
        locations.append("unrecognized_activity")
        continue

    target_dir = TARGET_DIR_BY_CATEGORY[category]

    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
        if mol is None:
            raise ValueError("RDKit could not parse a valid molecule from this SDF")

        preparer.prepare(mol)
        pdbqt_string = preparer.write_pdbqt_string()

        out_name = f"{pid}_{seq}.pdbqt".replace(" ", "_")
        out_path = target_dir / out_name
        with open(out_path, "w") as f:
            f.write(pdbqt_string)

        locations.append(category)
        n_converted += 1

    except Exception as e:
        print(f" ✗ Error preparing {pid} ({seq}): {e}")
        locations.append("conversion_error")
        n_conv_error += 1

merged["Location"] = locations
merged.to_csv(FINAL_CSV, index=False)

print(f"\nDone! Converted {n_converted} peptides to .pdbqt")
print("Breakdown by Location:")
print(merged["Location"].value_counts().to_string())
print(f"\nSkipped (no SDF available): {n_no_sdf}")
print(f"Unrecognized activity pattern: {n_unrecognized}")
print(f"Conversion errors: {n_conv_error}")
print(f"\nSummary CSV with Location column saved to {FINAL_CSV}")