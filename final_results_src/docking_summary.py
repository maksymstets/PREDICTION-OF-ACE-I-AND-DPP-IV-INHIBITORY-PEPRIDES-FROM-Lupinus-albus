"""
Docking + ADMET/GID summary, in the style of:
"Molecular docking analysis was performed on all ACE-I (n = X) and DPP-IV
(n = Y) inhibitory peptides... The binding energy... ranged from A to B
kcal/mol. Of them, Z% (n = W) were stable during in silico gastrointestinal
digestion and showed optimal ADMET properties. Their binding energies ranged
from C to D kcal/mol."

Computed separately for exhaustiveness=8 and exhaustiveness=32.



Requirements:
    pip install pandas
"""

import pandas as pd
from pathlib import Path

INPUT_FILE = Path("/home/maksym/lupinus/final_results_gds/ultimate_peptide_results_gds.csv")

df = pd.read_csv(INPUT_FILE)

# ---------------------------------------------------------------------------
# "Optimal ADMET" - admetlab3_final_results.csv was already filtered to only
# f20<0.3 & ROA<0.3 peptides before being merged into this file, so presence
# of these columns already implies a pass. Sanity-check that assumption too.
# ---------------------------------------------------------------------------
df["_optimal_admet"] = df["f20"].notna() & df["ROA"].notna()

sanity_fail = df[df["_optimal_admet"] & ((df["f20"] >= 0.3) | (df["ROA"] >= 0.3))]
if not sanity_fail.empty:
    print(f"WARNING: {len(sanity_fail)} row(s) have f20/ROA data but don't actually "
          f"satisfy f20<0.3 & ROA<0.3 - the 'already filtered' assumption may not hold "
          f"for this file. Re-check manually: {sanity_fail['Sequence'].tolist()}")

# ---------------------------------------------------------------------------
# GID stability - dual resistant (both pH conditions), coalescing native and
# _ADMETLab3 versions of the columns since either may be populated depending
# on which source file a given row originated from.
# ---------------------------------------------------------------------------
def resistant(row, native_col, lab3_col, label):
    val = row[native_col] if pd.notna(row[native_col]) else row[lab3_col]
    return val == label

df["_resistant_pH13"] = df.apply(
    lambda r: resistant(r, "GID resistance pH 1.3", "GID resistance pH 1.3_ADMETLab3", "GID (pH 1.3)"),
    axis=1
)
df["_resistant_pH2"] = df.apply(
    lambda r: resistant(r, "GID resistance ph>2", "GID resistance ph>2_ADMETLab3", "GID (pH > 2)"),
    axis=1
)
df["_gid_stable"] = df["_resistant_pH13"] & df["_resistant_pH2"]

df["_hit"] = df["_optimal_admet"] & df["_gid_stable"]

# ---------------------------------------------------------------------------
# Report per exhaustiveness level, per enzyme
# ---------------------------------------------------------------------------
EXHAUSTIVENESS_LEVELS = ["exhaustiveness_8", "exhaustiveness_32"]
ENZYMES = {
    "ACE-I": "ACE Affinity (kcal/mol)",
    "DPP-IV": "DPP-IV Affinity (kcal/mol)",
}

for ex in EXHAUSTIVENESS_LEVELS:
    print(f"\n{'='*70}\n{ex.upper()}\n{'='*70}")
    for enzyme_label, affinity_suffix in ENZYMES.items():
        col = f"{ex}_{affinity_suffix}"
        if col not in df.columns:
            print(f"  Column not found: {col}")
            continue

        docked = df[df[col].notna()]
        n_docked = len(docked)
        if n_docked == 0:
            print(f"\n{enzyme_label}: no docked peptides found for {ex}.")
            continue

        e_min, e_max = docked[col].min(), docked[col].max()

        hits = docked[docked["_hit"]]
        n_hits = len(hits)
        pct_hits = round(100 * n_hits / n_docked, 2) if n_docked else 0

        print(f"\n{enzyme_label} (n = {n_docked}):")
        print(f"  Binding energy range (all docked): {e_min:.1f} to {e_max:.1f} kcal/mol")

        if n_hits > 0:
            h_min, h_max = hits[col].min(), hits[col].max()
            print(f"  Stable in GI digestion + optimal ADMET: {pct_hits}% (n = {n_hits})")
            print(f"  Binding energy range (stable + optimal subset): {h_min:.1f} to {h_max:.1f} kcal/mol")
        else:
            print(f"  Stable in GI digestion + optimal ADMET: 0% (n = 0)")

print(f"\n{'='*70}")
