"""
Compare exhaustiveness=8 vs exhaustiveness=32 docking runs for ACE-I and
DPP-IV inhibitory peptides:
  1. Binding energy range per enzyme, per exhaustiveness level
  2. Per-peptide difference in binding energy between the two runs
  3. Whether the top candidates (best binders) differ between runs
  4. Comparison against reference drug binding energies:
       ACE-I:   Lisinopril  = -8.6 kcal/mol
       DPP-IV:  Saxagliptin = -8.4 kcal/mol, Vildagliptin = -8.8 kcal/mol

Requirements:
    pip install pandas
"""

import pandas as pd
from pathlib import Path

INPUT_FILE = Path("/home/maksym/lupinus/final_results_gds/ultimate_peptide_results_gds.csv")

TOP_N = 10  # how many top candidates to compare for overlap between runs

REFERENCE_DRUGS = {
    "ACE Affinity (kcal/mol)": {"Lisinopril": -8.6},
    "DPP-IV Affinity (kcal/mol)": {"Saxagliptin": -8.4, "Vildagliptin": -8.8},
}

df = pd.read_csv(INPUT_FILE)


def analyze_enzyme(df, affinity_suffix, enzyme_label):
    col8 = f"exhaustiveness_8_{affinity_suffix}"
    col32 = f"exhaustiveness_32_{affinity_suffix}"

    if col8 not in df.columns or col32 not in df.columns:
        print(f"  Columns not found for {enzyme_label}: {col8} / {col32}")
        return

    d8 = df[["Sequence", col8]].dropna().rename(columns={col8: "affinity_8"})
    d32 = df[["Sequence", col32]].dropna().rename(columns={col32: "affinity_32"})

    print(f"\n{'-'*60}\n{enzyme_label}\n{'-'*60}")
    print(f"Docked at exhaustiveness=8:  n = {len(d8)}, "
          f"range = {d8['affinity_8'].min():.1f} to {d8['affinity_8'].max():.1f} kcal/mol")
    print(f"Docked at exhaustiveness=32: n = {len(d32)}, "
          f"range = {d32['affinity_32'].min():.1f} to {d32['affinity_32'].max():.1f} kcal/mol")

    # ---- Per-peptide difference between runs (only peptides docked in both) ----
    both = d8.merge(d32, on="Sequence", how="inner")
    only_8 = set(d8["Sequence"]) - set(d32["Sequence"])
    only_32 = set(d32["Sequence"]) - set(d8["Sequence"])

    if only_8:
        print(f"  NOTE: {len(only_8)} peptide(s) docked at ex=8 only (missing from ex=32): "
              f"{sorted(only_8)}")
    if only_32:
        print(f"  NOTE: {len(only_32)} peptide(s) docked at ex=32 only (missing from ex=8): "
              f"{sorted(only_32)}")

    if not both.empty:
        both["delta"] = both["affinity_32"] - both["affinity_8"]
        print(f"\nPer-peptide difference (ex32 - ex8), n = {len(both)} peptides in both runs:")
        print(f"  Mean delta   : {both['delta'].mean():+.3f} kcal/mol")
        print(f"  Std delta    : {both['delta'].std():.3f} kcal/mol")
        print(f"  Min delta    : {both['delta'].min():+.3f} kcal/mol "
              f"({both.loc[both['delta'].idxmin(), 'Sequence']})")
        print(f"  Max delta    : {both['delta'].max():+.3f} kcal/mol "
              f"({both.loc[both['delta'].idxmax(), 'Sequence']})")
        n_unchanged = (both["delta"].abs() < 0.1).sum()
        print(f"  Peptides with |delta| < 0.1 kcal/mol (essentially unchanged): "
              f"{n_unchanged} / {len(both)} ({100*n_unchanged/len(both):.1f}%)")

    # ---- Top-N candidate overlap between runs ----
    top8 = set(d8.nsmallest(TOP_N, "affinity_8")["Sequence"])
    top32 = set(d32.nsmallest(TOP_N, "affinity_32")["Sequence"])
    overlap = top8 & top32
    only_top8 = top8 - top32
    only_top32 = top32 - top8

    print(f"\nTop {TOP_N} candidates (strongest binders):")
    print(f"  Overlap between ex=8 and ex=32: {len(overlap)} / {TOP_N} peptides: {sorted(overlap)}")
    if only_top8:
        print(f"  In top {TOP_N} at ex=8 only: {sorted(only_top8)}")
    if only_top32:
        print(f"  In top {TOP_N} at ex=32 only: {sorted(only_top32)}")

    # Best single candidate per run
    best8 = d8.loc[d8["affinity_8"].idxmin()]
    best32 = d32.loc[d32["affinity_32"].idxmin()]
    same_best = best8["Sequence"] == best32["Sequence"]
    print(f"\nBest single candidate at ex=8:  {best8['Sequence']} ({best8['affinity_8']:.1f} kcal/mol)")
    print(f"Best single candidate at ex=32: {best32['Sequence']} ({best32['affinity_32']:.1f} kcal/mol)")
    print(f"  Same peptide is #1 in both runs: {'YES' if same_best else 'NO'}")

    # ---- Comparison against reference drugs ----
    refs = REFERENCE_DRUGS.get(affinity_suffix, {})
    for drug, threshold in refs.items():
        n8 = (d8["affinity_8"] <= threshold).sum()
        n32 = (d32["affinity_32"] <= threshold).sum()
        print(f"\nPeptides matching/outperforming {drug} ({threshold} kcal/mol):")
        print(f"  ex=8 : {n8} / {len(d8)} ({100*n8/len(d8):.1f}%)")
        print(f"  ex=32: {n32} / {len(d32)} ({100*n32/len(d32):.1f}%)")


print("="*70)
print("DOCKING COMPARISON: exhaustiveness=8 vs exhaustiveness=32")
print("="*70)

analyze_enzyme(df, "ACE Affinity (kcal/mol)", "ACE-I")
analyze_enzyme(df, "DPP-IV Affinity (kcal/mol)", "DPP-IV")