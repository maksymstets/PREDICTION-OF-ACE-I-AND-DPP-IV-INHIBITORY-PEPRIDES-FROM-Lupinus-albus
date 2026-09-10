import pandas as pd
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
CSV_PATH = Path("/home/maksym/lupinus_analysis/ADMET_results/optimal_admet_activity_all_lengths.csv")

# Load dataset
if not CSV_PATH.exists():
    raise FileNotFoundError(f"Could not find the file: {CSV_PATH}")

df = pd.read_csv(CSV_PATH)

# Ensure the Sequence column is treated cleanly as uppercase strings
df['Sequence'] = df['Sequence'].fillna('').astype(str).str.upper()
total_peptides = len(df)

print(f"=========================================")
print(f"ANALYSIS SUMMARY (Total Peptides: {total_peptides})")
print(f"=========================================\n")

# 1. Lack Lysine (K)
lack_lysine = df[~df['Sequence'].str.contains('K', na=False)]
count_lack_k = len(lack_lysine)
pct_lack_k = (count_lack_k / total_peptides) * 100
print(f"1. Peptides lacking Lysine (K):")
print(f"   → Count: {count_lack_k} out of {total_peptides}")
print(f"   → Percentage: {pct_lack_k:.2f}%\n")

# 2a. Standard Aliphatic Definition (G, A, V, L, I, P)
aliphatic_standard = 'A|V|L|I|G|P'
has_aliphatic_std = df[df['Sequence'].str.contains(aliphatic_standard, na=False)]
count_aliph_std = len(has_aliphatic_std)
pct_aliph_std = (count_aliph_std / total_peptides) * 100
print(f"2a. Peptides with >=1 Aliphatic AA (Standard: A, V, L, I, G, P):")
print(f"   → Count: {count_aliph_std} out of {total_peptides}")
print(f"   → Percentage: {pct_aliph_std:.2f}%\n")

# 2b. Strict Explicit Examples Only (V, P, G, I)
aliphatic_explicit = 'V|P|G|I'
has_aliphatic_exp = df[df['Sequence'].str.contains(aliphatic_explicit, na=False)]
count_aliph_exp = len(has_aliphatic_exp)
pct_aliph_exp = (count_aliph_exp / total_peptides) * 100
print(f"2b. Peptides with >=1 Aliphatic AA (Strictly from examples: V, P, G, I):")
print(f"   → Count: {count_aliph_exp} out of {total_peptides}")
print(f"   → Percentage: {pct_aliph_exp:.2f}%\n")

# 3. Have Proline (P)
has_proline = df[df['Sequence'].str.contains('P', na=False)]
count_p = len(has_proline)
pct_p = (count_p / total_peptides) * 100
print(f"3. Peptides containing Proline (P):")
print(f"   → Count: {count_p} out of {total_peptides}")
print(f"   → Percentage: {pct_p:.2f}%\n")
print(f"=========================================")