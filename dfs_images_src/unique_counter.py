from pathlib import Path
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
CSV_PATH = Path(
    "/home/maksym/lupinus/ADMET_inputs/optimised_peptides_for_admet_metadata.csv"
)

if not CSV_PATH.exists():
  raise FileNotFoundError(f"Input file not found: {CSV_PATH}")

# ── Step 1: Load and Clean Data ───────────────────────────────────────────────
df = pd.read_csv(CSV_PATH, dtype=str).fillna("")

df["Sequence"] = df["Sequence"].astype(str).str.strip().str.upper()
df["Activity"] = df["Activity"].astype(str).str.strip()

# Exclude empty sequence rows
df = df[df["Sequence"] != ""]

# ── Step 2: Aggregate Activities per Unique Peptide Sequence ──────────────────
# Combines activities across rows if a sequence appears multiple times
df_unique = df.groupby("Sequence", as_index=False).agg(
    {"Activity": lambda x: "; ".join(set(v for v in x if v))}
)

# ── Step 3: Classify Activities ───────────────────────────────────────────────


def is_ace_inhibitor(act_str: str) -> bool:
  return "ace inhibitor" in act_str.lower()


def is_dpp4_inhibitor(act_str: str) -> bool:
  act_lower = act_str.lower()
  return (
      "dipeptidyl peptidase iv inhibitor" in act_lower
      or "dpp-iv" in act_lower
      or "dpp4" in act_lower
  )


df_unique["Is_ACE"] = df_unique["Activity"].apply(is_ace_inhibitor)
df_unique["Is_DPP4"] = df_unique["Activity"].apply(is_dpp4_inhibitor)
df_unique["Is_Dual"] = df_unique["Is_ACE"] & df_unique["Is_DPP4"]

# ── Step 4: Compute Metrics ───────────────────────────────────────────────────
total_unique = len(df_unique)
ace_count = int(df_unique["Is_ACE"].sum())
dpp4_count = int(df_unique["Is_DPP4"].sum())
dual_count = int(df_unique["Is_Dual"].sum())

# ── Step 5: Terminal Output ───────────────────────────────────────────────────
print("\n" + "=" * 55)
print("       UNIQUE PEPTIDES ACTIVITY ANALYSIS REPORT        ")
print("=" * 55)
print(f"Total Unique Peptides:                      {total_unique}")
print(f"ACE Inhibitor Peptides:                     {ace_count}")
print(f"Dipeptidyl Peptidase IV Inhibitor Peptides: {dpp4_count}")
print(f"Dual-Active Peptides (ACE + DPP-IV):        {dual_count}")
print("=" * 55 + "\n")