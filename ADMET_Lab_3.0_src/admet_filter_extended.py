"""
Reproduce the reference paper's "Table 2" style summary: for each parameter
(HIA, F20%, F30%, VD, T1/2, ROAT, Lipinski, GDS), count and percentage of
peptides meeting the "optimal" criterion, broken down by General / DPP-IV /
ACE-I / Dual activity class.

IMPORTANT: uses the RAW (unfiltered) ADMETLab3 results, not the already
-filtered admetlab3_final_results.csv - Table 2's denominators (e.g. "General
n=201") are the full evaluated peptide population, not a pre-filtered subset.
Running this against an already-filtered file will overstate every percentage.

Thresholds used (see accompanying explanation - adjust at the top if your
study defines these differently):
  - HIA, F20%, F30%, ROAT: optimal if probability < 0.3
    (ADMETLab's own "good" direction is low-probability for all four; the 0.3
    cutoff itself follows the convention from your quoted methods text)
  - VDss: optimal if 0.04 <= VDss <= 20 L/kg (ADMETLab's own documented range)
  - Lipinski: optimal if violations <= 1 (ADMETLab's own documented rule)
  - T1/2: NOT DEFINED by ADMETLab as a single cutoff - only categorical bins
    are documented (ultra-short <1h, short 1-4h, intermediate 4-8h, long >8h).
    No default is applied here; T1/2_OPTIMAL_MAX_HOURS below is left as None
    until you decide/confirm the paper's actual cutoff.
  - GDS: pulled from your own GID resistance columns, not from ADMETLab.
    Optimal if resistant at BOTH pH conditions.

Requirements:
    pip install pandas
"""

import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RAW_FILE = Path("/home/maksym/lupinus/ADMETLab_3.0_results/admetlab3_raw_results_gid_annotated.csv")
OUTPUT_FILE = Path("/home/maksym/lupinus/ADMETLab_3.0_results/table2_admet_frequency_summary.csv")

# Not needed anymore - RAW_FILE above already has GID resistance columns merged in.
GID_FILE = None

# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------
PROB_OPTIMAL_MAX = 0.3          # HIA, F20%, F30%, ROAT: optimal if value < this
VDSS_MIN, VDSS_MAX = 0.04, 20   # L/kg, per ADMETLab's own documented range
LIPINSKI_MAX_VIOLATIONS = 1     # per ADMETLab's own documented rule
T1_2_OPTIMAL_MAX_HOURS = None   # NOT DEFINED by ADMETLab - set this yourself if needed


def classify_activity(act_string):
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


df = pd.read_csv(RAW_FILE)
df["Activity_Class"] = df["Activity"].apply(classify_activity)

GID_COLS = ["GID resistance pH 1.3", "GID resistance ph>2"]
has_gid_cols = all(c in df.columns for c in GID_COLS)

if not has_gid_cols and GID_FILE is not None and GID_FILE.exists():
    gid_df = pd.read_csv(GID_FILE)
    missing = [c for c in GID_COLS if c not in gid_df.columns]
    if missing:
        print(f"WARNING: GID_FILE missing columns {missing} - skipping GDS row.")
    else:
        df = df.merge(gid_df[["Sequence"] + GID_COLS], on="Sequence", how="left")
        has_gid_cols = True

# ---------------------------------------------------------------------------
# Per-parameter "optimal" boolean columns
# ---------------------------------------------------------------------------
df["_opt_HIA"] = df["hia"] < PROB_OPTIMAL_MAX
df["_opt_F20"] = df["f20"] < PROB_OPTIMAL_MAX
df["_opt_F30"] = df["f30"] < PROB_OPTIMAL_MAX
df["_opt_VD"] = df["VDss"].between(VDSS_MIN, VDSS_MAX)
df["_opt_ROAT"] = df["ROA"] < PROB_OPTIMAL_MAX
df["_opt_LIPINSKI"] = df["Lipinski"] <= LIPINSKI_MAX_VIOLATIONS

if T1_2_OPTIMAL_MAX_HOURS is not None:
    df["_opt_T1_2"] = df["t0.5"] < T1_2_OPTIMAL_MAX_HOURS
else:
    df["_opt_T1_2"] = pd.NA  # left blank - no defined cutoff

if has_gid_cols:
    df["_opt_GDS"] = (df["GID resistance pH 1.3"] == True) & (df["GID resistance ph>2"] == True)
else:
    df["_opt_GDS"] = pd.NA

PARAM_COLS = {
    "HIA": "_opt_HIA",
    "F20%": "_opt_F20",
    "F30%": "_opt_F30",
    "VD": "_opt_VD",
    "T 1/2": "_opt_T1_2",
    "ROAT": "_opt_ROAT",
    "LIPINSKI": "_opt_LIPINSKI",
    "GDS": "_opt_GDS",
}

# ---------------------------------------------------------------------------
# Build the Table 2 style summary
# ---------------------------------------------------------------------------
groups = {
    "General": df,
    "DPP-IV": df[df["Activity_Class"].isin(["DPP", "Dual"])],
    "ACE-I": df[df["Activity_Class"].isin(["ACE", "Dual"])],
}

rows = []
header_row = {"Parameter": ""}
for label, subset in groups.items():
    header_row[f"{label} n_total"] = len(subset)
rows.append(header_row)

for display_name, col in PARAM_COLS.items():
    row = {"Parameter": display_name}
    for label, subset in groups.items():
        n_total = len(subset)
        if subset[col].isna().all():
            row[f"{label}_n"] = None
            row[f"{label}_%"] = None
            continue
        n_pass = int(subset[col].sum())
        pct = round(100 * n_pass / n_total, 2) if n_total else 0
        row[f"{label}_n"] = n_pass
        row[f"{label}_%"] = pct
    rows.append(row)

summary_df = pd.DataFrame(rows)
summary_df.to_csv(OUTPUT_FILE, index=False)

print(f"Saved Table 2 style summary to {OUTPUT_FILE}\n")
print(summary_df.to_string(index=False))

if T1_2_OPTIMAL_MAX_HOURS is None:
    print("\nNOTE: T 1/2 row is blank - ADMETLab does not define a single 'optimal' "
          "cutoff for half-life, only categorical bins. Set T1_2_OPTIMAL_MAX_HOURS "
          "at the top of this script once you've confirmed the paper's actual cutoff.")
if not has_gid_cols:
    print("NOTE: GDS row is blank - no GID resistance columns found in RAW_FILE, "
          "and no fallback GID_FILE was set.")