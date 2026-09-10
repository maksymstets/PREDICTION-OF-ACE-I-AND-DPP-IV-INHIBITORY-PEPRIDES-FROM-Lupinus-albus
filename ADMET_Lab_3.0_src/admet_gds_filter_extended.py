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
  - Lipinski: optimal if violations <= 1 (ADMETLab's own documented rule).
    Value is rounded to the nearest int before comparing, to guard against
    float artifacts (e.g. 0.999999997 or 1.00000001 from the API response)
    being misclassified at the boundary.
  - T1/2: optimal if t1/2 > 3 hours (your chosen cutoff for a sufficiently
    long, non-instantaneously-cleared half-life). ADMETLab itself does not
    define a single cutoff - only categorical bins (see prior explanation).
  - GID resistance: pulled from your own GID resistance columns, not from
    ADMETLab. Reported as three separate rows since resistance is tracked
    per pH condition, not as a single yes/no:
      - "GID (pH 1.3)"       : resistant at pH 1.3
      - "GID (pH > 2)"       : resistant at pH > 2
      - "GID Dual Resistant" : resistant at BOTH conditions

Requirements:
    pip install pandas
"""

import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RAW_FILE = Path("/home/maksym/lupinus/ADMETLab_3.0_results/admetlab3_raw_results_gid_annotated.csv")
OUTPUT_FILE = Path("/home/maksym/lupinus/ADMETLab_3.0_results/table_admet_gds_frequency_summary.csv")



# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------
PROB_OPTIMAL_MAX = 0.3          # HIA, F20%, F30%, ROAT: optimal if value < this
VDSS_MIN, VDSS_MAX = 0.04, 20   # L/kg, per ADMETLab's own documented range
LIPINSKI_MAX_VIOLATIONS = 1     # per ADMETLab's own documented rule (0-1 = acceptable)
T1_2_OPTIMAL_MIN_HOURS = 3      # optimal = "high" half-life, i.e. t1/2 > 3h

# GID resistance columns store categorical labels, not booleans - these are
# the exact string values that mark a peptide as resistant at each condition.
GID_PH13_COL = "GID resistance pH 1.3"
GID_PH2_COL = "GID resistance ph>2"
GID_PH13_RESISTANT_LABEL = "GID (pH 1.3)"
GID_PH2_RESISTANT_LABEL = "GID (pH > 2)"


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

missing_gid = [c for c in (GID_PH13_COL, GID_PH2_COL) if c not in df.columns]
if missing_gid:
    raise ValueError(
        f"Missing expected GID resistance columns {missing_gid} in {RAW_FILE}. "
        f"Available columns: {df.columns.tolist()}"
    )

# ---------------------------------------------------------------------------
# Per-parameter "optimal" boolean columns
# ---------------------------------------------------------------------------
df["_opt_HIA"] = df["hia"] < PROB_OPTIMAL_MAX
df["_opt_F20"] = df["f20"] < PROB_OPTIMAL_MAX
df["_opt_F30"] = df["f30"] < PROB_OPTIMAL_MAX
df["_opt_VD"] = df["VDss"].between(VDSS_MIN, VDSS_MAX)
df["_opt_ROAT"] = df["ROA"] < PROB_OPTIMAL_MAX
# Round before comparing - guards against float artifacts like 0.0000001
# from the API response being treated as a "violation" by a strict == 0 check.
df["_opt_LIPINSKI"] = df["Lipinski"].round().astype(int) <= LIPINSKI_MAX_VIOLATIONS

df["_opt_T1_2"] = df["t0.5"] > T1_2_OPTIMAL_MIN_HOURS

# Compare against the actual categorical labels, not True/False - this was the
# bug behind the previous all-zero GDS output.
df["_opt_GID_pH13"] = df[GID_PH13_COL] == GID_PH13_RESISTANT_LABEL
df["_opt_GID_pH2"] = df[GID_PH2_COL] == GID_PH2_RESISTANT_LABEL
df["_opt_GID_dual"] = df["_opt_GID_pH13"] & df["_opt_GID_pH2"]

PARAM_COLS = {
    "HIA": "_opt_HIA",
    "F20%": "_opt_F20",
    "F30%": "_opt_F30",
    "VD": "_opt_VD",
    "T 1/2": "_opt_T1_2",
    "ROAT": "_opt_ROAT",
    "LIPINSKI": "_opt_LIPINSKI",
    "GID (pH 1.3)": "_opt_GID_pH13",
    "GID (pH > 2)": "_opt_GID_pH2",
    "GID Dual Resistant": "_opt_GID_dual",
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
        n_pass = int(subset[col].sum())
        pct = round(100 * n_pass / n_total, 2) if n_total else 0
        row[f"{label}_n"] = n_pass
        row[f"{label}_%"] = pct
    rows.append(row)

summary_df = pd.DataFrame(rows)
summary_df.to_csv(OUTPUT_FILE, index=False)

print(f"Saved Table 2 style summary to {OUTPUT_FILE}\n")
print(summary_df.to_string(index=False))