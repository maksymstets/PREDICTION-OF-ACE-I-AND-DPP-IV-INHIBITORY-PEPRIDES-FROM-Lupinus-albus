from pathlib import Path
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_FILE = Path(
    "/home/maksym/lupinus/BIOPEP_results/Consolidated_A_B_coef.xlsx"
)
OUTPUT_FILE = Path(
    "/home/maksym/lupinus/BIOPEP_results/filterred_A_B_coef.xlsx"
)

TARGET_ACTIVITIES = ["ace inhibitor", "dipeptidyl peptidase iv inhibitor"]


def filter_biopep_ab_coef(input_path: Path, output_path: Path) -> pd.DataFrame:
  # Fallback path handler if running in 'lupinus_analysis' workspace
  if not input_path.exists():
    alt_path = Path(
        "/home/maksym/lupinus_analysis/BIOPEP_results/Consolidated_A_B_coef.xlsx"
    )
    if alt_path.exists():
      input_path = alt_path
    else:
      raise FileNotFoundError(f"Input file not found: {input_path}")

  print(f"Reading: {input_path}")
  df = pd.read_excel(input_path)

  # Validate required columns
  required_cols = ["ID of protein:", "NAME:", "Activity:", "A:", "B:"]
  missing_cols = [col for col in required_cols if col not in df.columns]
  if missing_cols:
    raise KeyError(f"Missing expected columns: {missing_cols}")

  # Clean string for matching
  activity_clean = df["Activity:"].fillna("").astype(str).str.strip().str.lower()

  # Filter for 'ACE inhibitor' or 'dipeptidyl peptidase IV inhibitor'
  mask = activity_clean.apply(
      lambda act: any(target in act for target in TARGET_ACTIVITIES)
  )
  df_filtered = df[mask].copy()

  # Ensure destination folder exists and save
  output_path.parent.mkdir(parents=True, exist_ok=True)
  df_filtered.to_excel(output_path, index=False)

  print(f"Filtered {len(df_filtered)} / {len(df)} rows.")
  print(f"Saved result to: {output_path}")

  return df_filtered


if __name__ == "__main__":
  filter_biopep_ab_coef(INPUT_FILE, OUTPUT_FILE)