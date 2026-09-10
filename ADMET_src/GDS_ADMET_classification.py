from pathlib import Path
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path("/home/maksym/lupinus/ADMET_results")
if not BASE_DIR.exists():
  BASE_DIR = Path("/home/maksym/lupinus_analysis/ADMET_results")

# Paths to inputs and output copies
CAT2_PATH = BASE_DIR / "classified_categories" / "Category2_Resistant_Intact.xlsx"
PREDICTIONS_PATH = BASE_DIR / "predictions.csv"
OPTIMAL_PATH = BASE_DIR / "optimal_admet_activity_all_lengths.csv"

PREDICTIONS_OUT = BASE_DIR / "predictions_gid_annotated.csv"
OPTIMAL_OUT = BASE_DIR / "optimal_admet_activity_all_lengths_gid_annotated.csv"


def load_category2_lookup(cat2_excel_path: Path) -> dict:
  """Builds a lookup map of peptide sequence -> GID_Mode string from Category 2."""
  if not cat2_excel_path.exists():
    raise FileNotFoundError(f"Category 2 Excel file not found: {cat2_excel_path}")

  df_cat2 = pd.read_excel(cat2_excel_path, dtype=str).fillna("")
  lookup = {}

  for _, row in df_cat2.iterrows():
    prot_name = str(row.get("Protein name", "")).strip().upper()
    seq_val = str(row.get("Sequence", "")).strip().upper()
    gid_mode = str(row.get("GID_Mode", "")).strip()

    if prot_name:
      lookup[prot_name] = gid_mode
    if seq_val:
      lookup[seq_val] = gid_mode

  return lookup


def annotate_csv_with_gid(
    csv_in: Path, csv_out: Path, resistance_map: dict
) -> pd.DataFrame:
  """Reads a CSV file, checks Category 2 resistance status for each peptide,"""
  if not csv_in.exists():
    raise FileNotFoundError(f"Input CSV not found: {csv_in}")

  df = pd.read_csv(csv_in, dtype=str)

  # Match sequence using 'Sequence' column or fallback to 'Protein name'
  match_col = (
      "Sequence"
      if "Sequence" in df.columns
      else ("Protein name" if "Protein name" in df.columns else None)
  )
  if not match_col:
    raise ValueError(f"No valid Sequence or Protein name column found in {csv_in}")

  res_ph13_list = []
  res_ph_gt2_list = []

  for _, row in df.iterrows():
    seq_key = str(row.get(match_col, "")).strip().upper()
    gid_mode = resistance_map.get(seq_key, "")

    # Check for pH 1.3 resistance
    if "1.3" in gid_mode or "1_3" in gid_mode:
      res_ph13_list.append("GID (pH 1.3)")
    else:
      res_ph13_list.append("Not Resistant")

    # Check for pH > 2 resistance
    if "pH > 2" in gid_mode or "gt_2" in gid_mode or (">" in gid_mode):
      res_ph_gt2_list.append("GID (pH > 2)")
    else:
      res_ph_gt2_list.append("Not Resistant")

  # Add requested columns
  df["GID resistance pH 1.3"] = res_ph13_list
  df["GID resistance ph>2"] = res_ph_gt2_list

  df.to_csv(csv_out, index=False)
  return df


# ── Execution ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
  print("Loading Category 2 GID resistance mappings...")
  lookup_map = load_category2_lookup(CAT2_PATH)
  print(f"Loaded {len(lookup_map)} intact resistant entries from Category 2.\n")

  # Annotate predictions.csv copy
  print(f"Processing: {PREDICTIONS_PATH.name}...")
  df_pred = annotate_csv_with_gid(PREDICTIONS_PATH, PREDICTIONS_OUT, lookup_map)
  print(f"  ✓ Saved copy to: {PREDICTIONS_OUT.name} ({len(df_pred)} rows)")

  # Annotate optimal_admet_activity_all_lengths.csv copy
  print(f"Processing: {OPTIMAL_PATH.name}...")
  df_opt = annotate_csv_with_gid(OPTIMAL_PATH, OPTIMAL_OUT, lookup_map)
  print(f"  ✓ Saved copy to: {OPTIMAL_OUT.name} ({len(df_opt)} rows)\n")

  # Summary of annotated results
  print("==================================================")
  print("               ANNOTATION SUMMARY                 ")
  print("==================================================")
  print(
      "predictions_gid_annotated.csv ->"
      f" pH 1.3 Resistant: {(df_pred['GID resistance pH 1.3'] == 'GID (pH 1.3)').sum()},"
      f" pH > 2 Resistant: {(df_pred['GID resistance ph>2'] == 'GID (pH > 2)').sum()}"
  )
  print(
      "optimal_admet_activity_all_lengths_gid_annotated.csv ->"
      f" pH 1.3 Resistant: {(df_opt['GID resistance pH 1.3'] == 'GID (pH 1.3)').sum()},"
      f" pH > 2 Resistant: {(df_opt['GID resistance ph>2'] == 'GID (pH > 2)').sum()}"
  )
  print("==================================================")