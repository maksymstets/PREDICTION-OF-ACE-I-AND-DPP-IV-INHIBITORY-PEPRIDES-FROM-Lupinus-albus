from pathlib import Path
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_FILE = Path("/home/maksym/lupinus/BIOPEP_results/hydrolysis_coef.xlsx")
OUTPUT_EXCEL = Path(
    "/home/maksym/lupinus/BIOPEP_results/filtered_hydrolysis_coef.xlsx"
)
OUTPUT_TXT = Path(
    "/home/maksym/lupinus/BIOPEP_results/hydrolysis_coef_min_max_summary.txt"
)

TARGET_ACTIVITIES = ["ace inhibitor", "dipeptidyl peptidase iv inhibitor"]
METRIC_COLS = ["Ae", "W", "Be", "V", "DHt[%]"]


def process_hydrolysis_coef(
    input_path: Path, excel_out: Path, txt_out: Path
) -> None:
  # Path fallback handler
  if not input_path.exists():
    alt_path = Path(
        "/home/maksym/lupinus_analysis/BIOPEP_results/hydrolysis_coef.xlsx"
    )
    if alt_path.exists():
      input_path = alt_path
    else:
      raise FileNotFoundError(f"Input file not found: {input_path}")

  print(f"Loading data from: {input_path}")
  df = pd.read_excel(input_path)

  # ── Step 1: Filter Activities ──────────────────────────────────────────────
  activity_clean = (
      df["Activity"].fillna("").astype(str).str.strip().str.lower()
  )
  mask = activity_clean.apply(
      lambda act: any(target in act for target in TARGET_ACTIVITIES)
  )
  df_filtered = df[mask].copy()

  # Ensure metric columns are numeric for accurate min/max calculations
  for col in METRIC_COLS:
    if col in df_filtered.columns:
      df_filtered[col] = pd.to_numeric(df_filtered[col], errors="coerce")

  # Save filtered dataset to Excel
  excel_out.parent.mkdir(parents=True, exist_ok=True)
  df_filtered.to_excel(excel_out, index=False)
  print(
      f"  ✓ Saved filtered dataset ({len(df_filtered)} rows) to: {excel_out}"
  )

  # ── Step 2: Calculate Min/Max Stats per Enzyme ID ──────────────────────────
  txt_lines = []
  txt_lines.append("=" * 85)
  txt_lines.append("        HYDROLYSIS COEFFICIENTS (Ae, W, Be, V, DHt[%]) SUMMARY")
  txt_lines.append("=" * 85 + "\n")

  for enzyme_id, group in df_filtered.groupby("Enzyme id"):
    txt_lines.append(f"ENZYME ID: {enzyme_id}")
    txt_lines.append("-" * 70)

    for metric in METRIC_COLS:
      if metric not in group.columns:
        continue

      valid = group.dropna(subset=[metric])
      if valid.empty:
        txt_lines.append(f"  {metric:8s} -> No numeric data available")
        continue

      # Maximum value and corresponding protein name(s)
      max_val = valid[metric].max()
      max_prots = valid[valid[metric] == max_val]["Protein Name"].unique()
      max_prots_str = "; ".join(str(p) for p in max_prots if pd.notna(p))

      # Minimum value and corresponding protein name(s)
      min_val = valid[metric].min()
      min_prots = valid[valid[metric] == min_val]["Protein Name"].unique()
      min_prots_str = "; ".join(str(p) for p in min_prots if pd.notna(p))

      txt_lines.append(f"  Metric: {metric}")
      txt_lines.append(
          f"    • Highest (Max): {max_val:<10.4f} | Protein(s): {max_prots_str}"
      )
      txt_lines.append(
          f"    • Lowest  (Min): {min_val:<10.4f} | Protein(s): {min_prots_str}"
      )

    txt_lines.append("\n" + "." * 70 + "\n")

  # Save summary text report
  with open(txt_out, "w", encoding="utf-8") as f:
    f.write("\n".join(txt_lines))

  print(f"  ✓ Saved min/max summary text report to: {txt_out}")


if __name__ == "__main__":
  process_hydrolysis_coef(INPUT_FILE, OUTPUT_EXCEL, OUTPUT_TXT)