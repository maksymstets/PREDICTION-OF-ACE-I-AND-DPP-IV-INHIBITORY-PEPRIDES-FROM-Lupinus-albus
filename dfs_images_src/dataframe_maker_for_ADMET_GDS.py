from pathlib import Path
import re
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_DIR = Path("/home/maksym/lupinus/ADMET_results/GID_ADMET_downloads_all")

DATA_COLUMNS = [
    "NO",
    "Peptide ID",
    "Sequence",
    "Location",
    "Name",
    "Function",
    "Activity",
    "Monoisotopic mass",
    "Chemical mass",
]


def parse_filename_metadata(filepath: Path) -> tuple[str, str]:
  """Extract clean (protein_name, enzyme_id) from filename via regex pattern matching."""
  stem = filepath.stem
  # Strip trailing output tags and duplicate run indicators
  stem = re.sub(
      r"(_Calculations|_ActiveFragments|\s*\(\d+\))+$",
      "",
      stem,
      flags=re.IGNORECASE,
  )

  # Regex rules ordered from specific to general
  rules = [
      (r"_GID_pH_1_?3$", "GID (pH 1.3)"),
      (r"_GID_pH_(gt_2|>_2|>2|___2)$", "GID (pH > 2)"),
      (r"_GID$", "GID"),
      (r"_pepsin_\(pH_?1\.3\)$", "13 | pepsin (pH 1.3) | EC 3.4.23.1"),
      (r"_pepsin_\(pH.*?\)$", "39 | pepsin (pH > 2) | 3.4.23.1"),
      (r"_chymotrypsin_\(A\)$", "11 | chymotrypsin (A) | EC 3.4.21.1"),
      (r"_Stem_bromelain$", "30 | Stem bromelain | EC 3.4.22.32"),
      (r"_subtilisin$", "41 | subtilisin | EC 3.4.21.62"),
      (r"_trypsin$", "12 | trypsin | EC 3.4.21.4"),
      (r"_papain$", "24 | papain | EC 3.4.22.2"),
      (r"_ficin$", "25 | ficin | EC 3.4.22.3"),
  ]

  for pattern, enzyme_id in rules:
    match = re.search(pattern, stem, flags=re.IGNORECASE)
    if match:
      protein_name = stem[: match.start()].strip("_")
      return protein_name, enzyme_id

  # Fallback for triple-underscore format (e.g. protein_11___chymotrypsin...)
  if "___" in stem:
    parts = stem.split("___")
    protein_name = re.sub(r"_\d+$", "", parts[0]).strip("_")
    enz_name = (
        parts[1].replace("_", " ").strip() if len(parts) > 1 else "unknown"
    )
    return protein_name, enz_name

  protein_name = re.sub(r"_\d+$", "", stem).strip("_")
  return protein_name, "unknown"


def parse_single_file(filepath: Path) -> tuple[pd.DataFrame, bool]:
  """Parse one BIOPEP ActiveFragments Excel file.

  Returns a tuple of (DataFrame, is_empty_flag). If file has no active sequence
  rows, returns a 1-row DataFrame preserving metadata with empty sequence
  fields.
  """
  protein_name, filename_enzyme_id = parse_filename_metadata(filepath)

  try:
    raw = pd.read_excel(filepath, header=None, dtype=str)
  except Exception as e:
    raise RuntimeError(f"Read error: {e}")

  # ── Extract Enzyme ID from metadata cells ──────────────────────────────────
  sheet_enzyme_id = None
  for row_idx in range(min(5, len(raw))):
    row_vals = raw.iloc[row_idx].dropna().tolist()
    row_str = " ".join(str(v) for v in row_vals).lower()
    if "enzyme" in row_str and "id" in row_str:
      for val in raw.iloc[row_idx]:
        if pd.notna(val) and "enzyme" not in str(val).lower():
          sheet_enzyme_id = str(val).strip()
          break
      break

  # Prioritize exact GID profile label from filename if present
  if "GID" in filepath.name:
    final_enzyme_id = filename_enzyme_id
  else:
    final_enzyme_id = (
        sheet_enzyme_id
        if sheet_enzyme_id and sheet_enzyme_id.lower() != "unknown"
        else filename_enzyme_id
    )

  # ── Find data header row containing "NO" ───────────────────────────────────
  header_row_idx = None
  for row_idx in range(len(raw)):
    first_cell = str(raw.iloc[row_idx, 0]).strip()
    if first_cell == "NO":
      header_row_idx = row_idx
      break

  # Extract data slice if header exists
  data_df = pd.DataFrame()
  if header_row_idx is not None:
    data_df = pd.read_excel(
        filepath, header=header_row_idx, dtype=str, keep_default_na=False
    )
    existing_cols = [c for c in DATA_COLUMNS if c in data_df.columns]
    data_df = data_df[existing_cols].copy()
    data_df.dropna(how="all", inplace=True)
    data_df.reset_index(drop=True, inplace=True)

  is_empty = len(data_df) == 0

  # ── Handle Empty Files: Generate single placeholder row ────────────────────
  if is_empty:
    empty_row = {col: "" for col in DATA_COLUMNS}
    data_df = pd.DataFrame([empty_row])

  # ── Inject Metadata Columns ────────────────────────────────────────────────
  data_df.insert(0, "Protein name", protein_name)
  data_df.insert(1, "Enzyme ID", final_enzyme_id)
  data_df.insert(2, "Source file", filepath.name)

  return data_df, is_empty


def natural_sort_key(path: Path) -> list:
  """Sort key that orders numbers numerically, not alphabetically."""
  return [
      int(c) if c.isdigit() else c.lower()
      for c in re.split(r"(\d+)", path.name)
  ]


def parse_all_files(input_dir: Path) -> pd.DataFrame:
  """Iterate over all *_ActiveFragments.xlsx files, concatenate, and report summary statistics."""
  files = sorted(input_dir.glob("*_ActiveFragments.xlsx"), key=natural_sort_key)
  if not files:
    raise FileNotFoundError(
        f"No *_ActiveFragments.xlsx files found in {input_dir}"
    )

  frames = []
  successful_files = 0
  empty_files_count = 0

  for fp in files:
    try:
      df, is_empty = parse_single_file(fp)
      frames.append(df)
      successful_files += 1

      if is_empty:
        empty_files_count += 1
        print(f"  ✓  {fp.name}  →  0 active fragments (placeholder added)")
      else:
        print(f"  ✓  {fp.name}  →  {len(df)} fragments parsed")

    except Exception as e:
      print(f"  ✗  {fp.name}  →  ERROR: {e}")

  if not frames:
    raise RuntimeError("No files were successfully parsed.")

  pct_empty = (
      (empty_files_count / successful_files) * 100 if successful_files > 0 else 0
  )

  print("\n==== FILE ANALYSIS SUMMARY ====")
  print(f"Total files found:          {len(files)}")
  print(f"Successfully parsed files:  {successful_files}")
  print(f"Files with 0 fragments:     {empty_files_count} ({pct_empty:.2f}%)")
  print(
      "Files with active data:    "
      f" {successful_files - empty_files_count} ({100.0 - pct_empty:.2f}%)"
  )
  print("===============================\n")

  return pd.concat(frames, ignore_index=True)


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
  print(f"Scanning: {INPUT_DIR}\n")
  combined_df = parse_all_files(INPUT_DIR)

  print(f"Final DataFrame shape: {combined_df.shape}")
  print(combined_df.head(10).to_string())

  out_path = Path(
      "/home/maksym/lupinus/ADMET_results/GID_ADMET_combined_active_fragments.xlsx"
  )
  out_path.parent.mkdir(parents=True, exist_ok=True)
  combined_df.to_excel(out_path, index=False)
  print(f"\nSaved → {out_path}")