import re
from pathlib import Path
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_DIR = Path(
    "/home/maksym/lupinus/BIOPEP_results/BIOPEP_Batch_Outputs_Combined"
)
OUTPUT_FILE = Path(
    "/home/maksym/lupinus/BIOPEP_results/combined_active_fragments.xlsx"
)

# Expected data columns in BIOPEP export tables
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

# Standardized Mapping Table for Filename Extraction
ENZYME_MAPPING = {
    "GID_pH_1_3": "GID (pH 1.3)",
    "GID_pH_gt_2": "GID (pH > 2)",
    "pepsin_(pH_1.3)": "13 | pepsin (pH 1.3) | EC 3.4.23.1",
    "pepsin_(pH___2)": "39 | pepsin (pH > 2) | 3.4.23.1",
    "chymotrypsin_(A)": "11 | chymotrypsin (A) | EC 3.4.21.1",
    "Stem_bromelain": "30 | Stem bromelain | EC 3.4.22.32",
    "subtilisin": "41 | subtilisin | EC 3.4.21.62",
    "trypsin": "12 | trypsin | EC 3.4.21.4",
    "papain": "24 | papain | EC 3.4.22.2",
    "ficin": "25 | ficin | EC 3.4.22.3",
}


def parse_filename_metadata(filepath: Path) -> tuple[str, str]:
  """Extract clean (protein_name, enzyme_id) from filename."""
  stem = filepath.stem
  stem = re.sub(r"_ActiveFragments$", "", stem)

  # Match longest enzyme suffix pattern first
  sorted_patterns = sorted(ENZYME_MAPPING.keys(), key=len, reverse=True)

  for pattern in sorted_patterns:
    if stem.endswith("_" + pattern):
      protein_name = stem[: -len("_" + pattern)]
      enzyme_id = ENZYME_MAPPING[pattern]
      return protein_name, enzyme_id

  return stem, "unknown"


def parse_single_file(filepath: Path) -> pd.DataFrame:
  """Parse one BIOPEP ActiveFragments Excel file and return a structured DataFrame."""
  raw = pd.read_excel(filepath, header=None, dtype=str)

  # Extract metadata from filename
  protein_name, filename_enzyme_id = parse_filename_metadata(filepath)

  # Check sheet header for single-enzyme runs
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

  # GID multi-enzyme runs must use filename metadata to override BIOPEP's first-enzyme default label
  if "GID" in filepath.name:
    enzyme_id = filename_enzyme_id
  else:
    enzyme_id = sheet_enzyme_id if sheet_enzyme_id else filename_enzyme_id

  # Locate table header row
  header_row_idx = None
  for row_idx in range(len(raw)):
    first_cell = str(raw.iloc[row_idx, 0]).strip()
    if first_cell == "NO":
      header_row_idx = row_idx
      break

  if header_row_idx is None:
    raise ValueError(f"Could not find data header row in {filepath.name}")

  # Read data table with NA protection
  data_df = pd.read_excel(
      filepath, header=header_row_idx, dtype=str, keep_default_na=False
  )

  # Keep expected columns
  existing_cols = [c for c in DATA_COLUMNS if c in data_df.columns]
  data_df = data_df[existing_cols].copy()

  # Clean rows
  data_df.dropna(how="all", inplace=True)
  data_df.reset_index(drop=True, inplace=True)

  # Insert metadata on left
  data_df.insert(0, "Protein name", protein_name)
  data_df.insert(1, "Enzyme ID", enzyme_id)
  data_df.insert(2, "Source file", filepath.name)

  return data_df


def natural_sort_key(path: Path) -> list:
  """Numerical natural sort key for sequence ordering."""
  return [
      int(c) if c.isdigit() else c.lower()
      for c in re.split(r"(\d+)", path.name)
  ]


def parse_all_files(input_dir: Path) -> pd.DataFrame:
  """Iterate over all *_ActiveFragments.xlsx files, parse, and combine."""
  files = sorted(input_dir.glob("*_ActiveFragments.xlsx"), key=natural_sort_key)

  if not files:
    raise FileNotFoundError(
        f"No *_ActiveFragments.xlsx files found in {input_dir}"
    )

  frames = []
  for fp in files:
    try:
      df_single = parse_single_file(fp)
      frames.append(df_single)
      print(f"  ✓  {fp.name:<55} → {len(df_single)} rows")
    except Exception as e:
      print(f"  ✗  {fp.name:<55} → ERROR: {e}")

  if not frames:
    raise RuntimeError("No files were successfully parsed.")

  return pd.concat(frames, ignore_index=True)


# ── Main Script ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
  print(f"Scanning directory: {INPUT_DIR}\n")
  combined_df = parse_all_files(INPUT_DIR)

  print(f"\nFinal DataFrame shape: {combined_df.shape}")
  print(combined_df.head(10).to_string())

  OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
  combined_df.to_excel(OUTPUT_FILE, index=False)

  print(f"\nSaved combined file → {OUTPUT_FILE}")