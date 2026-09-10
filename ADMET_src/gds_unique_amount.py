from pathlib import Path
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_FILE = Path(
    "/home/maksym/lupinus/ADMET_results/GID_ADMET_combined_active_fragments.xlsx"
)
OUTPUT_DIR = Path("/home/maksym/lupinus/ADMET_results/classified_categories")

# Define separate output file paths for each category
PATH_CAT1 = OUTPUT_DIR / "Category1_Cleaved_to_AA.xlsx"
PATH_CAT2 = OUTPUT_DIR / "Category2_Resistant_Intact.xlsx"
PATH_CAT3 = OUTPUT_DIR / "Category3_Smaller_Peptides.xlsx"


def detect_gid_mode(row: pd.Series) -> str:
  """Identifies the digestion pH profile from row metadata."""
  text = f"{row.get('Enzyme ID', '')} {row.get('Source file', '')}".lower()
  if "1.3" in text or "1_3" in text:
    return "GID (pH 1.3)"
  elif "2" in text or "gt_2" in text:
    return "GID (pH > 2)"
  return "Unknown GID"


def aggregate_peptides(df_sub: pd.DataFrame, group_key: str) -> pd.DataFrame:
  """Groups by key_col and combines unique entries for Activity, Name, Function,

  and other metadata separated by '; '.
  """

  if df_sub.empty:
    return df_sub.copy()

  def join_unique(series):
    vals = [
        str(v).strip()
        for v in series
        if pd.notna(v) and str(v).strip() != ""
    ]
    seen = set()
    unique_vals = []
    for v in vals:
      if v not in seen:
        seen.add(v)
        unique_vals.append(v)
    return "; ".join(unique_vals)

  agg_dict = {
      col: join_unique for col in df_sub.columns if col != group_key
  }
  return df_sub.groupby(group_key, as_index=False).agg(agg_dict)


def calculate_category_stats(
    df_raw_cat: pd.DataFrame, key_col: str, cat_name: str
) -> dict:
  """Calculates row counts and pH profile comparisons for a given category."""
  if df_raw_cat.empty:
    return {
        "Category": cat_name,
        "Total Rows": 0,
        "Total Unique": 0,
        "pH 1.3 Total": 0,
        "pH > 2 Total": 0,
        "pH 1.3 Only": 0,
        "pH > 2 Only": 0,
        "Both Profiles": 0,
    }

  seq_1_3 = set(
      df_raw_cat[df_raw_cat["GID_Mode"] == "GID (pH 1.3)"][key_col]
      .str.strip()
      .dropna()
      .unique()
  )
  seq_gt_2 = set(
      df_raw_cat[df_raw_cat["GID_Mode"] == "GID (pH > 2)"][key_col]
      .str.strip()
      .dropna()
      .unique()
  )

  seq_1_3.discard("")
  seq_gt_2.discard("")

  return {
      "Category": cat_name,
      "Total Rows": len(df_raw_cat),
      "Total Unique": len(seq_1_3 | seq_gt_2),
      "pH 1.3 Total": len(seq_1_3),
      "pH > 2 Total": len(seq_gt_2),
      "pH 1.3 Only": len(seq_1_3 - seq_gt_2),
      "pH > 2 Only": len(seq_gt_2 - seq_1_3),
      "Both Profiles": len(seq_1_3 & seq_gt_2),
  }


# ── Step 1: Load and Clean Raw Data ──────────────────────────────────────────
if not INPUT_FILE.exists():
  raise FileNotFoundError(f"Missing input file: {INPUT_FILE}")

print(f"Loading raw combined dataset from: {INPUT_FILE.name}")
df_raw = pd.read_excel(INPUT_FILE, dtype=str).fillna("")

for col in ["Sequence", "Protein name", "Activity", "Name", "Function"]:
  if col in df_raw.columns:
    df_raw[col] = df_raw[col].astype(str).str.strip()

df_raw["Sequence"] = df_raw["Sequence"].str.upper()
df_raw["Protein name"] = df_raw["Protein name"].str.upper()
df_raw["GID_Mode"] = df_raw.apply(detect_gid_mode, axis=1)

# ── Step 2: Categorize Raw Rows ──────────────────────────────────────────────

# Category 1: Cleaved to individual amino acids (Sequence empty, Protein name present)
mask_cat1 = (df_raw["Sequence"] == "") & (df_raw["Protein name"] != "")
df_cat1_raw = df_raw[mask_cat1].copy()

# Category 2: Resistant to GID (Sequence == Protein name)
mask_cat2 = (df_raw["Sequence"] != "") & (
    df_raw["Sequence"] == df_raw["Protein name"]
)
df_cat2_raw = df_raw[mask_cat2].copy()

# Category 3: Cleaved to smaller active fragments (Sequence != Protein name and non-empty)
mask_cat3 = (
    (df_raw["Sequence"] != "")
    & (df_raw["Protein name"] != "")
    & (df_raw["Sequence"] != df_raw["Protein name"])
)
df_cat3_raw = df_raw[mask_cat3].copy()

# ── Step 3: Build Actual Peptide DataFrames with Aggregated Activities ───────

df_cat1_final = aggregate_peptides(df_cat1_raw, group_key="Protein name")
df_cat2_final = aggregate_peptides(df_cat2_raw, group_key="Protein name")
df_cat3_final = aggregate_peptides(df_cat3_raw, group_key="Sequence")

# ── Step 4: Export to 3 Separate Excel Files ─────────────────────────────────
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df_cat1_final.to_excel(PATH_CAT1, index=False)
df_cat2_final.to_excel(PATH_CAT2, index=False)
df_cat3_final.to_excel(PATH_CAT3, index=False)

print(f"\nSaved separate Excel files to directory: {OUTPUT_DIR}")
print(f"  -> {PATH_CAT1.name}")
print(f"  -> {PATH_CAT2.name}")
print(f"  -> {PATH_CAT3.name}")

# ── Step 5: Terminal Report with Category Descriptions & Stats ───────────────
stats_cat1 = calculate_category_stats(
    df_cat1_raw, key_col="Protein name", cat_name="Category 1: Cleaved to AAs"
)
stats_cat2 = calculate_category_stats(
    df_cat2_raw, key_col="Protein name", cat_name="Category 2: GID Resistant"
)
stats_cat3 = calculate_category_stats(
    df_cat3_raw, key_col="Sequence", cat_name="Category 3: Smaller Peptides"
)

df_summary = pd.DataFrame([stats_cat1, stats_cat2, stats_cat3])

print("\n" + "=" * 85)
print(
    "                    GASTROINTESTINAL DIGESTION ANALYSIS REPORT                 "
)
print("=" * 85)

print("\nCATEGORY DESCRIPTIONS:")
print(
    "  * Category 1 (Fully Cleaved to AAs): Parent peptides completely"
    " degraded into"
)
print(
    "    individual amino acids during digestion (yielding no remaining active"
    " fragments)."
)
print(
    "  * Category 2 (Resistant / Intact): Parent peptides that remained 100%"
    " intact and"
)
print(
    "    unmodified after gastrointestinal digestion (Sequence == Protein"
    " name)."
)
print(
    "  * Category 3 (Smaller Peptides): New active peptide fragments generated"
    " by enzymatic"
)
print(
    "    hydrolysis of parent sequences (Sequence != Protein name).\n"
)

print("STATISTICAL COMPARISON (pH 1.3 vs pH > 2):")
print(df_summary.to_string(index=False))
print("=" * 85 + "\n")