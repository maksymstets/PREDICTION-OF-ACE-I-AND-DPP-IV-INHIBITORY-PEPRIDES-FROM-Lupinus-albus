import pandas as pd
from pathlib import Path

# ── 1. PATH CONFIGURATION ──────────────────────────────────────────────────────
INPUT_DIR  = Path("/home/maksym/lupinus/ADMET_results")
OUTPUT_DIR = Path("/home/maksym/lupinus/ADMET_results")

input_file  = INPUT_DIR / "predictions_gid_annotated.csv"
output_file = OUTPUT_DIR / "optimal_admet_gds.csv"

if not input_file.exists():
    raise FileNotFoundError(f"Could not find input file: {input_file}")

# Load full annotated dataset
df = pd.read_csv(input_file)
df.columns = df.columns.str.strip()

# ── 2. PREPROCESSING & DEDUPLICATION ───────────────────────────────────────────
df['Sequence'] = df['Sequence'].fillna('').astype(str).str.upper().str.strip()
df['Activity'] = df['Activity'].fillna('').astype(str).str.strip()
df['Protein name'] = df['Protein name'].fillna('').astype(str).str.strip()

# Deduplicate by unique Sequence
df = df.drop_duplicates(subset=['Sequence']).copy()

# Calculate Peptide Length
df['Length'] = df['Sequence'].apply(len)

# ── 3. ACTIVITY CLASSIFICATION ────────────────────────────────────────────────
def classify_activity(act_string):
    act_lower = act_string.lower()
    has_ace = 'ace' in act_lower
    has_dpp = 'dipeptidyl peptidase' in act_lower or 'dpp' in act_lower
    
    if has_ace and has_dpp:
        return 'Dual'
    elif has_ace:
        return 'ACE'
    elif has_dpp:
        return 'DPP'
    else:
        return 'Other'

df['Activity_Class'] = df['Activity'].apply(classify_activity)
starting_hits = df[df['Activity_Class'] != 'Other'].copy()

# ── 4. GDS RESISTANCE CLASSIFICATION (INFORMATIONAL ONLY) ──────────────────────
def determine_gid_resistance(row):
    r13 = str(row.get('GID resistance pH 1.3', '')).strip()
    r2  = str(row.get('GID resistance ph>2', '')).strip()
    
    is_13 = (r13.lower() not in ['not resistant', '', 'nan'])
    is_2  = (r2.lower()  not in ['not resistant', '', 'nan'])
    
    if is_13 and is_2:
        return 'Both (pH 1.3 & pH > 2)'
    elif is_13:
        return 'pH 1.3 only'
    elif is_2:
        return 'pH > 2 only'
    else:
        return 'Not Resistant'

starting_hits['GID_Resistance_Type'] = starting_hits.apply(determine_gid_resistance, axis=1)

# ── 5. ADMET FILTERING (BIOAVAILABILITY & CLINTOX ONLY) ───────────────────────
final_hits = starting_hits[
    (starting_hits['Bioavailability_Ma'] > 0.7) & 
    (starting_hits['ClinTox'] < 0.3)
].copy()

# ── 6. SELECTION OF OUTPUT COLUMNS & SAVING ───────────────────────────────────
columns_to_save = [
    'Sequence', 'Protein name', 'Length', 'Activity_Class', 'Activity', 
    'Enzyme ID', 'Peptide ID', 'smiles', 'molecular_weight', 'logP', 
    'Lipinski', 'HIA_Hou', 'VDss_Lombardo', 'Half_Life_Obach', 
    'Bioavailability_Ma', 'ClinTox',
    'LD50_Zhu',                            # Rat Toxicity (Informational)
    'GID_Resistance_Type',                # Concise GDS Resistance Summary
    'GID resistance pH 1.3',              # Raw GDS pH 1.3 column
    'GID resistance ph>2'                 # Raw GDS pH > 2 column
]

available_columns = [col for col in columns_to_save if col in final_hits.columns]
final_hits_clean = final_hits[available_columns].sort_values(by=['Activity_Class', 'Length'])

# ── 7. SUMMARY REPORTING ──────────────────────────────────────────────────────
print("==== OVERALL SUMMARY ====")
print(f"Total unique starting peptides:                  {len(df)}")
print(f"Target Activity starting peptides:             {len(starting_hits)}")
print(f"Passed Bioavailability (>0.7) & ClinTox (<0.3): {len(final_hits_clean)}\n")

print("==== GDS RESISTANCE BREAKDOWN (INFORMATIONAL) ====")
if 'GID_Resistance_Type' in final_hits_clean.columns:
    print(final_hits_clean['GID_Resistance_Type'].value_counts().to_string())

final_hits_clean.to_csv(output_file, index=False)
print(f"\nSaved successfully to: {output_file}")