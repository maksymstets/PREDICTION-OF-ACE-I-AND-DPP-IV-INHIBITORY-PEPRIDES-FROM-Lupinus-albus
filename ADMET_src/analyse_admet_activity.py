import pandas as pd
from pathlib import Path

# 1. Load the ADMET-AI predictions
INPUT_DIR  = Path("/home/maksym/lupinus/ADMET_results")
OUTPUT_DIR = Path("/home/maksym/lupinus/ADMET_results")
combined_path = INPUT_DIR / "predictions.csv"

# Load data
df = pd.read_csv(combined_path)

# 2. Format columns to prevent text errors (PRESERVING ORIGINAL CASE)
df['Sequence'] = df['Sequence'].fillna('').astype(str)
df['Activity'] = df['Activity'].fillna('').astype(str) 

# === NEW: Filter for UNIQUE sequences right at the beginning ===
df = df.drop_duplicates(subset=['Sequence']).copy()

# 3. Calculate Peptide Length 
df['Length'] = df['Sequence'].apply(len)

# 4. Classify ALL Unique Peptides by Activity Type (Moved up to process starting peptides)
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

# Apply the classification logic to the starting unique dataframe
df['Activity_Class'] = df['Activity'].apply(classify_activity)

# Filter unique starting peptides to keep only target activities (remove "Other")
starting_hits = df[df['Activity_Class'] != 'Other'].copy()

# 5. Apply the Optimal ADMET Filters (Bioavailability > 0.7 and Non-toxic < 0.3)
final_hits = starting_hits[
    (starting_hits['Bioavailability_Ma'] > 0.7) & 
    (starting_hits['ClinTox'] < 0.3)
].copy()

# 6. Select columns to display and save
columns_to_save = [
    'Sequence', 'Length', 'Activity_Class', 'Activity', 'smiles', 
    'Lipinski', 'HIA_Hou', 'VDss_Lombardo', 'Half_Life_Obach', 
    'Bioavailability_Ma', 'ClinTox'
]
final_hits_clean = final_hits[columns_to_save].sort_values(by=['Activity_Class', 'Length'])

# 7. Calculate Length Statistics on the Final Hits
if not final_hits_clean.empty:
    min_len = final_hits_clean['Length'].min()
    max_len = final_hits_clean['Length'].max()
    length_counts = final_hits_clean['Length'].value_counts().to_dict()
else:
    min_len, max_len = 0, 0
    length_counts = {}

# 8. Display Results Summary (Now perfectly mirrored for Starting vs Passed)
print("==== OVERALL SUMMARY ====")
print(f"Total unique starting peptides (All classes): {len(df)}")
print(f"Total unique starting peptides (Target Activities): {len(starting_hits)}")
print(f"Passed ADMET & Activity filters: {len(final_hits_clean)}\n")

print("==== STARTING PEPTIDES BY ACTIVITY CLASS ====")
print(starting_hits['Activity_Class'].value_counts().to_string())

print("\n==== HITS BY ACTIVITY CLASS (PASSED ADMET) ====")
print(final_hits_clean['Activity_Class'].value_counts().to_string())

print("\n==== PEPTIDE LENGTH STATISTICS ====")
print(f"Minimum peptide length: {min_len}")
print(f"Maximum peptide length: {max_len}\n")

print("==== LENGTH DISTRIBUTION (Step of 1) ====")
if min_len > 0:
    for i in range(min_len, max_len + 1):
        count = length_counts.get(i, 0)
        print(f"Length {i}: {count} peptides")
else:
    print("No peptides passed the filters.")

# 9. Save to a new CSV file
out_path = OUTPUT_DIR / "optimal_admet_activity_all_lengths.csv"
final_hits_clean.to_csv(out_path, index=False)
print(f"\nSaved successfully to: {out_path}")
