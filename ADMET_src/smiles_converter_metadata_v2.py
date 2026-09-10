import os

import pandas as pd
from pathlib import Path
from rdkit import Chem

# 1. Load Excel and FORCE Pandas to treat "NA" as text, not a missing value
INPUT_DIR  = Path("/home/maksym/lupinus/BIOPEP_results")
OUTPUT_DIR = Path("/home/maksym/lupinus/ADMET_inputs")
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
combined_path = INPUT_DIR / "combined_active_fragments.xlsx"
df = pd.read_excel(combined_path, dtype=str)

# 2. Filter for your target activities safely
is_ace = df['Activity'].str.contains("ACE inhibitor", case=False, na=False)
is_dpp = df['Activity'].str.contains("dipeptidyl peptidase IV inhibitor", case=False, na=False)
filtered_df = df[is_ace | is_dpp].copy()

print(f"Initial filter: Found {len(filtered_df)} active row entries.")

# 3. Group by Sequence and aggregate metadata to handle all duplication types
aggregation_rules = {
    'Protein name': lambda x: ", ".join(sorted(x.dropna().astype(str).unique())),
    'Activity': lambda x: ", ".join(sorted(x.dropna().astype(str).unique())),
    'Enzyme ID': lambda x: ", ".join(sorted(str(v) for v in x.dropna().unique())),
    'Peptide ID': 'first'  # Keeps the first ID as a reference
}

# This collapses all 3 duplicate scenarios perfectly
grouped_df = filtered_df.groupby('Sequence', as_index=False).agg(aggregation_rules)
print(f"Deduplication complete: {len(grouped_df)} unique peptide sequences retained.")

# 4. SMILES Converter
def sequence_to_smiles(seq):
    if not seq or pd.isna(seq):
        return None
    clean_seq = str(seq).strip().upper()
    mol = Chem.MolFromSequence(clean_seq)
    return Chem.MolToSmiles(mol) if mol else None

# 5. Apply and Export
grouped_df['smiles'] = grouped_df['Sequence'].apply(sequence_to_smiles)
grouped_df = grouped_df.dropna(subset=['smiles'])

output_file = "optimised_peptides_for_admet_metadata.csv"
out_path = OUTPUT_DIR / output_file
grouped_df.to_csv(out_path, index=False)

print(f"Success! Generated SMILES for active peptides and saved to '{out_path}'.")

