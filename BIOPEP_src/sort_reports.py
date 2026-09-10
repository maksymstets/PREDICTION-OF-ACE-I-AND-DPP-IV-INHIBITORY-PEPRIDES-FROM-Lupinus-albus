import os
import shutil
import glob
#I manually combined all calculation and activefragments from sequentilan and GID in one file for easier sorting. This script will move all calculation reports to a new folder for easier access and organization.
# ==================== CONFIGURATION ====================
# Source folder containing all the generated outputs
SOURCE_DIR = os.path.join(os.path.expanduser("~"), "lupinus", "BIOPEP_results", "BIOPEP_Batch_Outputs_Combined")

# New target folder specifically for your calculation reports
TARGET_DIR = os.path.join(os.path.expanduser("~"), "lupinus", "BIOPEP_results", "BIOPEP_Calculations_Outputs")
# =======================================================

# Create the new folder if it doesn't exist yet
if not os.path.exists(TARGET_DIR):
    os.makedirs(TARGET_DIR)
    print(f"-> Created new folder: {TARGET_DIR}")

# Search for any file ending with '_Calculations.xlsx'
search_pattern = os.path.join(SOURCE_DIR, "*_Calculations.xlsx")
calculation_files = glob.glob(search_pattern)

if not calculation_files:
    print(f"No calculation reports found in {SOURCE_DIR}.")
    exit()

print(f"Found {len(calculation_files)} calculation reports to organize.\n")

moved_counter = 0

# Loop through and move each calculation file
for file_path in calculation_files:
    file_name = os.path.basename(file_path)
    destination_path = os.path.join(TARGET_DIR, file_name)
    
    try:
        shutil.move(file_path, destination_path)
        print(f" [MOVED] {file_name}")
        moved_counter += 1
    except Exception as e:
        print(f" [ERROR] Could not move {file_name}: {e}")

print(f"\n*** Sorting complete! {moved_counter} files moved to: {TARGET_DIR} ***")
