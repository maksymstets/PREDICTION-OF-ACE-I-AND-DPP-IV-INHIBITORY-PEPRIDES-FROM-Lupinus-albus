import os
import csv

# Define paths matching your setup
RESULTS_DIR = "/home/maksym/lupinus/docking_results"
ACE_DIR = os.path.join(RESULTS_DIR, "ACE_results")
DPP_DIR = os.path.join(RESULTS_DIR, "DPP_results")
OUTPUT_DIR = "/home/maksym/lupinus/docking_results"

def parse_vina_log(file_path):
    """Parses a Vina log file and extracts the Mode 1 (best) affinity score."""
    try:
        with open(file_path, "r") as f:
            for line in f:
                parts = line.split()
                # Vina data rows look like: "1      -7.4      0.000      0.000"
                if len(parts) >= 2 and parts[0] == "1":
                    return float(parts[1])
    except Exception as e:
        print(f"Warning: Could not read {file_path}. Error: {e}")
    return None

def process_folder(folder_path, target_name):
    """Scans a folder for log files and compiles peptide affinities."""
    results = {}
    if not os.path.exists(folder_path):
        print(f"Directory not found: {folder_path}")
        return results

    for file_name in os.listdir(folder_path):
        if file_name.endswith("_log.txt"):
            peptide_name = file_name.replace("_log.txt", "")
            file_path = os.path.join(folder_path, file_name)
            affinity = parse_vina_log(file_path)
            if affinity is not None:
                results[peptide_name] = affinity
    return results

def save_ranked_csv(data_dict, output_path, target_label):
    """Sorts data from lowest to highest affinity and writes to a CSV file."""
    # Sort by value (lowest energy = strongest binding)
    sorted_data = sorted(data_dict.items(), key=lambda x: x[1])
    
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "Peptide_ID", f"{target_label} Affinity (kcal/mol)"])
        for rank, (peptide, affinity) in enumerate(sorted_data, 1):
            writer.writerow([rank, peptide, affinity])
    print(f"Saved ranked results to: {output_path}")

print("--- STARTING DOCKING DATA EXTRACTION ---")

# 1. Process ACE and DPP data
ace_data = process_folder(ACE_DIR, "ACE")
dpp_data = process_folder(DPP_DIR, "DPP")

print(f"Found {len(ace_data)} valid ACE log files.")
print(f"Found {len(dpp_data)} valid DPP log files.")

# 2. Export individual ranked files
save_ranked_csv(ace_data, os.path.join(OUTPUT_DIR, "ACE_ranked_results.csv"), "ACE")
save_ranked_csv(dpp_data, os.path.join(OUTPUT_DIR, "DPP_ranked_results.csv"), "DPP-IV")

# 3. Identify and rank Dual-Activity Peptides
dual_peptides = set(ace_data.keys()).intersection(set(dpp_data.keys()))
print(f"Identified {len(dual_peptides)} peptides with dual-activity testing.")

if dual_peptides:
    dual_results = []
    for pep in dual_peptides:
        # Calculate an average or simply list both for comparison
        avg_affinity = (ace_data[pep] + dpp_data[pep]) / 2
        dual_results.append((pep, ace_data[pep], dpp_data[pep], avg_affinity))
    
    # Sort dual results by the best average affinity
    dual_results.sort(key=lambda x: x[3])
    
    dual_csv_path = os.path.join(OUTPUT_DIR, "Dual_Activity_ranked_results.csv")
    with open(dual_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "Peptide_ID", "ACE Affinity (kcal/mol)", "DPP-IV Affinity (kcal/mol)", "Average Affinity"])
        for rank, (pep, ace, dpp, avg) in enumerate(dual_results, 1):
            writer.writerow([rank, pep, ace, dpp, round(avg, 2)])
    print(f"Saved compiled dual-activity rankings to: {dual_csv_path}")

print("--- EXTRACTION COMPLETE ---")
