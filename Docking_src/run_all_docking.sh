#!/bin/bash

# Enable nullglob so empty matches resolve to nothing instead of literal strings
shopt -s nullglob

PEPTIDES_DIR="/home/maksym/lupinus/Docking_inputs"
RESULTS_DIR="/home/maksym/lupinus/docking_results"
CONFIG_ACE="/home/maksym/lupinus/Docking_inputs/config_ACE.txt"
CONFIG_DPP="/home/maksym/lupinus/Docking_inputs/config_DPP.txt"

mkdir -p "$RESULTS_DIR/ACE_results"
mkdir -p "$RESULTS_DIR/DPP_results"

run_vina() {
    local ligand="$1"
    local config="$2"
    local out_dir="$3"
    local base_name=$(basename "$ligand" .pdbqt)

    echo "Running Vina for: $base_name"
    vina --config "$config" --ligand "$ligand" --out "$out_dir/${base_name}_out.pdbqt" > "$out_dir/${base_name}_log.txt"
}

echo "--- STARTING DOCKING PIPELINE ---"

# Function to process folder with check
process_folder() {
    local folder_path="$1"
    local config="$2"
    local out_dir="$3"
    local files=("$folder_path"/*.pdbqt)

    if [ ${#files[@]} -eq 0 ]; then
        echo "WARNING: No .pdbqt files found in $folder_path"
        return
    fi

    for pep in "${files[@]}"; do
        run_vina "$pep" "$config" "$out_dir"
    done
}

echo "Processing ACE_only folder..."
process_folder "$PEPTIDES_DIR/ACE_only_pdbqt" "$CONFIG_ACE" "$RESULTS_DIR/ACE_results"

echo "Processing DPP_only folder..."
process_folder "$PEPTIDES_DIR/DPP_only_pdbqt" "$CONFIG_DPP" "$RESULTS_DIR/DPP_results"

echo "Processing Dual_activity folder..."
dual_files=("$PEPTIDES_DIR/Dual"/*.pdbqt)
if [ ${#dual_files[@]} -eq 0 ]; then
    echo "WARNING: No .pdbqt files found in $PEPTIDES_DIR/Dual"
else
    for pep in "${dual_files[@]}"; do
        echo "--> Dual docking for $(basename "$pep" .pdbqt)"
        run_vina "$pep" "$CONFIG_ACE" "$RESULTS_DIR/ACE_results"
        run_vina "$pep" "$CONFIG_DPP" "$RESULTS_DIR/DPP_results"
    done
fi

echo "Check $RESULTS_DIR for your output poses and energy logs."
