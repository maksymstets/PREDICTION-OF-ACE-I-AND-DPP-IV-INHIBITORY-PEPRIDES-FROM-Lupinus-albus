import os
import glob
import subprocess
import numpy as np
from chimerax.core.commands import run

# Configuration
WORKDIR = "/home/maksym/lupinus/dock_res_2run/DPP_results"
RECEPTOR_PATH = "/home/maksym/lupinus/4PNZ.cif/4pnz.cif"
OUTPUT_DIR = "/home/maksym/lupinus/dock_res_2run/DPP_pdb"
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGET_PEPTIDES = [
    "WQ", "EYF", "HW", "WF", "WR",
    "NH", "AF", "YF", "HF", "AY"
]

POCKET_RESIDUES = {
    "S1":        [("TYR", 547), ("TYR", 631), ("VAL", 656), ("TRP", 659), ("TYR", 662), ("VAL", 711)],
    "Catalytic": [("SER", 630), ("ASP", 708), ("ASN", 710)],
    "S2":        [("GLU", 205), ("GLU", 206), ("ARG", 125)],
    "S2'":       [("VAL", 207), ("SER", 209), ("ARG", 358), ("PHE", 357)]
}

def min_distance(atoms1, atoms2):
    """Calculates minimum Euclidean distance between two sets of ChimeraX atoms."""
    if not atoms1 or not atoms2:
        return None
    coords1 = np.array([a.coord for a in atoms1])
    coords2 = np.array([a.coord for a in atoms2])
    dists = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=-1)
    return float(np.min(dists))

# Step 1: OpenBabel Conversion (.pdbqt -> _mod1.pdb)
converted_pdbs = {}

for pep in TARGET_PEPTIDES:
    matches = glob.glob(os.path.join(WORKDIR, f"*{pep}*.pdbqt"))
    if not matches:
        print(f"[Warning] No matching PDBQT file found for peptide sequence: {pep}")
        continue

    pdbqt_file = matches[0]
    out_pdb_name = f"{pep}_mod1.pdb"
    out_pdb_path = os.path.join(OUTPUT_DIR, out_pdb_name)

    cmd = ["obabel", pdbqt_file, "-O", out_pdb_path, "-f", "1", "-l", "1"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"[OpenBabel] Converted: {os.path.basename(pdbqt_file)} -> {out_pdb_name}")
        converted_pdbs[pep] = out_pdb_path
    except subprocess.CalledProcessError as e:
        print(f"[Error] Failed to convert {pep}: {e.stderr.decode('utf-8')}")


# Step 2: Load Receptor - must run BEFORE anything references `receptor`
run(session, "close all")

if os.path.isdir(RECEPTOR_PATH):
    found_files = glob.glob(os.path.join(RECEPTOR_PATH, "*.cif")) + glob.glob(os.path.join(RECEPTOR_PATH, "*.pdb"))
    if found_files:
        RECEPTOR_PATH = found_files[0]

receptor_models = run(session, f'open "{RECEPTOR_PATH}"')
receptor = receptor_models[0]

# Reference ligand (Omarigliptin: 2VH) and Sodium ion (NA) atom selections
vh2_atoms = [a for a in receptor.atoms if a.residue.name == "2VH"]
na_atoms = [a for a in receptor.atoms if a.residue.name == "NA"]

# Step 3: Iterate through converted peptides and measure distances
for pep, pdb_path in converted_pdbs.items():
    print("\n" + "=" * 60)
    print(f" ANALYSIS FOR PEPTIDE: {pep} ({os.path.basename(pdb_path)})")
    print("=" * 60)

    pep_models = run(session, f'open "{pdb_path}"')
    peptide = pep_models[0]
    pep_atoms = peptide.atoms

    print(f"{'Pocket':11s} {'Residue':10s} {f'{pep} dist (Å)':14s} {'2VH dist (Å)':14s} {'Chain'}")
    print("-" * 60)

    for pocket, residues in POCKET_RESIDUES.items():
        for resname, resnum in residues:
            # Check both chains - peptide could be docked near either monomer's pocket
            best_dist = None
            best_chain = None
            best_res_atoms = None
            for chain in ("A", "B"):
                res_atoms = [
                    a for a in receptor.atoms
                    if a.residue.number == resnum and a.residue.chain_id == chain
                ]
                if not res_atoms:
                    continue
                d = min_distance(res_atoms, pep_atoms)
                if d is not None and (best_dist is None or d < best_dist):
                    best_dist = d
                    best_chain = chain
                    best_res_atoms = res_atoms

            if best_dist is None:
                print(f"{pocket:11s} {resname}{resnum}: not found in either chain")
                continue

            vh2_chain_atoms = [a for a in vh2_atoms if a.residue.chain_id == best_chain]
            vh2_d = min_distance(best_res_atoms, vh2_chain_atoms)   # residue <-> 2VH, not peptide <-> 2VH

            pep_str = f"{best_dist:.2f}"
            vh2_str = f"{vh2_d:.2f}" if vh2_d is not None else "n/a"
            print(f"{pocket:11s} {resname}{resnum:<6d} {pep_str:14s} {vh2_str:14s} {best_chain}")
    # Sodium ion distance check (also checking both chains)
    if na_atoms:
        best_na_dist = None
        best_na_chain = None
        for chain in ("A", "B"):
            chain_na = [a for a in na_atoms if a.residue.chain_id == chain]
            d = min_distance(chain_na, pep_atoms)
            if d is not None and (best_na_dist is None or d < best_na_dist):
                best_na_dist = d
                best_na_chain = chain

        if best_na_dist is not None:
            vh2_na_d = min_distance(
                [a for a in vh2_atoms if a.residue.chain_id == best_na_chain],
                pep_atoms
            )
            vh2_na_str = f"{vh2_na_d:.2f}" if vh2_na_d is not None else "n/a"
            print(f"{'Na+':11s} {'ion':10s} {best_na_dist:.2f}          {vh2_na_str:14s} {best_na_chain}")

    run(session, f"close #{peptide.id_string}")
