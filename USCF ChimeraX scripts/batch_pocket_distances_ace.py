import os
import glob
import subprocess
import numpy as np
from chimerax.core.commands import run

# Configuration
WORKDIR = "/home/maksym/lupinus/dock_res_2run/ACE_results"
RECEPTOR_PATH = "/home/maksym/lupinus/1O86.cif/1o86.cif"

TARGET_PEPTIDES = [
    "VRF", "WMY", "VEF", "YPR", "TKY", 
    "IAY", "IVF", "VAF", "DF", "VVF"
]

POCKET_RESIDUES = {
    "S1":  [("ALA", 354), ("GLU", 384), ("TYR", 523)],
    "S2":  [("GLN", 281), ("HIS", 353), ("LYS", 511), ("HIS", 513), ("TYR", 520)],
    "S1'": [("GLN", 162)],
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
    # Match files containing the sequence (e.g., *VRF_out.pdbqt)
    matches = glob.glob(os.path.join(WORKDIR, f"*{pep}*.pdbqt"))
    if not matches:
        print(f"[Warning] No matching PDBQT file found for peptide sequence: {pep}")
        continue

    pdbqt_file = matches[0]
    out_pdb_name = f"{pep}_mod1.pdb"
    out_pdb_path = os.path.join(WORKDIR, out_pdb_name)

    # Convert pose 1 using OpenBabel
    cmd = ["obabel", pdbqt_file, "-O", out_pdb_path, "-f", "1", "-l", "1"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"[OpenBabel] Converted: {os.path.basename(pdbqt_file)} -> {out_pdb_name}")
        converted_pdbs[pep] = out_pdb_path
    except subprocess.CalledProcessError as e:
        print(f"[Error] Failed to convert {pep}: {e.stderr.decode('utf-8')}")

# Step 2: Load Receptor
run(session, "close all")
receptor_models = run(session, f'open "{RECEPTOR_PATH}"')
receptor = receptor_models[0]

lpr_atoms = [a for a in receptor.atoms if a.residue.name == "LPR"]
zn_atoms = [a for a in receptor.atoms if a.residue.name == "ZN"]

# Step 3: Iterate through converted peptides and measure distances
for pep, pdb_path in converted_pdbs.items():
    print("\n" + "=" * 58)
    print(f" ANALYSIS FOR PEPTIDE: {pep} ({os.path.basename(pdb_path)})")
    print("=" * 58)

    pep_models = run(session, f'open "{pdb_path}"')
    peptide = pep_models[0]
    pep_atoms = peptide.atoms

    print(f"{'Pocket':6s} {'Residue':10s} {f'{pep} dist (Å)':14s} {'LPR dist (Å)':14s}")
    print("-" * 52)

    # Pocket residue measurements
    for pocket, residues in POCKET_RESIDUES.items():
        for resname, resnum in residues:
            res_atoms = [
                a for a in receptor.atoms 
                if a.residue.number == resnum and a.residue.chain_id == "A"
            ]
            if not res_atoms:
                print(f"{pocket:6s} {resname}{resnum}: not found in receptor structure")
                continue

            pep_d = min_distance(res_atoms, pep_atoms)
            lpr_d = min_distance(res_atoms, lpr_atoms)

            pep_str = f"{pep_d:.2f}" if pep_d is not None else "n/a"
            lpr_str = f"{lpr_d:.2f}" if lpr_d is not None else "n/a"
            print(f"{pocket:6s} {resname}{resnum:<6d} {pep_str:14s} {lpr_str:14s}")

    # Zinc cofactor measurement
    if zn_atoms:
        pep_zn_d = min_distance(zn_atoms, pep_atoms)
        lpr_zn_d = min_distance(zn_atoms, lpr_atoms)
        pep_zn_str = f"{pep_zn_d:.2f}" if pep_zn_d is not None else "n/a"
        lpr_zn_str = f"{lpr_zn_d:.2f}" if lpr_zn_d is not None else "n/a"
        print(f"{'Zn':6s} {'cofactor':10s} {pep_zn_str:14s} {lpr_zn_str:14s}")
    else:
        print("Zn: not found in structure")

    # Unload current peptide model to maintain a clean workspace
    run(session, f"close #{peptide.id_string}")