This repo contains scripts for bulk BIOPEP digestion simulation, ADMET analysis via ADMETLab 3.0 and molecular docking via AutoDockVina.

The results are under review at Biomedical Engineering and Technology. 

## Folders and files description

| Directory       | Purpose                       | Key Contents                                            |
| --------------- | ------------------------------| ------------------------------------------------------- |
| ADMET_inputs    | Peptides for ADMET analysis   | optimised_peptides_for_admet_metadata.csv               |
| ADMET_results   | Results of an ADMET analysis  | predictions.csv, optimal_admet_activity_all_length.csv  |
| ADMET_src       | Scripts for ADMET analysis    | analyse_admet_activity.py, smiles_converter_metadata.py, ADMET-AI terminal commands|
| BIOPEP_src      | Scripts for BIOPEP automation | script_profiles.py, sort_reports.py, rename.py |
| Docking_src     | Preparatory scripts for docking | affinity.py, peptide_preparation.py|
| dfs_images_src | Dataframe- and image-making scripts | protein_protease_heatmap.py, dataframe_maker_for_active_fragments.py|
| images_results | Figures presented | peptide_table.png, protein_per_protease_heatmap.py
| inputs_fasta | Aminoacid sequences of lupin conglutins | All files here were used
| ACE ChimeraX commands and terminal meeko script | ChimeraX commands for the preparation for docking of a 3D structure of ACE and meeko script| adding polar hydrogens and charges and (2) removal of water molecules and crystallographic ligands|
| ACE_clean.pdb | ACE .pdb file prepared for docking | 
| ACE_prepared.pdbqt|
| DPP IV ChimeraX commands and terminal meeko script | ChimeraX commands for the preparation for docking of a 3D structure of DPP-IV and meeko script|
| DPP4_clean.pdb | DPP-IV .pdb file prepared for docking |
| DPP4_prepared.pdbqt
| README.md| Description| |
| config_ACE.txt| Configuration file for ACE docking| |
| config_DPP.txt| Configuration file for DPP-IV docking | |
| run_all_docking.sh| Bash script for docking
| uniprotkb_accession_A0A6A4QEF0_OR_access_2026_05_29.fasta| Original fasta file with lupin conglutin sequences| |

| Directory       | Purpose                       | Key Contents                                            |
| --------------- | ------------------------------| ------------------------------------------------------- |
| 1O86.cif 
| 4PNZ.cif
| ADMETLab_3.0_results
| ADMET_Lab_3.0_src
| ADMET_inputs
| ADMET_results
| ADMET_src
| BIOPEP_results
| BIOPEP_src
| Docking_inputs
| Docking_src
| USCF ChimeraX outputs
| USCF ChimeraX scripts
| dfs_images_src
| dock_inp_2run
| dock_res_2run
| docking_results
| final_results
| final_results_gds
| final_results_src
| images_results
| inputs_fasta | Aminoacid sequences of lupin conglutins | All files here were used. Original fasta file with lupin conglutin sequences - uniprotkb_accession_A0A6A4QEF0_OR_access_2026_05_29.fasta |
