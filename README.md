# This is a repository for scripts used for an article "IN SILICO PREDICTION OF ANTIDIABETIC AND ANTIHYPERTENSIVE PEPTIDES FROM LUPINUS ALBUS". 
The paper is under review at Biomedical Engineering and Technology. 
The repository contains all scripts used for this research, including intermediate and final versions as well as raw data. 




## Folders and files description


| Directory       | Purpose                       | Key Contents                                            |
| --------------- | ------------------------------| ------------------------------------------------------- |
| 1O86.cif | Crystallographic structure of an ACE-I| |
| 4PNZ.cif | Crystallographic structure of an DPP-IV| |
| ADMETLab_3.0_results | Csv files with an ADMET analysis| |
| ADMET_Lab_3.0_src | Scripts to perfrom an ADMET analysis and analyse data| |
| ADMET_inputs | Peptides for ADMET analysis   | optimised_peptides_for_admet_metadata.csv               |
| ADMET_results | Results of an ADMET analysis  | predictions.csv, optimal_admet_activity_all_length.csv , predictions_gid_annotated.csv, GID_ADMET_combined_active_fragments.xlsx |
| ADMET_src | Scripts for ADMET analysis    | analyse_admet_activity.py, smiles_converter_metadata_v2.py, admet_activity_gds_extended.py|
| BIOPEP_results | Biological activity profile and *in silico* hydrolysis results from BIOPEP_UWM| combined_active_fragments.xlsx, BIOPEP_Activity_Summary_Matrix.xlsx|
| BIOPEP_src | Scripts for BIOPEP automation | script_profiles.py, sort_reports.py, enzymes_action_full.py|
| Docking_inputs | 3D structures of peptides and ligands used for docking|  ACE_clean.pdb, ACE_prepared.pdbqt, DPP4_clean.pdb, DPP4_prepared.pdbqt, config_ACE.txt, config_DPP.txt |
| Docking_src | Preparatory and data analysis scripts for docking | sdf_to_pdbqt.py, run_all_docking.sh
| USCF ChimeraX outputs | Visualisation analysis, peptide location in an active center calculation| 
| USCF ChimeraX scripts | Visualisation analysis, peptide location in an active center calculation scripts| |
| dfs_images_src | Dataframe- and image-making scripts | protein_protease_heatmap.py, dataframe_maker_for_active_fragments.py, biological_activity_profile_report_and_figure.py, dataframe_maker_for_ADMET_GDS.py|
| dock_inp_2run | Config files for second run of docking (exhaustivenes 32)| |
| dock_res_2run | Results of second run of docking| |
| docking_results | first run of docking (exhaustivenes 8)| |
| final_results | Spreadsheets with all data| |
| final_results_gds |Spreadsheets with GDS data| |
| final_results_src |Scripts for combining all data together| |
| images_results | Figures | peptide_table.png, peptide_table_1.png, protein_per_proteas_heatmap.pdf|
| inputs_fasta | Aminoacid sequences of lupin conglutins | All files here were used. Original fasta file with lupin conglutin sequences uniprotkb_accession_A0A6A4QEF0_OR_access_2026_05_29.fasta| Original fasta file with lupin conglutin sequences |
