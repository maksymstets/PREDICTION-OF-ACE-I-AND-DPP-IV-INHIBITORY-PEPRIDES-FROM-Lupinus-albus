import os
import glob
import pandas as pd

def parse_ab_coef_file(file_path):
    """
    Reads a raw A_B_coef Excel file, skips the first row,
    and replaces the 'NAME:' column with the extracted filename prefix.
    """
    filename = os.path.basename(file_path)
    
    # 1. Extract protein name from filename (everything before '_A_B_coef')
    protein_name = filename.split("_Profile_Report")[0]
    
    # 2. Read the excel file, skipping the first row (index 0)
    df = pd.read_excel(file_path, skiprows=1)
    
    # Clean whitespace from column names to ensure clean targeting
    df.columns = [str(col).strip() for col in df.columns]
    
    # 3. Drop completely empty rows to keep the output clean
    df = df.dropna(how='all')
    
    # 4. Overwrite the 'NAME:' column with the parsed protein name
    # We check for both variations in case the raw file drops the colon
    if 'NAME:' in df.columns:
        df['NAME:'] = protein_name
    elif 'NAME' in df.columns:
        df['NAME'] = protein_name
    else:
        # Fallback: if the column is missing entirely, insert it
        df.insert(1, 'NAME:', protein_name)
        
    return df

def compile_all_data(folder_path):
    """
    Gathers all A_B_coef files in the directory and compiles them 
    into a single consolidated master Pandas DataFrame.
    """
    search_pattern = os.path.join(folder_path, "*_Profile_Report.xlsx")
    all_files = glob.glob(search_pattern)
    
    if not all_files:
        print(f"No files found in: {folder_path}")
        return pd.DataFrame()
        
    print(f"Found {len(all_files)} files. Parsing and consolidating...")
    
    dataframe_list = []
    for f in all_files:
        try:
            cleaned_df = parse_ab_coef_file(f)
            dataframe_list.append(cleaned_df)
        except Exception as e:
            print(f"  -> Error parsing file {os.path.basename(f)}: {e}")
            
    if not dataframe_list:
        return pd.DataFrame()
        
    # Stack all data vertically into a single consolidated master DataFrame
    master_df = pd.concat(dataframe_list, ignore_index=True)
    print("-> Master DataFrame successfully created!")
    return master_df

# ==================== RUN ENGINE ====================
if __name__ == "__main__":
    # Define Paths
    RAW_DATA_DIR = '/home/maksym/lupinus/BIOPEP_results/BIOPEP_A_B_coef_Outputs'
    OUTPUT_DIR = '/home/maksym/lupinus/BIOPEP_results'
    
    # Ensure output directory exists before saving
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 1. Generate consolidated dataframe in-memory
    master_df = compile_all_data(RAW_DATA_DIR)
    
    if master_df.empty:
        print("Process halted: No data loaded from files.")
    else:
        # 2. Save underlying Data Matrix as clean Excel File
        excel_out_path = os.path.join(OUTPUT_DIR, "Consolidated_A_B_coef.xlsx")
        master_df.to_excel(excel_out_path, index=False)
        print(f"Done → Consolidated Spreadsheet saved to: {excel_out_path}")
