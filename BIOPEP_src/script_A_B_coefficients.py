import os
import time
import glob
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException
from webdriver_manager.chrome import ChromeDriverManager

# ==================== CONFIGURATION ====================
# Clean URL for the Finding Profiles tool (removed the pre-submit GET parameters)
URL = "https://biochemia.uwm.edu.pl/biopep/finding_aby.php?ABYcalc=&sel_activity=*&but_list_pep.x=95&but_list_pep.y=14"

# 1. Path to the folder where you saved your FASTA files
FASTA_FOLDER = r"/home/maksym/lupinus/inputs_fasta" 

# Destination for organized outputs
DOWNLOAD_DIR = os.path.join(os.path.expanduser("~"), "lupinus", "BIOPEP_results", "BIOPEP_A_B_coef_Outputs")
# =======================================================

if not os.path.exists(DOWNLOAD_DIR):
    os.makedirs(DOWNLOAD_DIR)

# Chrome Driver Setup
chrome_options = webdriver.ChromeOptions()
chrome_options.add_experimental_option("prefs", {
    "download.default_directory": DOWNLOAD_DIR,
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing.enabled": True
})

driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=chrome_options)
wait = WebDriverWait(driver, 20)

def clean_filename(name):
    """Removes filesystem-unsafe characters from names."""
    for char in ['/', '\\', ':', '*', '?', '"', '<', '>', '|']:
        name = name.replace(char, '_')
    return name.replace(' ', '_')

def download_and_rename(click_action, download_dir, target_basename):
    """Triggers a download, waits for it to complete, and renames it safely."""
    initial_files = set(glob.glob(os.path.join(download_dir, "*")))
    
    # Click the Excel Export button
    click_action()
    
    start_time = time.time()
    timeout = 45  # Max seconds to wait per file download
    
    while time.time() - start_time < timeout:
        current_files = set(glob.glob(os.path.join(download_dir, "*")))
        new_files = current_files - initial_files
        
        # Filter out ongoing Chrome download temp files (.crdownload / .tmp)
        completed_files = [f for f in new_files if not f.endswith('.crdownload') and not f.endswith('.tmp')]
        
        if completed_files:
            time.sleep(1)  # Brief pause to let OS completely clear file locks
            downloaded_file = completed_files[0]
            ext = os.path.splitext(downloaded_file)[1]
            
            new_filename = f"{clean_filename(target_basename)}{ext}"
            new_path = os.path.join(download_dir, new_filename)
            
            # De-conflict in case file already exists
            if os.path.exists(new_path):
                base, extension = os.path.splitext(new_path)
                new_path = f"{base}_{int(time.time())}{extension}"
                
            try:
                os.rename(downloaded_file, new_path)
                print(f"    -> Saved: {os.path.basename(new_path)}")
                return
            except Exception:
                time.sleep(2)  # Secondary fallback wait if system is running slow
                os.rename(downloaded_file, new_path)
                print(f"    -> Saved: {os.path.basename(new_path)}")
                return
        time.sleep(0.5)
    print(f"    -> Error: Download timed out for {target_basename}")

def process_profile(fasta_sequence, target_name):
    """Loads the profile tool, pastes FASTA data, and exports the Excel report."""
    driver.get(URL)
    
    # Step 1: Input FASTA sequence
    print("  -> Pasting FASTA sequence...")
    fasta_input = wait.until(EC.presence_of_element_located((By.TAG_NAME, "textarea")))
    fasta_input.clear()
    fasta_input.send_keys(fasta_sequence)
    
    # Step 2: Click the 'Report' button
    print("  -> Locating the 'Report' button...")
    buttons = driver.find_elements(By.XPATH, "//input[@type='image' or @type='submit']")
    submit_btn = None

    for btn in buttons:
        name_attr = (btn.get_attribute("name") or "").lower()
        alt_attr = (btn.get_attribute("alt") or "").lower()
        src_attr = (btn.get_attribute("src") or "").lower()
        
        # Looking specifically for the report submission button
        if "prof" in name_attr or "find" in name_attr or "szukaj" in name_attr or "report" in alt_attr:
            submit_btn = btn
            break

    # Fallback if specific keywords aren't found
    if not submit_btn and buttons:
        submit_btn = buttons[-1]

    if not submit_btn:
        raise NoSuchElementException("Could not locate the 'Report' form submission button.")
        
    submit_btn.click()

    # Step 3: Click 'Export to Excel' on the next page
    print("  -> Waiting for the final report page to load...")
    wait.until(EC.presence_of_element_located((By.XPATH, "//input[@type='image' or @type='submit']")))
    time.sleep(1)
    
    final_buttons = driver.find_elements(By.XPATH, "//input[@type='image' or @type='submit']")
    excel_btn = None
    
    for btn in final_buttons:
        name_attr = (btn.get_attribute("name") or "").lower()
        alt_attr = (btn.get_attribute("alt") or "").lower()
        src_attr = (btn.get_attribute("src") or "").lower()
        value_attr = (btn.get_attribute("value") or "").lower()
        
        combined_metadata = f"{name_attr} {alt_attr} {src_attr} {value_attr}"
        if any(kw in combined_metadata for kw in ["excel", "xls", "export"]):
            excel_btn = btn
            print(f"  -> Found export button: '{btn.get_attribute('name')}'")
            break
            
    if not excel_btn and final_buttons:
        excel_btn = final_buttons[0]
        
    if not excel_btn:
        raise NoSuchElementException("Could not locate the Excel export button on the final summary page.")
        
    print("  -> Exporting to Excel...")
    download_and_rename(lambda: excel_btn.click(), DOWNLOAD_DIR, target_name)


# ==================== MAIN BATCH EXECUTION LOOP ====================
try:
    # Gather all FASTA, FA, and TXT files from targeted folder
    fasta_extensions = ["*.fasta", "*.fa", "*.txt"]
    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(glob.glob(os.path.join(FASTA_FOLDER, ext)))
        
    if not fasta_files:
        print(f"No valid FASTA files found in folder: {FASTA_FOLDER}")
        driver.quit()
        exit()

    print(f"Found {len(fasta_files)} files to process.")

    for file_index, fasta_path in enumerate(fasta_files, 1):
        fasta_filename = os.path.splitext(os.path.basename(fasta_path))[0]
        
        with open(fasta_path, "r") as f:
            fasta_sequence = f.read().strip()
            
        print(f"\n==================================================")
        print(f"[{file_index}/{len(fasta_files)}] PROCESSING FILE: {fasta_filename}")
        print(f"==================================================")
        
        # Format output filename
        report_name = f"{fasta_filename}_Profile_Report"
        
        # Process the profile tool
        process_profile(fasta_sequence, report_name)

    print(f"\n*** All tasks completed successfully! Outputs saved to: {DOWNLOAD_DIR} ***")

except Exception as e:
    print(f"\nA critical execution error occurred: {e}")
    input("\nPress Enter to exit and close the automated browser windows...")

finally:
    driver.quit()
