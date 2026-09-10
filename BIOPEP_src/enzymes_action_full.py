import glob
import os
import time
from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select, WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager

# ==================== CONFIGURATION ====================
URL = "https://biochemia.uwm.edu.pl/biopep/enzymes_act_for_your_seq.php?ABYcalc=&but_enzymes_seq.x=118&but_enzymes_seq.y=10"

# 1. Path to input FASTA directory
FASTA_FOLDER = r"/home/maksym/lupinus/inputs_fasta"

# 2. Sequential Enzymes List
ENZYMES_SEQUENTIAL = [
    "39 | pepsin (pH > 2) | 3.4.23.1",
    "12 | trypsin | EC 3.4.21.4",
    "11 | chymotrypsin (A) | EC 3.4.21.1",
    "24 | papain | EC 3.4.22.2",
    "25 | ficin | EC 3.4.22.3",
    "30 | Stem bromelain | EC 3.4.22.32",
    "41 | subtilisin | EC 3.4.21.62",
    "13 | pepsin (pH 1.3) | EC 3.4.23.1",
]

# 3. GID Profile 1: Standard Pepsin (pH > 2)
GID_PH_GT_2 = [
    "39 | pepsin (pH > 2) | 3.4.23.1",
    "12 | trypsin | EC 3.4.21.4",
    "11 | chymotrypsin (A) | EC 3.4.21.1",
]

# 4. GID Profile 2: High-Acid Pepsin (pH 1.3)
GID_PH_1_3 = [
    "13 | pepsin (pH 1.3) | EC 3.4.23.1",
    "12 | trypsin | EC 3.4.21.4",
    "11 | chymotrypsin (A) | EC 3.4.21.1",
]

# Target output directory
DOWNLOAD_DIR = os.path.join(
    os.path.expanduser("~"),
    "lupinus",
    "BIOPEP_results",
    "BIOPEP_Batch_Outputs_Combined",
)
# =======================================================

if not os.path.exists(DOWNLOAD_DIR):
    os.makedirs(DOWNLOAD_DIR)

# Selenium Chrome Preferences Setup
chrome_options = webdriver.ChromeOptions()
chrome_options.add_experimental_option(
    "prefs",
    {
        "download.default_directory": DOWNLOAD_DIR,
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True,
    },
)

driver = webdriver.Chrome(
    service=Service(ChromeDriverManager().install()), options=chrome_options
)
wait = WebDriverWait(driver, 20)


def clean_filename(name):
  """Removes filesystem-unsafe characters from target filenames."""
  for char in ["/", "\\", ":", "*", "?", '"', "<", ">", "|"]:
    name = name.replace(char, "_")
  return name.replace(" ", "_")


def download_and_rename(click_action, download_dir, target_basename):
  """Triggers file download and renames it cleanly to avoid collisions."""
  initial_files = set(glob.glob(os.path.join(download_dir, "*")))
  click_action()

  start_time = time.time()
  timeout = 45

  while time.time() - start_time < timeout:
    current_files = set(glob.glob(os.path.join(download_dir, "*")))
    new_files = current_files - initial_files

    completed_files = [
        f
        for f in new_files
        if not f.endswith(".crdownload") and not f.endswith(".tmp")
    ]

    if completed_files:
      time.sleep(1)
      downloaded_file = completed_files[0]
      ext = os.path.splitext(downloaded_file)[1]

      new_filename = f"{clean_filename(target_basename)}{ext}"
      new_path = os.path.join(download_dir, new_filename)

      if os.path.exists(new_path):
        base, extension = os.path.splitext(new_path)
        new_path = f"{base}_{int(time.time())}{extension}"

      try:
        os.rename(downloaded_file, new_path)
        print(f"    -> Saved: {os.path.basename(new_path)}")
        return
      except Exception:
        time.sleep(2)
        os.rename(downloaded_file, new_path)
        print(f"    -> Saved: {os.path.basename(new_path)}")
        return
    time.sleep(0.5)
  print(f"    -> Error: Download timed out for {target_basename}")


def submit_initial_form(fasta_sequence, enzymes):
  """Fills sequence data into textarea and maps enzyme selections into dropdown boxes (enz1, enz2, enz3...)."""
  driver.get(URL)

  fasta_input = wait.until(
      EC.presence_of_element_located((By.TAG_NAME, "textarea"))
  )
  fasta_input.clear()
  fasta_input.send_keys(fasta_sequence)

  for index, enz in enumerate(enzymes):
    box_name = f"enz{index + 1}"
    try:
      dropdown_element = driver.find_element(By.NAME, box_name)
      dropdown = Select(dropdown_element)

      target_enz = enz.strip().lower()
      enz_id = target_enz.split("|")[0].strip() if "|" in target_enz else None

      matched = False
      if enz_id:
        for option in dropdown.options:
          if option.get_attribute("value") == enz_id:
            dropdown.select_by_value(enz_id)
            matched = True
            break

      if not matched:
        for option in dropdown.options:
          option_text = option.text.strip()
          if target_enz in option_text.lower():
            dropdown.select_by_visible_text(option_text)
            matched = True
            break

      if not matched:
        print(
            f"  [Warning] Enzyme '{enz}' could not be matched in dropdown"
            f" '{box_name}'"
        )

    except NoSuchElementException:
      break

  buttons = driver.find_elements(
      By.XPATH, "//input[@type='image' or @type='submit']"
  )
  submit_btn = None

  for btn in buttons:
    name_attr = (btn.get_attribute("name") or "").lower()
    alt_attr = (btn.get_attribute("alt") or "").lower()
    src_attr = (btn.get_attribute("src") or "").lower()

    if (
        "enzymes" in name_attr
        or "report" in alt_attr
        or "report" in src_attr
        or "seq" in name_attr
    ):
      submit_btn = btn
      break

  if not submit_btn and buttons:
    submit_btn = buttons[-1]

  if not submit_btn:
    raise NoSuchElementException("Could not locate submission button.")

  submit_btn.click()


def process_report(fasta_sequence, enzymes, action_type, target_name):
  """Handles step navigation to pull Active Fragments and Calculations excel exports."""
  submit_initial_form(fasta_sequence, enzymes)

  time.sleep(1)
  initial_page_text = driver.find_element(By.TAG_NAME, "body").text.lower()
  if any(
      msg in initial_page_text
      for msg in ["no data", "no fragments", "brak danych", "nie znaleziono"]
  ):
    print(
        "     [Notice] No data available for this selection. Skipping"
        " download."
    )
    return

  print(f"  -> Processing intermediate page for: {action_type}...")
  wait.until(
      EC.presence_of_element_located(
          (By.XPATH, "//input[@type='image' or @type='submit']")
      )
  )

  buttons = driver.find_elements(
      By.XPATH, "//input[@type='image' or @type='submit']"
  )
  next_btn = None

  if action_type == "active_fragments":
    keywords = ["active", "akt", "szukaj", "search", "fragment", "sek"]
    preferred_index = 0
  else:
    keywords = ["calc", "oblicz", "param", "quant", "calcul"]
    preferred_index = 1

  for btn in buttons:
    combined_metadata = (
        f"{btn.get_attribute('name')} {btn.get_attribute('alt')}"
        f" {btn.get_attribute('src')} {btn.get_attribute('value')}".lower()
    )
    if any(kw in combined_metadata for kw in keywords):
      next_btn = btn
      break

  if not next_btn and len(buttons) > preferred_index:
    next_btn = buttons[preferred_index]
  if not next_btn and buttons:
    next_btn = buttons[0]

  next_btn.click()

  print("  -> Accessing final summary page to pull Excel export...")
  wait.until(
      EC.presence_of_element_located(
          (By.XPATH, "//input[@type='image' or @type='submit']")
      )
  )
  time.sleep(1)

  final_buttons = driver.find_elements(
      By.XPATH, "//input[@type='image' or @type='submit']"
  )
  excel_btn = None

  for btn in final_buttons:
    combined_metadata = (
        f"{btn.get_attribute('name')} {btn.get_attribute('alt')}"
        f" {btn.get_attribute('src')} {btn.get_attribute('value')}".lower()
    )
    if any(kw in combined_metadata for kw in ["excel", "xls", "export"]):
      excel_btn = btn
      break

  if not excel_btn:
    final_page_text = driver.find_element(By.TAG_NAME, "body").text.lower()
    if any(msg in final_page_text for msg in ["no data", "no fragments"]):
      print("     [Notice] No data to export. Skipping.")
      return
    else:
      raise NoSuchElementException("Excel export button missing.")

  download_and_rename(lambda: excel_btn.click(), DOWNLOAD_DIR, target_name)


# ==================== MAIN BATCH EXECUTION LOOP ====================
try:
  fasta_extensions = ["*.fasta", "*.fa", "*.txt"]
  fasta_files = []
  for ext in fasta_extensions:
    fasta_files.extend(glob.glob(os.path.join(FASTA_FOLDER, ext)))

  if not fasta_files:
    print(f"No valid FASTA files found in folder: {FASTA_FOLDER}")
    driver.quit()
    exit()

  print(f"Found {len(fasta_files)} FASTA files to process.")

  for file_index, fasta_path in enumerate(fasta_files, 1):
    fasta_filename = os.path.splitext(os.path.basename(fasta_path))[0]

    with open(fasta_path, "r") as f:
      fasta_sequence = f.read().strip()

    print(f"\n==================================================")
    print(
        f"[{file_index}/{len(fasta_files)}] PROCESSING FILE: {fasta_filename}"
    )
    print(f"==================================================")

    # --- TASK 1: Sequential Enzymes (One by One) ---
    print("\n--- Running Sequential Enzymes ---")
    for enz in ENZYMES_SEQUENTIAL:
      clean_enz_name = clean_filename(
          enz.split("|")[1].strip() if "|" in enz else enz
      )
      print(f"\n  Enzyme: {enz}")

      name_active = f"{fasta_filename}_{clean_enz_name}_ActiveFragments"
      process_report(fasta_sequence, [enz], "active_fragments", name_active)

      name_calc = f"{fasta_filename}_{clean_enz_name}_Calculations"
      process_report(fasta_sequence, [enz], "calculate", name_calc)

    # --- TASK 2: Simultaneous GID Profile 1 (Pepsin pH > 2) ---
    print("\n--- Running Simultaneous GID (Pepsin pH > 2) ---")
    name_gid1_active = f"{fasta_filename}_GID_pH_gt_2_ActiveFragments"
    process_report(
        fasta_sequence, GID_PH_GT_2, "active_fragments", name_gid1_active
    )

    name_gid1_calc = f"{fasta_filename}_GID_pH_gt_2_Calculations"
    process_report(fasta_sequence, GID_PH_GT_2, "calculate", name_gid1_calc)

    # --- TASK 3: Simultaneous GID Profile 2 (Pepsin pH 1.3) ---
    print("\n--- Running Simultaneous GID (Pepsin pH 1.3) ---")
    name_gid2_active = f"{fasta_filename}_GID_pH_1_3_ActiveFragments"
    process_report(
        fasta_sequence, GID_PH_1_3, "active_fragments", name_gid2_active
    )

    name_gid2_calc = f"{fasta_filename}_GID_pH_1_3_Calculations"
    process_report(fasta_sequence, GID_PH_1_3, "calculate", name_gid2_calc)

  print(
      f"\n*** All tasks completed successfully! Outputs saved to: {DOWNLOAD_DIR}"
      " ***"
  )

except Exception as e:
  print(f"\nA critical execution error occurred: {e}")
  input("\nPress Enter to exit and close browser drivers...")

finally:
  driver.quit()