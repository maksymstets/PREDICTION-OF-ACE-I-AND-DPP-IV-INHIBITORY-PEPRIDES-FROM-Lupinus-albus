import glob
import os
import time
import pandas as pd
from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select, WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager

# ==================== CONFIGURATION ====================
URL = "https://biochemia.uwm.edu.pl/biopep/enzymes_act_for_your_seq.php?ABYcalc=&but_enzymes_seq.x=141&but_enzymes_seq.y=15"

# 1. Path to your input CSV file containing optimized sequences
INPUT_CSV = r"/home/maksym/lupinus/ADMET_results/predictions.csv"

# 2. Gastrointestinal digestion profiles (Pepsin pH 1.3 vs Pepsin pH > 2)
DIGESTION_PROFILES = {
    "GID_pH_1_3": [
        "13 | pepsin (pH 1.3) | EC 3.4.23.1",
        "12 | trypsin | EC 3.4.21.4",
        "11 | chymotrypsin (A) | EC 3.4.21.1",
    ],
    "GID_pH_gt_2": [
        "39 | pepsin (pH > 2) | 3.4.23.1",
        "12 | trypsin | EC 3.4.21.4",
        "11 | chymotrypsin (A) | EC 3.4.21.1",
    ],
}

# Destination for organized Excel outputs
DOWNLOAD_DIR = os.path.join(
    os.path.expanduser("~"),
    "lupinus",
    "ADMET_results",
    "GID_ADMET_downloads_all",
)
# =======================================================

if not os.path.exists(DOWNLOAD_DIR):
  os.makedirs(DOWNLOAD_DIR)

# Chrome Driver Setup
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
  """Removes filesystem-unsafe characters from names."""
  for char in ["/", "\\", ":", "*", "?", '"', "<", ">", "|"]:
    name = name.replace(char, "_")
  return name.replace(" ", "_")


def download_and_rename(click_action, download_dir, target_basename):
  """Triggers a download, waits for it to complete, and renames it safely."""
  initial_files = set(glob.glob(os.path.join(download_dir, "*")))

  click_action()

  start_time = time.time()
  timeout = 45  # Max seconds to wait per file download

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


def submit_initial_form(peptide_sequence, enzymes):

  driver.get(URL)

  # Fill out the peptide sequence textarea
  fasta_input = wait.until(
      EC.presence_of_element_located((By.TAG_NAME, "textarea"))
  )
  fasta_input.clear()
  fasta_input.send_keys(peptide_sequence)

  # Loop over the necessary enzyme boxes (enz1, enz2, enz3...)
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
            print(
                f"  -> Box '{box_name}': Selected option ID '{enz_id}'"
                f" ({option.text.strip()})"
            )
            matched = True
            break

      if not matched:
        for option in dropdown.options:
          option_text = option.text.strip()
          if target_enz in option_text.lower():
            dropdown.select_by_visible_text(option_text)
            print(
                f"  -> Box '{box_name}': Selected by text '{option_text}'"
            )
            matched = True
            break

      if not matched:
        print(
            f"  [Warning] Enzyme '{enz}' could not be matched in box"
            f" '{box_name}'"
        )

    except NoSuchElementException:
      print(
          f"  [Notice] Dropdown box '{box_name}' does not exist on this page"
          " configuration. Skipping."
      )
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
    raise NoSuchElementException("Could not locate the submission button.")

  submit_btn.click()


def process_report(peptide_sequence, enzymes, action_type, target_name):
  """Navigates through intermediate pages to fetch the Excel sheet."""
  submit_initial_form(peptide_sequence, enzymes)

  time.sleep(1)
  initial_page_text = driver.find_element(By.TAG_NAME, "body").text.lower()
  if any(msg in initial_page_text for msg in ["no data", "no fragments"]):
    print(
        "     [Notice] No active fragments available for this sequence."
        " Skipping download."
    )
    return

  print(f"  -> Processing intermediate results for: {action_type}...")
  wait.until(
      EC.presence_of_element_located(
          (By.XPATH, "//input[@type='image' or @type='submit']")
      )
  )

  buttons = driver.find_elements(
      By.XPATH, "//input[@type='image' or @type='submit']"
  )
  next_btn = None

  keywords = ["active", "akt", "szukaj", "search", "fragment", "sek"]
  preferred_index = 0

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
  if not os.path.exists(INPUT_CSV):
    raise FileNotFoundError(f"Missing input database: {INPUT_CSV}")

  df_peptides = pd.read_csv(INPUT_CSV)
  print(f"Loaded {len(df_peptides)} optimized sequences from CSV for GID.")

  for file_index, row in df_peptides.iterrows():
    peptide_sequence = str(row["Sequence"]).strip()
    activity = row["Activity"]

    print("\n==================================================")
    print(
        f"[{file_index + 1}/{len(df_peptides)}] HYDROLYZING PEPTIDE:"
        f" {peptide_sequence}"
    )
    print(f"Target Activity: {activity}")
    print("==================================================")

    # Iterate through both digestion profiles sequentially
    for profile_name, enzymes in DIGESTION_PROFILES.items():
      print(f"\n--- Running Profile: {profile_name} ---")
      target_filename = f"{peptide_sequence}_{profile_name}_ActiveFragments"

      process_report(
          peptide_sequence, enzymes, "active_fragments", target_filename
      )
      time.sleep(1)

  print(
      "\n*** Extraction finished. Master data blocks written to:"
      f" {DOWNLOAD_DIR} ***"
  )

except Exception as e:
  print(f"\nAutomation error encountered: {e}")
  input("\nPress Enter to exit and safely wind down Selenium drivers...")

finally:
  driver.quit()
