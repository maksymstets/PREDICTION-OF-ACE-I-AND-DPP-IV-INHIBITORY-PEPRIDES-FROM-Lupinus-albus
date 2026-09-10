import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def parse_raw_profile(file_path):
    """
    Reads a raw BIOPEP Excel file, locates the 'Protein sequence:' landmark,
    skips the metadata, and builds a clean DataFrame with the 'Protein name' column.
    """
    filename = os.path.basename(file_path)
    
    # 1. Extract protein name from filename (everything before '_Profile_Report')
    protein_name = filename.split("_Profile_Report")[0]
    
    # 2. Read the raw excel file completely to locate the metadata anchor
    df_raw = pd.read_excel(file_path, header=None)
    
    anchor_row_idx = None
    for idx, row in df_raw.iterrows():
        # Combine all cell values in the row into a single string to inspect it
        row_text = " ".join([str(val) for val in row.values]).lower()
        
        if "protein sequence:" in row_text:
            anchor_row_idx = idx
            break
            
    if anchor_row_idx is None:
        # Fallback if the landmark isn't found for some reason
        print(f"  -> Warning: 'Protein sequence:' label not found in {filename}. Using default skip.")
        header_skip = 4
    else:
        # The sequence text is at (anchor_row_idx + 1)
        # The real table headers are exactly at (anchor_row_idx + 2)
        header_skip = anchor_row_idx + 1

    # 3. Reload the file cleanly, skipping all metadata rows up to the real table header
    df_clean = pd.read_excel(file_path, skiprows=header_skip)
    
    # Clean whitespace from column strings
    df_clean.columns = [str(col).strip() for col in df_clean.columns]
    
    # 4. Filter out any empty rows or stray text at the bottom of the export
    if len(df_clean) > 0:
        first_col = df_clean.columns[0]
        df_clean = df_clean.dropna(subset=[first_col])
        # Drop rows where the first column repeats header text or metadata accidentally
        df_clean = df_clean[~df_clean[first_col].astype(str).str.contains("protein|id of|name:", case=False, na=False)]

    # 5. Add the 'Protein name' column as the very first column on the left
    df_clean.insert(0, 'Protein name', protein_name)
    
    return df_clean

def load_all_profiles_to_df(folder_path):
    """
    Gathers all profile reports in the directory and compiles them 
    into a single consolidated master Pandas DataFrame.
    """
    search_pattern = os.path.join(folder_path, "*_Profile_Report.xlsx")
    all_files = glob.glob(search_pattern)
    
    if not all_files:
        print(f"No Profile Reports found in: {folder_path}")
        return pd.DataFrame()
        
    print(f"Found {len(all_files)} raw files. Parsing into memory using metadata anchor...")
    
    dataframe_list = []
    for f in all_files:
        try:
            cleaned_df = parse_raw_profile(f)
            dataframe_list.append(cleaned_df)
        except Exception as e:
            print(f"  -> Error parsing file {os.path.basename(f)}: {e}")
            
    if not dataframe_list:
        return pd.DataFrame()
        
    # Stack all data vertically into a single consolidated master DataFrame
    master_df = pd.concat(dataframe_list, ignore_index=True)
    print("-> Master DataFrame successfully created!")
    return master_df

def format_protein_label(raw_name):
    """
    Converts a raw protein/file identifier like 'cong_beta_1_1' or 'gama_cong_1'
    into a compact display label: 'beta'->'β', 'gamma'/'gama'->'γ', 'delta'->'δ',
    with any trailing numeric parts joined by underscore right after the symbol.
    e.g. cong_beta_1_1 -> 'β1_1', cong_beta_2 -> 'β2', cong_delta -> 'δ',
         gama_cong_1 -> 'γ1', legumin -> 'legumin' (no greek token -> unchanged).
    """
    GREEK_MAP = {"beta": "β", "gamma": "γ", "gama": "γ", "delta": "δ"}

    parts = raw_name.split("_")
    greek_symbol = None
    numbers = []

    for p in parts:
        low = p.lower()
        if low in GREEK_MAP:
            greek_symbol = GREEK_MAP[low]
        elif p.isdigit():
            numbers.append(p)

    if greek_symbol is None:
        return raw_name

    return greek_symbol + "_".join(numbers) if numbers else greek_symbol

def cell(ax, x, y, w, h, text,
         bg="#ffffff", fc="black", fs=7.0,
         bold=False, ha="center", va="center", border=True):
    """
    Helper function to paint an individual table cell with boundaries and text.
    """
    if bg != "none":
        rect = plt.Rectangle((x, y), w, h,
                             linewidth=0.5 if border else 0,
                             edgecolor="#999999" if border else "none",
                             facecolor=bg,
                             transform=ax.transAxes, clip_on=False)
        ax.add_patch(rect)
    
    # Add subtle left horizontal padding for left-aligned text to prevent border crowding
    tx = (x + w / 2) if ha == "center" else (x + 0.015)
    ax.text(tx, y + h / 2, text,
            ha=ha, va="center", fontsize=fs, color=fc,
            fontweight="bold" if bold else "normal",
            transform=ax.transAxes, clip_on=False,
            multialignment="center")

# ==================== RUN ENGINE ====================
if __name__ == "__main__":
    # Path to your raw folder containing the profile reports
    RAW_DATA_DIR = '/home/maksym/lupinus/BIOPEP_results/BIOPEP_Profiles_Outputs'
    
    # 1. Generate clean dataframe in-memory
    df = load_all_profiles_to_df(RAW_DATA_DIR)
    
    if df.empty:
        print("Process halted: No data loaded from files.")
    else:
        # ── 2. Build pivot: rows = Activity, cols = Protein name, values = count ─────
        pivot = (
            df.groupby(["Activity:", "Protein name"])
              .size()
              .unstack(fill_value=0)
              .reset_index()
        )
        pivot.columns.name = None
        pivot = pivot.rename(columns={"Activity:": "Activity"})

        protein_cols = [c for c in pivot.columns if c != "Activity"]
        pivot["Total General"] = pivot[protein_cols].sum(axis=1)

        # Display labels for headers (raw column keys stay unchanged for indexing)
        display_labels = {prot: format_protein_label(prot) for prot in protein_cols}

        # Replace 0 with "" for clean visual table display
        display_df = pivot.copy()
        for col in protein_cols + ["Total General"]:
            display_df[col] = display_df[col].apply(lambda x: "" if x == 0 else str(int(x)))

        # ── Footer totals: sum of peptides per protein (and grand total) ──────────────
        footer = {prot: str(int(pivot[prot].sum())) for prot in protein_cols}
        footer["Total General"] = str(int(pivot["Total General"].sum()))

        # ── 4. Dynamic Fitting using Matplotlib Renderer ──────────────────────────────
        n_rows      = len(display_df)
        n_proteins  = len(protein_cols)

        # 1. Define explicit font sizes for Activity cells (Decoupled from column width)
        fs_act_hdr = 10.0
        fs_act     = 9.5

        # 2. Measure actual rendered width via off-screen canvas attached to temp_fig
        temp_fig = plt.figure()
        temp_fig.canvas.draw()
        renderer = temp_fig.canvas.get_renderer()

        max_act_width_in = 0.0
        activity_texts = ["Activity", "Total general"] + display_df["Activity"].astype(str).tolist()

        for text_val in activity_texts:
            is_bold = text_val in ["Activity", "Total general"]
            fs = fs_act_hdr if is_bold else fs_act
            
            t = temp_fig.text(0, 0, text_val, fontsize=fs, fontweight="bold" if is_bold else "normal")
            bbox_width_in = t.get_window_extent(renderer).width / temp_fig.dpi
            if bbox_width_in > max_act_width_in:
                max_act_width_in = bbox_width_in

        plt.close(temp_fig)

        # 3. Add explicit horizontal padding
        ACT_PADDING_IN = 0.18
        act_inch = max_act_width_in + ACT_PADDING_IN

        # Header height calculation based on maximum line count
        max_lines = max([lbl.count("\n") + 1 for lbl in display_labels.values()]) if display_labels else 1
        hdr_inch  = max(0.5, 0.35 + (max_lines * 0.16))

        col_inch  = 0.65         
        tot_inch  = 0.75         
        row_inch  = 0.25         

        fig_w = act_inch + n_proteins * col_inch + tot_inch
        fig_h = hdr_inch + (n_rows + 1) * row_inch   # +1 for footer row

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.axis("off")

        # Column widths as fractions of figure width
        act_w  = act_inch / fig_w
        tot_w  = tot_inch / fig_w
        data_w = col_inch / fig_w

        # Row heights as fractions of figure height
        header_h1 = 0.35 / fig_h                        # Fixed space for top title banner
        header_h2 = (hdr_inch - 0.35) / fig_h           # Dynamic space for protein headers
        row_h     = row_inch / fig_h

        # Other font sizes scaled independently
        fs_banner  = min(12.0, fig_w * 0.8)
        fs_header  = min(9.0, col_inch * 72 * 0.22)   
        fs_data    = min(10.0, min(col_inch, row_inch) * 72 * 0.55)
        fs_total   = min(10.0, tot_inch * 72 * 0.22)

        y_top = 1.0 - header_h1
        y_h2  = y_top - header_h2

        HEADER_BG = "#dce6f1"
        STRIPE1   = "#f2f7fb"
        STRIPE2   = "#ffffff"

        # ── 5. Header row 1 – banner + "Activity" + "Total General" ──────────────────
        cell(ax, 0, y_h2, act_w, header_h1 + header_h2,
             "Activity", bg=HEADER_BG, bold=True, fs=fs_act_hdr)

        banner_w = data_w * n_proteins
        cell(ax, act_w, y_top, banner_w, header_h1,
             "Lupin Conglutin Proteins", bg=HEADER_BG, bold=True, fs=fs_banner)

        cell(ax, act_w + banner_w, y_h2, tot_w, header_h1 + header_h2,
             "Total\nGeneral", bg=HEADER_BG, bold=True, fs=fs_total)

        # ── 6. Header row 2 – protein names ──────────────────────────────────────────
        x_cur = act_w
        for prot in protein_cols:
            label = display_labels[prot]
            cell(ax, x_cur, y_h2, data_w, header_h2,
                 label, bg=HEADER_BG, bold=True, fs=fs_header)
            x_cur += data_w

        # ── 7. Data rows ──────────────────────────────────────────────────────────────
        for i, row in display_df.iterrows():
            y   = y_h2 - (i + 1) * row_h
            bg  = STRIPE1 if i % 2 == 0 else STRIPE2

            cell(ax, 0, y, act_w, row_h,
                 row["Activity"], bg=bg, ha="left", fs=fs_act)

            x_cur = act_w
            for prot in protein_cols:
                cell(ax, x_cur, y, data_w, row_h,
                     row[prot], bg=bg, fs=fs_data)
                x_cur += data_w

            cell(ax, act_w + banner_w, y, tot_w, row_h,
                 row["Total General"], bg=bg, bold=True, fs=fs_data)

        # ── 8. Footer row – "Total general" per protein ──────────────────────────────
        FOOTER_BG = "#dce6f1"
        y_footer = y_h2 - n_rows * row_h - row_h

        cell(ax, 0, y_footer, act_w, row_h,
             "Total general", bg=FOOTER_BG, bold=True, fs=fs_act, ha="left")

        x_cur = act_w
        for prot in protein_cols:
            cell(ax, x_cur, y_footer, data_w, row_h,
                 footer[prot], bg=FOOTER_BG, bold=True, fs=fs_data)
            x_cur += data_w

        cell(ax, act_w + banner_w, y_footer, tot_w, row_h,
             footer["Total General"], bg=FOOTER_BG, bold=True, fs=fs_data)

        # ── 9. Outer border box ───────────────────────────────────────────────────────
        table_h = header_h1 + header_h2 + n_rows * row_h + row_h
        outer = plt.Rectangle(
            (0, y_footer), 1, table_h,
            linewidth=1.2, edgecolor="#3a3a3a",
            facecolor="none", transform=ax.transAxes, clip_on=False)
        ax.add_patch(outer)

        # Save table image
        plt.tight_layout(pad=0)
        img_out_path = "/home/maksym/lupinus/images_results/peptide_table.png"
        plt.savefig(img_out_path, dpi=180, bbox_inches="tight", facecolor="white")
        plt.close()
        print(f"Done → Visual Table Image saved to: {img_out_path}")

        # ── 10. Save underlying Data Matrix as clean Excel File ─────────────────────
        excel_df = pivot.copy()
        
        footer_export = footer.copy()
        footer_export["Activity"] = "Total general"
        
        excel_df = pd.concat([excel_df, pd.DataFrame([footer_export])], ignore_index=True)
        
        excel_out_path = "/home/maksym/lupinus/BIOPEP_results/BIOPEP_Activity_Summary_Matrix.xlsx"
        excel_df.to_excel(excel_out_path, index=False)
        print(f"Done → Underlying Data Spreadsheet saved to: {excel_out_path}")