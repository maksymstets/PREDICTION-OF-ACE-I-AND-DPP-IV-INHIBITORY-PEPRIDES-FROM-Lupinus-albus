import matplotlib
matplotlib.use("Agg")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_DIR  = Path("/home/maksym/lupinus/BIOPEP_results")
OUTPUT_DIR = Path("/home/maksym/lupinus/images_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COL_PROTEIN  = "Protein name"
COL_ENZYME   = "Enzyme ID"
COL_ACTIVITY = "Activity"
COL_PEPTIDE  = "Peptide ID"

# ── Load DataFrame ────────────────────────────────────────────────────────────
combined_path = INPUT_DIR / "combined_active_fragments.xlsx"
df = pd.read_excel(combined_path, dtype=str)

# ── Protein Label Formatting ──────────────────────────────────────────────────
def format_protein_label(raw_name: str) -> str:
    """
    Converts a raw protein identifier like 'cong_beta_1_1' or 'gama_cong_1'
    into a compact display label: 'beta'->'β', 'gamma'/'gama'->'γ', 'delta'->'δ',
    with trailing numbers joined by underscore right after the symbol.
    """
    if pd.isna(raw_name):
        return raw_name
        
    GREEK_MAP = {"beta": "β", "gamma": "γ", "gama": "γ", "delta": "δ"}

    parts = str(raw_name).split("_")
    greek_symbol = None
    numbers = []

    for p in parts:
        low = p.lower()
        if low in GREEK_MAP:
            greek_symbol = GREEK_MAP[low]
        elif p.isdigit():
            numbers.append(p)

    if greek_symbol is None:
        return str(raw_name)

    return greek_symbol + "_".join(numbers) if numbers else greek_symbol

# ── Clean enzyme labels ───────────────────────────────────────────────────────
def short_enzyme_label(enzyme_id: str) -> str:
    if pd.isna(enzyme_id):
        return enzyme_id
    s = str(enzyme_id)
    if s == "GID":
        return "GID"
    
    parts = s.split("|")
    label = parts[1].strip() if len(parts) >= 2 else s.strip()
    
    if label:
        label = label[0].upper() + label[1:]
    
    label = label.replace("Ph ", "pH ").replace("(Ph", "(pH")
    return label

df["_protein_label"] = df[COL_PROTEIN].apply(format_protein_label)
df["_enzyme_label"]  = df[COL_ENZYME].apply(short_enzyme_label)

activities = df[COL_ACTIVITY].dropna().unique().tolist()

# ── Colormap ──────────────────────────────────────────────────────────────────
cmap = LinearSegmentedColormap.from_list(
    "white_red", ["#ffffff", "#fde8e8", "#f5a0a0", "#e03030", "#8b0000"]
)

# ── Build pivot ───────────────────────────────────────────────────────────────
def build_pivot(activity_df: pd.DataFrame) -> pd.DataFrame:
    pivot = (
        activity_df
        .groupby(["_protein_label", "_enzyme_label"])[COL_PEPTIDE]
        .nunique()
        .unstack(fill_value=0)
    )
    pivot.loc["Total General"] = pivot.sum()
    return pivot

# ── Plot one heatmap ──────────────────────────────────────────────────────────
def plot_heatmap(pivot: pd.DataFrame, activity_name: str, ax: plt.Axes,
                 cell_w_in: float, cell_h_in: float, dpi: int):
    proteins = pivot.index.tolist()
    enzymes  = pivot.columns.tolist()
    data     = pivot.values.astype(float)[::-1]   # flip: first protein at top

    n_rows, n_cols = data.shape

    data_only = data[1:]
    vmax = data_only.max() if data_only.size > 0 else 1
    norm = mcolors.Normalize(vmin=0, vmax=vmax)

    fontsize = max(4.0, cell_h_in * 72 * 0.60)

    for r in range(n_rows):
        is_total = (r == 0)
        for c in range(n_cols):
            val = data[r, c]
            color = "#e8eef4" if is_total else (cmap(norm(val)) if val > 0 else "#ffffff")
            ax.add_patch(plt.Rectangle(
                [c, r], 1, 1,
                facecolor=color, edgecolor="#cccccc", linewidth=0.4
            ))
            txt = str(int(val)) if val > 0 else ""
            ax.text(c + 0.5, r + 0.5, txt,
                    ha="center", va="center",
                    fontsize=fontsize,
                    fontweight="bold" if is_total else "normal",
                    color="#1a1a1a")

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)

    # x-axis on top
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(enzymes, rotation=40, ha="left", fontsize=fontsize)

    # y-axis: proteins top→bottom
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(proteins[::-1], fontsize=fontsize)
    for label in ax.get_yticklabels():
        if "Total" in label.get_text():
            label.set_fontweight("bold")

    ax.tick_params(length=0)
    ax.axhline(1, color="#888888", linewidth=0.8)
    
    title_pad = max(10.0, fontsize * 1.5)
    ax.set_title(f"Peptides — {activity_name}", fontsize=fontsize + 1,
                 fontweight="bold", pad=title_pad)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, aspect=25)
    cbar.set_label("Number of peptides", fontsize=fontsize - 1)
    cbar.ax.tick_params(labelsize=fontsize - 1.5)

# ── Cell and layout sizing ────────────────────────────────────────────────────
CELL_W   = 0.20   
CELL_H   = 0.14   
LABEL_W  = 1.2    
CBAR_W   = 0.55   
HEADER_H = 0.8    
FOOTER_H = 0.20   
DPI      = 150

out_path = OUTPUT_DIR / "protein_per_proteas_heatmap_short.pdf"
with PdfPages(out_path) as pdf:
    for idx, activity in enumerate(activities):
        subset = df[df[COL_ACTIVITY] == activity]
        pivot  = build_pivot(subset)
        n_r, n_c = pivot.shape

        fig_w = LABEL_W + n_c * CELL_W + CBAR_W
        fig_h = HEADER_H + n_r * CELL_H + FOOTER_H

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        plot_heatmap(pivot, activity, ax, CELL_W, CELL_H, DPI)

        left   = LABEL_W / fig_w
        right  = 1.0 - CBAR_W / fig_w
        bottom = FOOTER_H / fig_h
        top    = 1.0 - HEADER_H / fig_h
        fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)

        pdf.savefig(fig, dpi=DPI, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  [{idx+1}/{len(activities)}] {activity}")

print(f"\nSaved → {out_path}")
