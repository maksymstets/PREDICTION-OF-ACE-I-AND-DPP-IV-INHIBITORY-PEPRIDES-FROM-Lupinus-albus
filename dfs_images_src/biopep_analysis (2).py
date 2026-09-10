"""
BIOPEP-UWM coefficient analysis for Lupinus albus conglutins
=============================================================

Turns the long-format "enzyme action" export (Protein Name, Enzyme id, NO,
Activity, Ae, W, Be, V, DHt[%], Source file) into answers for a fixed set of
biological questions, rather than a table nobody reads.

READ THIS FIRST — directionality / caveats baked into the analysis below:
  - Ae, W  : higher = more/more-efficient release. No EC50 needed.
  - Be     : LOWER = more potent peptides released (same convention as B).
  - V      : CONFIRMED formula (from BIOPEP-UWM's own documentation):
                 V = BE / B
             where BE = potency-weighted release by this enzyme (Eq. 5) and
             B = whole-protein theoretical potency (Eq. 2). My earlier guess
             in this analysis (V = BE/B) was right; a subsequent "correction"
             based on misreading a screenshot's tiny decimal values was
             wrong and has been reverted. Interpretation, given the
             established convention that LOWER B/BE = more potent:
                 V < 1  ->  BE < B  ->  the peptides THIS enzyme actually
                            releases are, on average, MORE potent than the
                            protein's full theoretical fragment population.
                            The enzyme enriches for potency.
                 V > 1  ->  the released peptides are LESS potent on average
                            than the theoretical population — hydrolysis is
                            diluting potency, possibly missing the strongest
                            motifs or degrading them further.
                 V ~ 1  ->  released fragments track the theoretical potency
                            profile with no particular enrichment or dilution.
             V requires B (whole-protein value, from AB_FILE) alongside BE
             (already in your enzyme-action table) — see potency_enrichment()
             below for the paired computation.
  - DHt[%] : enzyme/protein property, NOT activity-specific. It repeats
             across every Activity row for the same Protein+Enzyme — the
             script de-duplicates before using it.

Dependencies: pandas, numpy, scipy, matplotlib
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# ======================================================================
# CONFIG — edit this block to match your real files/labels
# ======================================================================

DATA_FILE = Path("/home/maksym/biopep-automation/Excel tables results/hydrolysis_coef.xlsx")    # your compiled long-format table (csv or xlsx)
OUTPUT_DIR = Path("/home/maksym/biopep-automation/Excel tables results/analysis_output")  # where to write tables/plots
OUTPUT_DIR.mkdir(exist_ok=True)

# Exact strings as they appear in your Activity column
ACE_LABEL = "ACE inhibitor"
DPP_LABEL = "dipeptidyl peptidase IV inhibitor"
TARGET_ACTIVITIES = [ACE_LABEL, DPP_LABEL]

# ----------------------------------------------------------------------
# Enzyme ID constants — defined ONCE here, then reused everywhere below
# (ENZYME_LABELS, PEPSIN_PH13_ID, PEPSIN_PH2_ID, GID_ENZYME_ID) so a typo
# in one string can't silently desync a label from an actual filter/lookup.
# Edit the values on the right of "=" to match your exact "Enzyme id"
# column strings — everything else derives from these.
# ----------------------------------------------------------------------
TRYPSIN_ID       = "12|trypsin|EC 3.4.21.4"
CHYMOTRYPSIN_ID  = "11|chymotrypsin (A)|EC 3.4.21.1"
PAPAIN_ID        = "24|papain|EC 3.4.22.2"
FICIN_ID         = "25|ficin|EC 3.4.22.3"
BROMELAIN_ID     = "30|Stem bromelain|EC 3.4.22.32"
SUBTILISIN_ID    = "41|subtilisin|EC 3.4.21.62"
PEPSIN_PH2_ID    = "39|pepsin (pH > 2)|3.4.23.1"
PEPSIN_PH13_ID   = "13|pepsin (pH 1.3)|EC 3.4.23.1"   # set to None if you only ran pH > 2
GID_ID           = "GID"                       # set to None if you have no GID rows

# Enzyme id -> short label for readable tables/plots.
# Includes every single-enzyme condition PLUS both pepsins and GID, so
# anything that shows up in your raw data (e.g. in data_quality_report,
# which runs on the unfiltered df) renders with a readable name instead
# of the raw BIOPEP id string.
ENZYME_LABELS = {
    TRYPSIN_ID: "Trypsin",
    CHYMOTRYPSIN_ID: "Chymotrypsin",
    PAPAIN_ID: "Papain",
    FICIN_ID: "Ficin",
    BROMELAIN_ID: "Bromelain",
    SUBTILISIN_ID: "Subtilisin (Alcalase)",
    PEPSIN_PH2_ID: "Pepsin (pH>2)",
    PEPSIN_PH13_ID: "Pepsin (pH 1.3)",
    GID_ID: "GID",
}

# Exact Enzyme id string for the simultaneous pepsin+trypsin+chymotrypsin
# gastrointestinal-digestion run. GID is excluded from all single-enzyme
# ranking/comparison functions (it's a 3-enzyme mixture, not one enzyme —
# comparing it head-to-head with individual enzymes would conflate two
# different kinds of thing) but still gets a proper label wherever it
# appears outside those rankings (e.g. the raw data quality check).
GID_ENZYME_ID = GID_ID

# Optional second file: whole-protein A/B parameters from BIOPEP's
# "profile of potential biological activity" + "Calculation" tool.
# This is a DIFFERENT BIOPEP output from your enzyme-action table — A/B do
# not depend on enzyme. Expected columns: Protein Name, Activity, A, B.
AB_FILE = Path("/home/maksym/biopep-automation/Excel tables results/Consolidated_A_B_coef.xlsx")       # set to None if you don't have this yet

MIN_ENZYMES_FOR_CORR = 4   # need at least this many enzymes to trust a Spearman rho


# ======================================================================
# 1. LOAD & CLEAN
# ======================================================================
def load_long_table(path):
    path = Path(path)
    df = pd.read_csv(path) if path.suffix == ".csv" else pd.read_excel(path)
    df.columns = [c.strip() for c in df.columns]
    for c in ["Ae", "W", "Be", "V", "DHt[%]"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["Enzyme"] = df["Enzyme id"].map(ENZYME_LABELS).fillna(df["Enzyme id"])
    return df


def data_quality_report(df):
    print("=== Data quality check ===")
    for col in ["Ae", "W", "Be", "V"]:
        zero_frac = (df[col] == 0).mean()
        na_frac = df[col].isna().mean()
        print(f"  {col:<4s}: {zero_frac:6.1%} exact zeros, {na_frac:6.1%} missing")
    print("  If W/Be/V are ~100% zero even on high-Ae rows, spot-check a few")
    print("  cells directly in the BIOPEP web report before trusting stats below.\n")


def _load_ab_file(ab_path):
    """Load the whole-protein A/B table. Real BIOPEP export columns are
    ID, NAME, Activity, A, B — NAME holds the same protein identifiers
    (e.g. 'cong_beta_1_1') as 'Protein Name' in the enzyme-action table,
    so it's renamed here for joining."""
    ab = pd.read_csv(ab_path) if str(ab_path).endswith(".csv") else pd.read_excel(ab_path)
    ab.columns = [c.strip() for c in ab.columns]
    ab = ab.rename(columns={"NAME": "Protein Name"})
    for c in ["A", "B"]:
        if c in ab.columns:
            ab[c] = pd.to_numeric(ab[c], errors="coerce")
    return ab


def compare_conglutins_by_AB(ab_path, activity):
    if ab_path is None:
        print(f"[A/B comparison for '{activity}'] skipped — no AB_FILE configured.")
        print("  A and B come from BIOPEP's 'profile of potential biological activity'")
        print("  + 'Calculation' tool, a SEPARATE export from your enzyme-action table.")
        return None
    ab = _load_ab_file(ab_path)
    sub = ab[ab["Activity"] == activity].copy()
    if sub.empty:
        print(f"[A/B comparison for '{activity}'] no rows found for this activity.")
        return None
    sub = sub.sort_values("A", ascending=False)
    print(f"\n=== Q: Which conglutin is richest in '{activity}' motifs (pre-hydrolysis)? ===")
    print(sub[["Protein Name", "A", "B"]].to_string(index=False))
    print("  Highest A = most bioactive-motif-dense sequence.")
    print("  Lowest B  = strongest theoretical potency packed into the sequence.\n")
    return sub


def rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids=None):
    exclude_ids = exclude_ids or []
    sub = df[
        (df["Protein Name"] == protein)
        & (df["Activity"] == activity)
        & (~df["Enzyme id"].isin(exclude_ids))
    ].copy()
    if sub.empty:
        return None
    sub = sub[["Enzyme", "Ae", "W", "Be", "V", "DHt[%]"]].drop_duplicates()
    sub["rank_Ae"] = sub["Ae"].rank(ascending=False, method="min")
    sub["rank_Be"] = sub["Be"].rank(ascending=True, method="min")
    return sub.sort_values("rank_Ae")


def report_enzyme_rankings(df, proteins, activities, exclude_ids=None):
    for protein in proteins:
        for activity in activities:
            ranked = rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids)
            if ranked is None or ranked.empty:
                continue
            print(f"\n=== Q: Best enzyme for releasing '{activity}' peptides from {protein}? ===")
            print(ranked.to_string(index=False))
            best_ae = ranked.sort_values("rank_Ae").iloc[0]
            best_be = ranked.sort_values("rank_Be").iloc[0]
            if best_ae["Enzyme"] == best_be["Enzyme"]:
                print(f"  -> {best_ae['Enzyme']} wins on BOTH yield (Ae) and potency (Be).")
            else:
                print(f"  -> Highest yield (Ae): {best_ae['Enzyme']}")
                print(f"  -> Best potency (lowest Be): {best_be['Enzyme']}")
                print("     These disagree — see the AE/BE trade-off analysis below.")


def ae_be_tradeoff(df, protein, activity, exclude_ids=None):
    ranked = rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids)
    if ranked is None:
        return None
    valid = ranked.dropna(subset=["Ae", "Be"])
    if len(valid) < MIN_ENZYMES_FOR_CORR:
        print(f"[{protein} / {activity}] only {len(valid)} enzymes with both Ae & Be "
              f"(need >= {MIN_ENZYMES_FOR_CORR}) — Spearman skipped, reporting raw ranks only.")
        return {"protein": protein, "activity": activity, "n": len(valid),
                "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(valid["Ae"], valid["Be"])
    print(f"\n=== Q: Does enzyme yield (Ae) trade off against potency (Be) "
          f"for {protein} / {activity}? ===")
    print(f"  Spearman rho(Ae, Be) = {rho:.3f}  (p = {p:.3f}, n = {len(valid)} enzymes)")
    if rho > 0.3:
        print("  Positive: enzymes that release MORE peptides also tend to release "
              "MORE POTENT ones (rare, but no trade-off here).")
    elif rho < -0.3:
        print("  Negative: enzymes that release more peptides tend to release "
              "LESS potent ones (or vice versa) — a genuine yield/potency trade-off.")
    else:
        print("  Weak/no relationship: yield and potency vary roughly independently "
              "across enzymes for this conglutin+activity.")
    return {"protein": protein, "activity": activity, "n": len(valid), "rho": rho, "p": p}


def ae_be_tradeoff_all(df, proteins, activities, exclude_ids=None):
    rows = []
    for protein in proteins:
        for activity in activities:
            r = ae_be_tradeoff(df, protein, activity, exclude_ids)
            if r is not None:
                rows.append(r)
    return pd.DataFrame(rows)


def dht_vs_w(df, protein, activity, exclude_ids=None):
    ranked = rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids)
    if ranked is None:
        return None
    valid = ranked.dropna(subset=["DHt[%]", "W"])
    if len(valid) < MIN_ENZYMES_FOR_CORR:
        print(f"[{protein} / {activity}] only {len(valid)} enzymes with DHt & W "
              f"(need >= {MIN_ENZYMES_FOR_CORR}) — Spearman skipped.")
        return {"protein": protein, "activity": activity, "n": len(valid), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(valid["DHt[%]"], valid["W"])
    print(f"\n=== Q: Is high release efficiency (W) just brute-force hydrolysis (DHt), "
          f"for {protein} / {activity}? ===")
    print(f"  Spearman rho(DHt, W) = {rho:.3f}  (p = {p:.3f}, n = {len(valid)} enzymes)")
    if rho > 0.5:
        print("  Strong positive: efficiency tracks overall aggressiveness — not "
              "obviously 'smart' cleavage, just cutting more bonds overall.")
    else:
        print("  Weak/no correlation: some enzymes are efficient at releasing THIS "
              "activity's peptides without simply hydrolyzing more overall — worth "
              "highlighting as a specific-cleavage-site effect.")
    return {"protein": protein, "activity": activity, "n": len(valid), "rho": rho, "p": p}


def potency_enrichment(df, ab_path, protein, activity, exclude_ids=None):
    if ab_path is None:
        print(f"[Potency enrichment for {protein}/{activity}] skipped — "
              "no AB_FILE configured (need whole-protein B).")
        return None
    ab = _load_ab_file(ab_path)
    b_row = ab[(ab["Protein Name"] == protein) & (ab["Activity"] == activity)]
    if b_row.empty:
        print(f"[Potency enrichment for {protein}/{activity}] no B value found in AB_FILE.")
        return None
    b_val = b_row["B"].iloc[0]

    ranked = rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids)
    if ranked is None:
        return None
    ranked = ranked.copy()
    ranked["V_check"] = ranked["Be"] / b_val
    valid = ranked.dropna(subset=["Be"])
    if valid.empty:
        print(f"[Potency enrichment for {protein}/{activity}] no Be values available.")
        return None

    print(f"\n=== Q: Which enzymes enrich for potency vs. dilute it, for {protein}/{activity}? ===")
    print(f"  Whole-protein B = {b_val:.6g}")
    for _, row in valid.sort_values("V_check").iterrows():
        v = row["V_check"]
        if v < 0.9:
            verdict = "ENRICHES for potency (released peptides more potent than the whole-protein average)"
        elif v > 1.1:
            verdict = "DILUTES potency (released peptides weaker than the whole-protein average)"
        else:
            verdict = "roughly neutral (matches whole-protein average potency)"
        reported_v = row["V"] if pd.notna(row["V"]) else np.nan
        print(f"  {row['Enzyme']:<22s} Be={row['Be']:.4g}  V(reported)={reported_v:.4g}"
              f"  V(Be/B)={v:.4g}  -> {verdict}")
    return valid.sort_values("V_check")


def potency_enrichment_all(df, ab_path, proteins, activities, exclude_ids=None):
    for protein in proteins:
        for activity in activities:
            potency_enrichment(df, ab_path, protein, activity, exclude_ids)


def pepsin_ph_comparison(df, activity):
    if PEPSIN_PH13_ID is None:
        print(f"\n[pH comparison for '{activity}'] skipped — PEPSIN_PH13_ID not set.")
        print("  You need rows for BOTH pepsin conditions to compare them; set")
        print("  PEPSIN_PH13_ID once you've run the pH 1.3 hydrolysis in BIOPEP.")
        return None
    ph13 = df[(df["Enzyme id"] == PEPSIN_PH13_ID) & (df["Activity"] == activity)]
    ph2 = df[(df["Enzyme id"] == PEPSIN_PH2_ID) & (df["Activity"] == activity)]
    merged = ph13.merge(ph2, on="Protein Name", suffixes=("_pH1.3", "_pH>2"))
    if merged.empty:
        print(f"[pH comparison for '{activity}'] no overlapping proteins found.")
        return None
    print(f"\n=== Q: Does pepsin pH condition change what gets released, for '{activity}'? ===")
    for col in ["Ae", "W", "Be", "DHt[%]"]:
        a, b = merged[f"{col}_pH1.3"].dropna(), merged[f"{col}_pH>2"].dropna()
        n = min(len(a), len(b))
        if n < 6:
            print(f"  {col}: only {n} paired proteins — too few for a signed-rank test, "
                  "reporting means only.")
            print(f"    mean pH1.3 = {a.mean():.4g}, mean pH>2 = {b.mean():.4g}")
            continue
        stat, p = stats.wilcoxon(a.iloc[:n].values, b.iloc[:n].values)
        print(f"  {col}: Wilcoxon signed-rank p = {p:.3f}  "
              f"(mean pH1.3={a.mean():.4g}, mean pH>2={b.mean():.4g})")
    return merged


def cross_activity_comparison(df, proteins, exclude_ids=None):
    rows = []
    for protein in proteins:
        best = {}
        for activity in TARGET_ACTIVITIES:
            ranked = rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids)
            best[activity] = ranked["Ae"].max() if ranked is not None and not ranked.empty else np.nan
        rows.append({"Protein Name": protein, **best})
    out = pd.DataFrame(rows)
    print("\n=== Q: Is each conglutin a better ACE-I source or DPP-IV source "
          "(best achievable Ae across enzymes)? ===")
    print(out.to_string(index=False))
    return out


def enzyme_recommendation_table(df, proteins, activities, exclude_ids=None):
    all_rows = []
    for protein in proteins:
        for activity in activities:
            ranked = rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids)
            if ranked is None or ranked.empty:
                continue
            ranked = ranked.copy()
            ranked["composite_rank"] = ranked[["rank_Ae", "rank_Be"]].mean(axis=1)
            ranked["Protein Name"] = protein
            ranked["Activity"] = activity
            all_rows.append(ranked.sort_values("composite_rank"))
    if not all_rows:
        print("No data available to build recommendation table.")
        return None
    result = pd.concat(all_rows, ignore_index=True)
    out_path = OUTPUT_DIR / "enzyme_recommendation_table.csv"
    result.to_csv(out_path, index=False)
    print(f"\nFull enzyme recommendation table written to: {out_path}")
    return result


def enzyme_summary_boxplot(df, activity, exclude_ids=None):
    """One boxplot per enzyme, pooled across all proteins. No per-point
    labels needed at all — the enzyme identity IS the x-axis category.
    Answers 'which enzyme tends to yield more / more potent peptides,
    overall, across the whole species' without the label-overlap problem
    the original single scatter had at n~70-80 points.

    Built with plain matplotlib ax.boxplot() rather than pandas'
    DataFrame.boxplot(by=...) convenience wrapper — that wrapper interacts
    badly with sharex=True across two stacked subplots (it can silently
    double the internal tick-locator count), so grouping is done by hand
    here for reliability."""
    exclude_ids = exclude_ids or []
    sub = df[(df["Activity"] == activity) & (~df["Enzyme id"].isin(exclude_ids))].copy()
    if sub.empty:
        print(f"[enzyme summary boxplot for '{activity}'] no data.")
        return

    order = sub.groupby("Enzyme")["Ae"].median().sort_values().index.tolist()
    ae_groups = [sub.loc[sub["Enzyme"] == e, "Ae"].dropna().values for e in order]
    be_groups = [sub.loc[sub["Enzyme"] == e, "Be"].dropna().values for e in order]

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axes[0].boxplot(ae_groups, tick_labels=order, patch_artist=True)
    axes[0].set_ylabel("Ae (release frequency)")
    axes[0].set_title(f"Yield (Ae) across all conglutins — {activity}")

    axes[1].boxplot(be_groups, tick_labels=order, patch_artist=True)
    axes[1].set_ylabel("Be (release potency — LOWER = more potent)")
    axes[1].set_title(f"Potency (Be) across all conglutins — {activity}")
    plt.setp(axes[1].get_xticklabels(), rotation=30, ha="right")

    fig.tight_layout()
    fname = OUTPUT_DIR / f"enzyme_summary_{activity.replace(' ', '_')}.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"Saved: {fname}")
    print(f"  Each box pools all {sub['Protein Name'].nunique()} conglutins for that enzyme.")
    print("  An enzyme high in the top panel but also high (worse) in the bottom panel")
    print("  is a yield/potency trade-off visible in miniature.")


def quadrant_plot(df, proteins, activities, exclude_ids=None):
    """Small multiples: one panel PER PROTEIN, each with only its own
    ~7-8 enzymes labeled — replaces the original single crowded scatter
    where up to 10 proteins x 8 enzymes (~80 labeled points) were overlaid
    on one axes and became unreadable."""
    for activity in activities:
        n = len(proteins)
        ncols = 2
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4 * nrows), squeeze=False)
        for idx, protein in enumerate(proteins):
            ax = axes[idx // ncols][idx % ncols]
            ranked = rank_enzymes_for_protein_activity(df, protein, activity, exclude_ids)
            if ranked is None or ranked.empty:
                ax.set_visible(False)
                continue
            valid = ranked.dropna(subset=["Ae", "Be"])
            if valid.empty:
                ax.set_visible(False)
                continue
            ax.scatter(valid["Ae"], valid["Be"], alpha=0.8, s=60)
            for _, row in valid.iterrows():
                ax.annotate(row["Enzyme"], (row["Ae"], row["Be"]), fontsize=6.5,
                            alpha=0.8, xytext=(3, 3), textcoords="offset points")
            ax.set_title(protein, fontsize=10)
            ax.set_xlabel("Ae", fontsize=8)
            ax.set_ylabel("Be (lower=better)", fontsize=8)
            ax.tick_params(labelsize=7)
        for idx in range(n, nrows * ncols):
            axes[idx // ncols][idx % ncols].set_visible(False)
        fig.suptitle(f"Yield vs potency, one panel per conglutin — {activity}", fontsize=12)
        fig.tight_layout()
        fname = OUTPUT_DIR / f"quadrant_faceted_{activity.replace(' ', '_')}.png"
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved: {fname}")


def dht_vs_w_plot(df, proteins, activities, exclude_ids=None):
    """Degree of hydrolysis (DHt) vs release efficiency (W). Renamed from
    the original 'aggressiveness' framing — DHt is the standard term used
    throughout this script and the source literature, 'aggressiveness' was
    not. Colored by ENZYME (not protein), since the question this answers
    is enzyme-centric; a pooled linear trend + Spearman rho are annotated
    directly on the plot instead of leaving ~70-80 unexplained dots for
    the viewer to interpret unaided."""
    exclude_ids = exclude_ids or []
    for activity in activities:
        sub = df[(df["Activity"] == activity) & (~df["Enzyme id"].isin(exclude_ids))].copy()
        sub = sub.dropna(subset=["DHt[%]", "W"])
        if sub.empty:
            print(f"[DHt vs W plot for '{activity}'] no data.")
            continue

        fig, ax = plt.subplots(figsize=(7, 5.5))
        enzymes_present = sorted(sub["Enzyme"].unique())
        cmap = plt.get_cmap("tab10")
        for i, enzyme in enumerate(enzymes_present):
            e_sub = sub[sub["Enzyme"] == enzyme]
            ax.scatter(e_sub["DHt[%]"], e_sub["W"], label=enzyme,
                       color=cmap(i % 10), alpha=0.8, s=50)

        x, y = sub["DHt[%]"].values, sub["W"].values
        if len(x) >= 4:
            coeffs = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 100)
            ax.plot(xs, np.polyval(coeffs, xs), "k--", linewidth=1, alpha=0.6,
                    label="linear trend (all points pooled)")
            rho, p = stats.spearmanr(x, y)
            ax.text(0.03, 0.97, f"Pooled Spearman rho = {rho:.2f} (p = {p:.3g}, n = {len(x)})",
                    transform=ax.transAxes, fontsize=9, va="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

        ax.set_xlabel("Degree of hydrolysis, DHt [%]")
        ax.set_ylabel("W (bioactive-release efficiency)")
        ax.set_title(f"Degree of hydrolysis vs release efficiency — {activity}")
        ax.legend(fontsize=6.5, loc="lower right", ncol=1)
        fig.text(0.5, 0.005,
                  f"Each point = one conglutin x enzyme combination "
                  f"(n = {len(x)} across {sub['Protein Name'].nunique()} conglutins).",
                  ha="center", fontsize=8, style="italic")
        fig.tight_layout(rect=[0, 0.03, 1, 1])
        fname = OUTPUT_DIR / f"DHt_vs_W_{activity.replace(' ', '_')}.png"
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved: {fname}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    df = load_long_table(DATA_FILE)
    data_quality_report(df)

    proteins = sorted(df["Protein Name"].unique())
    exclude = [e for e in [GID_ENZYME_ID] if e is not None]

    print("############################################################")
    print("# Q1 — conglutin richness before hydrolysis (A, B)")
    print("############################################################")
    for activity in TARGET_ACTIVITIES:
        compare_conglutins_by_AB(AB_FILE, activity)

    print("\n############################################################")
    print("# Q2/Q3/Q4 — best enzyme per conglutin (Ae / W / Be)")
    print("############################################################")
    report_enzyme_rankings(df, proteins, TARGET_ACTIVITIES, exclude)

    print("\n############################################################")
    print("# Q5 — Ae vs Be trade-off (the correlation you asked for)")
    print("############################################################")
    tradeoff_df = ae_be_tradeoff_all(df, proteins, TARGET_ACTIVITIES, exclude)
    tradeoff_path = OUTPUT_DIR / "ae_be_tradeoff_spearman.csv"
    tradeoff_df.to_csv(tradeoff_path, index=False)
    print(f"\nSaved trade-off summary table to: {tradeoff_path}")

    print("\n############################################################")
    print("# Q6 — DHt vs W: is efficiency just brute-force hydrolysis?")
    print("############################################################")
    for protein in proteins:
        for activity in TARGET_ACTIVITIES:
            dht_vs_w(df, protein, activity, exclude)

    print("\n############################################################")
    print("# Q6b — potency enrichment vs dilution (V = Be/B)")
    print("############################################################")
    potency_enrichment_all(df, AB_FILE, proteins, TARGET_ACTIVITIES, exclude)

    print("\n############################################################")
    print("# Q7 — pepsin pH 1.3 vs pH > 2")
    print("############################################################")
    for activity in TARGET_ACTIVITIES:
        pepsin_ph_comparison(df, activity)

    print("\n############################################################")
    print("# Q8 — ACE-I source vs DPP-IV source, per conglutin")
    print("############################################################")
    cross_activity_comparison(df, proteins, exclude)

    print("\n############################################################")
    print("# Q9 — composite enzyme recommendation table + plots")
    print("############################################################")
    enzyme_recommendation_table(df, proteins, TARGET_ACTIVITIES, exclude)
    for activity in TARGET_ACTIVITIES:
        enzyme_summary_boxplot(df, activity, exclude)
    quadrant_plot(df, proteins, TARGET_ACTIVITIES, exclude)
    dht_vs_w_plot(df, proteins, TARGET_ACTIVITIES, exclude)

    print("\nDone. All tables/plots are in:", OUTPUT_DIR.resolve())
