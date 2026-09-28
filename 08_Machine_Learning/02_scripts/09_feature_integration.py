#!/usr/bin/env python3
# CELL 8 — FINAL FEATURE INTEGRATION (CURRENT FLOW)
# ============================================================
# Methods-ready explanation:
# 1️⃣ Integrates filter-based (Chi², MI) and wrapper-based (Boruta) features
# 2️⃣ STRICT set  = (Chi² ∩ MI) ∩ Boruta  → robust biomarkers
# 3️⃣ RELAXED set = (Chi² ∪ MI) ∪ Boruta  → exploratory ML features
# 4️⃣ Separate processing for chromosome, plasmid, combined datasets
# ============================================================

import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# -----------------------------
# Correct paths for current flow
# -----------------------------
FS_RESULTS = PROCESSED_DIR / "feature_selection" / "results"
VENN_DIR   = PROCESSED_DIR / "feature_selection/plots/venn_diagrams"
VENN_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Safe CSV loader
# -----------------------------
def load_genes(csv_path):
    if not csv_path.exists():
        raise FileNotFoundError(f"❌ Missing file: {csv_path}")
    df = pd.read_csv(csv_path)
    if df.empty:
        return set()
    return set(df["gene"].astype(str))

# -----------------------------
# Integration logic
# -----------------------------
def integrate_features(dataset):
    """
    Integrates Chi², MI and Boruta feature sets for a dataset
    """
    # ---- Load Chi² & MI ----
    chi2_genes = load_genes(FS_RESULTS / f"{dataset}_chi2_significant.csv")
    mi_genes   = load_genes(FS_RESULTS / f"{dataset}_mi_significant.csv")

    # ---- Load Boruta ----
    boruta_genes = load_genes(FS_RESULTS / f"{dataset}_boruta_selected.csv")

    # ---- Intermediate sets ----
    chi_mi_inter = chi2_genes & mi_genes
    chi_mi_union = chi2_genes | mi_genes

    # ---- Final feature sets ----
    strict_final  = chi_mi_inter & boruta_genes
    relaxed_final = chi_mi_union | boruta_genes

    # ---- Save final sets ----
    pd.DataFrame({"gene": sorted(strict_final)}).to_csv(
        FS_RESULTS / f"{dataset}_FINAL_STRICT_features.csv",
        index=False
    )

    pd.DataFrame({"gene": sorted(relaxed_final)}).to_csv(
        FS_RESULTS / f"{dataset}_FINAL_RELAXED_features.csv",
        index=False
    )

    # ---- Venn diagram (Boruta vs Chi²∩MI) ----
    plt.figure(figsize=(5, 5))
    venn2(
        [boruta_genes, chi_mi_inter],
        set_labels=("Boruta", "Chi² ∩ MI")
    )
    plt.title(f"{dataset.capitalize()} feature overlap")
    plt.tight_layout()
    plt.savefig(VENN_DIR / f"{dataset}_boruta_vs_chi2_mi.png", dpi=300)
    plt.close()

    # ---- Summary ----
    print(f"\n✔ {dataset.upper()}")
    print(f"  Chi² significant      : {len(chi2_genes)}")
    print(f"  MI significant        : {len(mi_genes)}")
    print(f"  Chi² ∩ MI             : {len(chi_mi_inter)}")
    print(f"  Boruta selected       : {len(boruta_genes)}")
    print(f"  FINAL STRICT features : {len(strict_final)}")
    print(f"  FINAL RELAXED features: {len(relaxed_final)}")

# -----------------------------
# Run integration for all datasets
# -----------------------------
for ds in ["chromosome", "plasmid", "combined"]:
    integrate_features(ds)

print("\n✅ CELL 8 COMPLETED SUCCESSFULLY")
print(f"📁 Final CSVs saved in: {FS_RESULTS}")
print(f"📊 Venn diagrams saved in: {VENN_DIR}")
