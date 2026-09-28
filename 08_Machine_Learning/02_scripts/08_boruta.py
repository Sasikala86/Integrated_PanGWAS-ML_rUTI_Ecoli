#!/usr/bin/env python3
# ============================================================
# ML_Accessory — Boruta all-relevant feature selection
# ============================================================
from pathlib import Path
# CELL 7 — BORUTA FEATURE SELECTION (Chromosome / Plasmid / Combined)
# =====================================================
# Comments / Methods:
# 1️⃣ Input: filtered ML matrices (constant genes removed, 1–99% prevalence, Cell 5)
# 2️⃣ Boruta: All-relevant feature selection using Random Forest
# 3️⃣ Output: Selected, Tentative, Rejected features
# 4️⃣ Saved CSVs in PROCESSED_DIR/feature_selection/results
# 5️⃣ Separate run for chromosome, plasmid, combined datasets
# 6️⃣ Later Venn diagram: Compare Boruta features with Chi² & MI (union/intersection)
# =====================================================

from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier
import pandas as pd

# Folder to save Boruta results
FS_RESULTS = PROCESSED_DIR / "feature_selection" / "results"
FS_RESULTS.mkdir(parents=True, exist_ok=True)

def run_boruta(df_filtered, dataset_name, n_estimators=100, random_state=42):
    """
    Perform Boruta all-relevant feature selection.
    Returns selected, tentative, rejected features and saves CSVs.
    """
    print(f"\n--- Running Boruta: {dataset_name} ---")
    
    gene_cols = df_filtered.columns.difference(["phenotype", "dataset"])
    X = df_filtered[gene_cols].astype(int).values
    y = df_filtered["phenotype"].astype(int).values  # ensure integer labels

    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=-1
    )
    
    boruta_selector = BorutaPy(
        estimator=rf,
        n_estimators='auto',
        random_state=random_state,
        verbose=2
    )
    
    boruta_selector.fit(X, y)
    
    # Selected / Tentative / Rejected features
    selected_features = gene_cols[boruta_selector.support_].tolist()
    tentative_features = gene_cols[boruta_selector.support_weak_].tolist()
    rejected_features = [g for g in gene_cols if g not in selected_features + tentative_features]
    
    print(f"✔ {dataset_name}: Selected={len(selected_features)}, Tentative={len(tentative_features)}, Rejected={len(rejected_features)}")
    
    # -------------------
    # Save CSVs
    # -------------------
    pd.DataFrame({"gene": selected_features}).to_csv(FS_RESULTS / f"{dataset_name}_boruta_selected.csv", index=False)
    pd.DataFrame({"gene": tentative_features}).to_csv(FS_RESULTS / f"{dataset_name}_boruta_tentative.csv", index=False)
    pd.DataFrame({"gene": rejected_features}).to_csv(FS_RESULTS / f"{dataset_name}_boruta_rejected.csv", index=False)
    
    return selected_features, tentative_features, rejected_features

# ====================================================
# 7a — Chromosome
# ====================================================
selected_chr, tentative_chr, rejected_chr = run_boruta(chrom_df_filtered, "chromosome")

# ====================================================
# 7b — Plasmid
# ====================================================
selected_pls, tentative_pls, rejected_pls = run_boruta(plasmid_df_filtered, "plasmid")

# ====================================================
# 7c — Combined
# ====================================================
selected_comb, tentative_comb, rejected_comb = run_boruta(combined_df_filtered, "combined")
