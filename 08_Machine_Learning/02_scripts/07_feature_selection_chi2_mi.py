#!/usr/bin/env python3
# ============================================================
# ML_Accessory — Chi² + Mutual Information feature screening
# ============================================================
from pathlib import Path
import pandas as pd
import numpy as np
import warnings
from sklearn.feature_selection import chi2, mutual_info_classif

warnings.filterwarnings("ignore", category=UserWarning)

BASE_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = BASE_DIR / "processed"
FS_RESULTS = PROCESSED_DIR / "feature_selection" / "results"
FS_PLOTS   = PROCESSED_DIR / "feature_selection/plots"
FS_RESULTS.mkdir(parents=True, exist_ok=True)
FS_PLOTS.mkdir(parents=True, exist_ok=True)

def run_chi2_mi(df_filtered, dataset_name):
    """
    Feature screening using Chi² and Mutual Information
    Input: ML matrices filtered for constant + 1-99% prevalence genes
    """
    # Gene columns
    gene_cols = df_filtered.columns.difference(["phenotype", "dataset"])

    X = df_filtered[gene_cols].astype(int).values
    y = df_filtered["phenotype"].astype(int).values

    # -------------------
    # Chi² (requires non-negative values)
    # -------------------
    chi_scores, chi_pvalues = chi2(X, y)

    chi_df = pd.DataFrame({
        "gene": gene_cols,
        "chi2_score": chi_scores,
        "p_value": chi_pvalues
    })

    # Significant genes at p <= 0.05
    chi_sig = chi_df[chi_df["p_value"] <= 0.05]
    chi_sig.to_csv(FS_RESULTS / f"{dataset_name}_chi2_significant.csv", index=False)

    # -------------------
    # Mutual Information
    # -------------------
    mi_scores = mutual_info_classif(X, y, random_state=42)
    mi_df = pd.DataFrame({
        "gene": gene_cols,
        "mi_score": mi_scores
    })

    # Significant genes with MI > 0
    mi_sig = mi_df[mi_df["mi_score"] > 0]
    mi_sig.to_csv(FS_RESULTS / f"{dataset_name}_mi_significant.csv", index=False)

    # -------------------
    # Intersection & Union
    # -------------------
    chi_set = set(chi_sig["gene"])
    mi_set  = set(mi_sig["gene"])

    inter = sorted(chi_set & mi_set)
    union = sorted(chi_set | mi_set)

    pd.DataFrame({"gene": inter}).to_csv(FS_RESULTS / f"{dataset_name}_chi2_mi_intersection.csv", index=False)
    pd.DataFrame({"gene": union}).to_csv(FS_RESULTS / f"{dataset_name}_chi2_mi_union.csv", index=False)

    print(
        f"✔ {dataset_name}: "
        f"Chi² sig={len(chi_sig)}, "
        f"MI sig={len(mi_sig)}, "
        f"Intersection={len(inter)}, "
        f"Union={len(union)}"
    )

    return chi_sig, mi_sig, inter, union


# -------------------
# Run for all datasets
# -------------------
chi2_chr, mi_chr, inter_chr, union_chr       = run_chi2_mi(chrom_df_filtered, "chromosome")
chi2_pls, mi_pls, inter_pls, union_pls       = run_chi2_mi(plasmid_df_filtered, "plasmid")
chi2_comb, mi_comb, inter_comb, union_comb   = run_chi2_mi(combined_df_filtered, "combined")
