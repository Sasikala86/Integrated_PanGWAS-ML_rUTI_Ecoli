#!/usr/bin/env python3
# ============================================================
# ML_Accessory — GENE FILTERING
# Constant genes + prevalence 1–99%
# ============================================================
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = BASE_DIR / "processed"

def filter_genes_ml(ml_df, dataset_name, min_prev=0.01, max_prev=0.99):
    gene_cols = ml_df.columns.difference(["phenotype", "dataset"])
    variance = ml_df[gene_cols].var()
    prevalence = ml_df[gene_cols].mean()
    constant_genes = variance[variance == 0].index.tolist()
    prev_genes = prevalence[
        (prevalence >= min_prev) & (prevalence <= max_prev)
    ].index.tolist()
    genes_to_keep = [g for g in prev_genes if g not in constant_genes]
    filtered_df = ml_df[genes_to_keep + ["phenotype", "dataset"]].copy()
    print(f"{dataset_name}: {len(genes_to_keep)} genes retained after filtering ({len(gene_cols)} original)")
    filtered_file = PROCESSED_DIR / f"{dataset_name}_ML_filtered_1-99.csv"
    filtered_df.to_csv(filtered_file, index=False)
    print(f"Saved: {filtered_file}")
    return filtered_df

for ds, fn in [
    ("chromosome", "ml_matrix_chromosomes.csv"),
    ("plasmid", "ml_matrix_plasmids.csv"),
    ("combined", "ml_matrix_combined.csv"),
]:
    df = pd.read_csv(PROCESSED_DIR / fn)
    filter_genes_ml(df, ds)
