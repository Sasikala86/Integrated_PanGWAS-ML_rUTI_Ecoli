#!/usr/bin/env python3
# ============================================================
# ML_Accessory — EDA
# ============================================================
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = BASE_DIR / "processed"
EDA_RESULTS = BASE_DIR / "results" / "EDA" / "results"
EDA_PLOTS = BASE_DIR / "results" / "EDA" / "plots"
EDA_RESULTS.mkdir(parents=True, exist_ok=True)
EDA_PLOTS.mkdir(parents=True, exist_ok=True)

def run_eda(ml_df, dataset_name):
    gene_cols = ml_df.columns.difference(["phenotype", "dataset"])
    prevalence = ml_df[gene_cols].mean()
    prevalence.to_csv(EDA_RESULTS / f"{dataset_name}_gene_prevalence.csv")
    variance = ml_df[gene_cols].var()
    variance.to_csv(EDA_RESULTS / f"{dataset_name}_gene_variance.csv")
    gene_counts = ml_df[gene_cols].sum(axis=1)
    gene_counts.to_csv(EDA_RESULTS / f"{dataset_name}_gene_count_per_sample.csv")
    summary = pd.DataFrame({
        "dataset": [dataset_name],
        "n_genes": [len(gene_cols)],
        "constant_genes": [(variance == 0).sum()],
        "mean_genes_per_sample": [gene_counts.mean()]
    })
    summary.to_csv(EDA_RESULTS / f"{dataset_name}_eda_summary.csv", index=False)
    return prevalence, variance, gene_counts

def plot_variance(variance, dataset_name):
    plt.figure(figsize=(6, 4))
    sns.histplot(variance, bins=50)
    plt.xlabel("Gene variance")
    plt.ylabel("Number of genes")
    plt.title(f"{dataset_name}: Gene variance distribution")
    plt.tight_layout()
    plt.savefig(EDA_PLOTS / f"{dataset_name}_gene_variance.png", dpi=600)
    plt.close()

def plot_gene_counts(gene_counts, dataset_name):
    plt.figure(figsize=(6, 4))
    sns.histplot(gene_counts, bins=40)
    plt.xlabel("Genes per isolate")
    plt.ylabel("Number of isolates")
    plt.title(f"{dataset_name}: Gene burden per isolate")
    plt.tight_layout()
    plt.savefig(EDA_PLOTS / f"{dataset_name}_gene_count_per_sample.png", dpi=600)
    plt.close()

files = {
    "chromosome": PROCESSED_DIR / "ml_matrix_chromosomes.csv",
    "plasmid": PROCESSED_DIR / "ml_matrix_plasmids.csv",
    "combined": PROCESSED_DIR / "ml_matrix_combined.csv",
}
for ds, path in files.items():
    df = pd.read_csv(path, index_col=0)
    prev, var, counts = run_eda(df, ds)
    plot_variance(var, ds)
    plot_gene_counts(counts, ds)

print("EDA completed.")
