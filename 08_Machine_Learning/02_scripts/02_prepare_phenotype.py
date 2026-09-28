#!/usr/bin/env python3
# ============================================================
# ML_Accessory — PHENOTYPE PREPARATION
# ============================================================
# Source logic preserved from the original workflow.
# rUTI -> 1; UTI -> 0.
# ============================================================
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = BASE_DIR / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

def load_pangenome(full_path):
    return pd.read_csv(full_path, sep="\t", index_col=0)

def create_phenotype(df, dataset_name):
    samples = df.columns
    phenotype = []

    for s in samples:
        if "rUTI" in s:
            phenotype.append(1)
        elif "UTI" in s:
            phenotype.append(0)
        else:
            raise ValueError(f"Phenotype not found in sample name: {s}")

    return pd.DataFrame({
        "sample": samples,
        "phenotype": phenotype,
        "dataset": dataset_name
    })

df_chr = load_pangenome(DATA_DIR / "gene_presence_absence_chromosomes.Rtab")
df_pls = load_pangenome(DATA_DIR / "gene_presence_absence_Plasmids.Rtab")
df_comb = load_pangenome(DATA_DIR / "gene_presence_absence_combined.Rtab")

pheno_chr = create_phenotype(df_chr, "chromosome")
pheno_pls = create_phenotype(df_pls, "plasmid")
pheno_comb = create_phenotype(df_comb, "combined")

pheno_chr.to_csv(PROCESSED_DIR / "phenotype_chromosome.csv", index=False)
pheno_pls.to_csv(PROCESSED_DIR / "phenotype_plasmid.csv", index=False)
pheno_comb.to_csv(PROCESSED_DIR / "phenotype_combined.csv", index=False)

print("Chromosome phenotype counts:\n", pheno_chr["phenotype"].value_counts())
print("Plasmid phenotype counts:\n", pheno_pls["phenotype"].value_counts())
print("Combined phenotype counts:\n", pheno_comb["phenotype"].value_counts())
print("\nPhenotype files saved in:", PROCESSED_DIR)
