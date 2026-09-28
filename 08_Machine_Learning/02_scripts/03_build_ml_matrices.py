#!/usr/bin/env python3
# ============================================================
# ML_Accessory — BUILD ML MATRICES
# ============================================================
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = BASE_DIR / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

def load_pangenome(full_path):
    return pd.read_csv(full_path, sep="\t", index_col=0)

def prepare_ml_matrix(df, pheno_df):
    X = df.T
    X = X.merge(pheno_df, left_index=True, right_on="sample")
    X.set_index("sample", inplace=True)
    return X

df_chr = load_pangenome(DATA_DIR / "gene_presence_absence_chromosomes.Rtab")
df_pls = load_pangenome(DATA_DIR / "gene_presence_absence_Plasmids.Rtab")
df_comb = load_pangenome(DATA_DIR / "gene_presence_absence_combined.Rtab")

pheno_chr = pd.read_csv(PROCESSED_DIR / "phenotype_chromosome.csv")
pheno_pls = pd.read_csv(PROCESSED_DIR / "phenotype_plasmid.csv")
pheno_comb = pd.read_csv(PROCESSED_DIR / "phenotype_combined.csv")

ml_chr = prepare_ml_matrix(df_chr, pheno_chr)
ml_pls = prepare_ml_matrix(df_pls, pheno_pls)
ml_comb = prepare_ml_matrix(df_comb, pheno_comb)

print("Chromosome ML shape:", ml_chr.shape)
print("Plasmid ML shape:", ml_pls.shape)
print("Combined ML shape:", ml_comb.shape)

ml_chr.to_csv(PROCESSED_DIR / "ml_matrix_chromosomes.csv")
ml_pls.to_csv(PROCESSED_DIR / "ml_matrix_plasmids.csv")
ml_comb.to_csv(PROCESSED_DIR / "ml_matrix_combined.csv")

print("ML matrices saved in:", PROCESSED_DIR)
