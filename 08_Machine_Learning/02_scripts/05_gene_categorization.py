#!/usr/bin/env python3
# ============================================================
# ML_Accessory — 4-TIER GENE CATEGORIZATION
# ============================================================
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
EDA_RESULTS = BASE_DIR / "results" / "EDA" / "results"
EDA_RESULTS.mkdir(parents=True, exist_ok=True)

def categorize_genes_4tier(prevalence, dataset_name):
    core = prevalence[prevalence >= 0.99].index.tolist()
    soft_core = prevalence[(prevalence >= 0.95) & (prevalence < 0.99)].index.tolist()
    shell = prevalence[(prevalence >= 0.15) & (prevalence < 0.95)].index.tolist()
    rare = prevalence[prevalence < 0.05].index.tolist()
    counts = pd.DataFrame({
        "dataset": [dataset_name],
        "total_genes": [len(prevalence)],
        "core_genes": [len(core)],
        "soft_core_genes": [len(soft_core)],
        "shell_genes": [len(shell)],
        "rare_genes": [len(rare)]
    })
    return counts, core, soft_core, shell, rare

all_counts = []
for ds in ["chromosome", "plasmid", "combined"]:
    prevalence = pd.read_csv(
        EDA_RESULTS / f"{ds}_gene_prevalence.csv", index_col=0
    ).iloc[:, 0]
    # Handle the normal Series CSV structure robustly.
    if prevalence.dtype == object:
        prevalence = pd.to_numeric(prevalence, errors="coerce")
    counts, core, soft, shell, rare = categorize_genes_4tier(prevalence, ds)
    all_counts.append(counts)
    for name, genes in [
        ("core", core), ("soft_core", soft), ("shell", shell), ("rare", rare)
    ]:
        pd.DataFrame({"gene": genes}).to_csv(
            EDA_RESULTS / f"{ds}_{name}_genes.csv", index=False
        )

pd.concat(all_counts, ignore_index=True).to_csv(
    EDA_RESULTS / "gene_category_4tier_counts.csv", index=False
)
print("4-tier gene categorization completed.")
