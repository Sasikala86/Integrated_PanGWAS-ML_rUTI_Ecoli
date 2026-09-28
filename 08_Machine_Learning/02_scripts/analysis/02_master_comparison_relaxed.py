#!/usr/bin/env python3
# FINAL MASTER BARPLOT
# Chromosome → Plasmid → Combined
# All Models (RELAXED by default)
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================
# PATH
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[2]
RESULTS_DIR = BASE_DIR / "results/Models"

FEATURE_TYPE = "RELAXED"   # Change to "STRICT" if needed

# ============================================================
# LOAD ALL METRICS (AUTO-DETECT)
# ============================================================

all_metrics = []

for metrics_file in RESULTS_DIR.rglob(f"*_{FEATURE_TYPE}_metrics.csv"):

    df = pd.read_csv(metrics_file)

    # Standardize model names
    if df["Model"].iloc[0] == "ExtraTrees":
        df["Model"] = "ExtraTrees"
    elif df["Model"].iloc[0] == "Logistic":
        df["Model"] = "Logistic"
    elif df["Model"].iloc[0] == "RF":
        df["Model"] = "RF"
    elif df["Model"].iloc[0] == "XGB":
        df["Model"] = "XGB"
    elif df["Model"].iloc[0] == "SVM":
        df["Model"] = "SVM"

    all_metrics.append(df)

metrics_df = pd.concat(all_metrics, ignore_index=True)

# ============================================================
# KEEP ONLY REQUIRED COLUMNS
# ============================================================

metrics_df = metrics_df[["Model", "Dataset", "Test_AUC"]]

# Order datasets
metrics_df["Dataset"] = pd.Categorical(
    metrics_df["Dataset"],
    categories=["chromosome", "plasmid", "combined"],
    ordered=True
)

# Order models (optional but cleaner)
metrics_df["Model"] = pd.Categorical(
    metrics_df["Model"],
    categories=["RF", "XGB", "ExtraTrees", "SVM", "Logistic"],
    ordered=True
)

# ============================================================
# PLOT
# ============================================================

plt.figure(figsize=(15,7))

ax = sns.barplot(
    data=metrics_df,
    x="Dataset",
    y="Test_AUC",
    hue="Model",
    errorbar=None
)

# Print AUC on bars
for container in ax.containers:
    ax.bar_label(container, fmt="%.3f", padding=3)

plt.title(f"Model Comparison Across Datasets ({FEATURE_TYPE})")
plt.ylabel("Test AUC")
plt.xlabel("Dataset")

plt.legend(title="Model", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()
plt.savefig(BASE_DIR / f"All_Models_{FEATURE_TYPE}_Comparison.png", dpi=300)
plt.close()

print("\n✔ Final comparison plot saved as:")
print(f"   All_Models_{FEATURE_TYPE}_Comparison.png")
# ============================================================
