#!/usr/bin/env python3
# MASTER BARPLOT – STRICT FEATURES ONLY
# Chromosome → Plasmid → Combined
# All 5 Models
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================
# PATH SETTINGS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[2]
RESULTS_DIR = BASE_DIR / "results/Models"

FEATURE_TYPE = "STRICT"

# ============================================================
# LOAD ALL STRICT METRICS FILES
# ============================================================

all_metrics = []

for metrics_file in RESULTS_DIR.rglob("*_STRICT_metrics.csv"):
    df = pd.read_csv(metrics_file)
    all_metrics.append(df)

if len(all_metrics) == 0:
    raise ValueError("❌ No STRICT metrics files found. Check folder path.")

metrics_df = pd.concat(all_metrics, ignore_index=True)

# ============================================================
# KEEP REQUIRED COLUMNS
# ============================================================

metrics_df = metrics_df[["Model", "Dataset", "Test_AUC"]]

# ============================================================
# ORDER DATASET
# ============================================================

metrics_df["Dataset"] = pd.Categorical(
    metrics_df["Dataset"],
    categories=["chromosome", "plasmid", "combined"],
    ordered=True
)

# ============================================================
# ORDER MODELS
# ============================================================

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

# Add AUC values on bars
for container in ax.containers:
    ax.bar_label(container, fmt="%.3f", padding=3)

plt.title("Model Comparison Across Datasets (STRICT Features)")
plt.ylabel("Test AUC")
plt.xlabel("Dataset")

plt.legend(title="Model", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()

output_file = BASE_DIR / "All_Models_STRICT_Comparison.png"
plt.savefig(output_file, dpi=300)
plt.close()

print("\n✔ STRICT comparison plot saved as:")
print(output_file)
# ============================================================
