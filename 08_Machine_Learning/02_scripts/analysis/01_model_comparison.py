#!/usr/bin/env python3
# MODEL PERFORMANCE COMPARISON
# Barplot 1 → Dataset Comparison
# Barplot 2 → Model Comparison (Best Dataset)
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

# ============================================================
# LOAD ALL METRICS
# ============================================================

all_metrics = []

models = ["RF", "XGB", "ExtraTrees", "SVM", "Logistic"]
datasets = ["chromosome", "plasmid", "combined"]
feature_sets = ["STRICT", "RELAXED"]

for model in models:
    for ds in datasets:
        for ft in feature_sets:

            file_path = RESULTS_DIR / model / ds / ft / f"{model}_{ds}_{ft}_metrics.csv"

            if file_path.exists():
                df = pd.read_csv(file_path)
                all_metrics.append(df)

metrics_df = pd.concat(all_metrics, ignore_index=True)

# ============================================================
# BARPLOT 1 → DATASET COMPARISON
# (Average Test AUC across all models)
# ============================================================

dataset_summary = (
    metrics_df
    .groupby("Dataset")[["Test_AUC", "CV_AUC_Mean"]]
    .mean()
    .reset_index()
)

plt.figure(figsize=(8,6))

sns.barplot(
    data=dataset_summary,
    x="Dataset",
    y="Test_AUC"
)

plt.title("Average Test AUC Across Models")
plt.ylabel("Mean Test AUC")
plt.xlabel("Dataset")
plt.tight_layout()

plt.savefig(BASE_DIR / "Dataset_Comparison_TestAUC.png", dpi=300)
plt.close()

print("\n✔ Barplot 1 saved: Dataset_Comparison_TestAUC.png")

# ============================================================
# FIND BEST DATASET
# ============================================================

best_dataset = dataset_summary.sort_values(
    by="Test_AUC", ascending=False
).iloc[0]["Dataset"]

print(f"\n🏆 Best Performing Dataset: {best_dataset.upper()}")

# ============================================================
# BARPLOT 2 → MODEL COMPARISON (BEST DATASET)
# ============================================================

best_df = metrics_df[metrics_df["Dataset"] == best_dataset]

plt.figure(figsize=(10,6))

sns.barplot(
    data=best_df,
    x="Model",
    y="Test_AUC",
    hue="Feature_Set"
)

plt.title(f"Model Comparison on {best_dataset.upper()} Dataset")
plt.ylabel("Test AUC")
plt.xlabel("Model")
plt.legend(title="Feature Set")

plt.tight_layout()
plt.savefig(BASE_DIR / f"{best_dataset}_Model_Comparison_TestAUC.png", dpi=300)
plt.close()

print(f"✔ Barplot 2 saved: {best_dataset}_Model_Comparison_TestAUC.png")

print("\n🎉 Model comparison completed successfully.")
