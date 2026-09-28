#!/usr/bin/env python3
# SHAP Analysis – XGB
# STRICT + COMBINED
# (Version-Proof Implementation)
# ============================================================

import os
import joblib
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

# ============================================================
# PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[2]

MODEL_PATH = f"{BASE_DIR}/results/Models/XGB/combined/STRICT/XGB_combined_STRICT_model.pkl"
FEATURE_LIST_PATH = f"{BASE_DIR}/processed/feature_selection/results/combined_FINAL_STRICT_features.csv"
FULL_MATRIX_PATH = f"{BASE_DIR}/processed/combined_ML_filtered_1-99.csv"
OUTPUT_DIR = f"{BASE_DIR}/results/SHAP/XGB/combined/STRICT"

TARGET_COLUMN = "phenotype"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# LOAD MODEL
# ============================================================

model = joblib.load(MODEL_PATH)
print("✅ Model loaded")

# ============================================================
# LOAD DATA
# ============================================================

full_df = pd.read_csv(FULL_MATRIX_PATH)
feature_list = pd.read_csv(FEATURE_LIST_PATH)["gene"].tolist()

X = full_df[feature_list]
y = full_df[TARGET_COLUMN].astype(int)

# ============================================================
# SAME TRAIN/TEST SPLIT
# ============================================================

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

print("Test set shape:", X_test.shape)

# ============================================================
# SHAP EXPLAINER (ROBUST METHOD)
# ============================================================

explainer = shap.Explainer(model.predict_proba, X_train)

shap_values = explainer(X_test)

# For binary classification → class 1
shap_array = shap_values.values[:, :, 1]

# ============================================================
# SAVE RAW SHAP VALUES
# ============================================================

shap_df = pd.DataFrame(shap_array, columns=X_test.columns)
shap_df.to_csv(
    f"{OUTPUT_DIR}/SHAP_values_XGB_combined_STRICT.csv",
    index=False
)

# ============================================================
# BEESWARM PLOT
# ============================================================

plt.figure()
shap.plots.beeswarm(shap_values[:, :, 1], show=False)
plt.savefig(
    f"{OUTPUT_DIR}/SHAP_summary_beeswarm_XGB_combined_STRICT.png",
    dpi=600,
    bbox_inches="tight"
)
plt.close()

# ============================================================
# BAR PLOT
# ============================================================

plt.figure()
shap.plots.bar(shap_values[:, :, 1], show=False)
plt.savefig(
    f"{OUTPUT_DIR}/SHAP_importance_bar_XGB_combined_STRICT.png",
    dpi=600,
    bbox_inches="tight"
)
plt.close()

# ============================================================
# MEAN ABS SHAP
# ============================================================

mean_abs_shap = np.abs(shap_array).mean(axis=0)

importance_df = pd.DataFrame({
    "Feature": X_test.columns,
    "Mean_Abs_SHAP": mean_abs_shap
}).sort_values(by="Mean_Abs_SHAP", ascending=False)

importance_df.to_csv(
    f"{OUTPUT_DIR}/SHAP_mean_abs_importance_XGB_combined_STRICT.csv",
    index=False
)

print("\n✅ SHAP completed successfully")
print("📂 Results saved in:", OUTPUT_DIR)
# ============================================================
