#!/usr/bin/env python3
# Calculate rUTI vs UTI Gene Proportions
# STRICT (56 genes)
# ============================================================

import pandas as pd

# ============================================================
# PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[2]

FULL_MATRIX_PATH = f"{BASE_DIR}/processed/combined_ML_filtered_1-99.csv"
FEATURE_LIST_PATH = f"{BASE_DIR}/processed/feature_selection/results/combined_FINAL_STRICT_features.csv"

TARGET_COLUMN = "phenotype"

# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(FULL_MATRIX_PATH)
feature_list = pd.read_csv(FEATURE_LIST_PATH)["gene"].tolist()

# Keep only phenotype + strict genes
df = df[[TARGET_COLUMN] + feature_list]

# ============================================================
# SPLIT GROUPS
# ============================================================

rUTI_df = df[df[TARGET_COLUMN] == 1]
UTI_df  = df[df[TARGET_COLUMN] == 0]

total_rUTI = len(rUTI_df)
total_UTI  = len(UTI_df)

print("Total rUTI:", total_rUTI)
print("Total UTI :", total_UTI)

# ============================================================
# CALCULATE PROPORTIONS
# ============================================================

results = []

for gene in feature_list:
    
    rUTI_present = rUTI_df[gene].sum()
    UTI_present  = UTI_df[gene].sum()
    
    rUTI_prop = rUTI_present / total_rUTI
    UTI_prop  = UTI_present / total_UTI
    
    difference = rUTI_prop - UTI_prop
    
    results.append([
        gene,
        rUTI_present,
        UTI_present,
        rUTI_prop,
        UTI_prop,
        difference
    ])

# ============================================================
# CREATE RESULT TABLE
# ============================================================

prop_df = pd.DataFrame(results, columns=[
    "Gene",
    "rUTI_count",
    "UTI_count",
    "rUTI_proportion",
    "UTI_proportion",
    "Difference_rUTI_minus_UTI"
])

prop_df = prop_df.sort_values(
    by="Difference_rUTI_minus_UTI",
    ascending=False
)

# ============================================================
# SAVE
# ============================================================

OUTPUT_PATH = f"{BASE_DIR}/results/SHAP/XGB/combined/STRICT/rUTI_vs_UTI_proportions_STRICT.csv"

prop_df.to_csv(OUTPUT_PATH, index=False)

print("\n✅ Proportion analysis completed")
print("Saved to:", OUTPUT_PATH)
