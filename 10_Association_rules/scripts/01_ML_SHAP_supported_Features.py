import pandas as pd

# ---------------------------
# 1. Load your CSV files
# ---------------------------
shap_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/SHAP/XGB/combined/STRICT/SHAP_mean_abs_importance_XGB_combined_STRICT.csv"
freq_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/SHAP/XGB/combined/STRICT/rUTI_vs_UTI_proportions_STRICT.csv"

shap_df = pd.read_csv(shap_file)
freq_df = pd.read_csv(freq_file)

# Rename 'Feature' column to 'Gene' for merging
shap_df = shap_df.rename(columns={"Feature": "Gene"})

# ---------------------------
# 2. Merge SHAP + frequency
# ---------------------------
merged_df = pd.merge(freq_df, shap_df, on="Gene", how="inner")

# ---------------------------
# 3. Set filtering thresholds
# ---------------------------
threshold_diff = 0.15
threshold_shap = 0.01

# ---------------------------
# 4. Filter genes
# ---------------------------
filtered_df = merged_df[
    (merged_df['Difference_rUTI_minus_UTI'] > threshold_diff) &
    (merged_df['Mean_Abs_SHAP'] > threshold_shap)
]

# Sort by Difference then SHAP importance
filtered_df = filtered_df.sort_values(
    by=['Difference_rUTI_minus_UTI', 'Mean_Abs_SHAP'],
    ascending=False
)

# ---------------------------
# 5. Save filtered genes
# ---------------------------
output_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/SHAP/XGB/combined/STRICT/filtered_rUTI_genes_STRICT.csv"

filtered_df.to_csv(output_file, index=False)

print(f"Filtered genes saved to: {output_file}")

print("Top filtered genes:")
print(
    filtered_df[
        ['Gene', 'Difference_rUTI_minus_UTI', 'Mean_Abs_SHAP']
    ].head(20)
)

# ---------------------------
# Load filtered genes
# ---------------------------
filtered_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/SHAP/XGB/combined/STRICT/filtered_rUTI_genes_STRICT.csv"

filtered_df = pd.read_csv(filtered_file)

# List of top SHAP features
shap_top_genes = [
    "tufB~~~tuf1",
    "gadB_1~~~gadB~~~gadA",
    "group_10659",
    "group_6673",
    "group_9085",
    "group_8601",
    "group_5784",
    "klcA_1~~~klcA~~~klcA_3~~~klcA_2~~~klcA_5~~~klcA_6",
    "group_4771",
    "group_6498",
    "group_852",
    "group_11659",
    "mazF~~~mazF_1",
    "group_8676",
    "group_6803",
    "flu_1~~~flu_2~~~flu~~~flu_3~~~flu_4",
    "btuB_3~~~cirA_1",
    "group_10205",
    "klcA_2~~~klcA_3~~~klcA_1~~~klcA~~~klcA_4",
    "group_5512",
    "group_5473",
    "group_9306",
    "group_487",
    "group_11078",
    "ldrD_2",
    "cbtA~~~cbtA_1~~~cbtA_2~~~cbtA_3",
    "group_10611",
    "group_11500",
    "group_10305",
    "group_4412",
    "papA~~~papA_2~~~papA_1",
    "klcA_5~~~klcA",
    "group_10976",
    "gadA~~~gadB_3~~~gadB_2~~~gadB",
    "group_7498",
    "group_3974",
    "group_12142",
    "cbeA_1~~~cbeA_2~~~cbeA_3~~~cbeA_4~~~cbeA~~~cbeA_5",
    "ldrD_1~~~ldrD_2~~~ldrD_3~~~ldrD",
    "group_10462",
    "group_6554",
    "higB2_1",
    "group_8337",
    "group_12144",
    "group_10489",
    "group_9905",
    "group_12145",
    "group_1721",
    "espC",
    "metG_1~~~metG_2~~~metG",
    "group_9743",
    "agaC_2",
    "glmU_3~~~glmU_1",
    "group_867",
    "group_9621",
    "group_7827",
]

# ---------------------------
# Find common genes
# ---------------------------
filtered_genes = set(filtered_df['Gene'])

common_genes = filtered_genes.intersection(shap_top_genes)

print(f"Number of common genes: {len(common_genes)}")

print("Common genes:")

for g in sorted(common_genes):
    print(g)
