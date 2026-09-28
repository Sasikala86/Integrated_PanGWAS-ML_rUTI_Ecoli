# New script for pyseer

import pandas as pd

# ==============================
# 1. Load Pyseer results
# ==============================
file_path = "/home/sastra/Sasikala/Pangenome/PanGWAS_Accessory/pyseer_accessory_results.txt"

df = pd.read_csv(file_path, sep="\t")

# ==============================
# 2. Clean data (remove bad model fits)
# ==============================
df_clean = df[
    df['notes'].isna() | (df['notes'] != 'bad-chisq')
].copy()

# ==============================
# 3. Convert numeric columns (safety)
# ==============================
numeric_cols = ["lrt-pvalue", "beta"]

for col in numeric_cols:
    df_clean[col] = pd.to_numeric(
        df_clean[col],
        errors='coerce'
    )

# Drop NA values after conversion
df_clean = df_clean.dropna(
    subset=numeric_cols
)

# ==============================
# 4. Statistical significance (Bonferroni)
# ==============================
bonferroni_threshold = 0.05 / len(df_clean)

df_sig = df_clean[
    df_clean["lrt-pvalue"] < bonferroni_threshold
].copy()

# ==============================
# 5. Filter rUTI-specific genes (IMPORTANT)
# ==============================
df_rUTI = df_sig[
    df_sig["beta"] > 1
].copy()

# ==============================
# 6. (Optional) Filter UTI-specific genes
# ==============================
df_UTI = df_sig[
    df_sig["beta"] < -1
].copy()

# ==============================
# 7. Sort by strongest effect
# ==============================
df_rUTI = df_rUTI.sort_values(
    by="beta",
    ascending=False
)

df_UTI = df_UTI.sort_values(
    by="beta"
)

# ==============================
# 8. Save outputs
# ==============================
rUTI_output = "/home/sastra/Sasikala/Pangenome/PanGWAS_Accessory/pyseer_rUTI_genes.csv"

UTI_output = "/home/sastra/Sasikala/Pangenome/PanGWAS_Accessory/pyseer_UTI_genes.csv"

df_rUTI.to_csv(
    rUTI_output,
    index=False
)

df_UTI.to_csv(
    UTI_output,
    index=False
)

# ==============================
# 9. Summary
# ==============================
print("========== PYSEER FILTER SUMMARY ==========")

print(f"Total variants: {len(df)}")
print(f"After cleaning: {len(df_clean)}")
print(f"Significant variants: {len(df_sig)}")
print(f"rUTI-specific genes: {len(df_rUTI)}")
print(f"UTI-specific genes: {len(df_UTI)}")

print("\nTop rUTI genes:")

display(df_rUTI.head())
