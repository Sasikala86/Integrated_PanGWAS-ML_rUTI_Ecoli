import pandas as pd

# ==============================
# 1. Load Scoary results
# ==============================
file_path = "/home/sastra/Sasikala/Pangenome/PanGWAS_Accessory/scoary_results/rUTI_28_02_2026_1639.results.csv"

df = pd.read_csv(file_path)

# ==============================
# 2. Convert numeric columns (safety)
# ==============================
numeric_cols = [
    "Sensitivity",
    "Specificity",
    "Odds_ratio",
    "Benjamini_H_p"
]

for col in numeric_cols:
    df[col] = pd.to_numeric(
        df[col],
        errors='coerce'
    )

# ==============================
# 3. Filter rUTI-enriched genes
# ==============================
filtered = df[
    (df["Benjamini_H_p"] <= 0.05) &
    (df["Odds_ratio"] > 1) &
    (df["Sensitivity"] >= 20) &
    (df["Specificity"] >= 60)
].copy()

# ==============================
# 4. Sort (most important first)
# ==============================
filtered = filtered.sort_values(
    by=["Benjamini_H_p", "Odds_ratio"],
    ascending=[True, False]
)

# ==============================
# 5. Save output
# ==============================
output_file = "/home/sastra/Sasikala/Pangenome/PanGWAS_Accessory/scoary_results/rUTI_filtered_genes.csv"

filtered.to_csv(
    output_file,
    index=False
)

# ==============================
# 6. Summary
# ==============================
print(f"Total genes in Scoary: {len(df)}")

print(
    f"Filtered rUTI-associated genes: {len(filtered)}"
)

print("\nTop genes:")

display(filtered.head())
