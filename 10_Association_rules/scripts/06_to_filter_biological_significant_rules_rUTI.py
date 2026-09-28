# ==========================================
# FINAL SCRIPT: FILTER BIOLOGICALLY SIGNIFICANT rUTI RULES
# ==========================================

import pandas as pd
import re
import os
import ast

# ==============================
# FILE PATHS
# ==============================
rules_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/rUTI_enriched_rules.csv"
output_dir = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/"
os.makedirs(output_dir, exist_ok=True)

# ==============================
# LOAD RULES
# ==============================
df = pd.read_csv(rules_file)

# ==============================
# CONVERT STRING → LIST
# ==============================
df["antecedents"] = df["antecedents"].apply(ast.literal_eval)
df["consequents"] = df["consequents"].apply(ast.literal_eval)

# ==============================
# FUNCTION: CLEAN GENE NAMES
# (collapse variants like higB2_1 → higB2)
# ==============================
def clean_gene(g):
    if g.startswith("group_"):
        return g
    return re.sub(r'_\d+$', '', g)

def clean_list(gene_list):
    return list(set([clean_gene(g) for g in gene_list if g != "rUTI" and g != "UTI"]))

df["antecedents_clean"] = df["antecedents"].apply(clean_list)
df["consequents_clean"] = df["consequents"].apply(clean_list)

# ==============================
# REMOVE TRIVIAL RULES
# (same gene in both sides)
# ==============================
df = df[
    df.apply(
        lambda row: len(set(row["antecedents_clean"]) & set(row["consequents_clean"])) == 0,
        axis=1
    )
]

# ==============================
# FILTER 1: rUTI in consequent
# ==============================
df = df[
    df["consequents"].apply(lambda x: "rUTI" in x)
]

# ==============================
# FILTER 2: STRONG RULES
# ==============================
df = df[
    (df["confidence"] >= 0.8) &
    (df["lift"] >= 1.5) &
    (df["support"] >= 0.05)
]

# ==============================
# FILTER 3: rUTI ENRICHMENT
# ==============================
df = df[
    (df["rUTI_prop"] >= 0.7) &
    (df["UTI_prop"] <= 0.3)
]

# ==============================
# FILTER 4: REMOVE WEAK SINGLE GENE RULES
# ==============================
df = df[
    df["antecedents_clean"].apply(lambda x: len(x) >= 2)
]

# ==============================
# SORT BEST RULES
# ==============================
df = df.sort_values(
    by=["lift", "confidence"],
    ascending=False
)

# ==============================
# SAVE FINAL RESULTS
# ==============================
df.to_csv(
    output_dir + "FINAL_biologically_significant_rUTI_rules.csv",
    index=False
)

# ==============================
# OPTIONAL: TOP RULES ONLY
# ==============================
top_df = df.head(30)

top_df.to_csv(
    output_dir + "TOP_30_rUTI_rules.csv",
    index=False
)

# ==============================
# SUMMARY
# ==============================
print("\n===== FINAL FILTER SUMMARY =====")
print("Total filtered rules:", len(df))

print("\nTop rules preview:")
print(df[["antecedents_clean", "confidence", "lift", "rUTI_prop"]].head(10))
