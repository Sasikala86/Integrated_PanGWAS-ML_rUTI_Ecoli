# ==========================================
# SEEDED APRIORI + rUTI vs UTI ENRICHMENT
# ==========================================

import pandas as pd
from mlxtend.frequent_patterns import apriori, association_rules

# ==============================
# FILE PATHS
# ==============================
rtab_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/gene_presence_absence_combined.Rtab"

gene_list_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/genes_present_in_atleast_2_methods.csv"

output_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/rUTI_enriched_rules.csv"

# ==============================
# LOAD DATA
# ==============================
df = pd.read_csv(
    rtab_file,
    sep="\t"
)

df = df.set_index(
    df.columns[0]
).T

# ==============================
# ADD PHENOTYPE COLUMNS
# ==============================
df["rUTI"] = (
    df.index.str.contains("rUTI")
    .astype(int)
)

df["UTI"] = (
    df.index.str.contains("UTI")
    .astype(int)
)

# IMPORTANT: avoid overlap
df["UTI"] = (
    (df["UTI"] == 1) &
    (df["rUTI"] == 0)
).astype(int)

print(
    df[["rUTI", "UTI"]].sum()
)

# ==============================
# LOAD SELECTED GENES
# ==============================
genes = pd.read_csv(
    gene_list_file
)["Gene"].tolist()

# ==============================
# MATCH GENES
# ==============================
def match_genes(genes, columns):

    matched = []

    for g in genes:

        for col in columns:

            if (
                col == g
                or
                col.startswith(g + "_")
            ):

                matched.append(col)

    return list(set(matched))

matched_genes = match_genes(
    genes,
    df.columns
)

print(
    "Matched genes:",
    len(matched_genes)
)

# ==============================
# SUBSET MATRIX
# ==============================
df_small = df[
    matched_genes + ["rUTI", "UTI"]
].astype(bool)

# ==============================
# APRIORI
# ==============================
freq_items = apriori(
    df_small,
    min_support=0.1,
    max_len=3,
    use_colnames=True
)

# ==============================
# RULES
# ==============================
rules = association_rules(
    freq_items,
    metric="confidence",
    min_threshold=0.6
)

rules["antecedents"] = rules[
    "antecedents"
].apply(list)

rules["consequents"] = rules[
    "consequents"
].apply(list)

# ==============================
# KEEP ONLY rUTI RULES
# ==============================
rules_rUTI = rules[
    rules["consequents"].apply(
        lambda x: "rUTI" in x
    )
].copy()

print(
    "Total rUTI rules:",
    len(rules_rUTI)
)

# ==============================
# CALCULATE PROPORTIONS
# ==============================
def rule_presence(row, df):

    genes = row["antecedents"]

    subset = df.copy()

    for g in genes:

        subset = subset[
            subset[g] == 1
        ]

    if len(subset) == 0:

        return 0, 0

    ruti_prop = subset["rUTI"].mean()

    uti_prop = subset["UTI"].mean()

    return (
        ruti_prop,
        uti_prop
    )

# Apply to all rules
rules_rUTI[
    ["rUTI_prop", "UTI_prop"]
] = rules_rUTI.apply(
    lambda row: pd.Series(
        rule_presence(row, df)
    ),
    axis=1
)

# ==============================
# FILTER: ENRICHED IN rUTI
# ==============================
final_rules = rules_rUTI[
    (rules_rUTI["rUTI_prop"] > 0.7) &
    (rules_rUTI["UTI_prop"] < 0.3)
]

# ==============================
# SORT
# ==============================
final_rules = final_rules.sort_values(
    by=[
        "lift",
        "confidence",
        "rUTI_prop"
    ],
    ascending=False
)

# ==============================
# SAVE
# ==============================
final_rules.to_csv(
    output_file,
    index=False
)

print(
    "Final enriched rules:",
    len(final_rules)
)

print(
    "Saved to:",
    output_file
)
