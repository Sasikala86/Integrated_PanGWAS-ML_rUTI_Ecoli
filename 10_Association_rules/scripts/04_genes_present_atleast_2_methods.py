# ==========================================
# FIXED SCRIPT: GENES ≥ 2 METHODS (APriori_Input)
# ==========================================

import pandas as pd
import re
import os

# ==============================
# FILE PATHS
# ==============================
ml_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/rUTI_enriched_genes_ML.csv"

pyseer_file = "/home/sastra/Sasikala/Pangenome/PanGWAS_Accessory/pyseer_rUTI_genes.csv"

scoary_file = "/home/sastra/Sasikala/Pangenome/PanGWAS_Accessory/scoary_results/rUTI_filtered_genes.csv"

output_dir = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/"

# ==============================
# ENSURE DIRECTORY EXISTS
# ==============================
os.makedirs(
    output_dir,
    exist_ok=True
)

print(
    f"\nSaving results to:\n{output_dir}"
)

# ==============================
# CLEAN FUNCTION
# ==============================
def clean_gene_list(series):

    genes = []

    for entry in series.dropna():

        parts = str(entry).split("~~~")

        for g in parts:

            g = g.strip()

            if g.startswith("group_"):

                genes.append(g)

            else:

                g = re.sub(
                    r'_\d+$',
                    '',
                    g
                )

                genes.append(g)

    return set(genes)

# ==============================
# LOAD FILES
# ==============================
ml_df = pd.read_csv(
    ml_file
)

pyseer_df = pd.read_csv(
    pyseer_file
)

scoary_df = pd.read_csv(
    scoary_file
)

print(
    "\nFiles loaded successfully!"
)

# ==============================
# DETECT COLUMNS
# ==============================
def find_gene_column(df):

    for col in df.columns:

        if (
            "gene" in col.lower()
            or
            "variant" in col.lower()
        ):

            return col

    return df.columns[0]

ml_col = find_gene_column(
    ml_df
)

pyseer_col = find_gene_column(
    pyseer_df
)

scoary_col = find_gene_column(
    scoary_df
)

print(
    "\nDetected columns:"
)

print(
    "ML:",
    ml_col
)

print(
    "Pyseer:",
    pyseer_col
)

print(
    "Scoary:",
    scoary_col
)

# ==============================
# CLEAN GENE SETS
# ==============================
ml_genes = clean_gene_list(
    ml_df[ml_col]
)

pyseer_genes = clean_gene_list(
    pyseer_df[pyseer_col]
)

scoary_genes = clean_gene_list(
    scoary_df[scoary_col]
)

# ==============================
# INTERSECTIONS
# ==============================
ml_pyseer = (
    ml_genes &
    pyseer_genes
)

ml_scoary = (
    ml_genes &
    scoary_genes
)

pyseer_scoary = (
    pyseer_genes &
    scoary_genes
)

genes_2plus = (
    ml_pyseer |
    ml_scoary |
    pyseer_scoary
)

# ==============================
# FILE PATHS FOR OUTPUT
# ==============================
file_main = os.path.join(
    output_dir,
    "genes_present_in_atleast_2_methods.csv"
)

file_ml_pyseer = os.path.join(
    output_dir,
    "ML_Pyseer_overlap.csv"
)

file_ml_scoary = os.path.join(
    output_dir,
    "ML_Scoary_overlap.csv"
)

file_pyseer_scoary = os.path.join(
    output_dir,
    "Pyseer_Scoary_overlap.csv"
)

# ==============================
# SAVE FILES
# ==============================
pd.DataFrame(
    {"Gene": sorted(genes_2plus)}
).to_csv(
    file_main,
    index=False
)

pd.DataFrame(
    {"Gene": sorted(ml_pyseer)}
).to_csv(
    file_ml_pyseer,
    index=False
)

pd.DataFrame(
    {"Gene": sorted(ml_scoary)}
).to_csv(
    file_ml_scoary,
    index=False
)

pd.DataFrame(
    {"Gene": sorted(pyseer_scoary)}
).to_csv(
    file_pyseer_scoary,
    index=False
)

# ==============================
# CONFIRMATION
# ==============================
print(
    "\n✅ FILES SAVED SUCCESSFULLY:"
)

print(
    file_main
)

print(
    file_ml_pyseer
)

print(
    file_ml_scoary
)

print(
    file_pyseer_scoary
)

# ==============================
# SUMMARY
# ==============================
print(
    "\n===== SUMMARY ====="
)

print(
    "ML genes:",
    len(ml_genes)
)

print(
    "Pyseer genes:",
    len(pyseer_genes)
)

print(
    "Scoary genes:",
    len(scoary_genes)
)

print(
    "\nFINAL genes (≥2 methods):",
    len(genes_2plus)
)
