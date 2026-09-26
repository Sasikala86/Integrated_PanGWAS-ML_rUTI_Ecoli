#!/usr/bin/env python3

# ============================================================
# PUBLICATION-READY PLASMIDOME ANALYSIS PIPELINE
# ============================================================
# Author: Sasikala
#
# PURPOSE:
#   - Extract plasmidome features from MOB-suite outputs
#   - Assign rUTI / UTI phenotype automatically
#   - Generate publication-quality figures
#   - Perform Fisher Exact Test + FDR correction
#   - Show ONLY significant labels
#   - Save 600 dpi PNG + JPEG outputs
#
# INPUT:
#   MOB-suite output directories containing:
#       mobtyper_results.txt
#       mge.report.txt (optional)
#
# OUTPUT:
#   plasmid_features_summary.csv
#   complete_plasmidome_summary.csv
#   feature-specific statistics CSV files
#   publication-quality PNG/JPEG figures
# ============================================================


import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


warnings.filterwarnings("ignore")


# ============================================================
# GLOBAL PLOT SETTINGS
# ============================================================

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titlesize"] = 18
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["legend.fontsize"] = 12
plt.rcParams["xtick.labelsize"] = 11
plt.rcParams["ytick.labelsize"] = 11
plt.rcParams["figure.dpi"] = 600


# ============================================================
# PATHS
# ============================================================

MOB_DIR = "/home/sastra/Sasikala/Pangenome/mob_output"

BASE_OUT = "/home/sastra/Sasikala/Pangenome/Final/plasmidome"

PLOT_DIR = os.path.join(BASE_OUT, "plots")
TABLE_DIR = os.path.join(BASE_OUT, "tables")
SUMMARY_DIR = os.path.join(BASE_OUT, "summary")

os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
os.makedirs(SUMMARY_DIR, exist_ok=True)

SUMMARY_CSV = os.path.join(
    SUMMARY_DIR,
    "plasmid_features_summary.csv"
)


# ============================================================
# COLORS
# ============================================================

COLOR_RUTI = "#ff7f0e"
COLOR_UTI = "#1f77b4"

PALETTE = {
    "rUTI": COLOR_RUTI,
    "UTI": COLOR_UTI
}


# ============================================================
# FEATURE EXTRACTION
# ============================================================

records = []


for sample_folder in sorted(os.listdir(MOB_DIR)):

    sample_path = os.path.join(MOB_DIR, sample_folder)

    if not os.path.isdir(sample_path):
        continue


    # --------------------------------------------------------
    # Assign phenotype from folder name
    # --------------------------------------------------------

    if sample_folder.startswith("rUTI"):
        phenotype = "rUTI"

    elif sample_folder.startswith("UTI"):
        phenotype = "UTI"

    else:
        phenotype = "unknown"


    mobtyper_file = os.path.join(
        sample_path,
        "mobtyper_results.txt"
    )

    if not os.path.exists(mobtyper_file):
        continue


    # --------------------------------------------------------
    # Read MOB-typer results
    # --------------------------------------------------------

    try:
        mob_df = pd.read_csv(
            mobtyper_file,
            sep="\t"
        )

    except Exception:
        continue


    # --------------------------------------------------------
    # Read MGE report if available
    # --------------------------------------------------------

    mge_file = os.path.join(
        sample_path,
        "mge.report.txt"
    )

    is_count_map = {}


    if os.path.exists(mge_file):

        try:

            mge_df = pd.read_csv(
                mge_file,
                sep="\t"
            )


            if "primary_cluster_id" in mge_df.columns:

                is_count_map = (
                    mge_df
                    .groupby("primary_cluster_id")
                    .size()
                    .to_dict()
                )

        except Exception:
            pass


    # --------------------------------------------------------
    # Extract plasmid-level features
    # --------------------------------------------------------

    for _, row in mob_df.iterrows():

        plasmid_id = row.get(
            "sample_id",
            ""
        )

        primary_cluster_id = row.get(
            "primary_cluster_id",
            ""
        )

        num_contigs = row.get(
            "num_contigs",
            0
        )

        circularity = row.get(
            "circularity_status",
            ""
        )

        rep_type = row.get(
            "rep_type(s)",
            ""
        )

        predicted_mobility = row.get(
            "predicted_mobility",
            ""
        )

        relaxase_type = row.get(
            "relaxase_type(s)",
            ""
        )

        mpf_type = row.get(
            "mpf_type",
            ""
        )

        orit_type = row.get(
            "orit_type(s)",
            ""
        )

        IS_count = is_count_map.get(
            primary_cluster_id,
            0
        )


        records.append({

            "sample_folder": sample_folder,

            "phenotype": phenotype,

            "plasmid_id": plasmid_id,

            "primary_cluster_id": primary_cluster_id,

            "num_contigs": num_contigs,

            "circularity": circularity,

            "rep_type": rep_type,

            "predicted_mobility": predicted_mobility,

            "relaxase_type": relaxase_type,

            "mpf_type": mpf_type,

            "orit_type": orit_type,

            "IS_count": IS_count

        })


# ============================================================
# CREATE PLASMID FEATURE SUMMARY
# ============================================================

df = pd.DataFrame(records)


if df.empty:

    raise RuntimeError(
        "No plasmid records were extracted from MOB-suite outputs."
    )


df["IS_count"] = pd.to_numeric(
    df["IS_count"],
    errors="coerce"
).fillna(0)


df["num_contigs"] = pd.to_numeric(
    df["num_contigs"],
    errors="coerce"
).fillna(0)


df.to_csv(
    SUMMARY_CSV,
    index=False
)


print(
    f"Plasmid feature summary saved to:\n{SUMMARY_CSV}"
)


# ============================================================
# EXPLODE COMMA-SEPARATED FEATURES
# ============================================================

def explode_column(df_input, column_name):

    temp = df_input.copy()

    temp[column_name] = (
        temp[column_name]
        .fillna("")
        .astype(str)
        .str.split(",")
    )

    temp = temp.explode(
        column_name
    )

    temp[column_name] = (
        temp[column_name]
        .astype(str)
        .str.strip()
    )

    temp = temp[
        (temp[column_name] != "") &
        (temp[column_name] != "-") &
        (temp[column_name] != "nan")
    ]

    return temp


df_rep = explode_column(
    df,
    "rep_type"
)

df_relax = explode_column(
    df,
    "relaxase_type"
)

df_mpf = explode_column(
    df,
    "mpf_type"
)


# ============================================================
# FISHER EXACT TEST + FDR
# ============================================================

def compute_fisher_statistics(
    df_full,
    df_exploded,
    feature_col
):

    results = []

    total_ruti = df_full[
        df_full["phenotype"] == "rUTI"
    ].shape[0]

    total_uti = df_full[
        df_full["phenotype"] == "UTI"
    ].shape[0]


    for feature in sorted(
        df_exploded[feature_col].unique()
    ):

        ruti_present = df_exploded[
            (df_exploded["phenotype"] == "rUTI") &
            (df_exploded[feature_col] == feature)
        ].shape[0]


        uti_present = df_exploded[
            (df_exploded["phenotype"] == "UTI") &
            (df_exploded[feature_col] == feature)
        ].shape[0]


        contingency = [

            [
                ruti_present,
                total_ruti - ruti_present
            ],

            [
                uti_present,
                total_uti - uti_present
            ]

        ]


        odds_ratio, pvalue = fisher_exact(
            contingency
        )


        results.append({

            "feature": feature,

            "rUTI_count": ruti_present,

            "UTI_count": uti_present,

            "odds_ratio": odds_ratio,

            "pvalue": pvalue

        })


    stats_df = pd.DataFrame(
        results
    )


    stats_df["FDR"] = multipletests(
        stats_df["pvalue"],
        method="fdr_bh"
    )[1]


    return stats_df


# ============================================================
# SIGNIFICANCE LABEL
# ============================================================

def significance_label(p):

    if p < 0.001:
        return "***"

    elif p < 0.01:
        return "**"

    elif p < 0.05:
        return "*"

    else:
        return ""


# ============================================================
# GENERIC BARPLOT FUNCTION
# ============================================================

def make_barplot(
    data,
    feature_col,
    title,
    outfile_prefix,
    top_n=None
):

    plot_data = data.copy()


    # --------------------------------------------------------
    # Count observations
    # --------------------------------------------------------

    count_df = (
        plot_data
        .groupby(
            [feature_col, "phenotype"]
        )
        .size()
        .reset_index(
            name="count"
        )
    )


    # --------------------------------------------------------
    # Feature ordering
    # --------------------------------------------------------

    feature_order = (
        count_df
        .groupby(feature_col)["count"]
        .sum()
        .sort_values(
            ascending=False
        )
        .index
        .tolist()
    )


    if top_n is not None:

        feature_order = feature_order[
            :top_n
        ]

        count_df = count_df[
            count_df[feature_col].isin(
                feature_order
            )
        ]


    # --------------------------------------------------------
    # Fisher statistics
    # --------------------------------------------------------

    stats_df = compute_fisher_statistics(
        data,
        data,
        feature_col
    )


    stats_file = os.path.join(
        TABLE_DIR,
        f"{outfile_prefix}_statistics.csv"
    )

    stats_df.to_csv(
        stats_file,
        index=False
    )


    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plt.figure(
        figsize=(12, 7)
    )


    ax = sns.barplot(
        data=count_df,
        x=feature_col,
        y="count",
        hue="phenotype",
        order=feature_order,
        palette=PALETTE
    )


    ax.set_title(
        title,
        fontweight="bold"
    )

    ax.set_xlabel(
        feature_col,
        fontweight="bold"
    )

    ax.set_ylabel(
        "Count",
        fontweight="bold"
    )


    plt.xticks(
        rotation=45,
        ha="right"
    )


    # --------------------------------------------------------
    # Significant labels based on FDR
    # --------------------------------------------------------

    significant_features = stats_df[
        stats_df["FDR"] < 0.05
    ]


    ymax = (
        count_df["count"].max()
        if not count_df.empty
        else 1
    )


    for _, row in significant_features.iterrows():

        feature = row["feature"]

        if feature not in feature_order:
            continue

        xpos = feature_order.index(
            feature
        )

        feature_max = count_df[
            count_df[feature_col] == feature
        ]["count"].max()


        label = significance_label(
            row["FDR"]
        )


        if label:

            ax.text(
                xpos,
                feature_max + ymax * 0.08,
                label,
                ha="center",
                va="bottom",
                fontweight="bold"
            )


    ax.set_ylim(
        0,
        ymax * 1.30
    )


    plt.tight_layout()


    # --------------------------------------------------------
    # Save 600 dpi PNG
    # --------------------------------------------------------

    png_file = os.path.join(
        PLOT_DIR,
        f"{outfile_prefix}.png"
    )


    plt.savefig(
        png_file,
        dpi=600,
        bbox_inches="tight"
    )


    # --------------------------------------------------------
    # Save 600 dpi JPEG
    # --------------------------------------------------------

    jpeg_file = os.path.join(
        PLOT_DIR,
        f"{outfile_prefix}.jpeg"
    )


    plt.savefig(
        jpeg_file,
        dpi=600,
        bbox_inches="tight"
    )


    plt.close()


# ============================================================
# GENERATE FEATURE DISTRIBUTION PLOTS
# ============================================================

make_barplot(
    df_rep,
    "rep_type",
    "Replicon Type Distribution",
    "replicon_distribution",
    top_n=15
)


make_barplot(
    df_relax,
    "relaxase_type",
    "Relaxase Type Distribution",
    "relaxase_distribution",
    top_n=10
)


make_barplot(
    df_mpf,
    "mpf_type",
    "MPF Type Distribution",
    "mpf_distribution",
    top_n=10
)


make_barplot(
    df,
    "predicted_mobility",
    "Predicted Mobility Distribution",
    "mobility_distribution"
)


# ============================================================
# IS COUNT BOXPLOT
# ============================================================

plt.figure(
    figsize=(8, 6)
)


sns.boxplot(
    data=df,
    x="phenotype",
    y="IS_count",
    hue="phenotype",
    palette=PALETTE,
    legend=False
)


plt.title(
    "Insertion Sequence Count",
    fontweight="bold"
)

plt.xlabel(
    "Phenotype",
    fontweight="bold"
)

plt.ylabel(
    "IS Count",
    fontweight="bold"
)


plt.tight_layout()


plt.savefig(
    os.path.join(
        PLOT_DIR,
        "IS_count_boxplot.png"
    ),
    dpi=600,
    bbox_inches="tight"
)


plt.savefig(
    os.path.join(
        PLOT_DIR,
        "IS_count_boxplot.jpeg"
    ),
    dpi=600,
    bbox_inches="tight"
)


plt.close()


# ============================================================
# NUMBER OF CONTIGS BOXPLOT
# ============================================================

plt.figure(
    figsize=(8, 6)
)


sns.boxplot(
    data=df,
    x="phenotype",
    y="num_contigs",
    hue="phenotype",
    palette=PALETTE,
    legend=False
)


plt.title(
    "Number of Contigs per Plasmid",
    fontweight="bold"
)

plt.xlabel(
    "Phenotype",
    fontweight="bold"
)

plt.ylabel(
    "Number of Contigs",
    fontweight="bold"
)


plt.tight_layout()


plt.savefig(
    os.path.join(
        PLOT_DIR,
        "num_contigs_boxplot.png"
    ),
    dpi=600,
    bbox_inches="tight"
)


plt.savefig(
    os.path.join(
        PLOT_DIR,
        "num_contigs_boxplot.jpeg"
    ),
    dpi=600,
    bbox_inches="tight"
)


plt.close()


# ============================================================
# COMPLETE PLASMIDOME SUMMARY
# ============================================================

complete_summary_file = os.path.join(
    TABLE_DIR,
    "complete_plasmidome_summary.csv"
)


df.to_csv(
    complete_summary_file,
    index=False
)


print("\n" + "=" * 70)
print("PLASMID-LEVEL PLASMIDOME ANALYSIS COMPLETED")
print("=" * 70)

print(
    f"Summary:\n{SUMMARY_CSV}"
)

print(
    f"Complete table:\n{complete_summary_file}"
)

print(
    f"Plots:\n{PLOT_DIR}"
)

print(
    f"Statistics:\n{TABLE_DIR}"
)

print("=" * 70)
