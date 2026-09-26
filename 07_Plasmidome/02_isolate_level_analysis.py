#!/usr/bin/env python3

# ============================================================
# ISOLATE-LEVEL PLASMIDOME ANALYSIS PIPELINE
# ============================================================
# Author: Sasikala
#
# PURPOSE:
#   1. Isolate-level plasmidome analysis
#   2. Presence/absence analysis per isolate
#   3. Fisher exact significance testing
#   4. FDR correction
#   5. Publication-quality figures
#   6. 600 dpi PNG + JPEG outputs
#
# IMPORTANT:
#   Unlike plasmid-level analysis:
#   Each isolate contributes ONLY ONCE
#   for a feature regardless of plasmid count.
#
# INPUT:
#   plasmid_features_summary.csv
#
# OUTPUT:
#   replicon_presence_absence.csv
#   relaxase_presence_absence.csv
#   mpf_presence_absence.csv
#   isolate-level statistics CSV files
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
# INPUT
# ============================================================

SUMMARY_CSV = (
    "/home/sastra/Sasikala/Pangenome/Final/"
    "plasmidome/summary/plasmid_features_summary.csv"
)


# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

BASE_OUT = (
    "/home/sastra/Sasikala/Pangenome/Final/"
    "plasmidome/isolate_level"
)


PLOT_DIR = os.path.join(
    BASE_OUT,
    "plots"
)


TABLE_DIR = os.path.join(
    BASE_OUT,
    "tables"
)


os.makedirs(
    PLOT_DIR,
    exist_ok=True
)


os.makedirs(
    TABLE_DIR,
    exist_ok=True
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
# LOAD PLASMID FEATURE SUMMARY
# ============================================================

df = pd.read_csv(
    SUMMARY_CSV
)


# ------------------------------------------------------------
# Keep only rUTI and UTI
# ------------------------------------------------------------

df = df[
    df["phenotype"].isin(
        ["rUTI", "UTI"]
    )
].copy()


# ============================================================
# EXPLODE COMMA-SEPARATED FEATURES
# ============================================================

def explode_column(
    df_input,
    column_name
):

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


# ============================================================
# CREATE FEATURE-SPECIFIC DATASETS
# ============================================================

rep_presence = explode_column(
    df,
    "rep_type"
)


relax_presence = explode_column(
    df,
    "relaxase_type"
)


mpf_presence = explode_column(
    df,
    "mpf_type"
)


# ============================================================
# CREATE ISOLATE-LEVEL PRESENCE/ABSENCE
# ============================================================

def make_presence_absence(
    df_input,
    feature_col
):

    presence = (

        df_input

        .groupby(
            [
                "sample_folder",
                "phenotype",
                feature_col
            ]
        )

        .size()

        .reset_index(
            name="present"
        )

    )


    # --------------------------------------------------------
    # Convert count to binary presence
    # --------------------------------------------------------

    presence["present"] = 1


    return presence


# ============================================================
# GENERATE PRESENCE TABLES
# ============================================================

rep_presence = make_presence_absence(
    rep_presence,
    "rep_type"
)


relax_presence = make_presence_absence(
    relax_presence,
    "relaxase_type"
)


mpf_presence = make_presence_absence(
    mpf_presence,
    "mpf_type"
)


# ============================================================
# SAVE PRESENCE/ABSENCE TABLES
# ============================================================

rep_presence.to_csv(
    os.path.join(
        TABLE_DIR,
        "replicon_presence_absence.csv"
    ),
    index=False
)


relax_presence.to_csv(
    os.path.join(
        TABLE_DIR,
        "relaxase_presence_absence.csv"
    ),
    index=False
)


mpf_presence.to_csv(
    os.path.join(
        TABLE_DIR,
        "mpf_presence_absence.csv"
    ),
    index=False
)


# ============================================================
# FISHER EXACT TEST — ISOLATE LEVEL
# ============================================================

def compute_isolate_fisher_statistics(
    feature_df,
    feature_col
):

    results = []


    # --------------------------------------------------------
    # Total number of unique isolates
    # --------------------------------------------------------

    total_ruti = feature_df[
        feature_df["phenotype"] == "rUTI"
    ]["sample_folder"].nunique()


    total_uti = feature_df[
        feature_df["phenotype"] == "UTI"
    ]["sample_folder"].nunique()


    # --------------------------------------------------------
    # Test each feature
    # --------------------------------------------------------

    for feature in sorted(
        feature_df[feature_col].unique()
    ):


        ruti_present = feature_df[
            (feature_df["phenotype"] == "rUTI") &
            (feature_df[feature_col] == feature)
        ]["sample_folder"].nunique()


        uti_present = feature_df[
            (feature_df["phenotype"] == "UTI") &
            (feature_df[feature_col] == feature)
        ]["sample_folder"].nunique()


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

            "rUTI_isolates": ruti_present,

            "UTI_isolates": uti_present,

            "odds_ratio": odds_ratio,

            "pvalue": pvalue

        })


    stats_df = pd.DataFrame(
        results
    )


    # --------------------------------------------------------
    # Benjamini-Hochberg FDR correction
    # --------------------------------------------------------

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
# ISOLATE-LEVEL BARPLOT
# ============================================================

def make_barplot(
    feature_df,
    feature_col,
    title,
    outfile_prefix,
    top_n=None
):


    # --------------------------------------------------------
    # Count unique isolates
    # --------------------------------------------------------

    count_df = (

        feature_df

        .groupby(
            [
                feature_col,
                "phenotype"
            ]
        )

        ["sample_folder"]

        .nunique()

        .reset_index(
            name="count"
        )

    )


    # --------------------------------------------------------
    # Feature ordering
    # --------------------------------------------------------

    feature_order = (

        count_df

        .groupby(
            feature_col
        )["count"]

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

    stats_df = compute_isolate_fisher_statistics(
        feature_df,
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
        "Number of isolates",
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
    # Save PNG
    # --------------------------------------------------------

    plt.savefig(
        os.path.join(
            PLOT_DIR,
            f"{outfile_prefix}.png"
        ),
        dpi=600,
        bbox_inches="tight"
    )


    # --------------------------------------------------------
    # Save JPEG
    # --------------------------------------------------------

    plt.savefig(
        os.path.join(
            PLOT_DIR,
            f"{outfile_prefix}.jpeg"
        ),
        dpi=600,
        bbox_inches="tight"
    )


    plt.close()


# ============================================================
# GENERATE ISOLATE-LEVEL PLOTS
# ============================================================

make_barplot(
    rep_presence,
    "rep_type",
    "Replicon Distribution (Isolate Level)",
    "replicon_distribution_isolate",
    top_n=15
)


make_barplot(
    relax_presence,
    "relaxase_type",
    "Relaxase Distribution (Isolate Level)",
    "relaxase_distribution_isolate",
    top_n=10
)


make_barplot(
    mpf_presence,
    "mpf_type",
    "MPF Distribution (Isolate Level)",
    "mpf_distribution_isolate",
    top_n=10
)


# ============================================================
# COMPLETION MESSAGE
# ============================================================

print("\n" + "=" * 70)
print("ISOLATE-LEVEL PLASMIDOME ANALYSIS COMPLETED")
print("=" * 70)

print(
    f"Input summary:\n{SUMMARY_CSV}"
)

print(
    f"Presence/absence tables:\n{TABLE_DIR}"
)

print(
    f"Statistics:\n{TABLE_DIR}"
)

print(
    f"Plots:\n{PLOT_DIR}"
)

print("=" * 70)
