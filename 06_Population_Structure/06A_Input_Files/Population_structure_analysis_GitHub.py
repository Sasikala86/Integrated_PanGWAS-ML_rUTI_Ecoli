# ==============================================================================
# PUBLICATION-QUALITY PCA + MDS ANALYSIS
# ==============================================================================
# INPUT:
# gene_presence_absence.Rtab
#
# OUTPUT:
# PCA_Plot_Publication.png
# PCA_Plot_Publication.pdf
# MDS_Jaccard_Publication.png
# MDS_Jaccard_Publication.pdf
#
# DESCRIPTION
# ------------------------------------------------------------------------------
# This script reproduces the R workflow:
#
# PCA:
#   prcomp()
#
# MDS:
#   dist(method="binary")
#   cmdscale()
#
# COLORS
# ------------------------------------------------------------------------------
# rUTI = Orange
# UTI  = Blue
#
# AUTHOR:
# Sasikala Project
# ==============================================================================


# ==============================================================================
# 1. IMPORT LIBRARIES
# ==============================================================================

import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from sklearn.decomposition import PCA

from scipy.stats import chi2
from scipy.spatial.distance import pdist, squareform


# ==============================================================================
# 2. GLOBAL PLOT SETTINGS
# ==============================================================================

plt.rcParams.update({

    "font.family": "DejaVu Sans",

    "font.weight": "bold",

    "axes.labelweight": "bold",

    "axes.titleweight": "bold",

    "axes.titlesize": 18,

    "axes.labelsize": 15,

    "xtick.labelsize": 12,

    "ytick.labelsize": 12,

    "legend.fontsize": 12,

    "figure.dpi": 600
})


# ==============================================================================
# REPOSITORY CONFIGURATION
# ==============================================================================
from pathlib import Path

REPO_ROOT = Path("../..").resolve()

input_file = REPO_ROOT / "04_Panaroo" / "results" / "combined" / "gene_presence_absence.Rtab"
OUTPUT_DIR = REPO_ROOT / "06_population_structure" / "results" / "pca_mds"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# ==============================================================================
# 4. COLORS
# ==============================================================================

PALETTE = {

    "rUTI": "#ff7f0e",   # ORANGE

    "UTI": "#1f77b4"     # BLUE
}


# ==============================================================================
# 5. LOAD GENE PRESENCE/ABSENCE MATRIX
# ==============================================================================

print("\nLoading Panaroo RTAB matrix...\n")

df = pd.read_csv(input_file, sep="\t")

# ----------------------------------------------------------
# FIRST COLUMN = GENE NAMES
# ----------------------------------------------------------

gene_names = df.iloc[:, 0]

# ----------------------------------------------------------
# REMAINING COLUMNS = GENOMES
# ----------------------------------------------------------

matrix = df.iloc[:, 1:]

# ----------------------------------------------------------
# GENOME NAMES
# ----------------------------------------------------------

genome_names = matrix.columns

print("Total genomes:", len(genome_names))

print("Total genes:", len(gene_names))


# ==============================================================================
# 6. CONVERT TO BINARY MATRIX
# ==============================================================================
# ROWS = GENOMES
# COLS = GENES
# ==============================================================================

binary_matrix = (matrix.T > 0).astype(int)

print("\nBinary matrix shape:")

print(binary_matrix.shape)


# ==============================================================================
# 7. DEFINE GROUPS
# ==============================================================================
# Genome names containing:
# "rUTI" --> rUTI
#
# otherwise --> UTI
# ==============================================================================

groups = np.where(

    genome_names.str.contains(
        "rUTI",
        case=False,
        na=False
    ),

    "rUTI",

    "UTI"
)

groups = pd.Series(groups, index=genome_names)

print("\nGroup distribution:")

print(groups.value_counts())


# ==============================================================================
# 8. CONFIDENCE ELLIPSE FUNCTION
# ==============================================================================

def confidence_ellipse(
    x,
    y,
    ax,
    edgecolor,
    facecolor,
    alpha=0.15
):

    if len(x) < 3:
        return

    cov = np.cov(x, y)

    vals, vecs = np.linalg.eigh(cov)

    order = vals.argsort()[::-1]

    vals = vals[order]

    vecs = vecs[:, order]

    theta = np.degrees(
        np.arctan2(*vecs[:, 0][::-1])
    )

    width, height = (
        2 * np.sqrt(vals * chi2.ppf(0.95, 2))
    )

    ellipse = Ellipse(

        xy=(np.mean(x), np.mean(y)),

        width=width,

        height=height,

        angle=theta,

        edgecolor=edgecolor,

        facecolor=facecolor,

        alpha=alpha,

        linewidth=2
    )

    ax.add_patch(ellipse)


# ==============================================================================
# 9. PCA ANALYSIS
# ==============================================================================

print("\nRunning PCA...\n")

pca = PCA(n_components=2)

pca_coords = pca.fit_transform(binary_matrix)

# ----------------------------------------------------------
# PCA DATAFRAME
# ----------------------------------------------------------

pca_df = pd.DataFrame({

    "PC1": pca_coords[:, 0],

    "PC2": pca_coords[:, 1],

    "Group": groups.values
})

# ----------------------------------------------------------
# VARIANCE EXPLAINED
# ----------------------------------------------------------

variance = pca.explained_variance_ratio_ * 100


# ==============================================================================
# 10. PCA PLOT
# ==============================================================================

fig, ax = plt.subplots(
    figsize=(8, 6),
    dpi=600
)

for grp in ["rUTI", "UTI"]:

    subset = pca_df[
        pca_df["Group"] == grp
    ]

    ax.scatter(

        subset["PC1"],

        subset["PC2"],

        s=90,

        color=PALETTE[grp],

        edgecolor="black",

        linewidth=0.5,

        alpha=0.85,

        label=grp
    )

    confidence_ellipse(

        subset["PC1"],

        subset["PC2"],

        ax=ax,

        edgecolor=PALETTE[grp],

        facecolor=PALETTE[grp],

        alpha=0.15
    )

# ----------------------------------------------------------
# LABELS
# ----------------------------------------------------------

ax.set_title(
    "PCA of Gene Presence/Absence",
    fontsize=18,
    fontweight="bold"
)

ax.set_xlabel(
    f"PC1 ({variance[0]:.1f}%)",
    fontsize=14,
    fontweight="bold"
)

ax.set_ylabel(
    f"PC2 ({variance[1]:.1f}%)",
    fontsize=14,
    fontweight="bold"
)

ax.legend(frameon=False)

ax.grid(alpha=0.2)

plt.tight_layout()

# ----------------------------------------------------------
# SAVE PCA
# ----------------------------------------------------------

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "PCA_Plot_Publication.png"
    ),

    dpi=600,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "PCA_Plot_Publication.pdf"
    ),

    bbox_inches="tight"
)

plt.close()

print("PCA plot saved.")


# ==============================================================================
# 11. JACCARD DISTANCE MATRIX
# ==============================================================================
# EXACT EQUIVALENT OF:
#
# R:
# dist(method="binary")
# ==============================================================================

print("\nCalculating Jaccard distance...\n")

binary_array = binary_matrix.astype(bool).to_numpy()

distance_matrix = squareform(
    pdist(binary_array, metric="jaccard")
)


# ==============================================================================
# 12. CLASSICAL MDS
# ==============================================================================
# EXACT EQUIVALENT OF:
#
# R:
# cmdscale()
#
# THIS IS TRUE CLASSICAL MDS / PCoA
# ==============================================================================

print("Running Classical MDS...\n")

# ----------------------------------------------------------
# NUMBER OF SAMPLES
# ----------------------------------------------------------

n = distance_matrix.shape[0]

# ----------------------------------------------------------
# DOUBLE CENTERING
# ----------------------------------------------------------

H = np.eye(n) - np.ones((n, n)) / n

B = -0.5 * H.dot(distance_matrix ** 2).dot(H)

# ----------------------------------------------------------
# EIGEN DECOMPOSITION
# ----------------------------------------------------------

eigvals, eigvecs = np.linalg.eigh(B)

# ----------------------------------------------------------
# SORT EIGENVALUES DESCENDING
# ----------------------------------------------------------

idx = np.argsort(eigvals)[::-1]

eigvals = eigvals[idx]

eigvecs = eigvecs[:, idx]

# ----------------------------------------------------------
# KEEP FIRST 2 DIMENSIONS
# ----------------------------------------------------------

coords = eigvecs[:, :2] * np.sqrt(eigvals[:2])

# ----------------------------------------------------------
# BUILD DATAFRAME
# ----------------------------------------------------------

mds_df = pd.DataFrame({

    "Dim1": coords[:, 0],

    "Dim2": coords[:, 1],

    "Group": groups.values
})


# ==============================================================================
# 13. MDS PLOT
# ==============================================================================

fig, ax = plt.subplots(
    figsize=(8, 6),
    dpi=600
)

for grp in ["rUTI", "UTI"]:

    subset = mds_df[
        mds_df["Group"] == grp
    ]

    ax.scatter(

        subset["Dim1"],

        subset["Dim2"],

        s=90,

        color=PALETTE[grp],

        edgecolor="black",

        linewidth=0.5,

        alpha=0.85,

        label=grp
    )

    confidence_ellipse(

        subset["Dim1"],

        subset["Dim2"],

        ax=ax,

        edgecolor=PALETTE[grp],

        facecolor=PALETTE[grp],

        alpha=0.15
    )

# ----------------------------------------------------------
# LABELS
# ----------------------------------------------------------

ax.set_title(
    "MDS Plot (Jaccard Distance)",
    fontsize=18,
    fontweight="bold"
)

ax.set_xlabel(
    "Dimension 1",
    fontsize=14,
    fontweight="bold"
)

ax.set_ylabel(
    "Dimension 2",
    fontsize=14,
    fontweight="bold"
)

ax.legend(frameon=False)

ax.grid(alpha=0.2)

plt.tight_layout()

# ----------------------------------------------------------
# SAVE MDS
# ----------------------------------------------------------

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "MDS_Jaccard_Publication.png"
    ),

    dpi=600,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "MDS_Jaccard_Publication.pdf"
    ),

    bbox_inches="tight"
)

plt.close()

print("MDS plot saved.")


# ==============================================================================
# DONE
# ==============================================================================

print("\n================================================")

print("ANALYSIS COMPLETED SUCCESSFULLY")

print("================================================")

print("\nOutput directory:")

print(OUTPUT_DIR)

# ==============================================================================
# PUBLICATION-QUALITY PANGENOME ANALYSIS PIPELINE
# ==============================================================================
# ANALYSES
# ------------------------------------------------------------------------------
# 1. Genetic Shrinkage
# 2. Pangenome Openness (Heaps' Law)
# 3. Core Genome Nucleotide Diversity
#
# INPUT:
# ------------------------------------------------------------------------------
# gene_presence_absence.Rtab
#
# OPTIONAL INPUT:
# ------------------------------------------------------------------------------
# core_gene_alignment.aln
# (required ONLY for nucleotide diversity)
#
# OUTPUT:
# ------------------------------------------------------------------------------
# Publication-quality PNG + PDF figures
# 600 dpi resolution
#
# COLORS:
# ------------------------------------------------------------------------------
# rUTI = Orange
# UTI  = Blue
#
# AUTHOR:
# Sasikala Project
# ==============================================================================


# ==============================================================================
# 1. IMPORT LIBRARIES
# ==============================================================================

import os
import random
import warnings

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import mannwhitneyu

from Bio import AlignIO

warnings.filterwarnings("ignore")


# ==============================================================================
# 2. GLOBAL SETTINGS
# ==============================================================================

plt.rcParams.update({

    "font.family": "DejaVu Sans",

    "font.weight": "bold",

    "axes.labelweight": "bold",

    "axes.titleweight": "bold",

    "axes.titlesize": 18,

    "axes.labelsize": 15,

    "xtick.labelsize": 12,

    "ytick.labelsize": 12,

    "legend.fontsize": 12,

    "figure.dpi": 600
})


# ==============================================================================
# 3. COLORS
# ==============================================================================

PALETTE = {

    "rUTI": "#ff7f0e",

    "UTI": "#1f77b4"
}


# ==============================================================================
# REPOSITORY CONFIGURATION
# ==============================================================================
from pathlib import Path

REPO_ROOT = Path("../..").resolve()

INPUT_FILE = REPO_ROOT / "04_Panaroo" / "results" / "combined" / "gene_presence_absence.Rtab"
CORE_ALIGNMENT = REPO_ROOT / "05_phylogeny" / "results" / "core_gene_alignment.aln"
OUTPUT_DIR = REPO_ROOT / "06_population_structure" / "results" / "pangenome_population_structure"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# ==============================================================================
# 5. LOAD PANAROO MATRIX
# ==============================================================================

print("\nLoading Panaroo matrix...\n")

df = pd.read_csv(INPUT_FILE, sep="\t")

gene_names = df.iloc[:, 0]

matrix = df.iloc[:, 1:]

genome_names = matrix.columns

# Binary matrix
binary_matrix = (matrix.T > 0).astype(int)

print("Binary matrix shape:")
print(binary_matrix.shape)


# ==============================================================================
# 6. DEFINE GROUPS
# ==============================================================================

groups = np.where(

    genome_names.str.contains(
        "rUTI",
        case=False,
        na=False
    ),

    "rUTI",

    "UTI"
)

groups = pd.Series(groups, index=genome_names)

print("\nGroup counts:")
print(groups.value_counts())


# ==============================================================================
# 7. GENETIC SHRINKAGE ANALYSIS
# ==============================================================================

print("\nRunning Genetic Shrinkage Analysis...\n")

gene_counts = binary_matrix.sum(axis=1)

shrinkage_df = pd.DataFrame({

    "Genome": genome_names,

    "Gene_Count": gene_counts.values,

    "Group": groups.values
})

# Statistics
ruti_counts = shrinkage_df[
    shrinkage_df["Group"] == "rUTI"
]["Gene_Count"]

uti_counts = shrinkage_df[
    shrinkage_df["Group"] == "UTI"
]["Gene_Count"]

stat, pvalue = mannwhitneyu(
    ruti_counts,
    uti_counts,
    alternative="two-sided"
)

print("rUTI mean genes:", round(ruti_counts.mean(), 2))
print("UTI mean genes:", round(uti_counts.mean(), 2))
print("P-value:", pvalue)


# ==============================================================================
# 8. GENETIC SHRINKAGE PLOT
# ==============================================================================

fig, ax = plt.subplots(
    figsize=(7, 6),
    dpi=600
)

positions = [1, 2]

data = [ruti_counts, uti_counts]

bp = ax.boxplot(
    data,
    patch_artist=True,
    widths=0.5
)

# Colors
for patch, grp in zip(bp["boxes"], ["rUTI", "UTI"]):

    patch.set_facecolor(PALETTE[grp])

    patch.set_alpha(0.7)

# Jitter points
for i, grp in enumerate(["rUTI", "UTI"]):

    subset = shrinkage_df[
        shrinkage_df["Group"] == grp
    ]

    x = np.random.normal(
        positions[i],
        0.05,
        size=len(subset)
    )

    ax.scatter(
        x,
        subset["Gene_Count"],
        color=PALETTE[grp],
        edgecolor="black",
        alpha=0.8,
        s=50
    )

# Labels
ax.set_xticks([1, 2])

ax.set_xticklabels(
    ["rUTI", "UTI"],
    fontweight="bold"
)

ax.set_ylabel(
    "Genes per Genome",
    fontweight="bold"
)

ax.set_title(
    "Genetic Shrinkage Analysis",
    fontweight="bold"
)

ax.text(
    1.5,
    max(gene_counts) * 1.02,
    f"P = {pvalue:.3e}",
    ha="center",
    fontsize=12,
    fontweight="bold"
)

plt.tight_layout()

# Save
plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Genetic_Shrinkage.png"
    ),
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Genetic_Shrinkage.pdf"
    ),
    bbox_inches="tight"
)

plt.close()

print("Genetic shrinkage plot saved.")


# ==============================================================================
# 9. PANGENOME OPENNESS ANALYSIS
# ==============================================================================
# HEAPS' LAW
# ==============================================================================

print("\nRunning Pangenome Openness Analysis...\n")


# ==============================================================================
# HEAPS FUNCTION
# ==============================================================================

def heaps_law(x, k, alpha):

    return k * (x ** alpha)


# ==============================================================================
# RAREFACTION FUNCTION
# ==============================================================================

def rarefaction_curve(
    subset_matrix,
    n_perm=100
):

    n_genomes = subset_matrix.shape[0]

    results = []

    for n in range(1, n_genomes + 1):

        pan_sizes = []

        for _ in range(n_perm):

            sampled = subset_matrix.sample(
                n=n,
                axis=0,
                replace=False
            )

            pan_size = (
                sampled.sum(axis=0) > 0
            ).sum()

            pan_sizes.append(pan_size)

        results.append({

            "Size": n,

            "Mean": np.mean(pan_sizes),

            "SD": np.std(pan_sizes)
        })

    return pd.DataFrame(results)


# ==============================================================================
# SPLIT GROUPS
# ==============================================================================

ruti_matrix = binary_matrix[
    groups == "rUTI"
]

uti_matrix = binary_matrix[
    groups == "UTI"
]

# Rarefaction
ruti_curve = rarefaction_curve(ruti_matrix)

uti_curve = rarefaction_curve(uti_matrix)

# ==============================================================================
# FIT HEAPS LAW
# ==============================================================================

ruti_params, _ = curve_fit(

    heaps_law,

    ruti_curve["Size"],

    ruti_curve["Mean"],

    maxfev=10000
)

uti_params, _ = curve_fit(

    heaps_law,

    uti_curve["Size"],

    uti_curve["Mean"],

    maxfev=10000
)

alpha_ruti = ruti_params[1]

alpha_uti = uti_params[1]

print("rUTI alpha:", round(alpha_ruti, 4))
print("UTI alpha:", round(alpha_uti, 4))


# ==============================================================================
# 10. OPENNESS PLOT
# ==============================================================================

fig, ax = plt.subplots(
    figsize=(8, 6),
    dpi=600
)

# rUTI
ax.plot(

    ruti_curve["Size"],

    ruti_curve["Mean"],

    color=PALETTE["rUTI"],

    linewidth=3,

    label=f"rUTI (α={alpha_ruti:.3f})"
)

ax.fill_between(

    ruti_curve["Size"],

    ruti_curve["Mean"] - ruti_curve["SD"],

    ruti_curve["Mean"] + ruti_curve["SD"],

    color=PALETTE["rUTI"],

    alpha=0.2
)

# UTI
ax.plot(

    uti_curve["Size"],

    uti_curve["Mean"],

    color=PALETTE["UTI"],

    linewidth=3,

    label=f"UTI (α={alpha_uti:.3f})"
)

ax.fill_between(

    uti_curve["Size"],

    uti_curve["Mean"] - uti_curve["SD"],

    uti_curve["Mean"] + uti_curve["SD"],

    color=PALETTE["UTI"],

    alpha=0.2
)

# Labels
ax.set_title(
    "Pangenome Openness",
    fontweight="bold"
)

ax.set_xlabel(
    "Number of Genomes",
    fontweight="bold"
)

ax.set_ylabel(
    "Pangenome Size",
    fontweight="bold"
)

ax.legend(frameon=False)

ax.grid(alpha=0.2)

plt.tight_layout()

# Save
plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Pangenome_Openness.png"
    ),
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Pangenome_Openness.pdf"
    ),
    bbox_inches="tight"
)

plt.close()

print("Pangenome openness plot saved.")


# ==============================================================================
# 11. CORE GENOME NUCLEOTIDE DIVERSITY
# ==============================================================================
# OPTIONAL
# Requires:
# core_gene_alignment.aln
# ==============================================================================

if os.path.exists(CORE_ALIGNMENT):

    print("\nRunning Core Genome Nucleotide Diversity...\n")

    alignment = AlignIO.read(
        CORE_ALIGNMENT,
        "fasta"
    )

    sequences = []

    names = []

    for record in alignment:

        sequences.append(str(record.seq))

        names.append(record.id)

    seq_df = pd.DataFrame({

        "Genome": names,

        "Sequence": sequences
    })

    seq_df["Group"] = np.where(

        seq_df["Genome"].str.contains(
            "rUTI",
            case=False,
            na=False
        ),

        "rUTI",

        "UTI"
    )

    # Simple pairwise diversity
    diversity = []

    for grp in ["rUTI", "UTI"]:

        subset = seq_df[
            seq_df["Group"] == grp
        ]

        seqs = subset["Sequence"].tolist()

        distances = []

        for i in range(len(seqs)):

            for j in range(i + 1, len(seqs)):

                s1 = seqs[i]
                s2 = seqs[j]

                mismatches = sum(
                    a != b
                    for a, b in zip(s1, s2)
                    if a != "-" and b != "-"
                )

                valid = sum(
                    a != "-" and b != "-"
                    for a, b in zip(s1, s2)
                )

                if valid > 0:

                    distances.append(
                        mismatches / valid
                    )

        diversity.append({

            "Group": grp,

            "Pi": np.mean(distances)
        })

    diversity_df = pd.DataFrame(diversity)

    # Plot
    fig, ax = plt.subplots(
        figsize=(6, 6),
        dpi=600
    )

    bars = ax.bar(

        diversity_df["Group"],

        diversity_df["Pi"],

        color=[
            PALETTE["rUTI"],
            PALETTE["UTI"]
        ],

        edgecolor="black"
    )

    ax.set_ylabel(
        "Nucleotide Diversity (π)",
        fontweight="bold"
    )

    ax.set_title(
        "Core Genome Nucleotide Diversity",
        fontweight="bold"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            OUTPUT_DIR,
            "Core_Genome_Nucleotide_Diversity.png"
        ),
        dpi=600,
        bbox_inches="tight"
    )

    plt.savefig(
        os.path.join(
            OUTPUT_DIR,
            "Core_Genome_Nucleotide_Diversity.pdf"
        ),
        bbox_inches="tight"
    )

    plt.close()

    print("Core genome diversity plot saved.")

else:

    print("\ncore_gene_alignment.aln NOT FOUND")
    print("Skipping nucleotide diversity analysis.")


# ==============================================================================
# 12. SAVE OPENNESS STATS
# ==============================================================================

stats_df = pd.DataFrame({

    "Group": ["rUTI", "UTI"],

    "Alpha": [
        alpha_ruti,
        alpha_uti
    ],

    "Status": [

        "Open" if alpha_ruti < 1 else "Closed",

        "Open" if alpha_uti < 1 else "Closed"
    ]
})

stats_df.to_csv(

    os.path.join(
        OUTPUT_DIR,
        "Openness_Statistics.csv"
    ),

    index=False
)

print("\nStatistics table saved.")


# ==============================================================================
# DONE
# ==============================================================================

print("\n================================================")
print("ALL ANALYSES COMPLETED SUCCESSFULLY")
print("================================================")

print("\nResults saved in:")
print(OUTPUT_DIR)

#!/usr/bin/env python3

# ==============================================================================
# PUBLICATION-QUALITY PANGENOME FIGURES
# ==============================================================================
#
# INPUT FILES
# ------------------------------------------------------------------------------
# 1. Openness_Stats.csv
# 2. Genetic_Shrinkage_Stats.csv
# 3. Genome_Gene_Counts.csv
# 4. Core_Genome_Diversity_Stats.csv
#
# OUTPUT
# ------------------------------------------------------------------------------
# High-resolution publication-quality figures:
#
# 1. Pangenome_Openness_Publication.png/pdf
# 2. Genetic_Shrinkage_Publication.png/pdf
# 3. Core_Genome_Diversity_Publication.png/pdf
#
# FEATURES
# ------------------------------------------------------------------------------
# - 600 dpi
# - vector PDF
# - rUTI = orange
# - UTI = blue
# - bold publication formatting
# - Nature/Frontiers-style aesthetics
#
# AUTHOR
# ------------------------------------------------------------------------------
# Sasikala Project
# ==============================================================================


# ==============================================================================
# 1. IMPORTS
# ==============================================================================

import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu


# ==============================================================================
# 2. GLOBAL SETTINGS
# ==============================================================================

plt.rcParams.update({

    "font.family": "DejaVu Sans",

    "font.weight": "bold",

    "axes.labelweight": "bold",

    "axes.titleweight": "bold",

    "axes.titlesize": 20,

    "axes.labelsize": 16,

    "xtick.labelsize": 13,

    "ytick.labelsize": 13,

    "legend.fontsize": 13,

    "figure.dpi": 600
})

sns.set_style("whitegrid")


# ==============================================================================
# 3. COLORS
# ==============================================================================

COLOR_RUTI = "#ff7f0e"
COLOR_UTI  = "#1f77b4"

PALETTE = {
    "rUTI": COLOR_RUTI,
    "UTI": COLOR_UTI
}


# ==============================================================================
# REPOSITORY CONFIGURATION
# ==============================================================================
from pathlib import Path

REPO_ROOT = Path("../..").resolve()

BASE_DIR = REPO_ROOT / "04_Panaroo" / "results"
OPENNESS_DIR = BASE_DIR / "openness"
CORE_DIR = BASE_DIR / "Core_genome_Diversity"
OUTPUT_DIR = REPO_ROOT / "06_population_structure" / "results" / "pangenome_population_structure"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# ==============================================================================
# 5. LOAD DATA
# ==============================================================================

# ------------------------------------------------------------------------------
# OPENNESS
# ------------------------------------------------------------------------------

openness_stats = pd.read_csv(
    os.path.join(
        OPENNESS_DIR,
        "Openness_Statistics.csv"
    )
)

# ------------------------------------------------------------------------------
# SHRINKAGE
# ------------------------------------------------------------------------------

shrinkage_stats = pd.read_csv(
    os.path.join(
        OPENNESS_DIR,
        "Genetic_Shrinkage_Stats.csv"
    )
)

# ------------------------------------------------------------------------------
# GENE COUNTS
# ------------------------------------------------------------------------------

gene_counts = pd.read_csv(
    os.path.join(
        OPENNESS_DIR,
        "Genome_Gene_Counts.csv"
    )
)

# ------------------------------------------------------------------------------
# CORE DIVERSITY
# ------------------------------------------------------------------------------

core_div = pd.read_csv(
    os.path.join(
        CORE_DIR,
        "Core_Genome_Diversity_Stats.csv"
    )
)


# ==============================================================================
# 6. PANGENOME OPENNESS PLOT
# ==============================================================================

print("\nGenerating openness plot...\n")

fig, ax = plt.subplots(
    figsize=(7, 6),
    dpi=600
)

bars = ax.bar(

    openness_stats["Group"],

    openness_stats["Alpha"],

    color=[
        PALETTE[g]
        for g in openness_stats["Group"]
    ],

    edgecolor="black",

    linewidth=1.5
)

# ------------------------------------------------------------------------------
# LABELS
# ------------------------------------------------------------------------------

for i, row in openness_stats.iterrows():

    ax.text(

        i,

        row["Alpha"] + 0.02,

        f"{row['Alpha']:.3f}",

        ha="center",

        fontsize=13,

        fontweight="bold"
    )

ax.set_ylim(0, 1)

ax.set_ylabel(
    "Pangenome Openness (Alpha)"
)

ax.set_title(
    "Pangenome Openness"
)

# ------------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------------

ax.grid(
    axis="y",
    alpha=0.25
)

plt.tight_layout()

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Pangenome_Openness_Publication.png"
    ),

    dpi=600,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Pangenome_Openness_Publication.pdf"
    ),

    bbox_inches="tight"
)

plt.close()

print("Openness plot saved.")


# ==============================================================================
# 7. GENETIC SHRINKAGE PLOT
# ==============================================================================

print("\nGenerating genetic shrinkage plot...\n")

fig, ax = plt.subplots(
    figsize=(7, 6),
    dpi=600
)

# ------------------------------------------------------------------------------
# BOXPLOT
# ------------------------------------------------------------------------------

sns.boxplot(

    data=gene_counts,

    x="Group",

    y="Gene_Count",

    palette=PALETTE,

    width=0.5,

    linewidth=1.5,

    ax=ax
)

# ------------------------------------------------------------------------------
# STRIP PLOT
# ------------------------------------------------------------------------------

sns.stripplot(

    data=gene_counts,

    x="Group",

    y="Gene_Count",

    color="black",

    size=4,

    alpha=0.5,

    ax=ax
)

# ------------------------------------------------------------------------------
# WILCOXON TEST
# ------------------------------------------------------------------------------

ruti = gene_counts[
    gene_counts["Group"] == "rUTI"
]["Gene_Count"]

uti = gene_counts[
    gene_counts["Group"] == "UTI"
]["Gene_Count"]

stat, pvalue = mannwhitneyu(

    ruti,

    uti,

    alternative="two-sided",

    method="asymptotic",

    use_continuity=True
)

# ------------------------------------------------------------------------------
# P-VALUE LABEL
# ------------------------------------------------------------------------------

ax.text(

    0.5,

    max(gene_counts["Gene_Count"]) + 50,

    f"Wilcoxon p = {pvalue:.4f}",

    ha="center",

    fontsize=13,

    fontweight="bold"
)

ax.set_title(
    "Genetic Shrinkage"
)

ax.set_ylabel(
    "Gene Count per Genome"
)

ax.set_xlabel("")

ax.grid(
    axis="y",
    alpha=0.25
)

plt.tight_layout()

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Genetic_Shrinkage_Publication.png"
    ),

    dpi=600,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Genetic_Shrinkage_Publication.pdf"
    ),

    bbox_inches="tight"
)

plt.close()

print("Genetic shrinkage plot saved.")


# ==============================================================================
# 8. CORE GENOME DIVERSITY
# ==============================================================================

print("\nGenerating core genome diversity plot...\n")

fig, ax = plt.subplots(
    figsize=(7, 6),
    dpi=600
)

bars = ax.bar(

    core_div["Group"],

    core_div["Nucleotide_Diversity_Pi"],

    color=[
        PALETTE[g]
        for g in core_div["Group"]
    ],

    edgecolor="black",

    linewidth=1.5
)

# ------------------------------------------------------------------------------
# LABELS
# ------------------------------------------------------------------------------

for i, row in core_div.iterrows():

    ax.text(

        i,

        row["Nucleotide_Diversity_Pi"] + 0.0003,

        f"{row['Nucleotide_Diversity_Pi']:.4f}",

        ha="center",

        fontsize=13,

        fontweight="bold"
    )

ax.set_ylabel(
    "Nucleotide Diversity (π)"
)

ax.set_title(
    "Core Genome Diversity"
)

ax.grid(
    axis="y",
    alpha=0.25
)

plt.tight_layout()

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Core_Genome_Diversity_Publication.png"
    ),

    dpi=600,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Core_Genome_Diversity_Publication.pdf"
    ),

    bbox_inches="tight"
)

plt.close()

print("Core genome diversity plot saved.")


# ==============================================================================
# DONE
# ==============================================================================

print("\n================================================")
print("PUBLICATION FIGURES GENERATED SUCCESSFULLY")
print("================================================")

print("\nOutput directory:")
print(OUTPUT_DIR)

# ==============================================================================
# PUBLICATION-QUALITY PHYLOGROUP + MLST ENRICHMENT ANALYSIS
# ==============================================================================
#
# FEATURES
# ------------------------------------------------------------------------------
# 1. Top 10 MLSTs shown
# 2. rUTI = Orange
# 3. UTI  = Blue
# 4. Statistical enrichment testing:
#       - Fisher Exact Test
#       - Chi-square automatically
# 5. Per-Phylogroup statistics
# 6. Per-ST statistics
# 7. Significant p-values added on barplots
# 8. Publication-quality figures
#
# OUTPUT
# ------------------------------------------------------------------------------
# 1. Phylogroup_Barplot.png/pdf
# 2. MLST_Barplot.png/pdf
# 3. Phylogroup_Statistics_Detailed.csv
# 4. MLST_Statistics_Detailed.csv
#
# ==============================================================================


# ==============================================================================
# 1. IMPORT LIBRARIES
# ==============================================================================

import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency


# ==============================================================================
# 2. GLOBAL SETTINGS
# ==============================================================================

plt.rcParams.update({

    "font.family": "DejaVu Sans",

    "font.weight": "bold",

    "axes.labelweight": "bold",

    "axes.titleweight": "bold",

    "axes.titlesize": 18,

    "axes.labelsize": 15,

    "xtick.labelsize": 12,

    "ytick.labelsize": 12,

    "legend.fontsize": 12,

    "figure.dpi": 600
})


# ==============================================================================
# 3. COLORS
# ==============================================================================

PALETTE = {

    "rUTI": "#ff7f0e",   # ORANGE

    "UTI": "#1f77b4"     # BLUE
}


# ==============================================================================
# REPOSITORY CONFIGURATION
# ==============================================================================
from pathlib import Path

REPO_ROOT = Path("../..").resolve()

phylogroup_file = REPO_ROOT / "06_population_structure" / "data" / "phylogroups.txt"
mlst_file = REPO_ROOT / "06_population_structure" / "data" / "rUTI_selected_STs.csv"
OUTPUT_DIR = REPO_ROOT / "06_population_structure" / "results" / "phylogroup_mlst"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# ==============================================================================
# 5. LOAD PHYLOGROUP DATA
# ==============================================================================

phylo = pd.read_csv(

    phylogroup_file,

    sep="\t",

    header=None,

    names=["Genome", "Phylogroup"]
)

phylo["Group"] = np.where(

    phylo["Genome"].str.contains(
        "rUTI",
        case=False,
        na=False
    ),

    "rUTI",

    "UTI"
)


# ==============================================================================
# 6. PHYLOGROUP STATISTICS
# ==============================================================================

print("\nRunning phylogroup statistics...\n")

phylo_stats = []

total_ruti = sum(phylo["Group"] == "rUTI")
total_uti  = sum(phylo["Group"] == "UTI")

for pg in sorted(phylo["Phylogroup"].unique()):

    a = len(phylo[
        (phylo["Phylogroup"] == pg) &
        (phylo["Group"] == "rUTI")
    ])

    b = len(phylo[
        (phylo["Phylogroup"] != pg) &
        (phylo["Group"] == "rUTI")
    ])

    c = len(phylo[
        (phylo["Phylogroup"] == pg) &
        (phylo["Group"] == "UTI")
    ])

    d = len(phylo[
        (phylo["Phylogroup"] != pg) &
        (phylo["Group"] == "UTI")
    ])

    contingency = [[a, b], [c, d]]

    oddsratio, pvalue = fisher_exact(contingency)

    phylo_stats.append({

        "Phylogroup": pg,

        "rUTI_Count": a,

        "UTI_Count": c,

        "Odds_Ratio": oddsratio,

        "P_value": pvalue,

        "Significant": "YES" if pvalue < 0.05 else "NO"
    })

phylo_stats_df = pd.DataFrame(phylo_stats)

phylo_stats_df.to_csv(

    os.path.join(
        OUTPUT_DIR,
        "Phylogroup_Statistics_Detailed.csv"
    ),

    index=False
)

print(phylo_stats_df)


# ==============================================================================
# 7. PHYLOGROUP BARPLOT
# ==============================================================================

phylo_plot_df = phylo.groupby(
    ["Phylogroup", "Group"]
).size().reset_index(name="Count")

fig, ax = plt.subplots(
    figsize=(10, 7),
    dpi=600
)

sns.barplot(

    data=phylo_plot_df,

    x="Phylogroup",

    y="Count",

    hue="Group",

    palette=PALETTE,

    ax=ax
)

# ----------------------------------------------------------
# ADD SIGNIFICANCE LABELS
# ----------------------------------------------------------

for i, row in phylo_stats_df.iterrows():

    if row["P_value"] < 0.05:

        max_count = phylo_plot_df[
            phylo_plot_df["Phylogroup"] == row["Phylogroup"]
        ]["Count"].max()

        ax.text(

            i,

            max_count + 1,

            f"p={row['P_value']:.3e}",

            ha="center",

            fontsize=10,

            fontweight="bold"
        )

ax.set_title(
    "Phylogroup Distribution",
    fontsize=20,
    fontweight="bold"
)

ax.set_xlabel(
    "Phylogroup",
    fontsize=15,
    fontweight="bold"
)

ax.set_ylabel(
    "Number of Isolates",
    fontsize=15,
    fontweight="bold"
)

ax.legend(frameon=False)

plt.tight_layout()

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Phylogroup_Barplot.png"
    ),

    dpi=1200,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "Phylogroup_Barplot.pdf"
    ),

    bbox_inches="tight"
)

plt.close()


# ==============================================================================
# 8. LOAD MLST DATA
# ==============================================================================

mlst = pd.read_csv(mlst_file)

mlst["Genome"] = mlst["Isolate"].str.replace(
    ".fasta",
    "",
    regex=False
)

mlst["Group"] = np.where(

    mlst["Genome"].str.contains(
        "rUTI",
        case=False,
        na=False
    ),

    "rUTI",

    "UTI"
)


# ==============================================================================
# 9. KEEP TOP 10 STs
# ==============================================================================

top10_sts = (
    mlst["ST"]
    .value_counts()
    .nlargest(10)
    .index
)

mlst = mlst[
    mlst["ST"].isin(top10_sts)
]


# ==============================================================================
# 10. MLST STATISTICS
# ==============================================================================

print("\nRunning MLST statistics...\n")

mlst_stats = []

total_ruti = sum(mlst["Group"] == "rUTI")
total_uti  = sum(mlst["Group"] == "UTI")

for st in sorted(mlst["ST"].unique()):

    a = len(mlst[
        (mlst["ST"] == st) &
        (mlst["Group"] == "rUTI")
    ])

    b = len(mlst[
        (mlst["ST"] != st) &
        (mlst["Group"] == "rUTI")
    ])

    c = len(mlst[
        (mlst["ST"] == st) &
        (mlst["Group"] == "UTI")
    ])

    d = len(mlst[
        (mlst["ST"] != st) &
        (mlst["Group"] == "UTI")
    ])

    contingency = [[a, b], [c, d]]

    oddsratio, pvalue = fisher_exact(contingency)

    mlst_stats.append({

        "ST": st,

        "rUTI_Count": a,

        "UTI_Count": c,

        "Odds_Ratio": oddsratio,

        "P_value": pvalue,

        "Significant": "YES" if pvalue < 0.05 else "NO"
    })

mlst_stats_df = pd.DataFrame(mlst_stats)

mlst_stats_df.to_csv(

    os.path.join(
        OUTPUT_DIR,
        "MLST_Statistics_Detailed.csv"
    ),

    index=False
)

print(mlst_stats_df)


# ==============================================================================
# 11. MLST BARPLOT
# ==============================================================================

mlst_plot_df = mlst.groupby(
    ["ST", "Group"]
).size().reset_index(name="Count")

fig, ax = plt.subplots(
    figsize=(12, 7),
    dpi=600
)

sns.barplot(

    data=mlst_plot_df,

    x="ST",

    y="Count",

    hue="Group",

    palette=PALETTE,

    ax=ax
)

# ----------------------------------------------------------
# SIGNIFICANCE LABELS
# ----------------------------------------------------------

for i, row in mlst_stats_df.iterrows():

    if row["P_value"] < 0.05:

        max_count = mlst_plot_df[
            mlst_plot_df["ST"] == row["ST"]
        ]["Count"].max()

        ax.text(

            i,

            max_count + 1,

            f"p={row['P_value']:.3e}",

            ha="center",

            fontsize=10,

            fontweight="bold"
        )

ax.set_title(
    "Top 10 MLST Distribution",
    fontsize=20,
    fontweight="bold"
)

ax.set_xlabel(
    "Sequence Type (ST)",
    fontsize=15,
    fontweight="bold"
)

ax.set_ylabel(
    "Number of Isolates",
    fontsize=15,
    fontweight="bold"
)

ax.legend(frameon=False)

plt.xticks(rotation=45)

plt.tight_layout()

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "MLST_Barplot.png"
    ),

    dpi=1200,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "MLST_Barplot.pdf"
    ),

    bbox_inches="tight"
)

plt.close()


# ==============================================================================
# 12. DONE
# ==============================================================================

print("\n====================================================")
print("PUBLICATION ANALYSIS COMPLETED SUCCESSFULLY")
print("====================================================")

print("\nOutput directory:")
print(OUTPUT_DIR)

# ==============================================================================
# PUBLICATION-QUALITY MLST ENRICHMENT ANALYSIS
# ==============================================================================
#
# INPUT FILE
# ------------------------------------------------------------------------------
# MLST TAB FILE FORMAT:
#
# Genome    Species    ST    adk    fumC ...
#
# Example:
#
# rUTI_ERR14226055.fasta    ecoli    421
# rUTI_ERR14226056.fasta    ecoli    69
#
#
# FEATURES
# ------------------------------------------------------------------------------
# 1. Extract Top 10 STs automatically
# 2. Compare rUTI vs UTI
# 3. Fisher exact test per ST
# 4. Significant STs highlighted on barplot
# 5. Publication-quality figures
# 6. rUTI = Orange
# 7. UTI  = Blue
#
#
# OUTPUT
# ------------------------------------------------------------------------------
# 1. MLST_Top10_Barplot.png
# 2. MLST_Top10_Barplot.pdf
# 3. MLST_Statistics_Detailed.csv
# 4. MLST_Count_Table.csv
#
# ==============================================================================


# ==============================================================================
# 1. IMPORT LIBRARIES
# ==============================================================================

import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import fisher_exact


# ==============================================================================
# 2. GLOBAL SETTINGS
# ==============================================================================

plt.rcParams.update({

    "font.family": "DejaVu Sans",

    "font.weight": "bold",

    "axes.labelweight": "bold",

    "axes.titleweight": "bold",

    "axes.titlesize": 18,

    "axes.labelsize": 15,

    "xtick.labelsize": 12,

    "ytick.labelsize": 12,

    "legend.fontsize": 12,

    "figure.dpi": 600
})


# ==============================================================================
# 3. COLORS
# ==============================================================================

PALETTE = {

    "rUTI": "#ff7f0e",   # ORANGE

    "UTI": "#1f77b4"     # BLUE
}


# ==============================================================================
# REPOSITORY CONFIGURATION
# ==============================================================================
from pathlib import Path

REPO_ROOT = Path("../..").resolve()

mlst_file = REPO_ROOT / "06_population_structure" / "data" / "mlst_results.txt"
OUTPUT_DIR = REPO_ROOT / "06_population_structure" / "results" / "phylogroup_mlst"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# ==============================================================================
# 5. LOAD MLST FILE
# ==============================================================================

print("\nLoading MLST data...\n")

mlst = pd.read_csv(

    mlst_file,

    sep="\t",

    header=None
)

# ----------------------------------------------------------
# RENAME COLUMNS
# ----------------------------------------------------------

mlst.columns = [

    "Genome",

    "Species",

    "ST",

    "adk",

    "fumC",

    "gyrB",

    "icd",

    "mdh",

    "purA",

    "recA"
]

# ----------------------------------------------------------
# CLEAN GENOME NAMES
# ----------------------------------------------------------

mlst["Genome"] = mlst["Genome"].str.replace(
    ".fasta",
    "",
    regex=False
)

# ----------------------------------------------------------
# DEFINE GROUPS
# ----------------------------------------------------------

mlst["Group"] = np.where(

    mlst["Genome"].str.contains(
        "rUTI",
        case=False,
        na=False
    ),

    "rUTI",

    "UTI"
)

print(mlst.head())


# ==============================================================================
# 6. TOP 10 STs
# ==============================================================================

print("\nSelecting Top 10 STs...\n")

top10_sts = (

    mlst["ST"]
    .value_counts()
    .nlargest(10)
    .index
)

mlst_top = mlst[
    mlst["ST"].isin(top10_sts)
]

print(top10_sts)


# ==============================================================================
# 7. COUNT TABLE
# ==============================================================================

count_table = pd.crosstab(

    mlst_top["ST"],

    mlst_top["Group"]
)

count_table.to_csv(

    os.path.join(
        OUTPUT_DIR,
        "MLST_Count_Table.csv"
    )
)

print("\nCount table:\n")
print(count_table)


# ==============================================================================
# 8. FISHER EXACT TEST FOR EACH ST
# ==============================================================================

print("\nRunning Fisher exact tests...\n")

stats_results = []

total_ruti = sum(mlst["Group"] == "rUTI")
total_uti  = sum(mlst["Group"] == "UTI")

for st in top10_sts:

    # ------------------------------------------------------
    # PRESENT
    # ------------------------------------------------------

    a = len(mlst[
        (mlst["ST"] == st) &
        (mlst["Group"] == "rUTI")
    ])

    c = len(mlst[
        (mlst["ST"] == st) &
        (mlst["Group"] == "UTI")
    ])

    # ------------------------------------------------------
    # ABSENT
    # ------------------------------------------------------

    b = total_ruti - a

    d = total_uti - c

    contingency = [

        [a, b],

        [c, d]
    ]

    oddsratio, pvalue = fisher_exact(contingency)

    # ------------------------------------------------------
    # ENRICHMENT
    # ------------------------------------------------------

    if a > c:

        enriched = "rUTI"

    elif c > a:

        enriched = "UTI"

    else:

        enriched = "Equal"

    stats_results.append({

        "ST": st,

        "rUTI_Count": a,

        "UTI_Count": c,

        "Enriched_In": enriched,

        "Odds_Ratio": oddsratio,

        "P_value": pvalue,

        "Significant": (
            "YES"
            if pvalue < 0.05
            else "NO"
        )
    })

stats_df = pd.DataFrame(stats_results)

# ----------------------------------------------------------
# SORT BY P VALUE
# ----------------------------------------------------------

stats_df = stats_df.sort_values(
    by="P_value"
)

# ----------------------------------------------------------
# SAVE TABLE
# ----------------------------------------------------------

stats_df.to_csv(

    os.path.join(
        OUTPUT_DIR,
        "MLST_Statistics_Detailed.csv"
    ),

    index=False
)

print(stats_df)


# ==============================================================================
# 9. BARPLOT DATA
# ==============================================================================

plot_df = mlst_top.groupby(
    ["ST", "Group"]
).size().reset_index(name="Count")


# ==============================================================================
# 10. PUBLICATION-QUALITY BARPLOT
# ==============================================================================

print("\nGenerating publication-quality plot...\n")

fig, ax = plt.subplots(
    figsize=(14, 8),
    dpi=600
)

sns.barplot(

    data=plot_df,

    x="ST",

    y="Count",

    hue="Group",

    palette=PALETTE,

    ax=ax
)

# ----------------------------------------------------------
# ADD SIGNIFICANCE LABELS
# ----------------------------------------------------------

for i, st in enumerate(top10_sts):

    row = stats_df[
        stats_df["ST"] == st
    ]

    if len(row) == 0:
        continue

    pval = row["P_value"].values[0]

    sig = row["Significant"].values[0]

    if sig == "YES":

        ymax = plot_df[
            plot_df["ST"] == st
        ]["Count"].max()

        ax.text(

            i,

            ymax + 1,

            f"p={pval:.2e}",

            ha="center",

            fontsize=10,

            fontweight="bold"
        )

# ----------------------------------------------------------
# LABELS
# ----------------------------------------------------------

ax.set_title(

    "Top 10 MLST Distribution",

    fontsize=22,

    fontweight="bold"
)

ax.set_xlabel(

    "Sequence Type (ST)",

    fontsize=16,

    fontweight="bold"
)

ax.set_ylabel(

    "Number of Isolates",

    fontsize=16,

    fontweight="bold"
)

ax.legend(

    title="Group",

    frameon=False
)

plt.xticks(rotation=45)

plt.tight_layout()


# ==============================================================================
# 11. SAVE FIGURES
# ==============================================================================

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "MLST_Top10_Barplot.png"
    ),

    dpi=1200,

    bbox_inches="tight"
)

plt.savefig(

    os.path.join(
        OUTPUT_DIR,
        "MLST_Top10_Barplot.pdf"
    ),

    bbox_inches="tight"
)

plt.close()


# ==============================================================================
# 12. DONE
# ==============================================================================

print("\n====================================================")
print("MLST ENRICHMENT ANALYSIS COMPLETED")
print("====================================================")

print("\nOutput directory:")
print(OUTPUT_DIR)

