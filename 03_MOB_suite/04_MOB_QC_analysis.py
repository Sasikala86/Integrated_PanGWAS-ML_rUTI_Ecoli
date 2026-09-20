# ============================================================
# MOB-SUITE CONTIG-LEVEL QC ANALYSIS
# FOR rUTI vs UTI E. coli GENOMES
# ============================================================
#
# INPUT ROOT:
# /home/sastra/Sasikala/Pangenome/mob_output/
#
# Expected structure:
#
# mob_output/
# ├── UTI_GCA_000776055/
# │   ├── contig_report.txt
# │   ├── mobtyper_results.txt
# │   ├── plasmid_*.fasta
# │   └── ...
# │
# ├── UTI_ERRxxxxxxx/
# │   └── contig_report.txt
# │
# ├── rUTI_GCA_xxxxx/
# │   └── contig_report.txt
# │
# └── rUTI_ERRxxxxxxx/
#     └── contig_report.txt
#
# Phenotype is determined ONLY from the genome folder/sample ID:
#
#   UTI_...   = UTI
#   rUTI_...  = rUTI
#
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np
import re
import warnings

warnings.filterwarnings("ignore")


# ============================================================
# 1. PATHS
# ============================================================

INPUT_DIR = Path(
    "/home/sastra/Sasikala/Pangenome/mob_output"
)

OUTPUT_DIR = INPUT_DIR / "MOB_QC_analysis"

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. BASIC SETTINGS
# ============================================================

CONTIG_REPORT_NAME = "contig_report.txt"


# ============================================================
# 3. PHENOTYPE FROM GENOME ID
# ============================================================

def get_phenotype(genome_id):

    """
    Determine phenotype from genome/sample ID.

    rUTI_XXXX  -> rUTI
    UTI_XXXX   -> UTI

    GCA/ERR accession type does not matter.
    """

    genome_id = str(genome_id).strip()

    if genome_id.startswith("rUTI_"):
        return "rUTI"

    elif genome_id.startswith("UTI_"):
        return "UTI"

    else:
        return "Unknown"


# ============================================================
# 4. CLEAN MOB-SUITE VALUES
# ============================================================

def clean_value(x):

    if pd.isna(x):
        return np.nan

    x = str(x).strip()

    missing_values = {
        "",
        "-",
        "--",
        "---",
        "NA",
        "N/A",
        "na",
        "n/a",
        "None",
        "none",
        "NULL",
        "null"
    }

    if x in missing_values:
        return np.nan

    return x


# ============================================================
# 5. MOLECULE TYPE
# ============================================================

def normalize_molecule_type(x):

    """
    Important:
    We use MOB-suite molecule_type as the compartment assignment.

    chromosome -> chromosome
    plasmid    -> plasmid

    Missing/ambiguous/other -> unclassified

    IMPORTANT:
    rep_type = '-' does NOT make a chromosome contig ambiguous.
    """

    if pd.isna(x):
        return "unclassified"

    x = str(x).strip().lower()

    if x == "chromosome":
        return "chromosome"

    if x == "plasmid":
        return "plasmid"

    if x in {
        "unclassified",
        "unknown",
        "ambiguous",
        "other",
        "",
        "-"
    }:
        return "unclassified"

    return "other"


# ============================================================
# 6. CIRCULARITY
# ============================================================

def classify_circularity(x):

    if pd.isna(x):
        return "unknown"

    x = str(x).strip().lower()

    if x in {
        "yes",
        "true",
        "circular",
        "complete"
    }:
        return "circular"

    if x in {
        "no",
        "false",
        "not_circular",
        "linear"
    }:
        return "not_circular"

    if "not tested" in x:
        return "not_tested"

    if x in {
        "",
        "-",
        "na",
        "n/a"
    }:
        return "unknown"

    return x


# ============================================================
# 7. MOBILITY
# ============================================================

def classify_mobility(x):

    if pd.isna(x):
        return "unknown"

    x = str(x).strip().lower()

    if x == "conjugative":
        return "conjugative"

    if x == "mobilizable":
        return "mobilizable"

    if x in {
        "non-mobilizable",
        "nonmobilizable",
        "non_mobilizable"
    }:
        return "non-mobilizable"

    return x


# ============================================================
# 8. BOOLEAN FIELD
# ============================================================

def has_value(x):

    return not pd.isna(
        clean_value(x)
    )


# ============================================================
# 9. FIND GENOME DIRECTORIES
# ============================================================

print("=" * 80)
print("MOB-SUITE CONTIG-LEVEL QC")
print("=" * 80)

print("\nInput directory:")
print(INPUT_DIR)

if not INPUT_DIR.exists():

    raise FileNotFoundError(
        f"\nInput directory does not exist:\n{INPUT_DIR}"
    )


# Find ONLY contig_report.txt files
report_files = sorted(
    INPUT_DIR.rglob(CONTIG_REPORT_NAME)
)


print(
    f"\ncontig_report.txt files found: "
    f"{len(report_files)}"
)


if len(report_files) == 0:

    raise FileNotFoundError(
        "\nNo contig_report.txt files were found."
    )


# ============================================================
# 10. CHECK GENOME IDs
# ============================================================

genome_records = []

for report_file in report_files:

    # Parent folder = genome ID
    genome_id = report_file.parent.name

    phenotype = get_phenotype(
        genome_id
    )

    genome_records.append({

        "genome_id": genome_id,

        "phenotype": phenotype,

        "contig_report":
            str(report_file),

        "genome_directory":
            str(report_file.parent)
    })


genome_records_df = pd.DataFrame(
    genome_records
)


# ============================================================
# 11. CHECK UNKNOWN IDs
# ============================================================

unknown_genomes = (
    genome_records_df[
        genome_records_df["phenotype"]
        == "Unknown"
    ]
)


print("\nGenome classification:")
print(
    genome_records_df["phenotype"]
    .value_counts(dropna=False)
)


if len(unknown_genomes) > 0:

    print("\nWARNING:")
    print(
        f"{len(unknown_genomes)} genome(s) "
        "could not be classified as UTI/rUTI."
    )

    print(
        unknown_genomes[
            ["genome_id", "contig_report"]
        ].to_string(index=False)
    )


# ============================================================
# 12. READ ALL CONTIG REPORTS
# ============================================================

all_reports = []
failed_reports = []

for _, genome_row in genome_records_df.iterrows():

    report_file = Path(
        genome_row["contig_report"]
    )

    genome_id = genome_row["genome_id"]

    phenotype = genome_row["phenotype"]

    try:

        df = pd.read_csv(
            report_file,
            sep="\t",
            dtype=str,
            comment="#"
        )

        # Remove accidental unnamed columns
        df = df.loc[
            :,
            ~df.columns.astype(str)
            .str.startswith("Unnamed")
        ]

        # Clean column names
        df.columns = [
            str(c).strip()
            for c in df.columns
        ]

        # Add our own metadata
        df["genome_id"] = genome_id
        df["phenotype"] = phenotype
        df["source_file"] = str(report_file)

        all_reports.append(df)

    except Exception as e:

        failed_reports.append({

            "genome_id": genome_id,

            "phenotype": phenotype,

            "file": str(report_file),

            "error": str(e)
        })


# ============================================================
# 13. COMBINE
# ============================================================

if len(all_reports) == 0:

    raise RuntimeError(
        "No contig reports could be successfully read."
    )


contigs = pd.concat(
    all_reports,
    ignore_index=True
)


print(
    f"\nTotal contig rows read: "
    f"{len(contigs):,}"
)


# ============================================================
# 14. REQUIRED COLUMNS
# ============================================================

required_columns = [

    "sample_id",
    "molecule_type",

    "primary_cluster_id",
    "secondary_cluster_id",

    "contig_id",

    "size",
    "gc",

    "circularity_status",

    "rep_type(s)",
    "rep_type_accession(s)",

    "relaxase_type(s)",
    "relaxase_type_accession(s)",

    "mpf_type",
    "mpf_type_accession(s)",

    "orit_type(s)",
    "orit_accession(s)",

    "predicted_mobility",

    "mash_nearest_neighbor",
    "mash_neighbor_distance",
    "mash_neighbor_identification",

    "repetitive_dna_id",
    "repetitive_dna_type",

    "filtering_reason"
]


missing_columns = [
    c for c in required_columns
    if c not in contigs.columns
]


if missing_columns:

    print("\nWARNING: These expected columns are missing:")

    for c in missing_columns:
        print("   ", c)

    print("\nAvailable columns:")

    for c in contigs.columns:
        print("   ", c)


# ============================================================
# 15. CLEAN ALL VALUES
# ============================================================

for col in contigs.columns:

    if contigs[col].dtype == "object":

        contigs[col] = (
            contigs[col]
            .apply(clean_value)
        )


# ============================================================
# 16. NORMALIZE MOLECULE TYPE
# ============================================================

contigs["molecule_class"] = (
    contigs["molecule_type"]
    .apply(normalize_molecule_type)
)


# ============================================================
# 17. SIZE
# ============================================================

if "size" in contigs.columns:

    contigs["size_bp"] = pd.to_numeric(
        contigs["size"],
        errors="coerce"
    )

else:

    contigs["size_bp"] = np.nan


# ============================================================
# 18. COMPARTMENT FLAGS
# ============================================================

contigs["is_chromosome"] = (
    contigs["molecule_class"]
    == "chromosome"
)

contigs["is_plasmid"] = (
    contigs["molecule_class"]
    == "plasmid"
)

contigs["is_unclassified"] = (
    contigs["molecule_class"]
    .isin([
        "unclassified",
        "other"
    ])
)


# ============================================================
# 19. PLASMID EVIDENCE FLAGS
# ============================================================

if "rep_type(s)" in contigs.columns:

    contigs["has_replicon"] = (
        contigs["rep_type(s)"]
        .apply(has_value)
    )

else:

    contigs["has_replicon"] = False


if "relaxase_type(s)" in contigs.columns:

    contigs["has_relaxase"] = (
        contigs["relaxase_type(s)"]
        .apply(has_value)
    )

else:

    contigs["has_relaxase"] = False


if "mpf_type" in contigs.columns:

    contigs["has_mpf"] = (
        contigs["mpf_type"]
        .apply(has_value)
    )

else:

    contigs["has_mpf"] = False


if "orit_type(s)" in contigs.columns:

    contigs["has_orit"] = (
        contigs["orit_type(s)"]
        .apply(has_value)
    )

else:

    contigs["has_orit"] = False


if "primary_cluster_id" in contigs.columns:

    contigs["has_primary_cluster"] = (
        contigs["primary_cluster_id"]
        .apply(has_value)
    )

else:

    contigs["has_primary_cluster"] = False


if "secondary_cluster_id" in contigs.columns:

    contigs["has_secondary_cluster"] = (
        contigs["secondary_cluster_id"]
        .apply(has_value)
    )

else:

    contigs["has_secondary_cluster"] = False


# ============================================================
# 20. CIRCULARITY
# ============================================================

if "circularity_status" in contigs.columns:

    contigs["circularity_class"] = (
        contigs["circularity_status"]
        .apply(classify_circularity)
    )

else:

    contigs["circularity_class"] = "unknown"


# ============================================================
# 21. MOBILITY
# ============================================================

if "predicted_mobility" in contigs.columns:

    contigs["mobility_class"] = (
        contigs["predicted_mobility"]
        .apply(classify_mobility)
    )

else:

    contigs["mobility_class"] = "unknown"


# ============================================================
# 22. CONTIG-LEVEL QC TABLE
# ============================================================

contig_qc = contigs[[
    "genome_id",
    "phenotype",
    "contig_id",
    "size_bp",
    "molecule_type",
    "molecule_class",
    "primary_cluster_id",
    "secondary_cluster_id",
    "circularity_status",
    "circularity_class",
    "rep_type(s)",
    "rep_type_accession(s)",
    "relaxase_type(s)",
    "relaxase_type_accession(s)",
    "mpf_type",
    "mpf_type_accession(s)",
    "orit_type(s)",
    "orit_accession(s)",
    "predicted_mobility",
    "mobility_class",
    "mash_nearest_neighbor",
    "mash_neighbor_distance",
    "mash_neighbor_identification",
    "filtering_reason"
]].copy()


# ============================================================
# 23. GENOME-LEVEL SUMMARY
# ============================================================

assembly_records = []

for genome_id, g in contigs.groupby(
    "genome_id",
    dropna=False
):

    phenotype = g["phenotype"].iloc[0]

    total = len(g)

    chromosome = (
        g["is_chromosome"].sum()
    )

    plasmid = (
        g["is_plasmid"].sum()
    )

    ambiguous = (
        g["is_unclassified"].sum()
    )

    assembly_records.append({

        "genome_id":
            genome_id,

        "phenotype":
            phenotype,

        "total_contigs":
            total,

        "chromosome_contigs":
            chromosome,

        "plasmid_contigs":
            plasmid,

        "unclassified_or_ambiguous_contigs":
            ambiguous,

        "percent_chromosome_contigs":
            100 * chromosome / total
            if total > 0 else np.nan,

        "percent_plasmid_contigs":
            100 * plasmid / total
            if total > 0 else np.nan,

        "percent_unclassified_contigs":
            100 * ambiguous / total
            if total > 0 else np.nan,

        "has_unclassified_contig":
            ambiguous > 0,

        "has_plasmid_contig":
            plasmid > 0,

        "replicon_supported_plasmid_contigs":
            (
                g["is_plasmid"]
                &
                g["has_replicon"]
            ).sum(),

        "relaxase_supported_plasmid_contigs":
            (
                g["is_plasmid"]
                &
                g["has_relaxase"]
            ).sum(),

        "oriT_supported_plasmid_contigs":
            (
                g["is_plasmid"]
                &
                g["has_orit"]
            ).sum(),

        "MPF_supported_plasmid_contigs":
            (
                g["is_plasmid"]
                &
                g["has_mpf"]
            ).sum(),

        "circular_plasmid_contigs":
            (
                g["is_plasmid"]
                &
                (
                    g["circularity_class"]
                    == "circular"
                )
            ).sum()
    })


assembly_summary = pd.DataFrame(
    assembly_records
)


# ============================================================
# 24. GROUP-LEVEL SUMMARY
# ============================================================

def group_summary(g, group_name):

    n_assemblies = (
        g["genome_id"]
        .nunique()
    )

    total_contigs = len(g)

    chromosome = (
        g["is_chromosome"].sum()
    )

    plasmid = (
        g["is_plasmid"].sum()
    )

    ambiguous = (
        g["is_unclassified"].sum()
    )

    genome_ids = (
        g["genome_id"]
        .dropna()
        .unique()
    )

    genome_table = assembly_summary[
        assembly_summary["genome_id"]
        .isin(genome_ids)
    ]

    assemblies_with_ambiguous = (
        genome_table[
            "has_unclassified_contig"
        ].sum()
    )

    assemblies_with_plasmid = (
        genome_table[
            "has_plasmid_contig"
        ].sum()
    )

    return {

        "group":
            group_name,

        "number_of_assemblies":
            n_assemblies,

        "total_contigs":
            total_contigs,

        "chromosome_contigs":
            chromosome,

        "plasmid_contigs":
            plasmid,

        "unclassified_or_ambiguous_contigs":
            ambiguous,

        "percent_chromosome_contigs":
            100 * chromosome / total_contigs
            if total_contigs else np.nan,

        "percent_plasmid_contigs":
            100 * plasmid / total_contigs
            if total_contigs else np.nan,

        "percent_unclassified_contigs":
            100 * ambiguous / total_contigs
            if total_contigs else np.nan,

        "assemblies_with_ambiguous_contig":
            assemblies_with_ambiguous,

        "percent_assemblies_with_ambiguous":
            100 * assemblies_with_ambiguous
            / n_assemblies
            if n_assemblies else np.nan,

        "assemblies_with_plasmid_contig":
            assemblies_with_plasmid,

        "percent_assemblies_with_plasmid":
            100 * assemblies_with_plasmid
            / n_assemblies
            if n_assemblies else np.nan,

        "plasmid_contigs_with_replicon":
            (
                g["is_plasmid"]
                &
                g["has_replicon"]
            ).sum(),

        "plasmid_contigs_with_relaxase":
            (
                g["is_plasmid"]
                &
                g["has_relaxase"]
            ).sum(),

        "plasmid_contigs_with_oriT":
            (
                g["is_plasmid"]
                &
                g["has_orit"]
            ).sum(),

        "plasmid_contigs_with_MPF":
            (
                g["is_plasmid"]
                &
                g["has_mpf"]
            ).sum(),

        "circular_plasmid_contigs":
            (
                g["is_plasmid"]
                &
                (
                    g["circularity_class"]
                    == "circular"
                )
            ).sum()
    }


summary_rows = []

for group in ["rUTI", "UTI"]:

    g = contigs[
        contigs["phenotype"]
        == group
    ]

    if len(g) > 0:

        summary_rows.append(
            group_summary(
                g,
                group
            )
        )


# Overall only for recognized genomes
recognized = contigs[
    contigs["phenotype"]
    .isin(["rUTI", "UTI"])
]

if len(recognized) > 0:

    summary_rows.append(
        group_summary(
            recognized,
            "Overall"
        )
    )


qc_summary = pd.DataFrame(
    summary_rows
)


# ============================================================
# 25. UNCLASSIFIED / AMBIGUOUS CONTIGS
# ============================================================

ambiguous_contigs = contigs[
    contigs["is_unclassified"]
].copy()


ambiguous_columns = [

    "genome_id",
    "phenotype",
    "contig_id",
    "size_bp",

    "molecule_type",
    "molecule_class",

    "primary_cluster_id",
    "secondary_cluster_id",

    "circularity_status",

    "rep_type(s)",
    "relaxase_type(s)",
    "mpf_type",
    "orit_type(s)",

    "predicted_mobility",

    "mash_nearest_neighbor",
    "mash_neighbor_distance",
    "mash_neighbor_identification",

    "filtering_reason"
]


ambiguous_columns = [
    c for c in ambiguous_columns
    if c in ambiguous_contigs.columns
]


ambiguous_contigs = ambiguous_contigs[
    ambiguous_columns
]


# ============================================================
# 26. PLASMID CONTIGS
# ============================================================

plasmid_contigs = contigs[
    contigs["is_plasmid"]
].copy()


# ============================================================
# 27. PLASMID EVIDENCE CATEGORY
# ============================================================

def evidence_category(row):

    evidence_count = sum([

        bool(row["has_replicon"]),

        bool(row["has_relaxase"]),

        bool(row["has_orit"]),

        bool(row["has_mpf"]),

        bool(row["has_primary_cluster"])
    ])

    if evidence_count >= 3:

        return "multiple_supporting_features"

    elif evidence_count >= 1:

        return "limited_supporting_features"

    else:

        return "no_detected_plasmid_markers"


plasmid_contigs[
    "MOB_suite_evidence_category"
] = plasmid_contigs.apply(
    evidence_category,
    axis=1
)


# ============================================================
# 28. PLASMID CONTIG EVIDENCE TABLE
# ============================================================

plasmid_evidence_columns = [

    "genome_id",
    "phenotype",

    "contig_id",
    "size_bp",

    "primary_cluster_id",
    "secondary_cluster_id",

    "circularity_status",
    "circularity_class",

    "rep_type(s)",
    "rep_type_accession(s)",

    "relaxase_type(s)",
    "relaxase_type_accession(s)",

    "mpf_type",
    "mpf_type_accession(s)",

    "orit_type(s)",
    "orit_accession(s)",

    "predicted_mobility",
    "mobility_class",

    "mash_nearest_neighbor",
    "mash_neighbor_distance",
    "mash_neighbor_identification",

    "filtering_reason",

    "has_replicon",
    "has_relaxase",
    "has_orit",
    "has_mpf",
    "has_primary_cluster",
    "has_secondary_cluster",

    "MOB_suite_evidence_category"
]


plasmid_evidence_columns = [
    c for c in plasmid_evidence_columns
    if c in plasmid_contigs.columns
]


plasmid_evidence = plasmid_contigs[
    plasmid_evidence_columns
].copy()


# ============================================================
# 29. PLASMID RECONSTRUCTION SUMMARY
# ============================================================

reconstruction_records = []


for (
    genome_id,
    cluster_id
), g in plasmid_contigs.groupby(
    [
        "genome_id",
        "primary_cluster_id"
    ],
    dropna=False
):

    phenotype = (
        g["phenotype"].iloc[0]
    )

    reconstruction_records.append({

        "genome_id":
            genome_id,

        "phenotype":
            phenotype,

        "primary_cluster_id":
            cluster_id,

        "num_contigs":
            len(g),

        "total_length_bp":
            g["size_bp"].sum(),

        "has_replicon":
            g["has_replicon"].any(),

        "has_relaxase":
            g["has_relaxase"].any(),

        "has_oriT":
            g["has_orit"].any(),

        "has_MPF":
            g["has_mpf"].any(),

        "has_circular_contig":
            (
                g["circularity_class"]
                == "circular"
            ).any(),

        "mobility_types":
            ";".join(
                sorted(
                    g["mobility_class"]
                    .dropna()
                    .unique()
                )
            ),

        "replicon_types":
            ";".join(
                sorted(
                    g["rep_type(s)"]
                    .dropna()
                    .astype(str)
                    .unique()
                )
            ),

        "relaxase_types":
            ";".join(
                sorted(
                    g["relaxase_type(s)"]
                    .dropna()
                    .astype(str)
                    .unique()
                )
            ),

        "MPF_types":
            ";".join(
                sorted(
                    g["mpf_type"]
                    .dropna()
                    .astype(str)
                    .unique()
                )
            ),

        "oriT_types":
            ";".join(
                sorted(
                    g["orit_type(s)"]
                    .dropna()
                    .astype(str)
                    .unique()
                )
            )
    })


reconstruction_summary = pd.DataFrame(
    reconstruction_records
)


# ============================================================
# 30. SAVE ALL RESULTS
# ============================================================

qc_summary.to_csv(
    OUTPUT_DIR /
    "01_MOB_QC_summary.csv",
    index=False
)

assembly_summary.to_csv(
    OUTPUT_DIR /
    "02_genome_level_QC.csv",
    index=False
)

contig_qc.to_csv(
    OUTPUT_DIR /
    "03_all_contigs_QC.csv",
    index=False
)

plasmid_evidence.to_csv(
    OUTPUT_DIR /
    "04_plasmid_contig_evidence.csv",
    index=False
)

reconstruction_summary.to_csv(
    OUTPUT_DIR /
    "05_plasmid_reconstruction_summary.csv",
    index=False
)

ambiguous_contigs.to_csv(
    OUTPUT_DIR /
    "06_unclassified_ambiguous_contigs.csv",
    index=False
)

genome_records_df.to_csv(
    OUTPUT_DIR /
    "07_genome_file_inventory.csv",
    index=False
)


if failed_reports:

    pd.DataFrame(
        failed_reports
    ).to_csv(
        OUTPUT_DIR /
        "08_failed_reports.csv",
        index=False
    )


# ============================================================
# 31. GENERATE HUMAN-READABLE REPORT
# ============================================================

report_lines = []

report_lines.append(
    "MOB-SUITE CONTIG-LEVEL QC REPORT"
)

report_lines.append(
    "=" * 70
)

report_lines.append(
    f"Input directory: {INPUT_DIR}"
)

report_lines.append(
    f"contig_report.txt files found: "
    f"{len(report_files)}"
)

report_lines.append(
    f"Reports successfully read: "
    f"{len(all_reports)}"
)

report_lines.append(
    f"Reports failed: "
    f"{len(failed_reports)}"
)

report_lines.append("")


report_lines.append(
    "GENOME CLASSIFICATION"
)

report_lines.append(
    "-" * 70
)

for group in [
    "rUTI",
    "UTI",
    "Unknown"
]:

    n = (
        genome_records_df[
            genome_records_df["phenotype"]
            == group
        ]["genome_id"]
        .nunique()
    )

    report_lines.append(
        f"{group}: {n}"
    )


report_lines.append("")


report_lines.append(
    "CONTIG-LEVEL QC"
)

report_lines.append(
    "-" * 70
)


for _, row in qc_summary.iterrows():

    report_lines.append("")

    report_lines.append(
        f"GROUP: {row['group']}"
    )

    report_lines.append(
        f"Assemblies: "
        f"{int(row['number_of_assemblies'])}"
    )

    report_lines.append(
        f"Total contigs: "
        f"{int(row['total_contigs'])}"
    )

    report_lines.append(
        f"Chromosome contigs: "
        f"{int(row['chromosome_contigs'])} "
        f"({row['percent_chromosome_contigs']:.2f}%)"
    )

    report_lines.append(
        f"Plasmid contigs: "
        f"{int(row['plasmid_contigs'])} "
        f"({row['percent_plasmid_contigs']:.2f}%)"
    )

    report_lines.append(
        f"Unclassified/ambiguous contigs: "
        f"{int(row['unclassified_or_ambiguous_contigs'])} "
        f"({row['percent_unclassified_contigs']:.2f}%)"
    )

    report_lines.append(
        f"Assemblies containing >=1 "
        f"unclassified/ambiguous contig: "
        f"{int(row['assemblies_with_ambiguous_contig'])} / "
        f"{int(row['number_of_assemblies'])} "
        f"({row['percent_assemblies_with_ambiguous']:.2f}%)"
    )

    report_lines.append(
        f"Assemblies containing >=1 plasmid contig: "
        f"{int(row['assemblies_with_plasmid_contig'])} / "
        f"{int(row['number_of_assemblies'])} "
        f"({row['percent_assemblies_with_plasmid']:.2f}%)"
    )

    report_lines.append(
        f"Plasmid contigs with replicon: "
        f"{int(row['plasmid_contigs_with_replicon'])}"
    )

    report_lines.append(
        f"Plasmid contigs with relaxase: "
        f"{int(row['plasmid_contigs_with_relaxase'])}"
    )

    report_lines.append(
        f"Plasmid contigs with oriT: "
        f"{int(row['plasmid_contigs_with_oriT'])}"
    )

    report_lines.append(
        f"Plasmid contigs with MPF: "
        f"{int(row['plasmid_contigs_with_MPF'])}"
    )

    report_lines.append(
        f"Circular plasmid contigs: "
        f"{int(row['circular_plasmid_contigs'])}"
    )


report_lines.append("")

report_lines.append(
    "IMPORTANT INTERPRETATION"
)

report_lines.append(
    "-" * 70
)

report_lines.append(
    "A contig was considered chromosome- or "
    "plasmid-assigned according to the MOB-suite "
    "molecule_type field."
)

report_lines.append(
    "Absence of a replicon was NOT treated as "
    "evidence of ambiguity when molecule_type "
    "was explicitly assigned."
)

report_lines.append(
    "The number of contigs within a plasmid "
    "reconstruction was NOT treated as the number "
    "of ambiguous contigs."
)


with open(
    OUTPUT_DIR /
    "09_MOB_QC_text_report.txt",
    "w"
) as f:

    f.write(
        "\n".join(report_lines)
    )


# ============================================================
# 32. FINAL CONSOLE OUTPUT
# ============================================================

print("\n")
print("=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)

print(
    f"\nResults directory:\n{OUTPUT_DIR}"
)

print("\n")
print("GROUP-LEVEL QC SUMMARY")
print("-" * 80)

print(
    qc_summary.to_string(
        index=False
    )
)


print("\n")
print("GENOME COUNTS")
print("-" * 80)

print(
    genome_records_df[
        "phenotype"
    ].value_counts(
        dropna=False
    ).to_string()
)


print("\n")
print("UNCLASSIFIED / AMBIGUOUS CONTIGS")
print("-" * 80)

print(
    f"Total unclassified/ambiguous contigs: "
    f"{len(ambiguous_contigs)}"
)

if len(ambiguous_contigs) > 0:

    print(
        ambiguous_contigs[
            [
                "genome_id",
                "phenotype",
                "contig_id",
                "size_bp",
                "molecule_type",
                "primary_cluster_id"
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


print("\n")
print("FILES GENERATED")
print("-" * 80)

for file in sorted(
    OUTPUT_DIR.iterdir()
):

    print(
        " ",
        file.name
    )


print("\nFinished successfully.")
