#!/usr/bin/env python3

# ============================================================
# NCBI ASSEMBLY QC SUMMARY
#
# Purpose:
#   Retrieve NCBI assembly statistics for an already-selected
#   genome list.
#
# IMPORTANT:
#
#   The genome list supplied to this script is ALREADY
#   SELECTED.
#
#   This script does NOT perform new genome selection.
#
#   No_of_contigs is taken specifically from:
#
#       NCBI assembly_summary_genbank.txt
#       -> contig_count
#
#   We DO NOT use:
#
#       assmstats-number-of-contigs
#
#   We DO NOT download a new assembly_summary file.
#
# FINAL OUTPUT:
#
#   Sample
#   assembly_accession
#   No_of_contigs
#   TotalLength
#   N50
#   L50
#   GC
#   Coverage
#   CheckM_Completeness
#   CheckM_Contamination
# ============================================================

from pathlib import Path
import argparse
import pandas as pd
import subprocess
import re


# ============================================================
# 1. COMMAND-LINE ARGUMENTS
# ============================================================

parser = argparse.ArgumentParser(
    description="Generate NCBI assembly QC summary."
)

parser.add_argument(
    "--genome-list",
    required=True,
    help="Text file containing already-selected genome IDs/accessions."
)

parser.add_argument(
    "--assembly-summary",
    required=True,
    help=(
        "Existing NCBI assembly_summary_genbank.txt "
        "(not downloaded by this script)."
    )
)

parser.add_argument(
    "--output-dir",
    required=True,
    help="Directory for output files."
)

args = parser.parse_args()


input_file = Path(
    args.genome_list
)

assembly_summary_file = Path(
    args.assembly_summary
)

output_dir = Path(
    args.output_dir
)


# ============================================================
# 2. CREATE OUTPUT DIRECTORY
# ============================================================

output_dir.mkdir(
    parents=True,
    exist_ok=True
)


print("=" * 75)
print("NCBI ASSEMBLY QC SUMMARY")
print("=" * 75)

print(
    f"\nGenome list:\n{input_file}"
)

print(
    f"\nExisting assembly summary:\n"
    f"{assembly_summary_file}"
)

print(
    f"\nOutput directory:\n{output_dir}"
)


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

if not input_file.exists():

    raise FileNotFoundError(
        f"Genome list not found:\n{input_file}"
    )


if not assembly_summary_file.exists():

    raise FileNotFoundError(
        "Existing NCBI assembly summary not found:\n"
        f"{assembly_summary_file}"
    )


# ============================================================
# 4. READ ALREADY-SELECTED GENOME LIST
# ============================================================

with open(
    input_file,
    "r"
) as f:

    raw_lines = [
        line.strip()
        for line in f
        if line.strip()
    ]


print(
    f"\nGenome-list entries: "
    f"{len(raw_lines)}"
)


# ============================================================
# 5. EXTRACT GCA/GCF ACCESSIONS
# ============================================================

def extract_accession(text):

    match = re.search(
        r"(GC[AF]_\d+\.\d+)",
        text
    )

    if match:

        return match.group(1)

    return None


records = []

for line in raw_lines:

    accession = extract_accession(
        line
    )

    if accession:

        records.append({

            "Sample": line,

            "assembly_accession": accession

        })


genomes = pd.DataFrame(
    records
)


if genomes.empty:

    raise ValueError(
        "No GCA/GCF accessions were found "
        "in the genome list."
    )


# Remove duplicate accessions
genomes = (
    genomes
    .drop_duplicates(
        subset="assembly_accession",
        keep="first"
    )
    .reset_index(drop=True)
)


print(
    f"Unique assembly accessions: "
    f"{len(genomes)}"
)


# ============================================================
# 6. NCBI ASSEMBLY SUMMARY COLUMN POSITIONS
#
# The existing file has NO HEADER ROW.
#
# Only columns 0 and 30 are needed:
#
#   0  = assembly_accession
#   30 = contig_count
#
# Therefore we do NOT load all 38 columns.
# ============================================================

assembly_usecols = [
    0,
    30
]

assembly_names = [
    "assembly_accession",
    "No_of_contigs"
]


# ============================================================
# 7. READ EXISTING 1.9 GB ASSEMBLY SUMMARY IN CHUNKS
# ============================================================

print("\n" + "=" * 75)
print("READING EXISTING NCBI ASSEMBLY SUMMARY")
print("=" * 75)

print(
    "\nThe large assembly_summary file will be processed "
    "in chunks."
)

print(
    "No new assembly_summary file will be downloaded."
)


requested_accessions = set(
    genomes[
        "assembly_accession"
    ]
)


selected_chunks = []

rows_read = 0


for chunk in pd.read_csv(

    assembly_summary_file,

    sep="\t",

    header=None,

    names=assembly_names,

    usecols=assembly_usecols,

    comment="#",

    dtype=str,

    chunksize=200000

):

    rows_read += len(chunk)

    selected = chunk[
        chunk["assembly_accession"].isin(
            requested_accessions
        )
    ]

    if not selected.empty:

        selected_chunks.append(
            selected.copy()
        )


print(
    f"\nRows scanned: "
    f"{rows_read:,}"
)


if selected_chunks:

    assembly_selected = pd.concat(
        selected_chunks,
        ignore_index=True
    )

else:

    assembly_selected = pd.DataFrame(
        columns=assembly_names
    )


print(
    f"Requested accessions: "
    f"{len(requested_accessions)}"
)

print(
    f"Found in assembly summary: "
    f"{len(assembly_selected)}"
)


# ============================================================
# 8. CHECK MISSING ACCESSIONS
# ============================================================

found_accessions = set(
    assembly_selected[
        "assembly_accession"
    ]
)


missing_accessions = sorted(
    requested_accessions -
    found_accessions
)


if missing_accessions:

    print(
        "\nWARNING: These accessions were not found "
        "in the existing assembly summary:"
    )

    for accession in missing_accessions:

        print(
            "   ",
            accession
        )

else:

    print(
        "\nAll requested accessions were found "
        "in the existing assembly summary."
    )


# ============================================================
# 9. REMOVE DUPLICATES
# ============================================================

assembly_selected = (
    assembly_selected
    .drop_duplicates(
        subset="assembly_accession",
        keep="first"
    )
)


# ============================================================
# 10. VERIFY No_of_contigs
# ============================================================

print(
    "\nNo_of_contigs is taken from "
    "NCBI assembly_summary: contig_count."
)


example = assembly_selected[
    assembly_selected[
        "assembly_accession"
    ] == "GCA_036422055.1"
]


if not example.empty:

    print(
        "\nExample verification:"
    )

    print(
        example.to_string(
            index=False
        )
    )


# ============================================================
# 11. CHECK NCBI DATASETS
# ============================================================

def command_exists(command):

    result = subprocess.run(

        [
            "bash",
            "-lc",
            f"command -v {command}"
        ],

        capture_output=True,

        text=True
    )

    return (
        result.returncode == 0
    )


if not command_exists(
    "datasets"
):

    raise RuntimeError(
        "NCBI 'datasets' command was not found."
    )


if not command_exists(
    "dataformat"
):

    raise RuntimeError(
        "NCBI 'dataformat' command was not found."
    )


print(
    "\nNCBI datasets and dataformat found."
)


# ============================================================
# 12. CREATE NCBI DATASETS JSON REPORT
#
# Only requested accessions are queried.
# ============================================================

json_file = (
    output_dir /
    "NCBI_assembly_data_report.jsonl"
)


print(
    "\nQuerying NCBI Datasets..."
)


datasets_cmd = [

    "datasets",

    "summary",

    "genome",

    "accession"

]


datasets_cmd.extend(
    genomes[
        "assembly_accession"
    ].tolist()
)


datasets_cmd.append(
    "--as-json-lines"
)


with open(
    json_file,
    "w"
) as outfile:

    result = subprocess.run(

        datasets_cmd,

        stdout=outfile,

        stderr=subprocess.PIPE,

        text=True
    )


if result.returncode != 0:

    print(
        "\nNCBI Datasets error:"
    )

    print(
        result.stderr
    )

    raise RuntimeError(
        "NCBI Datasets query failed."
    )


print(
    f"\nDatasets report created:\n"
    f"{json_file}"
)


# ============================================================
# 13. EXTRACT NCBI DATASETS STATISTICS
#
# IMPORTANT:
# We deliberately DO NOT request:
#
#   assmstats-number-of-contigs
#
# because No_of_contigs comes from the existing
# assembly_summary_genbank.txt contig_count field.
# ============================================================

datasets_tsv = (
    output_dir /
    "NCBI_Datasets_statistics.tsv"
)


fields = (

    "accession,"

    "assmstats-total-sequence-len,"

    "assmstats-contig-n50,"

    "assmstats-contig-l50,"

    "assmstats-gc-percent,"

    "assmstats-genome-coverage,"

    "checkm-completeness,"

    "checkm-contamination"

)


print(
    "\nExtracting NCBI Datasets statistics..."
)


dataformat_cmd = [

    "dataformat",

    "tsv",

    "genome",

    "--inputfile",

    str(json_file),

    "--fields",

    fields

]


with open(
    datasets_tsv,
    "w"
) as outfile:

    result = subprocess.run(

        dataformat_cmd,

        stdout=outfile,

        stderr=subprocess.PIPE,

        text=True
    )


if result.returncode != 0:

    print(
        "\ndataformat error:"
    )

    print(
        result.stderr
    )

    raise RuntimeError(
        "dataformat failed."
    )


print(
    f"\nDatasets statistics saved:\n"
    f"{datasets_tsv}"
)


# ============================================================
# 14. READ DATASETS STATISTICS
# ============================================================

df_stats = pd.read_csv(

    datasets_tsv,

    sep="\t",

    dtype=str

)


print(
    "\nDatasets columns:"
)

print(
    df_stats.columns.tolist()
)


# ============================================================
# 15. RENAME DATASETS COLUMNS
# ============================================================

rename_map = {

    "Assembly Accession":
        "assembly_accession",

    "Assembly Stats Total Sequence Length":
        "TotalLength",

    "Assembly Stats Contig N50":
        "N50",

    "Assembly Stats Contig L50":
        "L50",

    "Assembly Stats GC Percent":
        "GC",

    "Assembly Stats Genome Coverage":
        "Coverage",

    "CheckM completeness":
        "CheckM_Completeness",

    "CheckM contamination":
        "CheckM_Contamination"

}


df_stats = df_stats.rename(
    columns=rename_map
)


# ============================================================
# 16. CHECK REQUIRED DATASET COLUMNS
# ============================================================

required_columns = [

    "assembly_accession",

    "TotalLength",

    "N50",

    "L50",

    "GC",

    "Coverage",

    "CheckM_Completeness",

    "CheckM_Contamination"

]


missing_columns = [

    column

    for column in required_columns

    if column not in df_stats.columns

]


if missing_columns:

    print(
        "\nMissing columns:"
    )

    print(
        missing_columns
    )

    raise ValueError(
        "Required NCBI Datasets fields are missing."
    )


# ============================================================
# 17. MERGE
# ============================================================

print(
    "\nMerging assembly-summary contig_count "
    "with NCBI Datasets statistics..."
)


final = genomes.merge(

    assembly_selected,

    on="assembly_accession",

    how="left"

)


final = final.merge(

    df_stats,

    on="assembly_accession",

    how="left"

)


# ============================================================
# 18. NUMERIC CONVERSION
# ============================================================

numeric_columns = [

    "No_of_contigs",

    "TotalLength",

    "N50",

    "L50",

    "GC",

    "Coverage",

    "CheckM_Completeness",

    "CheckM_Contamination"

]


for column in numeric_columns:

    final[column] = pd.to_numeric(

        final[column],

        errors="coerce"

    )


# ============================================================
# 19. FINAL COLUMN ORDER
# ============================================================

final_columns = [

    "Sample",

    "assembly_accession",

    "No_of_contigs",

    "TotalLength",

    "N50",

    "L50",

    "GC",

    "Coverage",

    "CheckM_Completeness",

    "CheckM_Contamination"

]


final = final[
    final_columns
]


# ============================================================
# 20. RESULT CHECK
# ============================================================

print(
    "\n" + "=" * 75
)

print(
    "RESULT CHECK"
)

print(
    "=" * 75
)


print(
    f"\nNumber of genomes in input: "
    f"{len(genomes)}"
)

print(
    f"Number of genomes in final: "
    f"{len(final)}"
)


for column in final_columns[2:]:

    print(
        f"Missing {column}: "
        f"{final[column].isna().sum()}"
    )


# ============================================================
# 21. VERIFY <100 CONTIGS
#
# IMPORTANT:
# This is ONLY a verification.
#
# No genome is removed.
# ============================================================

print(
    "\n" + "=" * 75
)

print(
    "CONTIG COUNT CHECK"
)

print(
    "=" * 75
)


ge100 = final[
    final["No_of_contigs"] >= 100
]


if ge100.empty:

    print(
        "\nAll genomes have No_of_contigs <100."
    )

else:

    print(
        f"\nWARNING: {len(ge100)} genome(s) "
        "have No_of_contigs >=100."
    )

    print(
        ge100[
            [
                "Sample",
                "assembly_accession",
                "No_of_contigs"
            ]
        ].to_string(
            index=False
        )
    )


# ============================================================
# 22. SAVE CSV
# ============================================================

output_csv = (
    output_dir /
    "selected_genomes_NCBI_QC_summary.csv"
)


final.to_csv(
    output_csv,
    index=False
)


# ============================================================
# 23. SAVE EXCEL
# ============================================================

output_excel = (
    output_dir /
    "selected_genomes_NCBI_QC_summary.xlsx"
)


final.to_excel(
    output_excel,
    index=False
)


# ============================================================
# 24. SAVE CONTIG VERIFICATION
# ============================================================

verification = final[

    [
        "Sample",
        "assembly_accession",
        "No_of_contigs"
    ]

].copy()


verification["Contig_check"] = (

    verification[
        "No_of_contigs"
    ]

    .apply(

        lambda x:

        "PASS (<100)"

        if pd.notna(x) and x < 100

        else (

            "MISSING"

            if pd.isna(x)

            else "CHECK (>=100)"

        )

    )

)


verification_file = (
    output_dir /
    "contig_count_verification.csv"
)


verification.to_csv(

    verification_file,

    index=False

)


# ============================================================
# 25. DISPLAY FINAL TABLE
# ============================================================

print(
    "\n" + "=" * 75
)

print(
    "FINAL TABLE"
)

print(
    "=" * 75
)


print(
    final.head(20).to_string(
        index=False
    )
)


# ============================================================
# 26. OUTPUT FILES
# ============================================================

print(
    "\n" + "=" * 75
)

print(
    "FILES CREATED"
)

print(
    "=" * 75
)


print(
    f"\nCSV:\n{output_csv}"
)

print(
    f"\nExcel:\n{output_excel}"
)

print(
    f"\nContig verification:\n{verification_file}"
)

print(
    f"\nDatasets raw report:\n{json_file}"
)

print(
    f"\nDatasets statistics:\n{datasets_tsv}"
)


print(
    "\n" + "=" * 75
)

print(
    "DONE"
)

print(
    "=" * 75
)
