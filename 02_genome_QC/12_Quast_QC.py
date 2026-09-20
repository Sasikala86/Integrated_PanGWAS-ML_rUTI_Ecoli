#!/usr/bin/env python3

# ============================================================
# QUAST CONTIG COUNT EXTRACTION
#
# Purpose:
#   Extract "# contigs" from QUAST report.html and report.pdf
#   files.
#
# Genome ID:
#   Parent folder name containing the QUAST report.
#
# Output:
#   Genome_Contig_Counts.xlsx
#
# Columns:
#   Genome_ID
#   No_of_contigs
#   Source
#
# IMPORTANT:
#   This script ONLY extracts contig counts.
#   It does NOT perform genome selection/filtering.
# ============================================================

from pathlib import Path
import argparse
import pandas as pd
from bs4 import BeautifulSoup
import re
from pypdf import PdfReader


# ============================================================
# 1. COMMAND-LINE ARGUMENTS
# ============================================================

parser = argparse.ArgumentParser(
    description="Extract contig counts from QUAST HTML/PDF reports."
)

parser.add_argument(
    "--input-dir",
    required=True,
    help="Directory containing genome folders with QUAST reports."
)

parser.add_argument(
    "--output",
    default="Genome_Contig_Counts.xlsx",
    help="Output Excel filename."
)

args = parser.parse_args()

main_dir = Path(args.input_dir)
output_file = Path(args.output)


# ============================================================
# 2. CHECK INPUT DIRECTORY
# ============================================================

if not main_dir.exists():
    raise FileNotFoundError(
        f"Input directory not found:\n{main_dir}"
    )

if not main_dir.is_dir():
    raise NotADirectoryError(
        f"Input path is not a directory:\n{main_dir}"
    )


print("=" * 60)
print("QUAST REPORT SEARCH")
print("=" * 60)

print(f"Input directory: {main_dir}")
print()


# ============================================================
# 3. FIND ALL HTML AND PDF REPORTS
# ============================================================

html_reports = list(
    main_dir.rglob("report.html")
)

pdf_reports = list(
    main_dir.rglob("report.pdf")
)

print(
    f"HTML reports found : {len(html_reports)}"
)

print(
    f"PDF reports found  : {len(pdf_reports)}"
)

print()


# ============================================================
# 4. FUNCTION TO EXTRACT CONTIGS FROM HTML
# ============================================================

def extract_from_html(report):

    html = report.read_text(
        encoding="utf-8",
        errors="ignore"
    )

    soup = BeautifulSoup(
        html,
        "html.parser"
    )

    texts = [
        x.strip()
        for x in soup.stripped_strings
        if x.strip()
    ]

    for i, text in enumerate(texts):

        if text == "# contigs":

            for next_text in texts[i + 1:i + 5]:

                if next_text.isdigit():

                    return int(next_text)

    return None


# ============================================================
# 5. FUNCTION TO EXTRACT CONTIGS FROM PDF
# ============================================================

def extract_from_pdf(report):

    reader = PdfReader(
        str(report)
    )

    text = ""

    for page in reader.pages:

        page_text = page.extract_text()

        if page_text:

            text += "\n" + page_text

    # --------------------------------------------------------
    # Exact "# contigs" line
    #
    # Example:
    #
    # # contigs       65
    #
    # This avoids extracting:
    #
    # # contigs (>= 0 bp)
    # # contigs (>= 1000 bp)
    # --------------------------------------------------------

    pattern = (
        r"(?m)^\s*#\s*contigs\s+(\d+)\s*$"
    )

    match = re.search(
        pattern,
        text
    )

    if match:

        return int(
            match.group(1)
        )

    # --------------------------------------------------------
    # Alternative PDF formatting
    # --------------------------------------------------------

    pattern2 = (
        r"(?m)^\s*#\s*contigs\s*\n\s*(\d+)\s*$"
    )

    match2 = re.search(
        pattern2,
        text
    )

    if match2:

        return int(
            match2.group(1)
        )

    return None


# ============================================================
# 6. PROCESS REPORTS
# ============================================================

results = []
failed = []


# ============================================================
# 6A. PROCESS HTML REPORTS
# ============================================================

for report in html_reports:

    if report.parent == main_dir:

        print(
            f"Skipping main-directory report: {report}"
        )

        continue

    genome_id = report.parent.name

    try:

        contigs = extract_from_html(
            report
        )

        if contigs is not None:

            results.append({

                "Genome_ID": genome_id,

                "No_of_contigs": contigs,

                "Source": "HTML"

            })

            print(
                f"{genome_id:20s} "
                f"-> {contigs:4d} contigs "
                f"(HTML)"
            )

        else:

            failed.append({

                "Genome_ID": genome_id,

                "Report": str(report),

                "Source": "HTML",

                "Reason": "Could not find # contigs"

            })

            print(
                f"WARNING: Could not extract "
                f"# contigs -> {genome_id} (HTML)"
            )

    except Exception as e:

        failed.append({

            "Genome_ID": genome_id,

            "Report": str(report),

            "Source": "HTML",

            "Reason": str(e)

        })

        print(
            f"ERROR: {genome_id} (HTML) -> {e}"
        )


# ============================================================
# 6B. PROCESS PDF REPORTS
# ============================================================

for report in pdf_reports:

    if report.parent == main_dir:

        print(
            f"Skipping main-directory PDF: {report}"
        )

        continue

    genome_id = report.parent.name

    try:

        contigs = extract_from_pdf(
            report
        )

        if contigs is not None:

            results.append({

                "Genome_ID": genome_id,

                "No_of_contigs": contigs,

                "Source": "PDF"

            })

            print(
                f"{genome_id:20s} "
                f"-> {contigs:4d} contigs "
                f"(PDF)"
            )

        else:

            failed.append({

                "Genome_ID": genome_id,

                "Report": str(report),

                "Source": "PDF",

                "Reason": "Could not find # contigs"

            })

            print(
                f"WARNING: Could not extract "
                f"# contigs -> {genome_id} (PDF)"
            )

    except Exception as e:

        failed.append({

            "Genome_ID": genome_id,

            "Report": str(report),

            "Source": "PDF",

            "Reason": str(e)

        })

        print(
            f"ERROR: {genome_id} (PDF) -> {e}"
        )


# ============================================================
# 7. REMOVE DUPLICATES
#
# If both HTML and PDF exist for the same genome,
# keep the first successfully extracted result.
# ============================================================

if results:

    df = pd.DataFrame(
        results
    )

    df = (
        df
        .drop_duplicates(
            subset="Genome_ID",
            keep="first"
        )
        .sort_values("Genome_ID")
        .reset_index(drop=True)
    )

else:

    df = pd.DataFrame(
        columns=[
            "Genome_ID",
            "No_of_contigs",
            "Source"
        ]
    )


# ============================================================
# 8. FAILED REPORTS
# ============================================================

failed_df = pd.DataFrame(
    failed
)


# ============================================================
# 9. CREATE OUTPUT DIRECTORY
# ============================================================

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 10. SAVE EXCEL
# ============================================================

with pd.ExcelWriter(
    output_file,
    engine="openpyxl"
) as writer:

    df.to_excel(
        writer,
        sheet_name="Contig_Counts",
        index=False
    )

    if not failed_df.empty:

        failed_df.to_excel(
            writer,
            sheet_name="Not_Extracted",
            index=False
        )


# ============================================================
# 11. SUMMARY
# ============================================================

print()

print("=" * 60)
print("EXTRACTION COMPLETED")
print("=" * 60)

print(
    f"HTML reports found       : {len(html_reports)}"
)

print(
    f"PDF reports found        : {len(pdf_reports)}"
)

print(
    f"Unique genomes extracted : {len(df)}"
)

print(
    f"Reports not extracted    : {len(failed_df)}"
)

print()

print(
    "Output Excel:"
)

print(
    output_file
)

print("=" * 60)
