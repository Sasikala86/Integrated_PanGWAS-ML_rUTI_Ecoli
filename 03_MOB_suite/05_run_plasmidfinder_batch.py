from pathlib import Path
import subprocess
import json
import pandas as pd
import re
import shutil
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path("/home/sastra/Sasikala/Pangenome/mob_output")

PLASMID_FEATURES = ROOT / "plasmid_features_summary.csv"

OUTPUT = ROOT / "PlasmidFinder_validation"

JSON_DIR = OUTPUT / "json_results"
LOG_DIR = OUTPUT / "logs"

OUTPUT.mkdir(exist_ok=True)
JSON_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

PLASMIDFINDER_DB = Path(
    "/home/sastra/Sasikala/databases/plasmidfinder_db"
)

# PlasmidFinder settings
METHOD = "blastn"
MIN_COV = "0.6"
THRESHOLD = "0.9"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def phenotype_from_genome(genome_id):

    if genome_id.startswith("rUTI_"):
        return "rUTI"

    if genome_id.startswith("UTI_"):
        return "UTI"

    return "Unknown"


def clean_value(x):

    if pd.isna(x):
        return ""

    x = str(x).strip()

    if x in {"-", "nan", "None", "NA", "N/A", ""}:
        return ""

    return x


def normalize_replicon_name(x):

    x = clean_value(x)

    if not x:
        return ""

    # Remove common database formatting
    x = re.sub(r"_\d+$", "", x)

    # Some MOB-suite fields may contain multiple replicons
    return x


def extract_mob_replicons(value):

    value = clean_value(value)

    if not value:
        return []

    # MOB-suite commonly separates multiple replicons by /
    # but retain individual names for comparison.
    parts = re.split(r"[/;,]+", value)

    return [
        normalize_replicon_name(p)
        for p in parts
        if normalize_replicon_name(p)
    ]


def extract_plasmidfinder_json(json_file):

    result_rows = []

    with open(json_file, "r") as f:
        data = json.load(f)

    seq_regions = data.get("seq_regions", {})

    software_version = data.get(
        "software_version", ""
    )

    databases = data.get("databases", {})

    database_name = ""
    database_version = ""

    if databases:

        first_db = next(iter(databases.values()))

        database_name = first_db.get(
            "database_name", ""
        )

        database_version = first_db.get(
            "database_version", ""
        )

    if not seq_regions:

        result_rows.append({
            "plasmidfinder_positive": False,
            "plasmidfinder_replicon": "",
            "plasmidfinder_identity": "",
            "plasmidfinder_coverage": "",
            "plasmidfinder_ref_id": "",
            "plasmidfinder_ref_acc": "",
            "plasmidfinder_software_version": software_version,
            "plasmidfinder_database": database_name,
            "plasmidfinder_database_version": database_version
        })

        return result_rows

    for region in seq_regions.values():

        result_rows.append({
            "plasmidfinder_positive": True,
            "plasmidfinder_replicon": clean_value(
                region.get("name", "")
            ),
            "plasmidfinder_identity": region.get(
                "identity", ""
            ),
            "plasmidfinder_coverage": region.get(
                "coverage", ""
            ),
            "plasmidfinder_ref_id": clean_value(
                region.get("ref_id", "")
            ),
            "plasmidfinder_ref_acc": clean_value(
                region.get("ref_acc", "")
            ),
            "plasmidfinder_software_version": software_version,
            "plasmidfinder_database": database_name,
            "plasmidfinder_database_version": database_version
        })

    return result_rows


def compare_replicons(mob_replicon, pf_replicons):

    mob_list = extract_mob_replicons(mob_replicon)

    pf_list = [
        normalize_replicon_name(x)
        for x in pf_replicons
        if normalize_replicon_name(x)
    ]

    if not mob_list and not pf_list:
        return "both_negative"

    if mob_list and not pf_list:
        return "MOB-suite_only"

    if not mob_list and pf_list:
        return "PlasmidFinder_only"

    # Exact or partial overlap
    for mob in mob_list:
        for pf in pf_list:

            if mob.lower() == pf.lower():
                return "concordant"

            # Handle cases where names differ slightly
            if mob.lower() in pf.lower() or pf.lower() in mob.lower():
                return "concordant_partial"

    return "discordant"


# ============================================================
# LOAD MOB-SUITE PLASMID FEATURES
# ============================================================

print("=" * 70)
print("PLASMIDFINDER ORTHOGONAL VALIDATION")
print("=" * 70)

if not PLASMID_FEATURES.exists():

    raise FileNotFoundError(
        f"\nMissing file:\n{PLASMID_FEATURES}"
    )

mob = pd.read_csv(
    PLASMID_FEATURES,
    dtype=str
)

mob.columns = [
    str(c).strip()
    for c in mob.columns
]

print(
    f"\nMOB-suite plasmid reconstructions in summary: {len(mob)}"
)

# Check required columns
required = [
    "sample_folder",
    "plasmid_id"
]

missing = [
    c for c in required
    if c not in mob.columns
]

if missing:

    raise ValueError(
        f"Missing required columns in plasmid_features_summary.csv: {missing}"
    )

# ============================================================
# MAP MOB-SUITE RECORDS
# ============================================================

mob_lookup = {}

for _, row in mob.iterrows():

    key = (
        clean_value(row["sample_folder"]),
        clean_value(row["plasmid_id"])
    )

    mob_lookup[key] = row.to_dict()

# ============================================================
# FIND ALL PLASMID FASTA FILES
# ============================================================

plasmid_fastas = sorted(
    ROOT.glob(
        "**/plasmid_*.fasta"
    )
)

print(
    f"Plasmid FASTA files found: {len(plasmid_fastas)}"
)

# ============================================================
# RUN PLASMIDFINDER
# ============================================================

all_results = []
failed_results = []

processed = 0

for fasta in plasmid_fastas:

    genome_dir = fasta.parent

    genome_id = genome_dir.name

    phenotype = phenotype_from_genome(
        genome_id
    )

    # --------------------------------------------------------
    # Extract plasmid ID
    # Expected:
    # plasmid_AA279_rUTI_ERR4556422.fasta
    # --------------------------------------------------------

    filename = fasta.name

    plasmid_name = filename[
        len("plasmid_") :
        -len(".fasta")
    ]

    # Remove genome ID suffix
    suffix = "_" + genome_id

    if plasmid_name.endswith(suffix):

        plasmid_id = plasmid_name[
            :-len(suffix)
        ]

    else:

        # Fallback
        plasmid_id = plasmid_name.split(
            "_", 1
        )[0]

    key = (
        genome_id,
        plasmid_id
    )

    mob_row = mob_lookup.get(
        key,
        {}
    )

    mob_rep = clean_value(
        mob_row.get(
            "rep_type",
            mob_row.get(
                "rep_type(s)",
                ""
            )
        )
    )

    # --------------------------------------------------------
    # Output directory
    # --------------------------------------------------------

    safe_genome = re.sub(
        r"[^A-Za-z0-9_.-]",
        "_",
        genome_id
    )

    safe_plasmid = re.sub(
        r"[^A-Za-z0-9_.-]",
        "_",
        plasmid_id
    )

    result_dir = (
        JSON_DIR /
        safe_genome /
        safe_plasmid
    )

    result_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    json_file = (
        result_dir /
        f"{safe_plasmid}_plasmidfinder.json"
    )

    log_file = (
        LOG_DIR /
        f"{safe_genome}__{safe_plasmid}.log"
    )

    # --------------------------------------------------------
    # Run only if JSON does not already exist
    # --------------------------------------------------------

    if not json_file.exists():

        cmd = [
            "python",
            "-m",
            "plasmidfinder",

            "-i",
            str(fasta),

            "-o",
            str(result_dir),

            "-mp",
            METHOD,

            "-p",
            str(PLASMIDFINDER_DB),

            "-l",
            MIN_COV,

            "-t",
            THRESHOLD,

            "-j",
            str(json_file)
        ]

        try:

            with open(
                log_file,
                "w"
            ) as log:

                completed = subprocess.run(
                    cmd,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True
                )

            if completed.returncode != 0:

                failed_results.append({
                    "genome_id": genome_id,
                    "phenotype": phenotype,
                    "plasmid_id": plasmid_id,
                    "fasta": str(fasta),
                    "reason": (
                        f"PlasmidFinder exit code "
                        f"{completed.returncode}"
                    ),
                    "log": str(log_file)
                })

                continue

        except Exception as e:

            failed_results.append({
                "genome_id": genome_id,
                "phenotype": phenotype,
                "plasmid_id": plasmid_id,
                "fasta": str(fasta),
                "reason": str(e),
                "log": str(log_file)
            })

            continue

    # --------------------------------------------------------
    # Parse JSON
    # --------------------------------------------------------

    try:

        pf_hits = extract_plasmidfinder_json(
            json_file
        )

    except Exception as e:

        failed_results.append({
            "genome_id": genome_id,
            "phenotype": phenotype,
            "plasmid_id": plasmid_id,
            "fasta": str(fasta),
            "reason": (
                "JSON parsing error: "
                + str(e)
            ),
            "log": str(log_file)
        })

        continue

    # --------------------------------------------------------
    # Combine hits
    # --------------------------------------------------------

    pf_names = [
        x["plasmidfinder_replicon"]
        for x in pf_hits
        if x["plasmidfinder_positive"]
    ]

    pf_names_clean = [
        normalize_replicon_name(x)
        for x in pf_names
        if normalize_replicon_name(x)
    ]

    comparison = compare_replicons(
        mob_rep,
        pf_names_clean
    )

    # --------------------------------------------------------
    # One output row per PlasmidFinder hit
    # --------------------------------------------------------

    if pf_hits:

        for hit in pf_hits:

            result = {
                "genome_id": genome_id,
                "phenotype": phenotype,
                "plasmid_id": plasmid_id,
                "plasmid_fasta": str(fasta),

                "mob_rep_type": mob_rep,

                "mob_predicted_mobility":
                    clean_value(
                        mob_row.get(
                            "predicted_mobility",
                            ""
                        )
                    ),

                "mob_relaxase_type":
                    clean_value(
                        mob_row.get(
                            "relaxase_type",
                            ""
                        )
                    ),

                "mob_mpf_type":
                    clean_value(
                        mob_row.get(
                            "mpf_type",
                            ""
                        )
                    ),

                "mob_orit_type":
                    clean_value(
                        mob_row.get(
                            "orit_type",
                            ""
                        )
                    ),

                "plasmidfinder_positive":
                    hit[
                        "plasmidfinder_positive"
                    ],

                "plasmidfinder_replicon":
                    hit[
                        "plasmidfinder_replicon"
                    ],

                "plasmidfinder_identity":
                    hit[
                        "plasmidfinder_identity"
                    ],

                "plasmidfinder_coverage":
                    hit[
                        "plasmidfinder_coverage"
                    ],

                "plasmidfinder_ref_id":
                    hit[
                        "plasmidfinder_ref_id"
                    ],

                "plasmidfinder_ref_acc":
                    hit[
                        "plasmidfinder_ref_acc"
                    ],

                "plasmidfinder_software_version":
                    hit[
                        "plasmidfinder_software_version"
                    ],

                "plasmidfinder_database_version":
                    hit[
                        "plasmidfinder_database_version"
                    ],

                "replicon_comparison":
                    comparison,

                "mob_reconstruction_found":
                    bool(mob_row),

                "run_status":
                    "success"
            }

            all_results.append(result)

    else:

        all_results.append({

            "genome_id": genome_id,
            "phenotype": phenotype,
            "plasmid_id": plasmid_id,
            "plasmid_fasta": str(fasta),

            "mob_rep_type": mob_rep,

            "mob_predicted_mobility":
                clean_value(
                    mob_row.get(
                        "predicted_mobility",
                        ""
                    )
                ),

            "mob_relaxase_type":
                clean_value(
                    mob_row.get(
                        "relaxase_type",
                        ""
                    )
                ),

            "mob_mpf_type":
                clean_value(
                    mob_row.get(
                        "mpf_type",
                        ""
                    )
                ),

            "mob_orit_type":
                clean_value(
                    mob_row.get(
                        "orit_type",
                        ""
                    )
                ),

            "plasmidfinder_positive":
                False,

            "plasmidfinder_replicon":
                "",

            "plasmidfinder_identity":
                "",

            "plasmidfinder_coverage":
                "",

            "plasmidfinder_ref_id":
                "",

            "plasmidfinder_ref_acc":
                "",

            "plasmidfinder_software_version":
                "3.0.3",

            "plasmidfinder_database_version":
                "2.2.0",

            "replicon_comparison":
                comparison,

            "mob_reconstruction_found":
                bool(mob_row),

            "run_status":
                "success"
        })

    processed += 1

    if processed % 50 == 0:

        print(
            f"Processed {processed} plasmid FASTA files..."
        )

# ============================================================
# SAVE RAW RESULTS
# ============================================================

results_df = pd.DataFrame(
    all_results
)

failed_df = pd.DataFrame(
    failed_results
)

raw_results_file = (
    OUTPUT /
    "01_all_plasmidfinder_results.csv"
)

results_df.to_csv(
    raw_results_file,
    index=False
)

failed_file = (
    OUTPUT /
    "02_failed_plasmidfinder_runs.csv"
)

failed_df.to_csv(
    failed_file,
    index=False
)

# ============================================================
# CREATE ONE-ROW-PER-PLASMID SUMMARY
# ============================================================

if len(results_df) > 0:

    summary_rows = []

    group_cols = [
        "genome_id",
        "phenotype",
        "plasmid_id"
    ]

    for key, group in results_df.groupby(
        group_cols,
        dropna=False
    ):

        genome_id, phenotype, plasmid_id = key

        first = group.iloc[0]

        pf_positive = group[
            "plasmidfinder_positive"
        ].any()

        pf_replicons = sorted(
            set(
                group[
                    "plasmidfinder_replicon"
                ]
                .astype(str)
                .map(str.strip)
                .replace("", pd.NA)
                .dropna()
            )
        )

        pf_identities = pd.to_numeric(
            group[
                "plasmidfinder_identity"
            ],
            errors="coerce"
        ).dropna()

        pf_coverages = pd.to_numeric(
            group[
                "plasmidfinder_coverage"
            ],
            errors="coerce"
        ).dropna()

        summary_rows.append({

            "genome_id": genome_id,

            "phenotype": phenotype,

            "plasmid_id": plasmid_id,

            "mob_rep_type":
                first["mob_rep_type"],

            "plasmidfinder_positive":
                pf_positive,

            "plasmidfinder_replicons":
                ";".join(pf_replicons),

            "plasmidfinder_hit_count":
                len(pf_replicons),

            "max_identity":
                pf_identities.max()
                if len(pf_identities)
                else "",

            "max_coverage":
                pf_coverages.max()
                if len(pf_coverages)
                else "",

            "replicon_comparison":
                first["replicon_comparison"],

            "mob_predicted_mobility":
                first[
                    "mob_predicted_mobility"
                ],

            "mob_relaxase_type":
                first[
                    "mob_relaxase_type"
                ],

            "mob_mpf_type":
                first[
                    "mob_mpf_type"
                ],

            "mob_orit_type":
                first[
                    "mob_orit_type"
                ]
        })

    summary_df = pd.DataFrame(
        summary_rows
    )

else:

    summary_df = pd.DataFrame()

summary_file = (
    OUTPUT /
    "03_one_row_per_plasmid.csv"
)

summary_df.to_csv(
    summary_file,
    index=False
)

# ============================================================
# GROUP-LEVEL VALIDATION SUMMARY
# ============================================================

if len(summary_df) > 0:

    group_rows = []

    for phenotype in [
        "rUTI",
        "UTI",
        "Unknown"
    ]:

        sub = summary_df[
            summary_df["phenotype"] ==
            phenotype
        ]

        if len(sub) == 0:
            continue

        total = len(sub)

        pf_positive = (
            sub[
                "plasmidfinder_positive"
            ]
            .astype(bool)
            .sum()
        )

        group_rows.append({

            "phenotype":
                phenotype,

            "plasmid_reconstructions":
                total,

            "PlasmidFinder_positive":
                pf_positive,

            "PlasmidFinder_positive_percent":
                100 * pf_positive / total,

            "concordant":
                (
                    sub[
                        "replicon_comparison"
                    ] ==
                    "concordant"
                ).sum(),

            "concordant_partial":
                (
                    sub[
                        "replicon_comparison"
                    ] ==
                    "concordant_partial"
                ).sum(),

            "MOB_suite_only":
                (
                    sub[
                        "replicon_comparison"
                    ] ==
                    "MOB-suite_only"
                ).sum(),

            "PlasmidFinder_only":
                (
                    sub[
                        "replicon_comparison"
                    ] ==
                    "PlasmidFinder_only"
                ).sum(),

            "discordant":
                (
                    sub[
                        "replicon_comparison"
                    ] ==
                    "discordant"
                ).sum(),

            "both_negative":
                (
                    sub[
                        "replicon_comparison"
                    ] ==
                    "both_negative"
                ).sum()
        })

    group_df = pd.DataFrame(
        group_rows
    )

else:

    group_df = pd.DataFrame()

group_file = (
    OUTPUT /
    "04_group_level_validation.csv"
)

group_df.to_csv(
    group_file,
    index=False
)

# ============================================================
# PHENOTYPE COMPARISON OF REPLICON CONCORDANCE
# ============================================================

if len(summary_df) > 0:

    comparison_table = pd.crosstab(
        summary_df["phenotype"],
        summary_df["replicon_comparison"]
    )

    comparison_table.to_csv(
        OUTPUT /
        "05_replicon_concordance_by_phenotype.csv"
    )

# ============================================================
# PLASMIDFINDER POSITIVE HITS
# ============================================================

if len(results_df) > 0:

    positive_hits = results_df[
        results_df[
            "plasmidfinder_positive"
        ] == True
    ]

    positive_hits.to_csv(
        OUTPUT /
        "06_PlasmidFinder_positive_hits.csv",
        index=False
    )

# ============================================================
# MOB-SUITE ONLY REPLICONS
# ============================================================

if len(summary_df) > 0:

    mob_only = summary_df[
        summary_df[
            "replicon_comparison"
        ] == "MOB-suite_only"
    ]

    mob_only.to_csv(
        OUTPUT /
        "07_MOB_suite_only_replicons.csv",
        index=False
    )

# ============================================================
# DISCORDANT REPLICONS
# ============================================================

if len(summary_df) > 0:

    discordant = summary_df[
        summary_df[
            "replicon_comparison"
        ] == "discordant"
    ]

    discordant.to_csv(
        OUTPUT /
        "08_discordant_replicons.csv",
        index=False
    )

# ============================================================
# FINAL TEXT REPORT
# ============================================================

report_file = (
    OUTPUT /
    "09_PlasmidFinder_validation_report.txt"
)

with open(
    report_file,
    "w"
) as f:

    f.write(
        "PLASMIDFINDER ORTHOGONAL VALIDATION REPORT\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    )

    f.write(
        "PlasmidFinder version: 3.0.3\n"
    )

    f.write(
        "PlasmidFinder database: 2.2.0\n"
    )

    f.write(
        "Method: BLASTN\n"
    )

    f.write(
        f"Minimum coverage: {MIN_COV}\n"
    )

    f.write(
        f"Identity threshold: {THRESHOLD}\n\n"
    )

    f.write(
        f"Plasmid FASTA files found: "
        f"{len(plasmid_fastas)}\n"
    )

    f.write(
        f"Successfully processed: "
        f"{processed}\n"
    )

    f.write(
        f"Failed: "
        f"{len(failed_results)}\n\n"
    )

    if len(group_df) > 0:

        f.write(
            group_df.to_string(
                index=False
            )
        )

        f.write("\n")

print("\n" + "=" * 70)
print("BATCH ANALYSIS COMPLETE")
print("=" * 70)

print(
    f"\nPlasmid FASTA files found: {len(plasmid_fastas)}"
)

print(
    f"Successfully processed: {processed}"
)

print(
    f"Failed: {len(failed_results)}"
)

print(
    f"\nResults directory:\n{OUTPUT}"
)

print("\nMain files:")

print(
    "01_all_plasmidfinder_results.csv"
)

print(
    "02_failed_plasmidfinder_runs.csv"
)

print(
    "03_one_row_per_plasmid.csv"
)

print(
    "04_group_level_validation.csv"
)

print(
    "05_replicon_concordance_by_phenotype.csv"
)

print(
    "06_PlasmidFinder_positive_hits.csv"
)

print(
    "07_MOB_suite_only_replicons.csv"
)

print(
    "08_discordant_replicons.csv"
)

print(
    "09_PlasmidFinder_validation_report.txt"
)
