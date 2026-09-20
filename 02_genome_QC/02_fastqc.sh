#!/bin/bash

# Input folder containing all your FASTQ files
INPUT_DIR="../../01_data/fastq_files"

# Output folder for FastQC reports
OUTPUT_DIR="../results/fastqc_reports"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC on all FASTQ and FASTQ.GZ files
for FILE in "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fastq.gz; do

    if [ -f "$FILE" ]; then
        echo "Running FastQC on $FILE ..."
        fastqc "$FILE" -o "$OUTPUT_DIR" --threads 4
    fi

done

echo "FastQC analysis completed. Reports are in $OUTPUT_DIR"
