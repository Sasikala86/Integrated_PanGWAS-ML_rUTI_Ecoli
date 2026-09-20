#!/bin/bash

# Input folder containing all your trimmed FASTQ files
INPUT_DIR="../results/trimmed_fastq"

# Output folder for FastQC reports
OUTPUT_DIR="../results/posttrim_fastqc_reports"

# Number of threads
THREADS=4

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over all FASTQ and FASTQ.GZ paired files
for file in "$INPUT_DIR"/*_paired.fastq*; do

    if [ -f "$file" ]; then
        fastqc -t "$THREADS" "$file" -o "$OUTPUT_DIR"
        echo "FastQC completed for $file"
    fi

done

echo "All post-trimming FastQC analysis completed."
echo "Reports are in $OUTPUT_DIR"
