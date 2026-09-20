#!/bin/bash

# File containing list of SRA IDs
ID_LIST="SRR_ids.txt"

# Output folder for FASTQ files
OUTPUT_DIR="./fastq_files"
mkdir -p "$OUTPUT_DIR"

# Loop through each SRA ID in the list
while read -r SRR_ID; do

    # Skip empty lines
    if [ -z "$SRR_ID" ]; then
        continue
    fi

    echo "===== Processing $SRR_ID ====="

    # Download .sra file
    prefetch "$SRR_ID"

    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to prefetch $SRR_ID"
        continue
    fi

    # Convert .sra to FASTQ files
    fasterq-dump "$SRR_ID" --split-files -O "$OUTPUT_DIR"

    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to convert $SRR_ID to FASTQ"
        continue
    fi

    echo "$SRR_ID downloaded and converted successfully."

done < "$ID_LIST"
