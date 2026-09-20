#!/bin/bash

# ==============================
# Batch Unicycler Assembly Script
# ==============================

# File containing list of IDs
ID_LIST="PRJEB85317_5.txt"

# Directory containing your trimmed FASTQ files
FASTQ_DIR="../results/trimmed_fastq"

# Directory to store Unicycler assemblies
ASSEMBLY_DIR="../results/unicycler_output"

# Number of threads for Unicycler
THREADS=8

# Unicycler mode: choose from normal / conservative / bold
MODE="normal"

mkdir -p "$ASSEMBLY_DIR"

while read -r SAMPLE_ID; do

    echo "Running Unicycler for $SAMPLE_ID in $MODE mode..."

    # Input FASTQ files
    R1="$FASTQ_DIR/${SAMPLE_ID}_1_paired.fastq.gz"
    R2="$FASTQ_DIR/${SAMPLE_ID}_2_paired.fastq.gz"

    # Check if both FASTQs exist
    if [[ -f "$R1" && -f "$R2" ]]; then

        unicycler \
            -1 "$R1" \
            -2 "$R2" \
            -o "$ASSEMBLY_DIR/$SAMPLE_ID" \
            -t "$THREADS" \
            --mode "$MODE"

        echo "Finished $SAMPLE_ID. Output -> $ASSEMBLY_DIR/$SAMPLE_ID"

    else

        echo "FASTQ files not found for $SAMPLE_ID. Skipping..."

    fi

done < "$ID_LIST"

echo "All assemblies completed in $MODE mode."
