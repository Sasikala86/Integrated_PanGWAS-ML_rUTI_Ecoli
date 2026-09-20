#!/bin/bash

# Input file containing dataset IDs
ID_LIST="PRJEB85317_Quast.txt"

# Paths
ASSEMBLY_BASE="../results/unicycler_output"
OUTPUT_BASE="../results/quast_output"

# Number of threads for QUAST
THREADS=8

# Activate assembly environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate assembly_env

# Loop through each dataset ID
while read -r SAMPLE; do

    # Skip empty lines
    if [ -z "$SAMPLE" ]; then
        continue
    fi

    ASSEMBLY="$ASSEMBLY_BASE/$SAMPLE/assembly.fasta"
    OUTDIR="$OUTPUT_BASE/$SAMPLE"

    if [ -f "$ASSEMBLY" ]; then

        echo "Running QUAST for $SAMPLE..."

        mkdir -p "$OUTDIR"

        quast.py "$ASSEMBLY" \
            -o "$OUTDIR" \
            --threads "$THREADS"

    else

        echo "WARNING: Assembly not found for $SAMPLE ($ASSEMBLY)"

    fi

done < "$ID_LIST"

echo "All QUAST runs completed."
