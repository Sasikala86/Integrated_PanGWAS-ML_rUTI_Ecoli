#!/bin/bash

# ================================
# Panaroo workflow for bacterial chromosomes
# ================================

# Set paths
INPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/chromosomes/prokka_output_chromosomes"
OUTPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/chromosomes/panaroo_output_chromosomes"
VALID_GFF_LIST="$OUTPUT_DIR/valid_gffs.txt"
SKIPPED_REPORT="$OUTPUT_DIR/skipped_gffs.txt"

# Create output folder
mkdir -p "$OUTPUT_DIR"

# Clear previous lists if any
> "$VALID_GFF_LIST"
> "$SKIPPED_REPORT"

echo "Starting bacterial chromosome Panaroo workflow..."
echo "Checking GFF files for valid CDS features..."

# Loop through each subfolder and find .gff files
for gff in "$INPUT_DIR"/*/*.gff; do
    if [ -f "$gff" ]; then
        if grep -q "CDS" "$gff"; then
            echo "$gff" >> "$VALID_GFF_LIST"
        else
            echo "$gff" >> "$SKIPPED_REPORT"
        fi
    fi
done

echo "Valid GFFs listed in $VALID_GFF_LIST"
echo "Skipped GFFs report in $SKIPPED_REPORT"

# Check if there are any valid GFFs
if [ ! -s "$VALID_GFF_LIST" ]; then
    echo "No valid GFFs found. Exiting."
    exit 1
fi

# Run Panaroo
echo "Running Panaroo on valid bacterial chromosomes..."
panaroo \
-i $(cat "$VALID_GFF_LIST") \
-o "$OUTPUT_DIR" \
--clean-mode strict \
--remove-invalid-genes \
-t 16 \
--verbose

echo "Bacterial chromosome Panaroo workflow completed!"
echo "Check $OUTPUT_DIR for gene presence/absence results and graphs."
