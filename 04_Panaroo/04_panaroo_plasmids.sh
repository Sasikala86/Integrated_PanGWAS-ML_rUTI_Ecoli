#!/bin/bash
# panaroo_plasmids.sh
# One-step plasmid Panaroo workflow

set -e

PROKKA_DIR="/home/sastra/Sasikala/Pangenome/mob_output/plasmids/prokka_output_plasmids"
INPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/plasmids/panaroo_input_plasmids"
OUTPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/plasmids/panaroo_output_plasmids_cleaned"

VALID_GFF_LIST="$OUTPUT_DIR/valid_gffs.txt"
SKIPPED_REPORT="$OUTPUT_DIR/skipped_plasmids_report.txt"

mkdir -p "$OUTPUT_DIR"

echo "Starting plasmid Panaroo workflow..."

# Step 1: Clean .ffn FASTA files
echo "Cleaning .ffn FASTA files in $PROKKA_DIR ..."
for ffn in "$PROKKA_DIR"/*/*.ffn; do
    if [ -f "$ffn" ]; then
        # remove leading comment lines
        sed -i '/^>/!b;/^>/!d' "$ffn" 2>/dev/null || true
    fi
done

# Step 2: Check GFF files
echo "Checking GFF files in $INPUT_DIR for missing attributes..."
> "$VALID_GFF_LIST"
> "$SKIPPED_REPORT"

for gff in "$INPUT_DIR"/*.gff; do
    if [ ! -f "$gff" ]; then
        continue
    fi

    # Check for mandatory attributes (ID and gene)
    if grep -q "ID=" "$gff" && grep -q "gene=" "$gff"; then
        echo "$gff" >> "$VALID_GFF_LIST"
    else
        echo "$gff" >> "$SKIPPED_REPORT"
    fi
done

echo "Valid plasmids listed in $VALID_GFF_LIST"
echo "Skipped plasmids report in $SKIPPED_REPORT"

# Remove empty lines
sed -i '/^$/d' "$VALID_GFF_LIST"

# Step 3: Run Panaroo if valid GFFs exist
if [ -s "$VALID_GFF_LIST" ]; then
    echo "Running Panaroo..."
    GFF_LIST=$(tr '\n' ' ' < "$VALID_GFF_LIST")

    panaroo -i $GFF_LIST \
            -o "$OUTPUT_DIR" \
            -t 8 \
            --clean-mode strict \
            --remove-invalid-genes

    echo "Plasmid Panaroo workflow completed."
    echo "Check $OUTPUT_DIR for results and reports."
else
    echo "No valid GFFs found. Exiting."
fi

