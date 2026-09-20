#!/bin/bash

# Base directory containing MOBsuite output folders
BASE_DIR="/home/sastra/Sasikala/Pangenome/mob_output"

# Destination folder to store all plasmids
PLASMID_DIR="${BASE_DIR}/plasmids"

# Log file for missing plasmids
LOG_FILE="${BASE_DIR}/no_plasmid_datasets.txt"

# Create destination folder if not exists
mkdir -p "$PLASMID_DIR"

echo "🚀 Collecting plasmid FASTA files..."
echo "---------------------------------------"

# Loop through all subdirectories
for folder in "$BASE_DIR"/*/; do
    sample=$(basename "$folder")

    # Find plasmid files in the folder
    plasmid_files=("$folder"plasmid_*.fasta)

    if compgen -G "$folder"plasmid_*.fasta > /dev/null; then
        echo "📦 Copying plasmids from: $sample"
        for fasta in "${plasmid_files[@]}"; do
            [[ -e "$fasta" ]] || continue
            cp "$fasta" "$PLASMID_DIR/"
            echo "   ✅ Copied: $(basename "$fasta")"
        done
    else
        echo "⚠️  No plasmids found in: $sample"
        echo "$sample" >> "$LOG_FILE"
    fi
done

echo "---------------------------------------"
echo "🎯 All plasmids collected into: $PLASMID_DIR"
echo "🗒️  Datasets without plasmids are listed in: $LOG_FILE"
