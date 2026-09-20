
#!/bin/bash

# Base directory containing MOBsuite output folders
BASE_DIR="/home/sastra/Sasikala/Pangenome/mob_output"

# Log file to record folders without plasmids
LOG_FILE="${BASE_DIR}/no_plasmid_datasets.txt"

# Clear old log file
> "$LOG_FILE"

echo "🚀 Starting plasmid renaming process..."
echo "---------------------------------------"

# Loop through all subdirectories in mob_output
for folder in "$BASE_DIR"/*/; do
    sample=$(basename "$folder")
    plasmid_files=("$folder"plasmid_*.fasta)

    # Check if plasmid files exist
    if compgen -G "$folder"plasmid_*.fasta > /dev/null; then
        echo "🔹 Processing: $sample"
        for fasta in "${plasmid_files[@]}"; do
            [[ -e "$fasta" ]] || continue
            plasmid_name=$(basename "$fasta" .fasta)
            new_name="${folder}${plasmid_name}_${sample}.fasta"
            mv "$fasta" "$new_name"
            echo "   ✅ Renamed: $(basename "$fasta") → $(basename "$new_name")"
        done
    else
        echo "⚠️  No plasmid files found in: $sample"
        echo "$sample" >> "$LOG_FILE"
    fi
done

echo "---------------------------------------"
echo "🎯 Renaming process completed!"
echo "🗒️  Datasets without plasmids are listed in: $LOG_FILE"

