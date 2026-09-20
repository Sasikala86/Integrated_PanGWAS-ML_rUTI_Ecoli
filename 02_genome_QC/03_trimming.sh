#!/bin/bash

# Input and output directories
INPUT_DIR="../../01_data/fastq_files"
OUTPUT_DIR="../results/trimmed_fastq"

# Trimmomatic adapter file
ADAPTERS="/path/to/TruSeq3-PE.fa"

# Number of threads
THREADS=4

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Activate conda environment
source ~/anaconda3/bin/activate trimm_env

# Loop over all forward reads (_1.fastq or _1.fastq.gz)
for file in "$INPUT_DIR"/*_1.fastq*; do

    # Extract sample name
    sample=$(basename "$file" "_1.fastq")
    sample=${sample%.gz}

    # Input file paths
    input_forward="$INPUT_DIR/${sample}_1.fastq"
    input_reverse="$INPUT_DIR/${sample}_2.fastq"

    # Check if files are gzipped
    if [[ -f "${input_forward}.gz" ]]; then
        input_forward="${input_forward}.gz"
        input_reverse="${input_reverse}.gz"
    fi

    # Output file paths
    output_forward_paired="$OUTPUT_DIR/${sample}_1_paired.fastq.gz"
    output_forward_unpaired="$OUTPUT_DIR/${sample}_1_unpaired.fastq.gz"
    output_reverse_paired="$OUTPUT_DIR/${sample}_2_paired.fastq.gz"
    output_reverse_unpaired="$OUTPUT_DIR/${sample}_2_unpaired.fastq.gz"

    # Run Trimmomatic
    trimmomatic PE -threads "$THREADS" -phred33 \
        "$input_forward" "$input_reverse" \
        "$output_forward_paired" "$output_forward_unpaired" \
        "$output_reverse_paired" "$output_reverse_unpaired" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

    echo "Trimming completed for sample: $sample"

done

echo "All Trimmomatic processing completed."
