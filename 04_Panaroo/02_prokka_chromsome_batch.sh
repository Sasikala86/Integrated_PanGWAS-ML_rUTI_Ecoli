#!/bin/bash

INPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/chromosomes"  # or chromosomes
OUTPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/chromosomes/prokka_output_chromosomes"
THREADS=4

mkdir -p "$OUTPUT_DIR"

# Loop over all fasta files
for fasta in "$INPUT_DIR"/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    echo "🔹 Running Prokka for $sample ..."

    prokka --outdir "$OUTPUT_DIR/$sample" \
           --prefix "$sample" \
           --cpus "$THREADS" \
           "$fasta"

    echo "✅ Completed $sample"
done
