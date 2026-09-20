#!/bin/bash

# =========================
# MOB-suite batch processing
# =========================

# Input folder with genomes
INPUT_DIR="/home/sastra/Sasikala/Pangenome/FastANI/datasets"

# Output folder
OUTPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output"
mkdir -p "$OUTPUT_DIR"

# Prefixes for plasmids and chromosomes
PREFIX_PLASMID="p_"
PREFIX_CHROM="bc_"

# Threads per genome
THREADS=8

# Number of genomes to run in parallel
JOBS=3

# Load environment (if needed)
# conda activate mobsuite_env

cd "$INPUT_DIR" || exit

# Function to run MOB-suite on a single genome
run_mob() {
    genome="$1"
    name=$(basename "$genome" | cut -d. -f1)
    genome_out="$OUTPUT_DIR/$name"
    mkdir -p "$genome_out"

    echo "🟢 Running MOB-suite on $genome with $THREADS threads..."

    mob_recon -i "$genome" -o "$genome_out" -n $THREADS --force

    # Rename plasmid FASTA if exists
    if [ -f "$genome_out/plasmid.fasta" ]; then
        mv "$genome_out/plasmid.fasta" "$genome_out/${PREFIX_PLASMID}${name}.fasta"
    fi

    # Rename chromosome FASTA if exists
    if [ -f "$genome_out/chromosome.fasta" ]; then
        mv "$genome_out/chromosome.fasta" "$genome_out/${PREFIX_CHROM}${name}.fasta"
    fi

    echo "✅ Finished $genome"
}

export -f run_mob
export OUTPUT_DIR PREFIX_PLASMID PREFIX_CHROM THREADS

# Run genomes in parallel safely
find . -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fna" \) | \
    parallel -j $JOBS run_mob {}
    
echo "🎯 All genomes processed!"

