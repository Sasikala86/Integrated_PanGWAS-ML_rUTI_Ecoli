#!/bin/bash

INPUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/plasmids/merged_per_isolate_renamed"
OUT_DIR="/home/sastra/Sasikala/Pangenome/mob_output/plasmids/merged_per_isolate_renamed/prokka_plasmids"

mkdir -p "$OUT_DIR"

for fasta in "$INPUT_DIR"/*_merged.fasta; do
    isolate=$(basename "$fasta" _merged.fasta)
    locustag=$(echo "$isolate" | cut -c1-12)

    prokka \
        --outdir "$OUT_DIR/$isolate" \
        --prefix "$isolate" \
        --locustag "$locustag" \
        --kingdom Bacteria \
        --cpus 24 \
        --addgenes \
        "$fasta"
done

