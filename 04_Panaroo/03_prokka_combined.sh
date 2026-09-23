#!/bin/bash

input_dir="/home/sastra/Sasikala/Pangenome/FastANI/datasets"
output_dir="/home/sastra/Sasikala/Pangenome/Prokka"

mkdir -p "$output_dir"

total_cores=28
threads_per_genome=6
parallel_jobs=$(( total_cores / threads_per_genome ))

echo "⚡ Total cores: $total_cores"
echo "⚡ Threads per genome: $threads_per_genome"
echo "⚡ Max parallel genomes: $parallel_jobs"

if command -v parallel >/dev/null 2>&1; then
    echo "✅ GNU Parallel detected — running genomes in parallel"
    find "$input_dir" -maxdepth 1 -type f \( -iname "*.fna" -o -iname "*.fa" -o -iname "*.fasta" \) | \
    parallel -j "$parallel_jobs" --eta '
        fasta_file={}
        base=$(basename "$fasta_file")
        base=${base%.*}

        mkdir -p "'"$output_dir"'/$base"

        echo "=== Running Prokka on $base ==="
        prokka \
            --outdir "'"$output_dir"'/$base" \
            --prefix "$base" \
            --cpus '"$threads_per_genome"' \
            --force \
            "$fasta_file"

        echo "✅ Finished: $base"
        echo "-----------------------------------"
    '
else
    echo "⚠️ GNU Parallel not found — running genomes sequentially"
    for fasta in "$input_dir"/*.{fna,fa,fasta}; do
        [ -e "$fasta" ] || continue
        base=$(basename "$fasta")
        base=${base%.*}

        mkdir -p "$output_dir/$base"

        echo "=== Running Prokka on $base ==="
        prokka \
            --outdir "$output_dir/$base" \
            --prefix "$base" \
            --cpus "$threads_per_genome" \
            --force \
            "$fasta"

        echo "✅ Finished: $base"
        echo "-----------------------------------"
    done
fi

echo "🎯 All Prokka runs completed!"
echo "Results saved in: $output_dir"
