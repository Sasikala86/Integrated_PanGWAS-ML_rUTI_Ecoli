#!/bin/bash

# Input and output directories
UNICYCLER_DIR="../results/unicycler_output"
ASSEMBLY_DIR="../results/all_assemblies"
CHECKM_OUT="../results/checkm_output"

# Create directories if not exist
mkdir -p "$ASSEMBLY_DIR"
mkdir -p "$CHECKM_OUT"

# Copy and rename assemblies
for SAMPLE_DIR in "$UNICYCLER_DIR"/*; do

    SAMPLE=$(basename "$SAMPLE_DIR")

    if [ -f "$SAMPLE_DIR/assembly.fasta" ]; then

        cp "$SAMPLE_DIR/assembly.fasta" \
           "$ASSEMBLY_DIR/${SAMPLE}.fasta"

        echo "Copied $SAMPLE_DIR/assembly.fasta -> $ASSEMBLY_DIR/${SAMPLE}.fasta"

    else

        echo "No assembly.fasta found in $SAMPLE_DIR"

    fi

done

# Run CheckM lineage workflow
checkm lineage_wf \
    -x fasta \
    "$ASSEMBLY_DIR" \
    "$CHECKM_OUT" \
    -t 16
