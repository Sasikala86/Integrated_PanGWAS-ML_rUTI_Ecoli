#!/bin/bash

# Input folder containing FastQC results
FASTQC_DIR="../results/fastqc_results"

# Output folder for MultiQC report
MULTIQC_OUT="../results/multiqc_report"

mkdir -p "$MULTIQC_OUT"

echo "Running MultiQC..."

multiqc "$FASTQC_DIR" -o "$MULTIQC_OUT"

echo "MultiQC completed!"
