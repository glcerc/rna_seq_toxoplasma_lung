#!/bin/bash
# -----------------------------------------------------------------------
# SLURM Settings
# -----------------------------------------------------------------------
#SBATCH --job-name=multiqc
#SBATCH --time=00:30:00                 # 30 minutes
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M                     # 2GB RAM
#SBATCH --partition=pibu_el8
#SBATCH --output=logs/multiqc_before_%j.out
#SBATCH --error=logs/multiqc_before_%j.err

# ---
# Purpose: To collect all reports from different analyses (FastQC) and create a single summary report.
#
# Why is this NOT a job array?
# MultiQC is designed to run as a single job. It scans an entire directory (FASTQC_DIR)
# for log files, aggregates them, and outputs one single HTML report.
# We run FastQC (01) as an array (15 parallel jobs) but MultiQC (02) as one single job.
# ---

# -----------------------------------------------------------------------
# Directory and Variable Definitions
# -----------------------------------------------------------------------
# This is the correct container path you found:
CONTAINER="/containers/apptainer/multiqc-1.19.sif"
FASTQC_DIR="results/01_fastqc_output"
MULTIQC_DIR="results/01_multiqc"

# Create output directory
mkdir -p "$MULTIQC_DIR"

echo "--- Running MultiQC on FastQC reports ---"
echo "Input directory: $FASTQC_DIR"
echo "Output directory: $MULTIQC_DIR"

# Execute MultiQC
# -f: Force overwrite of existing reports
# -o: Output directory
# -n: Name of the output HTML file
# $FASTQC_DIR: The directory to scan for FastQC results
apptainer exec --bind /data/ "$CONTAINER" multiqc \
    -f \
    -o "$MULTIQC_DIR" \
    -n "multiqc_report_before_trim.html" \
    "$FASTQC_DIR"

echo "--- MultiQC Complete ---"
