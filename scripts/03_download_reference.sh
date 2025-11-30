#!/bin/bash

# --------------------------------------------
# SLURM JOB SETTINGS
# This job downloads the reference genome,
# verifies file integrity using Ensembl CHECKSUMS,
# and builds HISAT2 index inside Apptainer container.
# --------------------------------------------
#SBATCH --job-name=hisat2_index
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M
#SBATCH --partition=pibu_el8
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# --------------------------------------------
# PROJECT AND RESOURCE DIRECTORIES
# --------------------------------------------
PROJECT_DIR="/data/users/gercan/rna_seq_toxoplasma_lung"
HISAT2_SAMTOOLS_SIF="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
RESOURCE_DIR="$PROJECT_DIR/resources"

# Output prefix for HISAT2 index
INDEX_BASE="Mmusculus_idx"

# --------------------------------------------
# ENSEMBL FTP PATHS FOR MOUSE (GRCm39, release 110)
# --------------------------------------------
FTP_PATH="ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/"
GTF_FTP_PATH="ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/"

# File names we expect to download
FASTA_FILE="Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
GTF_FILE="Mus_musculus.GRCm39.110.gtf.gz"
CHECKSUMS_FILE="CHECKSUMS"

# --------------------------------------------
# CREATE RESOURCE FOLDER AND MOVE INTO IT
# --------------------------------------------
mkdir -p "$RESOURCE_DIR"
cd "$RESOURCE_DIR"

# --------------------------------------------
# STEP 1: DOWNLOAD FASTA, GTF, AND CHECKSUM FILES
# --------------------------------------------
echo "Downloading reference files..."

# Download FASTA, GTF, and Checksums
wget "$FTP_PATH/$FASTA_FILE" # Genome FASTA
wget "$FTP_PATH/$CHECKSUMS_FILE" # Checksum file for integrity verification
wget "$GTF_FTP_PATH/$GTF_FILE" # Annotation GTF

# --------------------------------------------
# STEP 2: VERIFY FASTA FILE INTEGRITY
# Using Ensembl CHECKSUMS to ensure correct download
# --------------------------------------------
if [ -f "$FASTA_FILE" ] && [ -f "$CHECKSUMS_FILE" ]; then
    echo "Verifying $FASTA_FILE integrity..."
    
    # Extract expected checksum values for FASTA file
    EXPECTED_SUM=$(grep "$FASTA_FILE" "$CHECKSUMS_FILE" | awk '{print $1}')
    EXPECTED_BLOCKS=$(grep "$FASTA_FILE" "$CHECKSUMS_FILE" | awk '{print $2}')
    
    # Calculate checksum for the downloaded file
    CALCULATED_SUM_LINE=$(sum "$FASTA_FILE")
    CALCULATED_SUM=$(echo "$CALCULATED_SUM_LINE" | awk '{print $1}')
    CALCULATED_BLOCKS=$(echo "$CALCULATED_SUM_LINE" | awk '{print $2}')
    
    # Compare expected vs calculated checksum
    if [ "$CALCULATED_SUM" = "$EXPECTED_SUM" ] && [ "$CALCULATED_BLOCKS" = "$EXPECTED_BLOCKS" ]; then
        echo "Checksum passed for $FASTA_FILE. File is intact."
    else
        echo "Checksum FAILED for $FASTA_FILE. Expected sum: $EXPECTED_SUM, Calculated sum: $CALCULATED_SUM."
        exit 1
    fi
else
    echo "Download FAILED. FASTA or CHECKSUMS file not found."
    exit 1
fi

echo "Building Hisat2 index for $FASTA_FILE..."

# --------------------------------------------
# STEP 3: BUILD HISAT2 INDEX USING APPTAINER
# --------------------------------------------
apptainer exec --bind /data/ "$HISAT2_SAMTOOLS_SIF" \
    hisat2-build "$RESOURCE_DIR/$FASTA_FILE" "$RESOURCE_DIR/$INDEX_BASE"

echo "Hisat2 index built successfully."
