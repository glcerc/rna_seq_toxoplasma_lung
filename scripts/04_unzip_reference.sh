#!/bin/bash
# --------------------------------------------
# SLURM SETTINGS
# This job checks the resources folder and
# unzips the mouse reference FASTA file (.gz)
# using gunzip inside the Hisat2/Samtools container.
# --------------------------------------------
#SBATCH --job-name=gunzip_fasta
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/%u/rna_seq_toxoplasma_lung/logs/unzip_%j.out
#SBATCH --error=/data/users/%u/rna_seq_toxoplasma_lung/logs/unzip_%j.err

# --------------------------------------------
# PROJECT AND RESOURCE PATHS
# Using $USER makes the script reusable for any user
# --------------------------------------------
PROJECT_DIR="/data/users/$USER/rna_seq_toxoplasma_lung"
RES="$PROJECT_DIR/resources"

# Apptainer/Singularity container that includes gunzip (via samtools image)
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# --------------------------------------------
# MOVE INTO RESOURCES FOLDER
# FASTA .gz file must be inside this directory
# --------------------------------------------
cd "$RES"

# --------------------------------------------
# STEP 1: UNZIP FASTA FILE
# -k keeps the .gz file after extraction
# This avoids accidental deletion of compressed FASTA
# --------------------------------------------
echo "--- Unzipping FASTA file if needed ---"
apptainer exec --bind /data "$CONTAINER" \
    gunzip -k Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# --------------------------------------------
# FINAL MESSAGE
# --------------------------------------------
echo "--- Unzip complete ---"

