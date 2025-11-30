#!/bin/bash
# -----------------------------------------------------------
# SLURM SETTINGS
# This job builds a HISAT2 index from the decompressed
# reference genome FASTA file. The index will be stored in
# results/02_alignment_index.
# -----------------------------------------------------------
#SBATCH --job-name=hisat2_index
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M
#SBATCH --partition=pibu_el8
#SBATCH --output=logs/index_%j.out
#SBATCH --error=logs/index_%j.err

# -----------------------------------------------------------
# PATHS AND DIRECTORIES
# -----------------------------------------------------------
PROJECT="/data/users/$USER/rna_seq_toxoplasma_lung"
RES="$PROJECT/resources"
INDEX_OUT="$PROJECT/results/02_alignment_index"

# Container that includes HISAT2 and samtools
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Reference FASTA (must be uncompressed)
FA="$RES/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# -----------------------------------------------------------
# SAFETY CHECK: Ensure FASTA file exists
# HISAT2 cannot index .gz files; decompressed FASTA required.
# -----------------------------------------------------------
if [ ! -f "$FA" ]; then
    echo "ERROR: FASTA file does not exist: $FA"
    echo "Please run the download/decompress script first!"
    exit 1
fi

# -----------------------------------------------------------
# PREPARE OUTPUT DIRECTORY
# -----------------------------------------------------------
mkdir -p "$INDEX_OUT"

echo "Cleaning old index..."
rm -f "$INDEX_OUT"/*.ht2 # Remove any previous index to avoid mixing files

# -----------------------------------------------------------
# BUILD HISAT2 INDEX
# The prefix determines index file names: prefix.{1..8}.ht2
# -----------------------------------------------------------
echo "Building HISAT2 index..."
apptainer exec --bind /data "$CONTAINER" \
    hisat2-build "$FA" "$INDEX_OUT/Mus_musculus_GRCm39"

# -----------------------------------------------------------
# FINAL MESSAGE + LIST RESULTS
# -----------------------------------------------------------
echo "DONE. Index is ready in:"
ls -lh "$INDEX_OUT"

