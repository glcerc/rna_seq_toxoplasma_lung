#!/bin/bash
# -----------------------------------------------------------
# SLURM SETTINGS
# Gene-level read counting using FeatureCounts (Subread).
# This script processes ALL sorted BAM files at once and
# generates a single count matrix (gene_counts.txt).
# -----------------------------------------------------------
#SBATCH --job-name=featureCounts_lung
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=2000M
#SBATCH --partition=pibu_el8
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# ===============================
# Paths
# ===============================
PROJECT_DIR="/data/users/gercan/rna_seq_toxoplasma_lung"

BAM_DIR="$PROJECT_DIR/results/03_alignment_bam"  # Sorted BAM files
GTF_FILE="$PROJECT_DIR/resources/Mus_musculus.GRCm39.110.gtf.gz" # Annotation
OUT_DIR="$PROJECT_DIR/results/04_featureCounts"  # Output folder

# Subread container (contains featureCounts)
CONTAINER="/containers/apptainer/subread_2.0.6.sif"

mkdir -p "$OUT_DIR"

echo "=== Running FeatureCounts ==="
echo "Using GTF: $GTF_FILE"
echo "Input BAM directory: $BAM_DIR"

# -----------------------------------------------------------
# FEATURECOUNTS COMMAND
# -T: number of threads (use SLURM_CPUS_PER_TASK)
# -p: paired-end mode (count fragments instead of reads)
# -s 2: strandedness (2 = reverse stranded, for TruSeq RNA)
# -a: GTF annotation
# -o: output count matrix
# -----------------------------------------------------------
apptainer exec --bind /data/ "$CONTAINER" featureCounts \
    -T $SLURM_CPUS_PER_TASK \
    -p \
    -s 2 \
    -a "$GTF_FILE" \
    -o "$OUT_DIR/gene_counts.txt" \
    "$BAM_DIR"/*.sorted.bam

echo "=== FeatureCounts Completed ==="

