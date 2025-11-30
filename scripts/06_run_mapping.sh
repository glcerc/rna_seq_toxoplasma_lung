#!/bin/bash
# -----------------------------------------------------------
# SLURM SETTINGS
# HISAT2 paired-end RNA-seq alignment (array job)
# Each array index maps one sample from lung_samples.txt
# Produces sorted BAM + BAM index for downstream quantification.
# -----------------------------------------------------------
#SBATCH --job-name=hisat2_mapping
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=25000M
#SBATCH --partition=pibu_el8
#SBATCH --array=0-14
#SBATCH --output=/data/users/%u/rna_seq_toxoplasma_lung/logs/mapping_%A_%a.out
#SBATCH --error=/data/users/%u/rna_seq_toxoplasma_lung/logs/mapping_%A_%a.err

# -----------------------------------------------------------
# PATHS AND DIRECTORIES
# -----------------------------------------------------------
PROJECT="/data/users/$USER/rna_seq_toxoplasma_lung"

RAW="$PROJECT/data/raw_data"  # Input FASTQ folder
INDEX="$PROJECT/results/02_alignment_index/Mus_musculus_GRCm39" # HISAT2 index prefix
OUT="$PROJECT/results/03_alignment_bam"  # Output BAM folder
LIST="$PROJECT/scripts/lung_samples.txt"   # Sample names list
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p "$OUT"

# -----------------------------------------------------------
# LOAD SAMPLE NAME USING TASK ID
# Each line in lung_samples.txt corresponds to one sample
# -----------------------------------------------------------
mapfile -t SAMPLES < "$LIST"
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

R1="$RAW/${SAMPLE}_1.fastq.gz"
R2="$RAW/${SAMPLE}_2.fastq.gz"

echo "=== Mapping $SAMPLE ==="

# -----------------------------------------------------------
# STEP 1: HISAT2 MAPPING
# -x: index prefix
# -1/-2: paired-end reads
# --rna-strandness RF: library type (Illumina TruSeq)
# -p 4: use 4 threads
# -S: write SAM output
# HISAT2 stderr (alignment summary) redirected to summary.txt
# -----------------------------------------------------------
apptainer exec --bind /data "$CONTAINER" \
    hisat2 -x "$INDEX" \
    -1 "$R1" -2 "$R2" \
    --rna-strandness RF \
    -p 4 \
    -S "$OUT/${SAMPLE}.sam" \
    2> "$OUT/${SAMPLE}_hisat2_summary.txt"

# -----------------------------------------------------------
# STEP 2: CONVERT SAM → BAM
# Using -hbS to produce compressed BAM directly
# -----------------------------------------------------------
apptainer exec --bind /data "$CONTAINER" \
    samtools view -hbS "$OUT/${SAMPLE}.sam" \
    > "$OUT/${SAMPLE}.bam"

# -----------------------------------------------------------
# STEP 3: SORT BAM
# -@ threads
# -m memory per thread
# Sorting required for indexing + featureCounts
# -----------------------------------------------------------
apptainer exec --bind /data "$CONTAINER" \
    samtools sort -@ 4 -m 2G \
    -o "$OUT/${SAMPLE}.sorted.bam" \
    "$OUT/${SAMPLE}.bam"

# -----------------------------------------------------------
# STEP 4: INDEX SORTED BAM
# Creates .bai file used by IGV and featureCounts
# -----------------------------------------------------------
apptainer exec --bind /data "$CONTAINER" \
    samtools index "$OUT/${SAMPLE}.sorted.bam"

# -----------------------------------------------------------
# CLEANUP INTERMEDIATE FILES
# -----------------------------------------------------------
rm "$OUT/${SAMPLE}.sam" "$OUT/${SAMPLE}.bam"

echo "=== DONE: $SAMPLE ==="

