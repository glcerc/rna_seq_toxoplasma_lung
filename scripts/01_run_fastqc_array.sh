#!/bin/bash
# -----------------------------------------------------------------------
# SLURM Settings (Using JOB ARRAY for 15 samples)
# -----------------------------------------------------------------------
#SBATCH --job-name=fastqc_array_lung    # This indicates that it is a Bash script. It reports the job name, duration, and amount of CPU and RAM to the cluster manager (SLURM).
#SBATCH --time=01:00:00                 # 1 hour max time per sample
#SBATCH --cpus-per-task=1               # 1 CPU
#SBATCH --mem=1000M                     # 1000MB RAM
#SBATCH --partition=pibu_el8
#SBATCH --output=logs/fastqc_array_%A_%a.out  # Unique output file per job
#SBATCH --error=logs/fastqc_array_%A_%a.err
#SBATCH --array=1-15                    # Starts 15 independent jobs Job Array Command: Tells Slurm to create a job array that will run this script from 1 to 15 (15 times).

# -----------------------------------------------------------------------
# Directory and Variable Definitions 
# -----------------------------------------------------------------------
FASTQC_OUT_DIR="results/01_fastqc_output"

#FastQC Container Path--> Defines the full path to the Apptainer/Singularity container where the FastQC tool is located.
CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"


#Path to text file containing list of 15 file paths to process. This is done manually (or with a separate script) beforerunning the script. It lists the full paths of all R1 files and saves them as text in the file lung_samples.txt.
#The file contains text (TXT), not data.
SAMPLE_LIST="scripts/lung_samples.txt" 


# Create output directory
mkdir -p "$FASTQC_OUT_DIR"

# -----------------------------------------------------------------------
# CORE LOGIC: Get the specific FASTQ file for this array job ID
# -----------------------------------------------------------------------
# Reads the corresponding file path based on the job array index (SLURM_ARRAY_TASK_ID)
R1_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < "$SAMPLE_LIST")

# Define the paired R2 file path (Using the _1.fastq.gz -> _2.fastq.gz naming)
R2_FILE="${R1_FILE/_1.fastq.gz/_2.fastq.gz}"

# Log the current array task ID and the specific R1/R2 files being processed.
echo "--- FastQC Analysis for Array Task ID: $SLURM_ARRAY_TASK_ID ---"
echo "Processing R1: $R1_FILE"
echo "Processing R2: $R2_FILE"

# Execute FastQC inside the Apptainer container
#GIVING INPUT: The FastQC command is run. The paths of both the dynamically retrieved files R1and R2 are given as input.FastQC reads these two files and generates quality reports.
apptainer exec --bind /data/ "$CONTAINER" fastqc -o "$FASTQC_OUT_DIR" "$R1_FILE" "$R2_FILE" 

#If an error occurs while the FastQC program is running, it catches it and writes to the error log which file is causingthe problem (if block).
if [ $? -ne 0 ]; then
    echo "Error: FastQC failed for \$(basename \$R1_FILE)" >&2
fi
echo "--- FastQC Analysis Complete ---"
