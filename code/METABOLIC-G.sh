#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=60G
#SBATCH --array=0-39             # Adjust to number of genome folders minus 1
#SBATCH --job-name=metabolic_g-mags
#SBATCH --output=%A_%a-%x.out
#SBATCH --mail-user=esmee@u.northwestern.edu
#SBATCH --mail-type=END,FAIL

# Parent directory containing genome folders
PARENT_DIR="/projects/p32449/maca_mags_metabolic/data/2025-09-16_mags-in-folders"
OUT_DIR="/projects/b1042/OsburnLab/2025-09-16_metabolic_out"
mkdir -p "$OUT_DIR"

# Activate METABOLIC environment
module purge all
module load python-miniconda3
source activate /projects/p32449/goop_stirrers/miniconda3/envs/METABOLIC_v4.0

# METABOLIC-G executable
MB="/projects/p32449/goop_stirrers/METABOLIC_2025-09-02/METABOLIC/METABOLIC-G.pl"

# -------------------------------
# List all genome folders
# -------------------------------
FOLDERS=($(find "$PARENT_DIR" -mindepth 1 -maxdepth 1 -type d))
if [ ${#FOLDERS[@]} -eq 0 ]; then
    echo "No genome folders found in $PARENT_DIR"
    exit 1
fi

# Select folder for this array task
GENOME_FOLDER="${FOLDERS[$SLURM_ARRAY_TASK_ID]}"
BASE_NAME=$(basename "$GENOME_FOLDER")

# Output folder (METABOLIC will create subfolders inside)
GENOME_OUT="$OUT_DIR/$BASE_NAME"
mkdir -p "$GENOME_OUT"

# -------------------------------
# Run METABOLIC-G on the folder
# -------------------------------
perl "$MB" -in-gn "$GENOME_FOLDER" -o "$GENOME_OUT" -t $SLURM_CPUS_PER_TASK

echo "METABOLIC-G finished for $BASE_NAME"

