#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH --array=0-39             # Adjust to number of genome folders minus 1
#SBATCH --job-name=gtdbtk-classify-wf_MaCa21
#SBATCH --output=%A_%a-%x.out
#SBATCH --mail-user=esmee@u.northwestern.edu
#SBATCH --mail-type=END,FAIL

# set directories
in_dir="/projects/p32449/maca_mags_metabolic/data/2025-09-16_mags-in-folders"
out_parent_dir="/scratch/jhr1326/2025-09-22_gtdbk_out"

mkdir -p "$out_parent_dir" 

# Activate METABOLIC environment
module purge all
module load python-miniconda3
source activate /projects/p32449/goop_stirrers/miniconda3/envs/gtdbtk-2.5.2 # may need to change this 

# List all genome folders
folders=($(find "$in_dir" -mindepth 1 -maxdepth 1 -type d))
if [ ${#folders[@]} -eq 0 ]; then
    echo "No genome folders found in $in_dir"
    exit 1
fi

# Select folder for this array task
genome_folder="${folders[$SLURM_ARRAY_TASK_ID]}"
base_name=$(basename "$genome_folder")

# Set outdir for each array task
out_dir="${out_parent_dir}/${base_name}"

mkdir -p "$out_dir"

# Run gtdbtk classify
gtdbtk classify_wf \
    --genome_dir "$genome_folder" \
    --out_dir "$out_dir" \
    --extension "fasta" \
    --cpus $SLURM_CPUS_PER_TASK 

