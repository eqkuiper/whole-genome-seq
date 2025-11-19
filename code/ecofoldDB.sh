#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --array=0-5 # adjust to number of genomes
#SBATCH --mem=10G
#SBATCH --job-name=EcoFoldDB
#SBATCH --mail-user=esmee@u.northwestern.edu # change to your email
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=%A_%a-%x.out

# USER INPUTS
# tsv with genome fp and genome id as cols
faa_dir=/projects/p32449/isolate_genomes/data/prokka
out_dir=/projects/p32449/isolate_genomes/data/EcoFoldDB
#############

module load cuda   # If your cluster requires a CUDA module
echo "Using GPU: $CUDA_VISIBLE_DEVICES"
nvidia-smi

mkdir -p "${out_dir}/tmp"

# Generate genome list (all .faa files in prokka sample dirs)
find ${faa_dir} -name "*.faa" > "${out_dir}/tmp/genomes.txt"

# Pick genome for this SLURM_ARRAY_TASK_ID
genome=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${out_dir}/tmp/genomes.txt")

# Add Foldseek/EcoFoldDB to PATH
export PATH=/projects/p31618/software/EcoFoldDB/foldseek/bin/:$PATH

# Move to EcoFoldDB directory
cd /projects/p31618/software/EcoFoldDB

# Get the filename without path
genome_name=$(basename "$genome" .faa)

# Define output directory
genome_dir="${out_dir}/annotated/${genome_name}"

# Run annotation
./EcoFoldDB-annotate.sh \
    --EcoFoldDB_dir /projects/p31618/software/EcoFoldDB/EcoFoldDB_v2.0 \
    --gpu 1 \
    --ProstT5_dir /projects/p31618/software/EcoFoldDB/ProstT5_dir \
    -o "${genome_dir}" \
    "$genome"

rm ${out_dir}/tmp



