#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 148:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --array=0-5  # change to number of metagenomes
#SBATCH --mem=150GB
#SBATCH --mail-user=esmee@u.northwestern.edu # change to your email
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name="spades-assembly"
#SBATCH --output=/scratch/jhr1326/%A_%a.out       

### USER INPUTS
metagenome_list="/projects/p32449/maca_mags_metabolic/data/mags_to_annotate_assemblies.txt" # list of metagenomes
work_dir="/projects/p32449/maca_mags_metabolic/data/2025-10-07_maca_metaG" # directory where inputs and outputs live
data_in="${work_dir}/01_trimmomatic_out" # fp with trimmomatic output
out_dir="/scratch/jhr1326" # directory with enough space for big old outputs
###############

# define metagenome
IFS=$'\n' read -d '' -r -a input_args < "${metagenome_list}"
metagenome=${input_args[$SLURM_ARRAY_TASK_ID]}

# make output directory for each genome
data_out="${out_dir}/02_assembled-spades_10Nov2025/${metagenome}"
mkdir -p $data_out

# submit an assembly job for each metagenome from trimmomatic output
echo "Assembling ${metagenome} using SPAdes ..."

module purge all
module load python-miniconda3
source activate /projects/p31618/software/spades/4.2.0

# double check to make the data required is present
if [[ ! -f ${data_in}/${metagenome}_R1_paired.fastq ]] || [[ ! -f ${data_in}/${metagenome}_R2_paired.fastq ]]; then
    echo "Error: input files for ${metagenome} not found!"
    exit 1
fi

# log start time
start_time=$(date +%s)
echo "Assembly of ${metagenome} started at: $(date)"

# run SPAdes
spades.py \
-t $SLURM_NTASKS_PER_NODE \
-k 15,21,33,55,77 \
--only-assembler \
--meta \
-1 ${data_in}/${metagenome}_R1_paired.fastq \
-2 ${data_in}/${metagenome}_R2_paired.fastq \
-o $data_out

# record end time
end_time=$(date +%s)
echo "Assembly of ${metagenome} finished at: $(date)"

# calculate total runtime in seconds
total_time=$((end_time - start_time))

# convert seconds to hours:minutes:seconds
hours=$((total_time / 3600))
minutes=$(((total_time % 3600) / 60))
seconds=$((total_time % 60))

echo "Total runtime for ${metagenome}: ${hours}h ${minutes}m ${seconds}s"

echo "Assembly complete! Have a great $(date +%A)!"