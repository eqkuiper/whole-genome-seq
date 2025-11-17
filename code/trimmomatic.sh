#!/bin/bash
#SBATCH -A p32449               				  				
#SBATCH -p short            				  					 
#SBATCH -t 04:00:00            				      					
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-20 # change to number of metagenomes
#SBATCH --mem=8GB
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=esmee@u.northwestern.edu # change to your email
#SBATCH --mail-type=END 
#SBATCH --job-name="trimmomatic"
#SBATCH --output=%A_%a.out      

module purge all
module load python-miniconda3
source activate /projects/p31618/software_2025/trimmomatic

# list of metagenomes 
# naming scheme should match input fastq.gz files. for example, given the following files:
    # Metagenome_1_S1_R1_001.fastq.gz Metagenome_1_S1_R2_002.fastq.gz
# the corresponding metagenome name will be:
    # Metagenome_1_S1
IFS=$'\n' read -d '' -r -a input_args < /projects/p32449/maca_mags_metabolic/data/mags_to_annotate_assemblies.txt
metagenome=${input_args[$SLURM_ARRAY_TASK_ID]}

# input data folder contains paired-end fastq.gz files
DATA_IN=/projects/p31618/nu-seq/Osburn02_12.10.2021_metagenomes/reads
# output directory of choice
DATA_OUT=/projects/p32449/maca_mags_metabolic/data/2025-10-07_maca_metaG
mkdir -p $DATA_OUT

printf "\n | trimmomatic version: `trimmomatic -version`\n\n"

cd $DATA_IN

printf "
| Input metagenomic reads:
`ls -f $DATA_IN | grep fastq.gz | sort -u`
"

printf "
 ----------------------------------\n   ->trimming adapters off paired-end reads from $f\n"

trimmomatic PE \
-threads $SLURM_NTASKS \
${metagenome}_R1_001.fastq.gz \
${metagenome}_R2_001.fastq.gz  \
$DATA_OUT/${metagenome}_R1_paired.fastq \
$DATA_OUT/${metagenome}_R1_unpaired.fastq \
$DATA_OUT/${metagenome}_R2_paired.fastq \
$DATA_OUT/${metagenome}_R2_unpaired.fastq \
ILLUMINACLIP:/projects/p31618/databases/adapters.fa:2:40:15 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:35


printf "\n
fastq input directory: $DATA_IN
data output directory: $DATA_OUT
--------------------------------
job $SLURM_JOB_ID completed.
 | node          | `hostname`
 | date          | `date`
 | mem per core  | 2GB
 | cores         | $SLURM_NTASKS
"
