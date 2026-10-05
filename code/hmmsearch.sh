#!/bin/bash
#SBATCH --account=p32449
#SBATCH --partition=short
#SBATCH --job-name=hmmsearch_mags
#SBATCH --output=logs/hmmsearch_%A_%a.out
#SBATCH --error=logs/hmmsearch_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-7        # set upper bound to number of MAGs; %20 = max 20 at once

set -euo pipefail
mkdir -p logs

module purge
module load mamba/23.1.0
source activate /projects/p31618/software/hmmer-3.4   # needs hmmsearch AND prodigal; see note below

# ---- EDIT THESE --------------------------------------------------------
MAG_DIR="/projects/p32449/isolate_genomes/data/spades"     # directory of fasta files
MAG_EXT="fasta"                                      # fa, fasta, or fna
MAG_LIST="/projects/p32449/isolate_genomes/data/isolate_list.txt"                           # created once, see below

HMM_FILE="/projects/p32449/maca_take2/hmm_profiles/degredation_ko.hmm"
OUT_ROOT="data/hmmsearch/degradation"
# -------------------------------------------------------------------------

CPUS=${SLURM_CPUS_PER_TASK}

# Pick this task's MAG from the list
MAG_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MAG_LIST")
[ -n "$MAG_FILE" ] || { echo "No MAG for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
#MAG=$(basename "$MAG_FILE" ".${MAG_EXT}")
MAG=$(basename "$(dirname "$MAG_FILE")")

OUT_DIR="${OUT_ROOT}/${MAG}"
mkdir -p "$OUT_DIR"
FAA="${OUT_DIR}/${MAG}.faa"

echo "Task ${SLURM_ARRAY_TASK_ID}: ${MAG} on $(hostname) at $(date)"

# Predict proteins (skip if already done, so reruns are cheap)
if [ ! -s "$FAA" ]; then
  prodigal -i "$MAG_FILE" -a "$FAA" -p meta -q -o /dev/null
fi

hmmsearch \
  --cpu "$CPUS" \
  -E 1e-10 \
  --domtblout "${OUT_DIR}/${MAG}.domtblout.tsv" \
  --tblout    "${OUT_DIR}/${MAG}.tblout.tsv" \
  -o          "${OUT_DIR}/${MAG}.hmmsearch.out" \
  "$HMM_FILE" "$FAA"

echo "Done ${MAG} at $(date)"