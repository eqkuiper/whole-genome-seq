#!/usr/bin/env bash
# stitch_domtblout.sh
#
# Recursively finds all *.domtblout.tsv files under a base directory,
# concatenates them, and adds a "subdir" column recording which parent
# directory each row came from.
#
# Usage:
#   ./stitch_domtblout.sh <base_dir> <output.tsv>

set -euo pipefail
BASE_DIR="$1"
OUTFILE="$2"

# Header: subdir + the 22 hmmsearch domtblout fields + description
printf 'subdir\ttarget_name\ttarget_accession\ttlen\tquery_name\tquery_accession\tqlen\tfull_seq_evalue\tfull_seq_score\tfull_seq_bias\tdom_num\tdom_of\tc_evalue\ti_evalue\tdom_score\tdom_bias\thmm_from\thmm_to\tali_from\tali_to\tenv_from\tenv_to\tacc\ttarget_description\n' > "$OUTFILE"

find "$BASE_DIR" -type f -iname "*.domtblout.tsv" | while read -r f; do
    subdir=$(basename "$(dirname "$f")")
    # grep exits 1 when a file has no non-comment lines (no domain hits).
    # With `set -o pipefail` that would otherwise kill the whole script,
    # so we swallow that specific "no matches" case with `|| true`.
    { grep -v '^#' "$f" || true; } | grep -v '^[[:space:]]*$' | awk -v OFS='\t' -v s="$subdir" '
    {
        desc = ""
        for (i = 23; i <= NF; i++) desc = (desc=="") ? $i : desc" "$i
        printf "%s", s
        for (i = 1; i <= 22; i++) printf "\t%s", $i
        printf "\t%s\n", desc
    }'
done >> "$OUTFILE"

echo "Done. $(($(wc -l < "$OUTFILE") - 1)) data rows written to $OUTFILE" >&2