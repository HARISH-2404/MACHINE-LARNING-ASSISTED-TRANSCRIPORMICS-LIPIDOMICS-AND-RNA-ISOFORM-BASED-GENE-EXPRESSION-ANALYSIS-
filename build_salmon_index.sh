#!/usr/bin/env bash
# Build Salmon index from a transcript FASTA matching the annotation used for tx2gene.
set -euo pipefail
if [[ $# -lt 2 ]]; then echo "Usage: bash scripts/build_salmon_index.sh transcripts.fa index_dir [kmer_length]"; exit 2; fi
FASTA="$1"; INDEX="$2"; K="${3:-31}"
command -v salmon >/dev/null || { echo 'Install Salmon first'; exit 1; }
[[ -s "$FASTA" ]] || { echo "Missing transcript FASTA: $FASTA"; exit 1; }
salmon index -t "$FASTA" -i "$INDEX" -k "$K"
echo "Index created at $INDEX. Keep FASTA and annotation versions recorded."