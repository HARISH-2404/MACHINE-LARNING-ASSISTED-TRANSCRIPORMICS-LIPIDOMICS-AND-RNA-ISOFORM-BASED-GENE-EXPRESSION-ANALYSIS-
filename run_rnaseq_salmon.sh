#!/usr/bin/env bash
# RNA-seq read QC, optional trimming, and Salmon transcript quantification.
# Requires: fastqc, multiqc, fastp (optional), salmon; paired-end FASTQ files.
# Input samples.tsv columns: sample_id, read1, read2 (paths relative to project root or absolute).
set -euo pipefail
if [[ $# -lt 3 ]]; then
  echo "Usage: bash scripts/run_rnaseq_salmon.sh samples.tsv transcriptome_index output_dir [--trim]" >&2; exit 2
fi
SAMPLES="$1"; INDEX="$2"; OUT="$3"; TRIM="${4:-}"
mkdir -p "$OUT/qc/raw" "$OUT/qc/trimmed" "$OUT/quant"
command -v fastqc >/dev/null || { echo 'fastqc not found'; exit 1; }
command -v multiqc >/dev/null || { echo 'multiqc not found'; exit 1; }
command -v salmon >/dev/null || { echo 'salmon not found'; exit 1; }
if [[ "$TRIM" == "--trim" ]]; then command -v fastp >/dev/null || { echo 'fastp not found'; exit 1; }; fi
# Skip header; TSV must have sample_id, read1, read2 columns.
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample r1 r2; do
  [[ -n "$sample" && -f "$r1" && -f "$r2" ]] || { echo "Invalid row or missing FASTQ for $sample" >&2; exit 1; }
  fastqc "$r1" "$r2" --outdir "$OUT/qc/raw"
  q1="$r1"; q2="$r2"
  if [[ "$TRIM" == "--trim" ]]; then
    mkdir -p "$OUT/trimmed"
    fastp -i "$r1" -I "$r2" -o "$OUT/trimmed/${sample}_R1.fastq.gz" -O "$OUT/trimmed/${sample}_R2.fastq.gz" \
      --html "$OUT/qc/trimmed/${sample}.html" --json "$OUT/qc/trimmed/${sample}.json"
    q1="$OUT/trimmed/${sample}_R1.fastq.gz"; q2="$OUT/trimmed/${sample}_R2.fastq.gz"
  fi
  salmon quant -i "$INDEX" -l A -1 "$q1" -2 "$q2" --validateMappings \
    -p "${THREADS:-4}" -o "$OUT/quant/$sample"
done
multiqc "$OUT/qc" "$OUT" --outdir "$OUT/qc/multiqc" || true
echo "Finished. Review QC reports and confirm quantification before downstream analysis."