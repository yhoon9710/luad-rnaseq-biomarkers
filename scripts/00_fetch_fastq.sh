#!/usr/bin/env bash
# Download paired-end FASTQ for each run in config/samples.tsv and subsample
# to N read pairs so the pipeline runs on a laptop. Same seed for R1/R2 keeps pairs in sync.
set -euo pipefail

SAMPLES=${1:-config/samples.tsv}
N=${N_READS:-2000000}
OUT=data/raw
TMP=data/tmp
mkdir -p "$OUT" "$TMP"

grep -v '^#' "$SAMPLES" | tail -n +2 | while IFS=$'\t' read -r sample condition patient srr; do
  [[ -z "${sample:-}" ]] && continue
  if [[ -s "$OUT/${sample}_R1.fastq.gz" ]]; then
    echo "[skip] $sample"; continue
  fi
  echo "[fetch] $sample <- $srr"
  prefetch "$srr" -O "$TMP"
  fasterq-dump --split-files --threads 4 -O "$TMP" "$TMP/$srr"
  seqtk sample -s42 "$TMP/${srr}_1.fastq" "$N" | gzip > "$OUT/${sample}_R1.fastq.gz"
  seqtk sample -s42 "$TMP/${srr}_2.fastq" "$N" | gzip > "$OUT/${sample}_R2.fastq.gz"
  rm -rf "$TMP/${srr}"*
done

# Transcriptome reference for Salmon (once)
mkdir -p resources
if [[ ! -s resources/gencode.transcripts.fa.gz ]]; then
  wget -O resources/gencode.transcripts.fa.gz \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
fi
echo "done"
