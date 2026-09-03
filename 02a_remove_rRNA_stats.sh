#!/bin/bash

set -euo pipefail

read -rp "Directory with raw files (input FASTQ): " RAW_DIR
read -rp "Directory with clean files (after fastp): " CLEAN_DIR
read -rp "Directory with files after rRNA removal: " RRNA_DIR
read -rp "Output directory: " OUTDIR
read -rp "Output file name [default: rrna_removal_stats.tsv]: " OUTFILE

OUTFILE=${OUTFILE:-rrna_removal_stats.tsv}

mkdir -p "$OUTDIR"

OUTFILE_PATH="${OUTDIR}/${OUTFILE}"
LOGFILE="${OUTDIR}/missing_files.log"

echo -e "Sample\tReads_input\tReads_clean\tReads_rRNAfree\tPercent_remaining" > "$OUTFILE_PATH"
echo -e "Missing files log (generated on $(date))" > "$LOGFILE"

for dir in "$RAW_DIR" "$CLEAN_DIR" "$RRNA_DIR"; do
  [[ -d "$dir" ]] || { echo "Directory does not exist: $dir"; exit 1; }
done

mapfile -t R1_CLEAN_LIST < <(find "$CLEAN_DIR" -type f -iname "*_R1_clean.fastq.gz" | sort)

if [[ "${#R1_CLEAN_LIST[@]}" -eq 0 ]]; then
  echo "No *_R1_clean.fastq.gz files found in directory: $CLEAN_DIR"
  exit 1
fi

for R1 in "${R1_CLEAN_LIST[@]}"; do
  SAMPLE=$(basename "$R1" | sed -E 's/_R1_clean\.f.*//')
  R2="$CLEAN_DIR/${SAMPLE}_R2_clean.fastq.gz"
  OUT1="$RRNA_DIR/${SAMPLE}.rRNAfree_R1.fastq.gz"
  OUT2="$RRNA_DIR/${SAMPLE}.rRNAfree_R2.fastq.gz"

  RAW1=$(find "$RAW_DIR" -type f -iname "${SAMPLE}_R1*.f*q*" | head -n 1 || true)
  RAW2=$(find "$RAW_DIR" -type f -iname "${SAMPLE}_R2*.f*q*" | head -n 1 || true)

  MISSING=0
  for f in "$RAW1" "$RAW2" "$R1" "$R2" "$OUT1" "$OUT2"; do
    [[ -n "$f" && -f "$f" ]] || { echo "No file: $f" >> "$LOGFILE"; MISSING=1; }
  done

  if [[ "$MISSING" -eq 1 ]]; then
    echo "⚠️  $SAMPLE – missing files, skipping"
    continue
  fi

  READS_INPUT=$( (zcat "$RAW1"; zcat "$RAW2") | wc -l )
  READS_INPUT=$((READS_INPUT / 4))

  READS_CLEAN=$( (zcat "$R1"; zcat "$R2") | wc -l )
  READS_CLEAN=$((READS_CLEAN / 4))

  READS_RRNA=$( (zcat "$OUT1"; zcat "$OUT2") | wc -l )
  READS_RRNA=$((READS_RRNA / 4))

  if [[ "$READS_CLEAN" -eq 0 ]]; then
    echo "⚠️  $SAMPLE – READS_CLEAN = 0, skipping" | tee -a "$LOGFILE"
    continue
  fi

  PERCENT=$(awk -v a="$READS_CLEAN" -v b="$READS_RRNA" 'BEGIN {printf "%.2f", (b/a)*100}')

  echo -e "$SAMPLE\t$READS_INPUT\t$READS_CLEAN\t$READS_RRNA\t$PERCENT" >> "$OUTFILE_PATH"
  echo "✅ $SAMPLE: $READS_INPUT → $READS_CLEAN → $READS_RRNA ($PERCENT%)"
done

echo -e "\n📁 Statistics saved to: $OUTFILE_PATH"
echo "📄 Missing files logged in: $LOGFILE"
