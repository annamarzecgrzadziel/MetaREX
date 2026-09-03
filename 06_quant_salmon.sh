#!/bin/bash

set -euo pipefail

echo "Activating Conda environment: metatrascriptomics_base"
eval "$(command conda 'shell.bash' 'hook' 2>/dev/null)"
conda activate metatrascriptomics_base

read -rp "Directory containing CDS files (*.cds): " CDS_DIR
read -rp "Directory containing FASTQ files (*.rRNAfree_R*.fq.gz / *.fastq.gz): " FQ_DIR
read -rp "Output directory: " OUT_DIR
read -rp "Threads per sample [default: 8]: " THREADS
read -rp "Maximum number of parallel samples [default: 4]: " MAX_PARALLEL

THREADS=${THREADS:-8}
MAX_PARALLEL=${MAX_PARALLEL:-4}

[[ -d "$CDS_DIR" ]] || { echo "Directory does not exist: $CDS_DIR"; exit 1; }
[[ -d "$FQ_DIR" ]] || { echo "Directory does not exist: $FQ_DIR"; exit 1; }

mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR/logs"

mapfile -t CDS_FILES < <(find "$CDS_DIR" -type f | grep -E '\.cds(\.fa)?$' | sort)

if [[ "${#CDS_FILES[@]}" -eq 0 ]]; then
  echo "No *.cds or *.cds.fa files found in: $CDS_DIR"
  exit 1
fi

find_fastq() {
  local dir="$1"
  local pattern1="$2"
  local pattern2="$3"

  find "$dir" -type f \( -name "$pattern1" -o -name "$pattern2" \) | head -n 1 || true
}

run_sample() {
  local CDS="$1"

  local SAMPLE
  SAMPLE=$(basename "$CDS" | sed -E 's/\.contigs.*$//')

  local R1 R2
  R1=$(find_fastq "$FQ_DIR" "${SAMPLE}.rRNAfree_R1.fq.gz" "${SAMPLE}.rRNAfree_R1.fastq.gz")
  R2=$(find_fastq "$FQ_DIR" "${SAMPLE}.rRNAfree_R2.fq.gz" "${SAMPLE}.rRNAfree_R2.fastq.gz")

  if [[ -z "$R1" || -z "$R2" || ! -f "$R1" || ! -f "$R2" ]]; then
    echo "R1/R2 files not found for $SAMPLE – skipping"
    echo "$SAMPLE" >> "$OUT_DIR/logs/missing_fastq_samples.txt"
    return 0
  fi

  local WORKDIR="$OUT_DIR/$SAMPLE"
  mkdir -p "$WORKDIR"

  echo "[$SAMPLE] Start quant"

  {
    echo "=== SAMPLE: $SAMPLE ==="
    echo "CDS: $CDS"
    echo "R1:  $R1"
    echo "R2:  $R2"
    echo
    echo "[$(date '+%F %T')] Building index..."
  } > "$WORKDIR/run.log"

  if ! salmon index -t "$CDS" -i "$WORKDIR/salmon_index" \
      > "$WORKDIR/index.log" 2>&1; then
    echo "[$SAMPLE] Error during salmon index"
    echo "FAILED_INDEX" > "$WORKDIR/status.txt"
    return 1
  fi

  echo "[$(date '+%F %T')] Running quantification..." >> "$WORKDIR/run.log"

  if ! salmon quant \
      -i "$WORKDIR/salmon_index" \
      -l A \
      -1 "$R1" -2 "$R2" \
      -p "$THREADS" \
      -o "$WORKDIR" \
      > "$WORKDIR/quant.log" 2>&1; then
    echo "[$SAMPLE] Error during salmon quant"
    echo "FAILED_QUANT" > "$WORKDIR/status.txt"
    return 1
  fi

  if [[ -f "$WORKDIR/quant.sf" ]]; then
    echo "OK" > "$WORKDIR/status.txt"
    echo "[$SAMPLE] Completed"
  else
    echo "FAILED_NO_QUANT_SF" > "$WORKDIR/status.txt"
    echo "[$SAMPLE] quant.sf nie powstał"
    return 1
  fi
}

running_jobs=0

for CDS in "${CDS_FILES[@]}"; do
  run_sample "$CDS" &
  ((running_jobs+=1))

  if (( running_jobs >= MAX_PARALLEL )); then
    wait -n || true
    ((running_jobs-=1))
  fi
done

wait || true

OK_COUNT=0
FAIL_COUNT=0

for status_file in "$OUT_DIR"/*/status.txt; do
  [[ -f "$status_file" ]] || continue
  status=$(cat "$status_file")
  if [[ "$status" == "OK" ]]; then
    ((OK_COUNT+=1))
  else
    ((FAIL_COUNT+=1))
  fi
done

echo
echo "Processing completed."
echo "Successful: $OK_COUNT"
echo "Failed:     $FAIL_COUNT"

if [[ -f "$OUT_DIR/logs/missing_fastq_samples.txt" ]]; then
  echo "Samples without FASTQ files were saved to: $OUT_DIR/logs/missing_fastq_samples.txt"
fi
