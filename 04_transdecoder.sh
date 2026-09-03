#!/bin/bash

echo "Activating Conda environment: metatrascriptomics_base"
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate metatrascriptomics_base

set -euo pipefail

read -rp "Directory containing *.fa files: " INDIR
read -rp "Output directory: " OUTDIR
read -rp "Maximum number of parallel jobs [default: 4]: " MAX_PARALLEL
MAX_PARALLEL=${MAX_PARALLEL:-4}

mkdir -p "$OUTDIR"

mapfile -t CONTIGS < <(find "$INDIR" -type f -name "*.fa" | sort)

for FILE in "${CONTIGS[@]}"; do
  SAMPLE=$(basename "$FILE" | sed 's/\.fa$//')
  WORKDIR="$OUTDIR/$SAMPLE"
  mkdir -p "$WORKDIR"
  cp "$FILE" "$WORKDIR/$SAMPLE.transcripts.fa"

  (
    echo "[$SAMPLE] Start"
    cd "$WORKDIR"
    TransDecoder.LongOrfs -t "$SAMPLE.transcripts.fa"
    TransDecoder.Predict  -t "$SAMPLE.transcripts.fa"
    echo "[$SAMPLE] Finished"
  ) &

  while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do
    sleep 1
  done
done

wait
echo -e "\n All samples have been processed."
