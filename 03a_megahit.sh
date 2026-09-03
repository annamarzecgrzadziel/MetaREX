#!/bin/bash

echo "Activating Conda environment: metagenomics_base"
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate metagenomics_base

set -euo pipefail

read -rp "Directory containing rRNA-free files: " INDIR
read -rp "Output directory: " OUTDIR
read -rp "Number of threads [default: 16]: " THREADS
read -rp "List of k-mers (e.g. 21,41,61) [default: 21,41,61,81,101,121]: " KLIST

THREADS=${THREADS:-16}
KLIST=${KLIST:-21,41,61,81,101,121}

mkdir -p "$OUTDIR"

echo "🔍 Szukam par plików .rRNAfree_R1/R2.fastq.gz w: $INDIR"
mapfile -t R1_FILES < <(find "$INDIR" -type f -iname "*.rRNAfree_R1.fastq.gz" | sort)

for R1 in "${R1_FILES[@]}"; do
  R2="${R1/_R1/_R2}"
  SAMPLE=$(basename "$R1" | sed 's/.rRNAfree_R1.fastq.gz//')
  OUT="$OUTDIR/$SAMPLE"

  if [[ ! -f "$R2" ]]; then
    echo "Paired file not found for $SAMPLE – skipping"
    continue
  fi

  echo "Assembling: $SAMPLE"
  megahit -1 "$R1" -2 "$R2" -o "$OUT" -t "$THREADS" \
    --k-list "$KLIST" \
    --min-contig-len 300

   echo "Assembly completed: $OUT/final.contigs.fa"
done

echo -e "\n Transcriptome assembly completed."
