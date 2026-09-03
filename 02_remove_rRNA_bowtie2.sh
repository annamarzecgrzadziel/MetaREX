#!/bin/bash

read -rp "Enter the directory containing fastp-processed files (clean.fq): " FQDIR
read -rp "Enter the path to the Bowtie2 index (without .bt2): " BT2_IDX
read -rp "How many CPU cores should be used [default: 8]? " THREADS
read -rp "Enter the output directory [default: ./rRNA_removed]: " OUTDIR

THREADS=${THREADS:-8}
OUTDIR=${OUTDIR:-./rRNA_removed}
LOGDIR="$OUTDIR/bowtie2_logs"

mkdir -p "$OUTDIR" "$LOGDIR"

echo "Searching for R1 files in: $FQDIR"

mapfile -t R1_FILES < <(find "$FQDIR" -type f -iname "*_R1*.f*q*" | sort)

if [ ${#R1_FILES[@]} -eq 0 ]; then
  echo "No *_R1*.f*q(.gz) files found in $FQDIR"
  exit 1
fi

echo "Found ${#R1_FILES[@]} R1 files"

for R1 in "${R1_FILES[@]}"; do
  R2="${R1/_R1/_R2}"

  if [[ ! -f "$R2" ]]; then
    echo "Paired file not found for $R1 → skipping"
    continue
  fi

  BASENAME=$(basename "$R1")
  SAMPLE=$(echo "$BASENAME" | sed -E 's/_R1.*//')

  echo "Sample: $SAMPLE"

  OUT1="$OUTDIR/${SAMPLE}.rRNAfree_R1.fastq.gz"
  OUT2="$OUTDIR/${SAMPLE}.rRNAfree_R2.fastq.gz"
  TMP_OUT1="$OUTDIR/${SAMPLE}.rRNAfree.1.fastq.gz"
  TMP_OUT2="$OUTDIR/${SAMPLE}.rRNAfree.2.fastq.gz"
  LOG="$LOGDIR/${SAMPLE}.log"

  bowtie2 -x "$BT2_IDX" \
  -1 "$R1" -2 "$R2" \
  --end-to-end --very-sensitive \
  -L 31 -N 0 --score-min L,0,-0.6 \
  --threads "$THREADS" \
  --un-conc-gz "$OUTDIR/${SAMPLE}.rRNAfree.fastq.gz" \
  -S /dev/null \
  2> "$LOG"

  if [[ -f "$TMP_OUT1" && -f "$TMP_OUT2" ]]; then
    mv "$TMP_OUT1" "$OUT1"
    mv "$TMP_OUT2" "$OUT2"
    echo "$SAMPLE: saved $OUT1 + $OUT2"
  else
    echo "$SAMPLE: output files not found (all reads = rRNA?) – skipping"
    continue
  fi
done

echo -e "\n Processing of all samples completed."
