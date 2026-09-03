#!/bin/bash

set -u

read -rp "📂 Directory containing FASTQ files (*.R1/R2.fastq.gz): " FQDIR
read -rp "📦 Path to Bowtie2 index (without .bt2): " BT2_IDX
read -rp "📁 Output directory: " OUTDIR
read -rp "⚙️ Number of threads [default: 8]: " THREADS
read -rp "📄 Output file name [default: rrna_content_summary.tsv]: " OUTFILE

THREADS=${THREADS:-8}
OUTFILE=${OUTFILE:-rrna_content_summary.tsv}

mkdir -p "$OUTDIR"/logs

OUTFILE_PATH="${OUTDIR}/${OUTFILE}"

echo -e "Sample\tTotal_Reads\trRNA_Aligned\tmRNA_Unaligned\tPercent_rRNA\tPercent_mRNA" > "$OUTFILE_PATH"

echo "🔎 Searching for R1+R2 pairs in: $FQDIR"
mapfile -t R1_LIST < <(find "$FQDIR" -type f -iname "*_R1*.f*q.gz" | sort)

for R1 in "${R1_LIST[@]}"; do
  R2="${R1/_R1/_R2}"

  if [[ ! -f "$R2" ]]; then
    echo "Paired file not found for $R1 – skipping"
    continue
  fi

  SAMPLE=$(basename "$R1" | sed -E 's/_R1.*//')
  echo "Sample: $SAMPLE"

  LOG="${OUTDIR}/logs/bowtie2_${SAMPLE}.log"

  bowtie2 -x "$BT2_IDX" \
    -1 "$R1" -2 "$R2" \
    --end-to-end --very-sensitive -L 31 -N 0 --score-min L,0,-0.6 \
    -p "$THREADS" \
    -S /dev/null \
    2> "$LOG"

  TOTAL=$(grep "were paired" "$LOG" | awk '{print $1}')
  ALIGNED=$(grep "aligned exactly 1 time" "$LOG" | awk '{s+=$1} END {print s}')
  ALIGNED2=$(grep "aligned >1 times" "$LOG" | awk '{s+=$1} END {print s}')
  UNALIGNED=$(grep "aligned 0 times" "$LOG" | awk '{print $1}')

  ALIGNED_TOTAL=$((ALIGNED + ALIGNED2))

  if [[ -z "$TOTAL" || "$TOTAL" -eq 0 ]]; then
    echo "No data found for sample $SAMPLE – skipping"
    continue
  fi

  PCT_RRNA=$(awk -v a=$ALIGNED_TOTAL -v b=$TOTAL 'BEGIN {printf "%.2f", (a/b)*100}')
  PCT_MRNA=$(awk -v a=$UNALIGNED -v b=$TOTAL 'BEGIN {printf "%.2f", (a/b)*100}')

  echo -e "$SAMPLE\t$TOTAL\t$ALIGNED_TOTAL\t$UNALIGNED\t$PCT_RRNA\t$PCT_MRNA" >> "$OUTFILE_PATH"

  echo "✅ $SAMPLE → rRNA: $PCT_RRNA%, mRNA: $PCT_MRNA%"
done

echo -e "\n📁 Results saved to: $OUTFILE_PATH"
echo "📂 Logs: ${OUTDIR}/logs/"
