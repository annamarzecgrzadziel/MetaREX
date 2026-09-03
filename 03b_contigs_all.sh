#!/bin/bash

set -euo pipefail

read -rp "Source directory (MEGAHIT output): " SRC_DIR
read -rp "Output directory (contigs_all): " OUT_DIR

if [[ ! -d "$SRC_DIR" ]]; then
    echo "Source directory does not exist: $SRC_DIR"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "Searching for final.contigs.fa files..."

COUNT=0

find "$SRC_DIR" -type f -name "final.contigs.fa" | while read -r file; do
    sample=$(basename "$(dirname "$file")")
    cp "$file" "$OUT_DIR/${sample}.contigs.fa"
    echo "Copied: $sample"
    COUNT=$((COUNT + 1))
done

echo
echo "All files saved in: $OUT_DIR"
echo "Number of files: $COUNT"
