#!/bin/bash

set -euo pipefail

read -rp "Source directory (TransDecoder): " SRC_DIR
read -rp "Output directory (all_cds): " OUT_DIR

[[ -d "$SRC_DIR" ]] || { echo "Directory does not exist: $SRC_DIR"; exit 1; }


mkdir -p "$OUT_DIR"

echo "Searching for *.transdecoder.cds files..."

mapfile -t FILES < <(find "$SRC_DIR" -type f -name "*.transdecoder.cds")

if [[ "${#FILES[@]}" -eq 0 ]]; then
    echo "No *.transdecoder.cds files found"
    exit 1
fi


for file in "${FILES[@]}"; do
    sample=$(basename "$(dirname "$file")")
    cp "$file" "$OUT_DIR/${sample}.cds.fa"
    echo "✅ $sample"
done

echo
echo "📁 All files saved in: $OUT_DIR"
echo "🔢 Number of files: ${#FILES[@]}"
