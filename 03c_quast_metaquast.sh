#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /data/conda_envs/nanopore_assembly || exit 1

echo "Select analysis mode:"
echo "1) QUAST (for genomes)"
echo "2) MetaQUAST (for metagenomes)"
read -p "Select 1 or 2: " MODE

read -p "Enter the directory containing contigs (.fa or .fasta in subfolders): " INPUT_DIR

read -p "Enter the output directory for the report: " OUTDIR
mkdir -p "$OUTDIR"

THREADS=32
FILES=()
LABELS=()

while IFS= read -r -d '' FILE; do
  SAMPLE=$(basename "${FILE%.*}")
  FILES+=("$FILE")
  LABELS+=("$SAMPLE")
done < <(find "$INPUT_DIR" -type f \( -iname "*.fa" -o -iname "*.fasta" \) -print0)

if [ ${#FILES[@]} -eq 0 ]; then
  echo "No .fa or .fasta files found in directory $INPUT_DIR"
  exit 1
fi

if [[ "$MODE" == "1" ]]; then
  echo "Running QUAST..."
  quast.py -t $THREADS --min-contig 200 -o "$OUTDIR" \
    --labels "$(IFS=,; echo "${LABELS[*]}")" \
    "${FILES[@]}"
elif [[ "$MODE" == "2" ]]; then
  echo "Running MetaQUAST..."
  metaquast.py -t $THREADS --min-contig 200 -o "$OUTDIR" \
    --labels "$(IFS=,; echo "${LABELS[*]}")" \
    "${FILES[@]}"
else
  echo "Invalid mode selection. Aborting."
  exit 1
fi

conda activate /data/conda_envs/bioinfo_base || echo "Failed to activate the bioinfo_base environment — MultiQC may not work."
echo "Running MultiQC on $OUTDIR"
multiqc "$OUTDIR" -n multiqc_quast.html -o "$OUTDIR"

echo "QUAST/MetaQUAST report: $OUTDIR/report.html"
echo "MultiQC report: $OUTDIR/multiqc_quast.html"

