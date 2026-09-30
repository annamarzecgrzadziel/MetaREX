#!/usr/bin/env bash

set -euo pipefail

echo "CARD/RGI analysis of all assembled proteins"

echo "Activating conda environment: card"
eval "$(conda shell.bash hook)"
conda activate card

if ! command -v rgi >/dev/null 2>&1; then
    echo "ERROR: rgi was not found in the card environment." >&2
    exit 1
fi

read -rp "Directory containing per-sample assemblies (03_assembly_megahit): " ASM_DIR
read -rp "CARD output directory [15_CARD]: " OUTDIR
read -rp "Number of threads [16]: " THREADS

OUTDIR=${OUTDIR:-15_CARD_all_proteins}
THREADS=${THREADS:-16}

if [[ ! -d "$ASM_DIR" ]]; then
    echo "ERROR: Assembly directory does not exist: $ASM_DIR" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
ASM_DIR=$(realpath "$ASM_DIR")
OUTDIR=$(realpath "$OUTDIR")

mapfile -t CONTIG_FILES < <(
    find "$ASM_DIR" -mindepth 2 -maxdepth 2 -type f -name "final.contigs.fa" -print | sort
)

if [[ "${#CONTIG_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No per-sample final.contigs.fa files were found under $ASM_DIR" >&2
    exit 1
fi

printf 'sample\tassembly\tprotein_file\trgi_output\n' > "$OUTDIR/CARD_manifest.tsv"

for ASM_FASTA in "${CONTIG_FILES[@]}"; do
    SAMPLE=$(basename "$(dirname "$ASM_FASTA")")
    SAMPLE_OUT="$OUTDIR/$SAMPLE"
    mkdir -p "$SAMPLE_OUT"

    ALL_CONTIGS="$SAMPLE_OUT/${SAMPLE}_all_contigs.fasta"
    PROTEIN_FASTA="$SAMPLE_OUT/${SAMPLE}_all_contigs.fasta.transdecoder.pep"
    RGI_PREFIX="card_all_${SAMPLE}"

    echo "================================================="
    echo "Processing sample: $SAMPLE"
    echo "Assembly: $ASM_FASTA"
    echo "No TPM pre-filter is applied. All assembled contigs are translated."
    echo "================================================="

    cp "$ASM_FASTA" "$ALL_CONTIGS"

    cd "$SAMPLE_OUT"
    conda activate metatrascriptomics_base
    TransDecoder.LongOrfs -t "$(basename "$ALL_CONTIGS")"
    TransDecoder.Predict -t "$(basename "$ALL_CONTIGS")" --no_refine_starts

    if [[ ! -s "$PROTEIN_FASTA" ]]; then
        echo "ERROR: TransDecoder produced no protein FASTA for $SAMPLE" >&2
        exit 1
    fi

    conda activate card
    rgi \
        -i "$(basename "$PROTEIN_FASTA")" \
        -o "$RGI_PREFIX" \
        -t protein \
        -n "$THREADS" \
        -a BLAST \
        -e YES

    cd - >/dev/null

    printf '%s\t%s\t%s\t%s\n' \
        "$SAMPLE" \
        "$ASM_FASTA" \
        "$PROTEIN_FASTA" \
        "$SAMPLE_OUT/${RGI_PREFIX}.txt" \
        >> "$OUTDIR/CARD_manifest.tsv"

    echo "Completed sample: $SAMPLE"
done

echo "CARD/RGI analysis of all assembled proteins completed."
echo "Results: $OUTDIR"
echo "Manifest: $OUTDIR/CARD_all_proteins_manifest.tsv"
