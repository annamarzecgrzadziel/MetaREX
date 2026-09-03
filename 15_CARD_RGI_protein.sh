#!/bin/bash

set -euo pipefail

echo "CARD pipeline PER SAMPLE started"

echo "Activating conda environment: card"
eval "$(conda shell.bash hook)"
conda activate card

if ! command -v rgi &> /dev/null; then
  echo " RGI not found in PATH"
  exit 1
fi

read -rp " Enter the salmon_TPM_all.tsv file: " TPM_FILE
read -rp " Enter the directory containing per-sample assemblies (03_assembly_megahit): " ASM_DIR
read -rp " Enter the CARD output directory [15_CARD]: " OUTDIR
read -rp " Number of threads [32]: " THREADS

OUTDIR=${OUTDIR:-15_CARD}
THREADS=${THREADS:-32}

mkdir -p "${OUTDIR}"

SAMPLES=$(head -n 1 "${TPM_FILE}" | tr '\t' '\n' | tail -n +2)

echo "${SAMPLES}"

for SAMPLE in ${SAMPLES}; do

  echo "================================================="
  echo " Processing sample: ${SAMPLE}"
  echo "================================================="

  SAMPLE_OUT="${OUTDIR}/${SAMPLE}"
  mkdir -p "${SAMPLE_OUT}"

  echo " Filtering transcripts with TPM > 1"

  awk -v sample="${SAMPLE}" '
    BEGIN { FS=OFS="\t" }
    NR==1 {
      for (i=1;i<=NF;i++) if ($i==sample) col=i
      next
    }
    $col > 1 { print $1 }
  ' "${TPM_FILE}" > "${SAMPLE_OUT}/${SAMPLE}_TPM_gt1.p1_ids.txt"

  N_IDS=$(wc -l < "${SAMPLE_OUT}/${SAMPLE}_TPM_gt1.p1_ids.txt")
  echo " Transcripts TPM > 1: ${N_IDS}"

  if [[ "${N_IDS}" -eq 0 ]]; then
    echo " No expressed transcripts for ${SAMPLE}, skipping"
    continue
  fi

  sed 's/\.p[0-9]\+$//' \
    "${SAMPLE_OUT}/${SAMPLE}_TPM_gt1.p1_ids.txt" \
    | sort -u \
    > "${SAMPLE_OUT}/${SAMPLE}_TPM_gt1_contigs.txt"

  ASM_FASTA="${ASM_DIR}/${SAMPLE}/final.contigs.fa"

  if [[ ! -f "${ASM_FASTA}" ]]; then
    echo " Assembly FASTA not found: ${ASM_FASTA}"
    continue
  fi

  echo " Extracting FASTA for expressed contigs"

  seqtk subseq \
    "${ASM_FASTA}" \
    "${SAMPLE_OUT}/${SAMPLE}_TPM_gt1_contigs.txt" \
    > "${SAMPLE_OUT}/${SAMPLE}_TPM_gt1.fasta"

 echo " Running TransDecoder"

echo "🔧 Activating conda environment: metatrascriptomics_base"
conda activate metatrascriptomics_base

cd "${SAMPLE_OUT}"

TransDecoder.LongOrfs \
  -t "${SAMPLE}_TPM_gt1.fasta"

TransDecoder.Predict \
  -t "${SAMPLE}_TPM_gt1.fasta"

cd - > /dev/null

conda activate card


 echo " Running CARD (RGI)"

  cd "${SAMPLE_OUT}"

  rgi \
    -i "${SAMPLE}_TPM_gt1.fasta.transdecoder.pep" \
    -o "card_amr_${SAMPLE}" \
    -t protein \
    -n "${THREADS}" \
    -e loose

  cd - > /dev/null

  echo " Sample ${SAMPLE} finished"

done

echo " CARD pipeline PER SAMPLE DONE"
