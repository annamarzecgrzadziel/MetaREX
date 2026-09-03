#!/usr/bin/env bash
set -eo pipefail

ENV_PATH="/data/conda_envs/eggnog_mapper_v2"
export EGGNOG_DATA_DIR="/home/bio/data/bazy/eggnog"

echo "Activating Conda environment: $ENV_PATH"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate "$ENV_PATH"

if ! command -v emapper.py &>/dev/null; then
  echo "emapper.py is not available!"
  exit 1
fi

VERSION=$(emapper.py --version)
echo "eggNOG-mapper: $VERSION"
echo "Database location: $EGGNOG_DATA_DIR"
echo

read -rp "Directory containing TransDecoder results (*.transdecoder.pep): " PEP_DIR
read -rp "Output directory (e.g. 08_eggnog_annotation): " OUT_DIR

CPU_DEFAULT=32
read -rp "  Number of CPU cores [default: ${CPU_DEFAULT}]: " CPU
CPU=${CPU:-$CPU_DEFAULT}

if [[ ! -f "$EGGNOG_DATA_DIR/eggnog.db" || ! -f "$EGGNOG_DATA_DIR/eggnog_proteins.dmnd" ]]; then
  echo " ERROR: Key database files (eggnog.db or .dmnd) are missing in $EGGNOG_DATA_DIR!"
  exit 1
fi

mkdir -p "$OUT_DIR"

ALL_PEP="$OUT_DIR/all_cds.pep"

mapfile -t PEP_FILES < <(find "$PEP_DIR" -type f -name "*.transdecoder.pep")

if [[ "${#PEP_FILES[@]}" -eq 0 ]]; then
  echo " No *.transdecoder.pep files found"
  exit 1
fi

cat "${PEP_FILES[@]}" > "$ALL_PEP"
PEP_COUNT=$(grep -c "^>" "$ALL_PEP")

echo -e "\n Running eggNOG-mapper2..."

mkdir -p "$OUT_DIR/tmp"

emapper.py \
  -i "$ALL_PEP" \
  --itype proteins \
  -o eggnog_cds \
  --cpu "$CPU" \
  --output_dir "$OUT_DIR" \
  --temp_dir "$OUT_DIR/tmp" \
  --go_evidence non-electronic \
  --override

echo -e "\n PROCESS COMPLETED SUCCESSFULLY"
echo " Results: $OUT_DIR/eggnog_cds.emapper.annotations"
