#!/bin/bash


echo "Activating Conda environment: metagenomics_base"
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate metagenomics_base

echo "Enter the full path to the directory containing *_R1_*.fastq.gz files:"
read -r INPUT_DIR

echo "Enter the path to the output directory (FASTQ + reports):"
read -r OUTPUT_DIR

mkdir -p "$OUTPUT_DIR/fastq_clean"
mkdir -p "$OUTPUT_DIR/reports"

echo -e "\n Select data cleaning strategy (trimming):"
echo "----------------------------------------------------------"
echo "1. Adapter removal only"
echo "   Fast and safe – for high-quality data"
echo "   --detect_adapter_for_pe"
echo
echo "2. Adapters + end trimming (cut_front + cut_tail)"
echo "   For data with low-quality read ends"
echo "   --detect_adapter_for_pe --cut_front --cut_tail"
echo
echo "3. Adapters + quality trimming (Phred + length)"
echo "   When you want control over read quality and length"
echo "   --detect_adapter_for_pe --qualified_quality_phred --length_required"
echo
echo "4. Adapters + ends + quality trimming (full trimming)"
echo "   For challenging datasets or publication-quality processing"
echo "   --detect_adapter_for_pe --cut_front --cut_tail --qualified_quality_phred --length_required"
echo
echo "5. Custom fastp parameters"
echo "   For advanced users"
echo "----------------------------------------------------------"

read -p "Select an option (1-5): " STRATEGY

case $STRATEGY in
  1)
    TRIM_OPTS="--detect_adapter_for_pe"
    ;;
  2)
    TRIM_OPTS="--detect_adapter_for_pe --cut_front --cut_tail"
    ;;
  3)
    read -p "Minimum Phred quality [default: 20]: " PHRED
    read -p "Minimum read length [default: 50]: " MINLEN
    PHRED=${PHRED:-20}
    MINLEN=${MINLEN:-50}
    TRIM_OPTS="--detect_adapter_for_pe --qualified_quality_phred $PHRED --length_required $MINLEN"
    ;;
  4)
    read -p "Minimum Phred quality [default: 20]: " PHRED
    read -p "Minimum read length [default: 50]: " MINLEN
    PHRED=${PHRED:-20}
    MINLEN=${MINLEN:-50}
    TRIM_OPTS="--detect_adapter_for_pe --cut_front --cut_tail --qualified_quality_phred $PHRED --length_required $MINLEN"
    ;;
  5)
     echo "Enter all fastp parameters manually:"
    read -r TRIM_OPTS
    ;;
  *)
    echo "Invalid selection"
    exit 1
    ;;
esac

echo "Running Fastp..."

for fq1 in "$INPUT_DIR"/*_R1.fastq.gz; do
    fq2="${fq1/_R1/_R2}"
    [[ ! -f "$fq2" ]] && echo "Paired file not found for $fq1 – skipping" && continue

    sample=$(basename "$fq1" | sed 's/_R1.*.fastq.gz//')
    echo "🔬 Sample: $sample"

    fastp \
        -i "$fq1" -I "$fq2" \
        -o "$OUTPUT_DIR/fastq_clean/${sample}_R1_clean.fastq.gz" \
        -O "$OUTPUT_DIR/fastq_clean/${sample}_R2_clean.fastq.gz" \
        -h "$OUTPUT_DIR/reports/${sample}_fastp.html" \
        -j "$OUTPUT_DIR/reports/${sample}_fastp.json" \
        -w 8 \
        $TRIM_OPTS
done

echo "Generating MultiQC report..."
multiqc "$OUTPUT_DIR/reports" -o "$OUTPUT_DIR/multiqc"

echo -e "\n✅ All done!"
echo "📦 Cleaned FASTQ files: $OUTPUT_DIR/fastq_clean/"
echo "📄 Fastp reports:       $OUTPUT_DIR/reports/"
echo "📊 MultiQC report:      $OUTPUT_DIR/multiqc/"

