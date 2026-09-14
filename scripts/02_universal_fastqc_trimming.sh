#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="../data/raw_fastq"
TRIM_DIR="../data/trimmed_fastq"
QC_DIR="../data/qc_reports"
THREADS=4

mkdir -p "${TRIM_DIR}" "${QC_DIR}"

echo "[+] Running FastQC on raw files..."
fastqc -t "${THREADS}" "${INPUT_DIR}"/*.fastq.gz -o "${QC_DIR}"

for R1 in "${INPUT_DIR}"/*_1.fastq.gz; do
    if [[ -f "${R1}" ]]; then
        BASE=$(basename "${R1}" _1.fastq.gz)
        R2="${INPUT_DIR}/${BASE}_2.fastq.gz"
        
        echo "[+] Trimming paired-end library: ${BASE}"
        fastp -i "${R1}" -I "${R2}" \
              -o "${TRIM_DIR}/${BASE}_1.trimmed.fastq.gz" \
              -O "${TRIM_DIR}/${BASE}_2.trimmed.fastq.gz" \
              --thread "${THREADS}" --html "${QC_DIR}/${BASE}_fastp.html"
    fi
done

for SE in "${INPUT_DIR}"/*.fastq.gz; do
    if [[ -f "${SE}" ]] && [[ "${SE}" != *"_1"* ]] && [[ "${SE}" != *"_2"* ]]; then
        BASE=$(basename "${SE}" .fastq.gz)
        echo "[+] Trimming single-end library: ${BASE}"
        fastp -i "${SE}" -o "${TRIM_DIR}/${BASE}.trimmed.fastq.gz" \
              --thread "${THREADS}" --html "${QC_DIR}/${BASE}_fastp.html"
    fi
done

echo "[✓] Quality control and trimming completed."
