#!/usr/bin/env bash
set -euo pipefail

THREADS=4
OUT_DIR="../data/raw_fastq"
TEMP_DIR="../data/sra_cache"
ACCESSION_FILE="../data/accessions.txt"

mkdir -p "${OUT_DIR}" "${TEMP_DIR}"

command -v prefetch >/dev/null 2>&1 || { echo "[-] prefetch missing. Install SRA Toolkit."; exit 1; }
command -v fasterq-dump >/dev/null 2>&1 || { echo "[-] fasterq-dump missing. Install SRA Toolkit."; exit 1; }

mapfile -t ACCESSIONS < <(grep -v '^[[:space:]]*$' "${ACCESSION_FILE}" | grep -v '^[[:space:]]*#')

for SRA_ID in "${ACCESSIONS[@]}"; do
    SRA_ID=$(echo "${SRA_ID}" | tr -d '\r' | xargs)
    echo "[+] Processing SRA run: ${SRA_ID}"
    
    prefetch --output-directory "${TEMP_DIR}" --max-size 100G "${SRA_ID}"
    
    SRA_FILE="${TEMP_DIR}/${SRA_ID}/${SRA_ID}.sra"
    if [[ ! -f "${SRA_FILE}" ]]; then SRA_FILE="${TEMP_DIR}/${SRA_ID}.sra"; fi

    fasterq-dump "${SRA_FILE}" --outdir "${OUT_DIR}" --temp "${TEMP_DIR}" --threads "${THREADS}" --split-3

    if command -v pigz >/dev/null 2>&1; then
        pigz -f -p "${THREADS}" "${OUT_DIR}/${SRA_ID}"*.fastq
    else
        gzip -f "${OUT_DIR}/${SRA_ID}"*.fastq
    fi

    rm -rf "${TEMP_DIR}/${SRA_ID}" "${TEMP_DIR}/${SRA_ID}.sra"
done

rm -rf "${TEMP_DIR}"
echo "[✓] SRA download and conversion step completed successfully."
