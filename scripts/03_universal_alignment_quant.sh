#!/usr/bin/env bash
set -euo pipefail

TRIM_DIR="../data/trimmed_fastq"
ALIGN_DIR="../data/aligned_bam"
COUNTS_DIR="../data/counts"
THREADS=4

GENOME_INDEX="../reference/genome_index"
GTF_ANNOTATION="../reference/annotation.gtf"

mkdir -p "${ALIGN_DIR}" "${COUNTS_DIR}"

for R1 in "${TRIM_DIR}"/*_1.trimmed.fastq.gz; do
    if [[ -f "${R1}" ]]; then
        BASE=$(basename "${R1}" _1.trimmed.fastq.gz)
        R2="${TRIM_DIR}/${BASE}_2.trimmed.fastq.gz"
        
        echo "[+] Aligning paired-end samples: ${BASE}"
        hisat2 -p "${THREADS}" -x "${GENOME_INDEX}" -1 "${R1}" -2 "${R2}" | \
            samtools view -bS - | samtools sort -o "${ALIGN_DIR}/${BASE}.sorted.bam"
        samtools index "${ALIGN_DIR}/${BASE}.sorted.bam"
    fi
done

echo "[+] Quantifying gene expression levels..."
featureCounts -T "${THREADS}" -p -t exon -g gene_id \
    -a "${GTF_ANNOTATION}" \
    -o "${COUNTS_DIR}/gene_counts_matrix.txt" \
    "${ALIGN_DIR}"/*.sorted.bam

echo "[✓] Alignment and quantification completed successfully."
