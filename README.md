# Conserved Pan-Stress Transcriptomic Regulators and Interactome Hub Genes in *Arabidopsis thaliana*

An end-to-end computational pipeline for multi-condition RNA-seq meta-analysis, cross-stress consensus DEG discovery, functional enrichment, and protein interactome hub prioritization across 6 biotic and abiotic stress conditions.

---

## 1. Overview

Plants under natural conditions frequently encounter concurrent environmental stresses. This repository contains an automated, reproducible workflow used to identify the core pan-stress transcriptomic signature of *Arabidopsis thaliana* across six distinct experimental regimes:
* **Abiotic Stresses:** Heat Shock (37°C), Cold (4°C), Drought / Dehydration, Salinity (150 mM NaCl)
* **Biotic Stresses:** Bacterial pathogen (*Pseudomonas syringae* pv. *tomato* DC3000) and Viral pathogen (*Turnip Crinkle Virus*, TCV)

Using individual negative binomial differential expression modeling, random-effects meta-analysis (DerSimonian-Laird), STRING interactome topology analysis, and cytoHubba Maximal Clique Centrality (MCC), this workflow isolates the conserved multi-organelle proteostasis network required for universal stress resilience.

---

## 2. Repository Structure

```text
Arabidopsis_Stress_MetaAnalysis/
├── scripts/
│   ├── 01_universal_sra_download.sh         # Parallel SRA prefetch & fasterq-dump conversion
│   ├── 02_universal_fastqc_trimming.sh       # FastQC inspection & fastp adapter trimming
│   ├── 03_universal_alignment_quant.sh       # HISAT2 alignment & featureCounts quantification
│   ├── 04_universal_meta_analysis.R          # DESeq2, random-effects meta-analysis & intersection
│   ├── 05_universal_go_kegg_enrichment.R     # clusterProfiler GO (BP/MF/CC) & KEGG pathways
│   └── generate_figures.R                    # Publication-grade visual suite (Heatmaps, Forest plots)
├── data/
│   ├── accessions.txt                        # Target NCBI SRA run accessions
│   └── raw_fastq/                            # Raw FASTQ storage (ignored by git)
├── analysis/
│   ├── meta_analysis/                        # Pooled effect sizes & FDR outputs
│   └── Arabidopsis_Network/                  # High-res publication visuals & hub tables
├── outputs/
│   └── functional_enrichment/                # GO & KEGG tabular results
└── README.md
