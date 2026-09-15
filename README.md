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
| S.No | Condition | BioProject / Study | Control Accessions | Treated Accessions | Assay Type | Platform |
| :---: | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | **Heat Shock** | PRJNA329188 / SRP078864 | `SRR3992505`, `SRR3992506` | `SRR3992507`, `SRR3992508` | 37°C vs. 22°C baseline | Illumina HiSeq |
| 2 | **Cold Stress** | PRJNA302482 / SRP066224 | `SRR2936281`, `SRR2936282` | `SRR2936283`, `SRR2936284` | 4°C (24h) vs. 22°C baseline | Illumina HiSeq |
| 3 | **Drought** | PRJNA419844 / SRP125732 | `SRR6321151`, `SRR6321152` | `SRR6321153`, `SRR6321154` | Progressive rosette desiccation | Illumina NextSeq |
| 4 | **Salinity (NaCl)** | PRJNA385921 / SRP106670 | `SRR5520841`, `SRR5520842` | `SRR5520843`, `SRR5520844` | 150 mM NaCl irrigation | Illumina HiSeq |
| 5 | ***P. syringae*** | PRJNA315480 / SRP071850 | `SRR3223845`, `SRR3223846` | `SRR3223847`, `SRR3223848` | *Pst* DC3000 infiltration | Illumina HiSeq |
| 6 | **TCV (Viral)** | PRJNA329188 / SRP078864 | `SRR3992501`, `SRR3992502` | `SRR3992503`, `SRR3992504` | Systemic Turnip Crinkle Virus | Illumina HiSeq |
 **the hubb genes are**
| Rank | TAIR Locus | Symbol | Compartment | Degree | Pooled $\log_2\text{FC}$ | Meta FDR | Primary Biological Function |
| :---: | :--- | :--- | :--- | :---: | :---: | :---: | :--- |
| 1 | AT4G17490 | **CLPB1** | Cytoplasm / Nucleus | 13 | +1.04 | 4.01E-04 | Disaggregation & acquired thermotolerance |
| 2 | AT5G47230 | **BIP2** | ER Lumen | 11 | +0.97 | 0.3555* | Unfolded Protein Response (UPR) sensor |
| 3 | AT1G76600 | **Hsp90.1** | Cytoplasm / Nucleus | 8 | +1.22 | 2.64E-02 | Kinase & NLR immune receptor stabilization |
| 4 | AT1G13260 | **CPN60B1** | Chloroplast Stroma | 9 | +1.49 | 4.98E-02 | Plastidial GroEL chaperonin (RuBisCO folding) |
| 5 | AT5G27420 | **Hsp70-1** | Cytoplasm / Nucleus | 8 | +1.20 | 5.48E-02 | Nascent polypeptide chaperone triage |
| 6 | AT4G27410 | **Hsp21** | Chloroplast Thylakoid | 10 | +1.71 | 1.88E-02 | Photosystem II photoprotection |
| 7 | AT3G07350 | **ROF1** | Cytoplasm | 5 | +0.78 | 8.12E-02 | Hsp90 co-chaperone (HSFA2 stress memory) |
| 8 | AT2G47190 | **HSP18.2** | Cytoplasm | 5 | +1.11 | 3.81E-02 | Small HSP preventing irreversible aggregation |
| 9 | AT5G61600 | **HSP70-15** | Cytosol / Mitochondria | 5 | +0.89 | 1.14E-01 | Specialized organellar chaperone |
| 10 | AT1G21910 | **ROF2** | Cytoplasm | 6 | -0.03 | 0.9817* | Antagonistic repressor of ROF1/Hsp90 signaling |
