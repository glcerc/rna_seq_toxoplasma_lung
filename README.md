# RNA-seq Analysis of Interferon Signaling During Toxoplasma gondii Infection
This repository contains an RNA-seq analysis workflow investigating how Type I and Type II interferon receptor signaling shapes the transcriptional response to Toxoplasma gondii infection in mouse lung tissue.
The analysis is based on publicly available bulk RNA-seq data and was performed as part of an RNA-seq course project.

## Project Overview
Biological question: How does the absence of interferon receptor signaling (IFNAR/IFNGR) affect gene expression during T. gondii infection in lung tissue?

Approach:
Bulk RNA-seq analysis

Genotype × infection interaction model

Functional enrichment analysis of differentially expressed genes

Organism: Mus musculus
Tissue: Lung
Data source: GEO (GSE119855)

🧪 Experimental Design
| Genotype  | Condition  | Samples |
|-----------|------------|---------|
| WT        | Control    | 3       |
| WT        | Infection  | 3       |
| DoubleKO  | Control    | 3       |
| DoubleKO  | Infection  | 6       |


Total: 15 samples, paired-end RNA-seq.

## Analysis Workflow

### Quality Control

  FastQC
  
  MultiQC

### Read Alignment

  HISAT2 (RF strandedness)
  
  Reference: Mus musculus GRCm39 (Ensembl)

### Read Quantification

  FeatureCounts (gene-level counts)

### Exploratory Data Analysis

  VST normalization
  
  PCA
  
  Sample distance and correlation heatmaps

### Differential Expression Analysis

  DESeq2
  
  Design: ~ genotype + condition + genotype:condition

### Functional Enrichment

  GO Biological Process enrichment
  
  clusterProfiler

## Key Results
🔹 Differential Expression

570 genes showed genotype-dependent infection responses

258 genes: stronger response in DoubleKO

312 genes: stronger response in WT

🔹 Functional Enrichment

Enriched GO terms include:

Interferon signaling pathways

Antigen processing and presentation

Inflammatory and immune response processes

## Figures

### Figure 1. PCA of VST-transformed gene expression

Principal component analysis (PCA) of variance-stabilized gene expression values showing clustering of samples based on genotype and infection status. Each point represents one sample.

<img width="664" height="664" alt="basicpca" src="https://github.com/user-attachments/assets/319d23e2-b2eb-428f-9f09-05fc528b65b5" />

### Figure 2. Volcano plot of genotype × infection interaction

Volcano plot displaying differential gene expression results for the genotype × infection interaction term. Each point represents a gene, plotted by log2 fold change and adjusted p-value.

<img width="530" height="454" alt="genotype×infection_interaction" src="https://github.com/user-attachments/assets/98a424f8-9e12-4411-971d-44cd5ba3b473" />

### Figure 3. Normalized expression of selected interferon-stimulated genes

Boxplots showing DESeq2-normalized counts for Gbp5, Ifit1, and Isg15 across experimental groups. These genes display genotype-dependent infection responses.

<img width="530" height="454" alt="interferon-stimulated_genes" src="https://github.com/user-attachments/assets/97c36676-3e66-4848-ab7b-6b6a1b4cef02" />

### Figure 4. GO Biological Process enrichment analysis

Dot plot of significantly enriched Gene Ontology Biological Process terms among genes showing genotype-dependent infection responses. Dot size indicates gene count and color represents adjusted p-value.

<img width="865" height="547" alt="GO_Biological_Process_enrichment" src="https://github.com/user-attachments/assets/e4a110ea-639a-4830-963a-4db0c682528e" />

All figures are generated from the analysis scripts in this repository.
## Notes and Limitations

This study uses bulk RNA-seq, therefore cell type–specific effects cannot be resolved.
No trimming was applied, as quality control and downstream alignment metrics did not indicate a need for it.
Results are based on transcriptome-level measurements and were not experimentally validated.

## Acknowledgements
This project was carried out as part of the RNA-seq course at the University of Bern.
