# LUAD Transcriptomic Analysis

## Overview

This project presents a transcriptomic analysis of lung adenocarcinoma (LUAD) using publicly available gene expression data.

The main goal is to identify differentially expressed genes between lung adenocarcinoma tumor samples and normal lung tissue, followed by protein-protein interaction analysis and identification of potential hub genes.

## Dataset

- **Dataset:** GSE43458
- **Organism:** Homo sapiens
- **Platform:** GPL6244
- **Data source:** NCBI Gene Expression Omnibus (GEO)

The dataset contains gene expression profiles from lung adenocarcinoma tumor and normal lung tissue samples.

## Analysis Workflow

The project follows the workflow below:

1. Dataset retrieval and preprocessing
2. Sample annotation and group assignment
3. Quality control
4. Differential expression analysis using `limma`
5. Identification of upregulated and downregulated genes
6. Visualization of differential expression results
7. Protein-protein interaction (PPI) network analysis
8. Network analysis using STRING and Cytoscape
9. Hub gene identification using `cytoHubba`
10. Functional enrichment analysis
11. Survival analysis
12. Drug-gene interaction analysis

## Quality Control

Sample quality and expression distributions were evaluated using:

- Boxplot
- Density plot
- Principal Component Analysis (PCA)

These analyses were used to examine the distribution of expression values and assess overall sample structure before differential expression analysis.

## Differential Expression Analysis

Differential expression analysis was performed using the `limma` package in R.

The analysis was used to compare gene expression between tumor and normal samples and identify differentially expressed genes.

The following visualizations were generated:

- Volcano plot
- MA plot
- Heatmap
- Mean-variance trend
- Top upregulated genes
- Top downregulated genes

## PPI Network Analysis

Differentially expressed genes were used as input for protein-protein interaction analysis.

### STRING

The STRING database was used to construct a high-confidence PPI network.

Parameters included:

- Organism: *Homo sapiens*
- Network type: Full STRING network
- Minimum interaction score: 0.700
- FDR stringency: Medium (5%)

### Cytoscape

The STRING interaction network was imported into Cytoscape for network visualization and further topology analysis.

### cytoHubba

The cytoHubba plugin was used to identify highly ranked hub genes using the Maximum Clique Centrality (MCC) algorithm.

The top-ranked genes will be evaluated for downstream analyses.

## Downstream Analysis

Selected hub genes will be further investigated using:

- Functional enrichment analysis
- UALCAN
- Kaplan-Meier survival analysis using KM Plotter
- Drug-gene interaction analysis using DGIdb

## Project Status

**Current status:** In progress

Completed:

- Dataset retrieval
- Sample annotation
- Quality control
- Differential expression analysis
- DEG visualization
- Initial PPI network analysis
- STRING network construction

In progress:

- Cytoscape network analysis
- Hub gene identification
- Functional enrichment analysis
- Survival analysis
- Drug-gene interaction analysis

## Tools and Resources

- R
- GEOquery
- limma
- STRING
- Cytoscape
- cytoHubba
- Enrichr
- UALCAN
- KM Plotter
- DGIdb

## Repository Structure

```text
LUAD-Transcriptomics_Project/
│
├── README.md
├── LUAD_limma.Rmd
├── R/
├── Data/
├── Images/
└── Results/
