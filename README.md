# LUAD Transcriptomic Analysis

## Overview

This project presents a transcriptomic analysis of lung adenocarcinoma (LUAD) using publicly available gene expression data.

The main goal is to identify differentially expressed genes between lung adenocarcinoma tumor samples and normal lung tissue, followed by protein-protein interaction (PPI) network analysis, hub gene identification, functional enrichment analysis, survival analysis, and drug-gene interaction analysis.

## Dataset

* **Dataset:** GSE43458
* **Organism:** *Homo sapiens*
* **Platform:** GPL6244
* **Data source:** NCBI Gene Expression Omnibus (GEO)

The dataset contains gene expression profiles from lung adenocarcinoma tumor and normal lung tissue samples.

## Analysis Workflow

The project follows the workflow below:

1. Dataset retrieval and preprocessing
2. Sample annotation and group assignment
3. Quality control
4. Differential expression analysis using `limma`
5. Identification of upregulated and downregulated genes
6. Differential expression visualization
7. Protein-protein interaction (PPI) network analysis using STRING
8. Network visualization and analysis using Cytoscape
9. Hub gene identification using cytoHubba and the MCC algorithm
10. Functional enrichment analysis
11. Survival analysis
12. Drug-gene interaction analysis

## Quality Control

Sample quality and expression distributions were evaluated using:

* Boxplot
* Density plot
* Principal Component Analysis (PCA)

These analyses were used to examine the distribution of expression values and assess overall sample structure before differential expression analysis.

## Differential Expression Analysis

Differential expression analysis was performed using the `limma` package in R.

The analysis compared gene expression between tumor and normal samples to identify differentially expressed genes.

The following visualizations were generated:

* Volcano plot
* MA plot
* Heatmap
* Mean-variance trend
* Top upregulated genes
* Top downregulated genes

Upregulated and downregulated gene lists were subsequently used for downstream PPI and functional enrichment analyses.

## PPI Network Analysis

Differentially expressed genes were used as input for protein-protein interaction analysis.

### STRING

The STRING database was used to construct PPI networks for upregulated and downregulated genes.

Parameters included:

* **Organism:** *Homo sapiens*
* **Network:** Full STRING network
* **Minimum interaction score:** 0.700
* **FDR stringency:** Medium (5%)

The resulting STRING networks and cleaned gene lists are included in the repository.

### Cytoscape

STRING PPI networks were imported into Cytoscape for network visualization and topology analysis.

Network layouts were generated for the upregulated and downregulated gene networks.

### cytoHubba

The cytoHubba plugin in Cytoscape was used to identify highly ranked hub genes using the **Maximum Clique Centrality (MCC)** algorithm.

The top 10 and top 20 hub genes were identified from the analyzed network.

## Functional Enrichment Analysis

Functional enrichment analysis was performed for upregulated and downregulated gene sets.

The analysis included:

* Gene Ontology (GO)

  * Biological Process (BP)
  * Cellular Component (CC)
  * Molecular Function (MF)
* KEGG pathways
* Reactome pathways
* WikiPathways

The resulting enrichment tables are included in the repository.

## Downstream Analysis

Selected hub genes will be further investigated using:

* UALCAN
* Kaplan-Meier survival analysis using KM Plotter
* Drug-gene interaction analysis using DGIdb

## Project Status

**Current status:** In progress

### Completed

* Dataset retrieval
* Sample annotation and group assignment
* Quality control
* Differential expression analysis
* DEG visualization
* Upregulated and downregulated gene identification
* STRING PPI network construction
* Cytoscape network visualization
* cytoHubba hub gene analysis using MCC
* Top 10 and Top 20 hub gene identification
* STRING functional enrichment analysis for upregulated genes

### In Progress

* Functional enrichment analysis for downregulated genes
* Survival analysis of selected hub genes
* Drug-gene interaction analysis
* Final interpretation and integration of results

## Tools and Resources

* R
* GEOquery
* limma
* STRING
* Cytoscape
* cytoHubba
* Enrichr
* UALCAN
* KM Plotter
* DGIdb

## Repository Structure

```text

```
