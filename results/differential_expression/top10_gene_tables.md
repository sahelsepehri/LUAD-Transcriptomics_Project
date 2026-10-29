# Differential Expression: Top Gene Tables

Selection criteria: adjusted P-value < 0.05 and |log2FC| >= 1, consistent with the thresholds used for differential-expression classification.

Genes are ranked by adjusted P-value, with absolute log2FC used to break ties. When multiple features map to the same gene symbol, only the feature with the smallest adjusted P-value is retained. Symbols found in both directions are excluded pending review.

## Top 10 Upregulated Genes

| Gene | log2FC | P-value | Adjusted P-value |
| --- | --- | --- | --- |
| SPP1 | 4.340 | 3.08e-29 | 4.45e-26 |
| GOLM1 | 1.500 | 1.40e-19 | 1.76e-17 |
| PYCR1 | 1.131 | 2.28e-19 | 2.72e-17 |
| KIF26B | 1.150 | 1.08e-18 | 1.15e-16 |
| ABCC3 | 1.857 | 2.07e-18 | 2.08e-16 |
| SPTBN2 | 1.114 | 5.61e-18 | 5.19e-16 |
| HIST1H2BM | 1.644 | 1.03e-17 | 9.06e-16 |
| B3GNT3 | 1.269 | 1.08e-17 | 9.44e-16 |
| NME1-NME2///NME1 | 1.089 | 1.29e-17 | 1.11e-15 |
| TERC | 1.655 | 6.81e-17 | 5.19e-15 |

## Top 10 Downregulated Genes

| Gene | log2FC | P-value | Adjusted P-value |
| --- | --- | --- | --- |
| SERTM1 | -1.818 | 7.27e-40 | 2.42e-35 |
| ITLN2 | -1.840 | 6.16e-38 | 1.02e-33 |
| GRK5 | -1.662 | 2.91e-37 | 3.23e-33 |
| C10orf67 | -1.690 | 1.00e-36 | 8.31e-33 |
| ANGPTL7 | -1.230 | 2.17e-34 | 1.45e-30 |
| FGFBP2 | -1.686 | 3.98e-32 | 1.89e-28 |
| PTPN21 | -1.660 | 1.03e-31 | 3.42e-28 |
| SASH1 | -1.609 | 9.56e-32 | 3.42e-28 |
| ADRA1A | -1.013 | 8.31e-32 | 3.42e-28 |
| GRIA1 | -2.071 | 1.78e-31 | 5.38e-28 |

## Interpretation

Positive log2FC indicates higher expression in Tumor than Normal samples; negative log2FC indicates lower expression in Tumor samples. These genes are candidates for downstream enrichment and protein interaction analyses.
