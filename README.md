# Imaging-transcriptomics

This repository contains the code to:
a) Calculate_t-map.py
Computes region-wise marker statistics (t-values) using linear regression to adjust for covariates such as age, sex, and education.

b) PLS.py
Performs partial least squares (PLS) regression between the regional t-map and gene expression data from the Allen Human Brain Atlas (AHBA).

c) Analysis_of_AD-related_gene_from_Genecards.py
Analyzes the association between PLS-derived gene weights and AD-related genes obtained from GeneCards (https://www.genecards.org).

d) Specificity_analysis.R
Assesses gene specificity from three perspectives: cell type, tissue region, and developmental stage. Cell-type datasets are sourced from the Allen Cell Types Database (https://celltypes.brain-map.org/), DroNC-seq, and the Karolinska Institute. Gene expression data with hippocampal subfield and entorhinal cortex annotations were obtained from GSE278723 (https://www.ncbi.nlm.nih.gov/). Gene expression data across different developmental stages of the human brain were sourced from the BrainSpan Atlas (https://www.brainspan.org/static/download.html).
