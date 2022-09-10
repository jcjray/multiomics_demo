Date: 2022.09.07
Title: PBMC 3K Multiome Demo

This is a demonstration project showing how to use a multiomic single cell dataset to assign cell states or types to human PBMCs. It uses public data
and open source bioinformatics packages in R and Python.

The dataset used is 10X PBMC 3K multiome granulocyte sorted v1 chemistry processed on Cell Ranger ARC 2.0.0 available at the following URL:
https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0

In this demo, we process the scATAC data in R using the package ArchR, export the results for use in Python, and process the gene expression
data using Scanpy in Python. We perform quality control filtering on the data and then integrate the modes together, comparing the signals for
markers of immune cell lineages detected by single nuclear gene expression and by scATAC.

This analysis shows that ATAC-seq is a more robust signal of cell state than gene expression. However, gene expression plays a critical role in ensuring
the quality of the data by indicating low-quality cells that may be missed by ATAC-seq alone.

The 10X data should be in the subdirectory ./pbmc_3k/

The same results can be produced by running the Snakemake pipeline.
