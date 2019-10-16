# Gerson_single-cell
The analysis workflow is presented here in support of Babcock & Ocasio, et al. manuscript "Single cell RNA-seq shows cellular heterogeneity and lineage expansion in a mouse model of SHH-driven medulloblastoma support resistance to SHH inhibitor therapy" 

The complete analysis workflow used in this manuscript has been made available as an R notebook and included within this repository.
Instructions for complete reproduction of the figures and data presented in the author's manuscript may be found within Single-Cell_Workflow.Rmd. Some repetitive tasks in this workflow were stored for convenience as functions within SingleCell.Utilities.R.

The raw data files used in this analysis can be accessed via NCBI's Gene Expression Omnibus (GEO) at accession GSE129730.
The workflow presented here begins with the Digital Expression Matrices (DGEs), a matrix in the format of genes X barcodes, where values are UMI-corrected read counts. These DGEs may be downloaded from GEO as a supplementary file, or may be produced from the provided fastq files following a standard single-cell pre-processing workflow (we recommend the Drop-Seq Tools, found here: https://github.com/broadinstitute/Drop-seq/releases)

To facilitate exploration of our data by the scientific community we have published a web-based visualization application in support of the data, available at: http://gershon-lab.med.unc.edu/single-cell/
The app allows users to plot gene expression over top of a t-SNE plot for all novel data sets presented in the manuscript. Our hope is that users will verify our findings, as well as explore genes of particular interest to them.

prep-genes.R contains the scripts used to split and organize the expression data (DGE) into smaller, memory-friendly matrices capable of being loaded into R on our web server.
app.R contains a single-file shiny application.
