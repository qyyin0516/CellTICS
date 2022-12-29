# CellTICS: an interpretable neural network for cell type identification based on single-cell RNA-seq data
CellTICS is a biologically interpretable neural network for (sub-) cell type identification based on single-cell RNA-seq data. It prioritizes marker genes with cell-type specific high and low expression score, uses a hierarchy of biological pathways for neural network construction, and applies a two-stage strategy to predict cell types and sub-cell types. CellTICS corresponds to the following paper:

Yin, Q., Chen, L.. CellTICS: an interpretable neural network for cell type identification based on single-cell RNA-seq data, under review.

## Dependencies
CellTICS is built with Python 3 (>= 3.9.2) with the following packages:

* numpy >= 1.23.1
* pandas >= 1.2.3
* keras >= 2.9.0
* tensorflow >= 2.9.0
* networkx >= 2.6.3
* biomart >= 0.9.2

## Installation
Clone the github repository and enter CellTICS directory with

    $ git clone https://github.com/qyyin0516/CellTICS.git
    $ cd CellTICS
  
However, `CellTICS/reactome/Ensembl2Reactome_All_Levels.txt` and `CellTICS/example_data/example_data.zip` are stored with Git LFS because they are larger than 25MB. Please download the two files directly via the github page. After downloading them, please put `Ensembl2Reactome_All_Levels.txt` to `CellTICS/reactome`. Then, after unzipping `example_data.zip`, please put the unzipped folder `example_data` to `CellTICS`. Sorry for any inconvenience! 

## Usage
