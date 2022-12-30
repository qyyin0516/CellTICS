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
Placeholder.

### Options
The following options are not null and should always be specified.

For `code/get_gene_pathways.py`:

    -ensembl_pathway_relation       the relationship of all genes and pathways, given as 'reactome/Ensembl2Reactome_All_Levels.txt'

For `code/main.py`:

    -dataset_name                   name of dataset
    -reference_data_path            path of reference scRNA-seq data
    -query_data_path                path of query scRNA-seq data    
    -reference_label_path           path of reference scRNA-seq label
    -query_label_path               path of query scRNA-seq label 
    -pathway_names                  pathway information, given as 'reactome/ReactomePathways.txt'
    -pathway_relation               the hierarchy among the pathways, given as 'reactome/ReactomePathwaysRelation.txt'
    

### An example
Run the following codes:

        $ python -u code/get_gene_pathways.py -ensembl_pathway_relation 'reactome/Ensembl2Reactome_All_Levels.txt'
        $ python -u code/main.py -dataset_name 'L5MB'\
                                 -reference_data_path 'example_data/L5MB_rdata.csv'\
                                 -query_data_path 'example_data/L5MB_qdata.csv'\
                                 -reference_label_path 'example_data/L5MB_rlabel.csv'\
                                 -query_label_path 'example_data/L5MB_qlabel.csv'\
                                 -pathway_names 'reactome/ReactomePathways.txt'\
                                 -pathway_relation 'reactome/ReactomePathwaysRelation.txt'\

The outputs, containing predicted labels and important pathways, are in folder `CellTICS/L5MB_results`. The name before the underline (L5MB here) is the name of dataset.

To get the ACC and macro F1 score of the prediction, run the following codes:
        
        $ python -u code/evaluate.py -true_label_path 'example_data/L5MB_qlabel.csv'\
                                     -prediction_label_path 'L5MB_results/pred_y.csv'\
