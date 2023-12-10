# CellTICS: an explainable neural network for cell-type identification and interpretation based on single-cell RNA-seq data
CellTICS is a biologically interpretable neural network for (sub-) cell-type identification and interpretation based on single-cell RNA-seq data. It prioritizes marker genes with cell-type specific high and low expression score, uses a hierarchy of biological pathways for neural network construction, and applies a two-stage strategy to predict cell types and sub-cell types. CellTICS corresponds to the following paper:

Yin, Q., Chen, L.. CellTICS: an explainable neural network for cell-type identification and interpretation based on single-cell RNA-seq data, Brief Bioinform, 2024, 25(1): bbad449. https://doi.org/10.1093/bib/bbad449

scIAE corresponds to the following paper:

Yin, Q., Wang, Y., Guan, J., Ji, G.. scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data, 

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
The input of CellTICS are reference scRNA-seq data, reference label, and query data. Reference data and query data should be a gene-by-cell matrix. Reference label should be a two-column matrix representing cell type and sub-cell type of each cell. Dataset name should be specified. Pathway information, pathway hierarchy and relationship of genes and pathways should be also specified, while they are all in folder `reactome`. Then, after running CellTICS, the outputs can be obtained. One file is the predicted labels of the query data, a two-column matrix representing cell type and sub-cell type of each cell. Another files are important pathways for each cell type and sub-cell type. We also offer an evaluating function to get the ACC and macro F1 score of the prediction.

### Options
The following options should always be specified.

For `code/main.py`:

    -dataset_name                   name of dataset
    -reference_data_path            path of reference scRNA-seq data
    -query_data_path                path of query scRNA-seq data    
    -reference_label_path           path of reference scRNA-seq label
    -ensembl_pathway_relation       the relationship of all genes and pathways, given as 'reactome/Ensembl2Reactome_All_Levels.txt'
    -pathway_names                  pathway information, given as 'reactome/ReactomePathways.txt'
    -pathway_relation               the hierarchy among the pathways, given as 'reactome/ReactomePathwaysRelation.txt'

For `code/evaluate.py`:

    -true_label_path                path of true labels of query data
    -prediction_label_path          path of predicted labels of query data
    
The following options have default settings and do not have to be specified.
    
For `code/main.py`:
    
    -print_information              if the information of training procedure is printed, default: True
    -ensembl                        if genes are represented as Ensembl ID, default: True
    -species                        species which the dataset is from, default: 'mouse'
    -normalization                  if the data is normalized, default: True
    -marker                         if marker genes are extracted, default: True
    -highly_expressed_threshold     threshold for highly expressed genes, default: 0.95
    -lowly_expressed_threshold      threshold for lowly expressed genes, default: 0.9
    -n_hidden_layer                 number of hidden layers of the neural network, default: 5
    -epoch_ctp                      epochs for training the neural network for cell type, default: 10
    -epoch_subctp                   epochs for training the neural network for sub-cell type, default: 10
    -learning_rate                  learning rate, default: 0.001
    -batch_size                     batch size, default: 32
    -l2_regularization              L2 regularization parameter, default: 0.0001
    -print_cost                     if the cost of each epoch is printed, default: False
    -pathway_importance_threshold   threshold to measure the pathway is important for the cell type, default: 0.1


### An example
Run the following codes:

        $ python -u code/main.py -dataset_name 'L5MB'\
                                 -reference_data_path 'example_data/L5MB_rdata.csv'\
                                 -query_data_path 'example_data/L5MB_qdata.csv'\
                                 -reference_label_path 'example_data/L5MB_rlabel.csv'\
                                 -ensembl_pathway_relation 'reactome/Ensembl2Reactome_All_Levels.txt'\
                                 -pathway_names 'reactome/ReactomePathways.txt'\
                                 -pathway_relation 'reactome/ReactomePathwaysRelation.txt'\
                                 
The outputs, containing predicted labels and important pathways, are in folder `L5MB_results`. The name before the underline (L5MB here) is the name of dataset.

To get the ACC and macro F1 score of the prediction, run the following codes:
        
        $ python -u code/evaluate.py -true_label_path 'example_data/L5MB_qlabel.csv'\
                                     -prediction_label_path 'L5MB_results/pred_y.csv'\
