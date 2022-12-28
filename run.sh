#!/usr/bin/env pythoget_gene_pathways.pyn
# -*- coding: encoding -*- 
#coding=gbk
module load python
module load cuda/11.2
module load cudnn/8.1

python -u code/get_gene_pathways.py -ensembl_pathway_relation 'reactome/Ensembl2Reactome_All_Levels.txt'

python -u code/main.py -dataset_name 'L5MB'\
                       -reference_data_path 'example_data/L5MB_rdata.csv'\
                       -query_data_path 'example_data/L5MB_qdata.csv'\
                       -reference_label_path 'example_data/L5MB_rlabel.csv'\
                       -query_label_path 'example_data/L5MB_qlabel.csv'\
                       -pathway_names 'reactome/ReactomePathways.txt'\
                       -pathway_genes 'ReactomeGenes.csv'\
                       -pathway_relation 'reactome/ReactomePathwaysRelation.txt'\

python -u code/evaluate.py -true_label_path 'example_data/L5MB_qlabel.csv'\
                           -prediction_label_path 'L5MB_results/pred_y.csv'\

