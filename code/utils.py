import numpy as np
import pandas as pd
import re
import biomart


def symbol_to_ensembl(dt, species='mouse'):
    server = biomart.BiomartServer('http://useast.ensembl.org/biomart')
    if species == 'mouse':
        mart = server.datasets['mmusculus_gene_ensembl']
    elif species == 'human':
        mart = server.datasets['hsapiens_gene_ensembl']
    elif species == 'rat':
        mart = server.datasets['rnorvegicus_gene_ensembl']

    attributes = ['external_gene_name', 'ensembl_gene_id']
    response = mart.search({'attributes': attributes})
    data = response.raw.data.decode('ascii')
    transformation = {}
    for line in data.splitlines():
        line = line.split('\t')
        gene_symbol = line[0]
        ensembl = line[1]
        transformation[gene_symbol] = ensembl

    gene_dt = list(dt.columns)
    gene_db = list(transformation.keys())
    gene_intersect = list(set(gene_dt) & set(gene_db))
    ensembl = []
    for i in gene_intersect:
        ensembl.append(transformation[i])
    dt = dt.loc[:, gene_intersect]
    dt.columns = ensembl
    return dt


def ctp_subctp_relation(train_y):
    ctp_unique = train_y.iloc[:, 0].unique()
    ctp_subctp = {}
    for l in ctp_unique:
        corrsponding_cells = train_y[(train_y["celltype"] == l)]
        subctp_unique = corrsponding_cells.iloc[:, 1].unique()
        ctp_subctp[l] = subctp_unique
    return ctp_subctp


def softmax(x):
    x = np.asarray(x)
    x_col_max = x.max(axis=0)
    x_col_max = x_col_max.reshape([1,x.shape[1]])
    x = x - x_col_max
    x_exp = np.exp(x)
    x_exp_col_sum = x_exp.sum(axis=0).reshape([1,x.shape[1]])
    softmax = x_exp / x_exp_col_sum
    return softmax


def one_hot_coding(train_y, datatype):
    Y_train_df = pd.get_dummies(train_y[datatype])
    Y_train = np.array(Y_train_df).T
    return Y_train


def get_prediction(output, train_y, test_x, datatype):
    output_softmax = output
    y_pred = pd.DataFrame(data=0, index=test_x.columns, columns=[datatype])
    for i in range(output_softmax.shape[1]):
        value = softmax(output[:, i].reshape(output[:, i].shape[0], 1))
        output_softmax[:, i] = value.reshape(value.shape[0],)
        y_pred.iloc[i, 0] = train_y.columns.values[output_softmax[:, i].argmax()]
    return y_pred


def get_pathway_importance(train_y, activation_output, datatype, bigctp=None, pathways_bigctp=None, thr=0.1):
    ctp_unique = list(activation_output[2][2].index.values)
    if (len(ctp_unique) == 1) and (datatype == 'subcelltype'):
        pathways_final = pathways_bigctp[pathways_bigctp['celltype'] == bigctp]
        pathways_final.loc[:, 'celltype'] = ctp_unique[0]
        return pathways_final
    else:
        pathway_importance = {}
        for output_layer in range(2, len(activation_output) + 2):
            pathway_importance_layer = {}
            for j in range(1, output_layer):
                activation_value = activation_output[output_layer][j].T
                activation_average = pd.DataFrame(data=0, index=ctp_unique, columns=activation_value.columns)
                for k in range(len(ctp_unique)):
                    activation_ctp = activation_value[train_y[datatype] == ctp_unique[k]]
                    activation_average.iloc[k, :] = activation_ctp.mean()
                pathway_importance_layer[j] = activation_average
            pathway_importance[output_layer] = pathway_importance_layer
        pathways = {}
        for output_layer in range(2, len(activation_output) + 2):
            pathway_ctp_value = pd.DataFrame(data=None, columns=['pathway', 'celltype', 'value'])
            count = 0
            for j in range(1, output_layer):
                pathway_value = pathway_importance[output_layer][j]
                for i in range(pathway_value.shape[1]):
                    pathway_ordered = pathway_value.iloc[:, i].sort_values(ascending=False)
                    if (pathway_ordered[0] - pathway_ordered[1]) >= thr:
                        pathway_ctp_value.loc[count] = [pathway_value.columns[i], pathway_ordered.index[0],
                                                        pathway_ordered[0] - pathway_ordered[1]]
                        count = count + 1
            pathways[output_layer] = pathway_ctp_value
            pathways[output_layer] = pathways[output_layer].sort_values(by='value', ascending=False)
            for j in range(pathways[output_layer].shape[0]):
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy1') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy2') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy3') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy4') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy5') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy6') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy7') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy8') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
                pathways[output_layer].iloc[j, 0] = re.sub(re.escape('_copy9') + '$', '',
                                                           pathways[output_layer].iloc[j, 0])
            pathways[output_layer] = pathways[output_layer].drop_duplicates(subset='pathway', keep='first')
        if len(activation_output) == 1:
            pathways_final = pathways[2]
        else:
            pathways_final = pathways[2]
            for output_layer in range(3, len(activation_output) + 2):
                pathways_final = pd.concat([pathways_final, pathways[output_layer]], axis=0)
        pathways_final = pathways_final.sort_values(by='value', ascending=False)
        pathways_final = pathways_final.drop_duplicates(subset='pathway', keep='first')
        return pathways_final
