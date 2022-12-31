import os
import argparse
from utils import *
from gene_expression import *
from pathway_hierarchy import *
from neural_network import *


parser = argparse.ArgumentParser()
parser.add_argument('-dataset_name', type=str)
parser.add_argument('-reference_data_path', type=str)
parser.add_argument('-query_data_path', type=str)
parser.add_argument('-reference_label_path', type=str)
parser.add_argument('-ensembl', type=bool, default=True)
parser.add_argument('-highly_expressed_threshold', type=float, default=0.95)
parser.add_argument('-lowly_expressed_threshold', type=float, default=0.9)
parser.add_argument('-normalization', type=bool, default=True)
parser.add_argument('-marker', type=bool, default=True)
parser.add_argument('-print_information', type=bool, default=True)

parser.add_argument('-species', type=str, default='mouse')
parser.add_argument('-ensembl_pathway_relation', type=str)
parser.add_argument('-pathway_names', type=str)
parser.add_argument('-pathway_relation', type=str)
parser.add_argument('-n_hidden_layer', type=int, default=5)

parser.add_argument('-epoch_ctp', type=int, default=10)
parser.add_argument('-epoch_subctp', type=int, default=50)
parser.add_argument('-learning_rate', type=float, default=0.001)
parser.add_argument('-batch_size', type=int, default=32)
parser.add_argument('-l2_regularization', type=float, default=0.0001)
parser.add_argument('-print_cost', type=bool, default=False)
parser.add_argument('-pathway_importance_threshold', type=float, default=0.1)


def main():
    # get gene expression
    args = parser.parse_args()
    dataset_name = args.dataset_name
    if not os.path.exists(dataset_name + '_results'):
        os.mkdir(dataset_name + '_results')
    rdata = pd.read_csv(args.reference_data_path, index_col=0)
    qdata = pd.read_csv(args.query_data_path, index_col=0)
    rlabel = pd.read_csv(args.reference_label_path)
    if not args.ensembl:
        rdata = symbol_to_ensembl(rdata, species=args.species)
        qdata = symbol_to_ensembl(qdata, species=args.species)

    train_x, test_x, train_y = get_expression(rdata,
                                              qdata,
                                              rlabel,
                                              thrh=args.highly_expressed_threshold,
                                              thrl=args.lowly_expressed_threshold,
                                              normalization=args.normalization,
                                              marker=args.marker)
    ctp_subctp = ctp_subctp_relation(train_y)
    pathway_genes = get_gene_pathways(args.ensembl_pathway_relation)

    # big cell type prediction
    if args.print_information:
        print("Current prediction is for big cell type.")
    masking, layers_node, train_x, test_x = get_masking(args.pathway_names,
                                                        pathway_genes,
                                                        args.pathway_relation,
                                                        train_x,
                                                        test_x,
                                                        train_y,
                                                        datatype='celltype',
                                                        species=args.species,
                                                        n_hidden=args.n_hidden_layer)
    pred_y_df = pd.DataFrame(data=0, index=test_x.columns, columns=list(range(2, len(masking) + 2)))
    activation_output = {}
    for output_layer in range(2, len(masking) + 2):
        if args.print_information:
            print("Current sub-neural network has " + str(output_layer - 1) + " hidden layers.")
        output_train, output_test = model(np.array(train_x),
                                          one_hot_coding(train_y, 'celltype'),
                                          np.array(test_x),
                                          layers_node,
                                          masking,
                                          output_layer,
                                          num_epochs=args.epoch_ctp,
                                          learning_rate=args.learning_rate,
                                          minibatch_size=args.batch_size,
                                          gamma=args.l2_regularization,
                                          print_cost=args.print_cost)
        for j in range(len(output_train)):
            ctp_sort = layers_node[0]
            ctp_sort.sort()
            if j != output_layer - 1:
                output_train[j + 1] = pd.DataFrame(data=output_train[j + 1],
                                                   index=layers_node[len(layers_node) - 2 - j],
                                                   columns=train_x.columns)
            else:
                output_train[j + 1] = pd.DataFrame(data=output_train[j + 1],
                                                   index=ctp_sort,
                                                   columns=train_x.columns)
        activation_output[output_layer] = output_train
        pred_y_df[output_layer] = get_prediction(output_test[output_layer],
                                                 pd.get_dummies(train_y['celltype']),
                                                 test_x,
                                                 datatype='celltype')
    pred_y = pd.DataFrame(data=0, index=test_x.columns, columns=['celltype'])
    pred_y['celltype'] = pred_y_df.T.describe().T['top']
    pred_y.insert(1, 'subcelltype', 0)
    if args.print_information:
        print("Big cell type is predicted!")
    pathways_bigctp = get_pathway_importance(train_y,
                                             activation_output,
                                             datatype='celltype',
                                             thr=args.pathway_importance_threshold)
    pathways_bigctp.to_csv(dataset_name + '_results' + "/pathways_bigctp.csv")

    # sub-cell type prediction
    pathways_subctp = {}
    for l in list(ctp_subctp.keys()):
        if args.print_information:
            print("Current prediction is for sub-cell type of " + str(l) + '.')
        train_y_sub = train_y[(train_y["celltype"] == l)]
        train_x_sub = train_x.loc[:, train_y_sub.index.tolist()]
        test_x_sub = test_x.loc[:, pred_y[(pred_y["celltype"] == l)].index.tolist()]
        if test_x_sub.shape[1] == 0:
            continue
        masking_sub, layers_node_sub, train_x_sub, test_x_sub = get_masking(args.pathway_names,
                                                                            pathway_genes,
                                                                            args.pathway_relation,
                                                                            train_x_sub,
                                                                            test_x_sub,
                                                                            train_y_sub,
                                                                            datatype='subcelltype',
                                                                            species=args.species,
                                                                            n_hidden=args.n_hidden_layer)
        pred_y_df_sub = pd.DataFrame(data=0, index=test_x_sub.columns, columns=list(range(2, len(masking_sub) + 2)))
        activation_output_sub = {}
        for output_layer_sub in range(2, len(masking_sub) + 2):
            if args.print_information:
                print("Current sub-neural network has " + str(output_layer_sub - 1) + " hidden layers.")
            output_train_sub, output_test_sub = model(np.array(train_x_sub),
                                                      one_hot_coding(train_y_sub, 'subcelltype'),
                                                      np.array(test_x_sub),
                                                      layers_node_sub,
                                                      masking_sub,
                                                      output_layer_sub,
                                                      num_epochs=args.epoch_subctp,
                                                      learning_rate=args.learning_rate,
                                                      minibatch_size=args.batch_size,
                                                      gamma=args.l2_regularization,
                                                      print_cost=args.print_cost)
            for j in range(len(output_train_sub)):
                ctp_sort_sub = layers_node_sub[0]
                ctp_sort_sub.sort()
                if j != output_layer_sub - 1:
                    output_train_sub[j + 1] = pd.DataFrame(data=output_train_sub[j + 1],
                                                           index=layers_node_sub[len(layers_node_sub) - 2 - j],
                                                           columns=train_x_sub.columns)
                else:
                    output_train_sub[j + 1] = pd.DataFrame(data=output_train_sub[j + 1], index=ctp_sort_sub,
                                                           columns=train_x_sub.columns)
            activation_output_sub[output_layer_sub] = output_train_sub
            pred_y_df_sub[output_layer_sub] = get_prediction(output_test_sub[output_layer_sub],
                                                             pd.get_dummies(train_y_sub['subcelltype']),
                                                             test_x_sub,
                                                             'subcelltype')

        pred_y.loc[pred_y_df_sub.index.tolist(), 'subcelltype'] = pred_y_df_sub.T.describe().T['top']
        pathways_subctp[l] = get_pathway_importance(train_y_sub,
                                                    activation_output_sub,
                                                    datatype='subcelltype',
                                                    bigctp=l,
                                                    pathways_bigctp=pathways_bigctp,
                                                    thr=args.pathway_importance_threshold)
        pathways_subctp[l].to_csv(dataset_name + '_results' + '/pathways_subctp_' + str(l) + '.csv')
    if args.print_information:
        print("Sub-cell type is predicted!")
    pred_y.to_csv(dataset_name + '_results' + '/pred_y.csv')


if __name__ == "__main__":
    main()
