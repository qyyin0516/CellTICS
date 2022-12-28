import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-species', type=str, default='mouse')
parser.add_argument('-ensembl_pathway_relation', type=str)


def get_gene_pathways(input_file, species='mouse'):   # From "Ensembl2Reactome_All_Levels"(all species) to "ReactomeGenes"(mouse)
    dt = pd.read_table(input_file, header=None)
    if species == 'mouse':
        ensembl = dt[(dt[0] >= 'ENSMUSG00000000000') & (dt[0] <= 'ENSMUST99999999999')]
    elif species == 'human':
        ensembl = dt[((dt[0] >= 'ENSG00000000000000') & (dt[0] <= 'ENSG00099999999999'))
                     | ((dt[0] >= 'ENSP00000000000000') & (dt[0] <= 'ENSP00099999999999'))
                     | ((dt[0] >= 'ENST00000000000000') & (dt[0] <= 'ENST00099999999999'))]
    elif species == 'rat':
        ensembl = dt[(dt[0] >= 'ENSRNOG00000000000') & (dt[0] <= 'ENSRNOT99999999999')]

    ensembl = ensembl.iloc[:, 0:2]
    ensembl.columns = ['gene', 'group']
    ensembl = pd.DataFrame(ensembl, columns=['group', 'gene'])
    ensembl.index = range(0, ensembl.shape[0])

    if species == 'mouse':
        for i in range(0, ensembl.shape[0]):
            ensembl.iloc[i, 1] = 'ENSMUSG' + ensembl.iloc[i, 1][7:18]
    elif species == 'human':
        for i in range(0, ensembl.shape[0]):
            ensembl.iloc[i, 1] = 'ENSG000' + ensembl.iloc[i, 1][7:18]
    elif species == 'rat':
        for i in range(0, ensembl.shape[0]):
            ensembl.iloc[i, 1] = 'ENSRNOG' + ensembl.iloc[i, 1][7:18]
    ensembl = ensembl.drop_duplicates()
    ensembl.index = range(0, ensembl.shape[0])
    ensembl.to_csv("ReactomeGenes.csv")


def main():
    args = parser.parse_args()
    get_gene_pathways(args.ensembl_pathway_relation, args.species)


if __name__ == "__main__":
    main()
