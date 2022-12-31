import numpy as np
import pandas as pd
import networkx as nx


def get_gene_pathways(input_file, species='mouse'):
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
    return ensembl


class Reactome():

    def __init__(self, pathway_names, pathway_genes, relations_file_name, species):
        self.pathway_names = self.load_names(pathway_names)
        self.pathway_genes = pathway_genes
        self.hierarchy = self.load_hierarchy(relations_file_name)
        self.species = species

    def load_names(self, pathway_names):
        filename = pathway_names
        df = pd.read_csv(filename, sep='\t')
        df.columns = ['reactome_id', 'pathway_name', 'species']
        return df

    def load_hierarchy(self, relations_file_name):
        filename = relations_file_name
        df = pd.read_table(filename, header=None)
        df.columns = ['child', 'parent']
        return df


class ReactomeNetwork():

    def __init__(self, pathway_names, pathway_genes, relations_file_name, species):
        self.reactome = Reactome(pathway_names, pathway_genes, relations_file_name, species)  # low level access to reactome pathways and genes
        self.netx = self.get_reactome_networkx()

    def get_reactome_networkx(self):
        if hasattr(self, 'netx'):
            return self.netx
        hierarchy = self.reactome.hierarchy
        # filter hierarchy to have mouse pathways only
        if self.reactome.species == 'mouse':
            abbr = 'MMU'
        elif self.reactome.species == 'human':
            abbr = 'HSA'
        elif self.reactome.species == 'rat':
            abbr = 'RNO'
        species_hierarchy = hierarchy[hierarchy['child'].str.contains(abbr)]
        net = nx.from_pandas_edgelist(species_hierarchy, 'child', 'parent', create_using=nx.DiGraph())
        net.name = 'reactome'

        # add root node
        roots = [n for n, d in net.in_degree() if d == 0]
        root_node = 'root'
        edges = [(root_node, n) for n in roots]
        net.add_edges_from(edges)

        return net


def add_duplicated_edges(G, node, n_levels):
    edges = []
    source = node
    for l in range(n_levels):
        target = node + '_copy' + str(l + 1)
        edge = (source, target)
        source = target
        edges.append(edge)

    G.add_edges_from(edges)
    return G, target


def add_gene_edges(G, node, pathways, genes_df):
    genes = []
    if type(pathways) == str:
        genes = genes + list(genes_df[genes_df['group'] == pathways]['gene'])
    else:
        for i in range(len(pathways)):
            genes = genes + list(genes_df[genes_df['group'] == pathways[i]]['gene'])
    genes = list(set(genes))
    edges = []
    source = node
    for target in genes:
        edge = (source, target)
        edges.append(edge)
    G.add_edges_from(edges)
    return G


def get_nodes_at_level(net, distance):
    nodes = set(nx.ego_graph(net, 'root', radius=distance))     # get all nodes within distance around the query node
    if distance >= 1.:
        nodes -= set(nx.ego_graph(net, 'root', radius=distance - 1))     # remove nodes that are not at the specified distance but closer
    return list(nodes)


def gene_mapping(gene, df):
    inter_gene = list(set(gene) & set(df['gene']))
    genedict = {}
    genelist = [df.iloc[0, 0]]
    genedict[df.iloc[0, 1]] = genelist
    for i in range(1, df.shape[0]):
        if df.iloc[i, 1] == df.iloc[i - 1, 1]:
            genelist.append(df.iloc[i, 0])
        else:
            genedict[df.iloc[i - 1, 1]] = genelist
            genelist = [df.iloc[i, 0]]
    mappingdf = pd.DataFrame(data=None, columns=['group', 'gene'])
    for j in range(len(inter_gene)):
        mappingdf_iter = {'group': genedict[inter_gene[j]],
                  'gene': [inter_gene[j]] * len(genedict[inter_gene[j]])}
        mappingdf_iter = pd.DataFrame(mappingdf_iter)
        mappingdf = pd.concat([mappingdf, mappingdf_iter])
    return mappingdf


def get_masking(pathway_names, pathway_genes, relations_file_name, train_x, test_x, train_y, datatype, species='mouse', n_hidden=5):
    reactome_net = ReactomeNetwork(pathway_names, pathway_genes, relations_file_name, species)
    genes_df = reactome_net.reactome.pathway_genes
    genes_df = gene_mapping(train_x.index.tolist(), genes_df)

    original_network = reactome_net.netx
    original_terminal_nodes = [n for n, d in original_network.out_degree() if d == 0]
    in_genes_df = [False for x in range(0, len(original_terminal_nodes))]
    while in_genes_df.count(False) > 0:
        in_genes_df = [False for x in range(0, len(original_terminal_nodes))]
        for i in range(len(original_terminal_nodes)):
            if original_terminal_nodes[i] in genes_df['group'].to_list():
                in_genes_df[i] = True
        for i in range(len(original_terminal_nodes)):
            if in_genes_df[i] == False:
                original_network.remove_node(original_terminal_nodes[i])
        original_terminal_nodes = [n for n, d in original_network.out_degree() if d == 0]

    sub_graph = nx.ego_graph(original_network, 'root', radius=n_hidden) #subgraph of neighbors centered at node "root" <= a given radius (n_level).
    sub_terminal_nodes = [n for n, d in sub_graph.out_degree() if d == 0]
    for node in sub_terminal_nodes:
        distance = len(nx.shortest_path(sub_graph, source='root', target=node)) #len of distance: num of nodes in the shortest path
        if (distance == n_hidden + 1) & (node not in original_terminal_nodes):
            part_graph = nx.ego_graph(original_network, node, radius=100)
            corresponding_terminal_nodes = [n for n, d in part_graph.out_degree() if d == 0]
            sub_graph = add_gene_edges(sub_graph, node, corresponding_terminal_nodes, genes_df)
        elif (distance == n_hidden + 1) & (node in original_terminal_nodes):
            corresponding_terminal_nodes = node
            sub_graph = add_gene_edges(sub_graph, node, corresponding_terminal_nodes, genes_df)
        elif distance <= n_hidden:
            diff = n_hidden - distance + 1
            sub_graph, copy_node = add_duplicated_edges(sub_graph, node, diff)
            sub_graph = add_gene_edges(sub_graph, copy_node, node, genes_df)
    final_network = sub_graph

    n_level = n_hidden + 2
    layers_node = {}
    for i in range(n_level):
        nodes = get_nodes_at_level(final_network, i)
        layers_node[i] = nodes
    layers_node[0] = list(set(train_y[datatype]))

    layers_relationship = []
    for i in range(n_hidden):
        nodes = get_nodes_at_level(final_network, i + 1)
        dict = {}
        for n in nodes:
            next = list(final_network.successors(n))
            dict[n] = next
        layers_relationship.append(dict)

    masking = {}
    masking[0] = pd.DataFrame(data=1, index=layers_node[0], columns=layers_node[1])
    for i in range(1, n_level-1):
        masking[i] = pd.DataFrame(data=0, index=layers_node[i], columns=layers_node[i + 1])
        for n in layers_node[i]:
            for node in layers_relationship[i - 1][n]:
                if node in layers_node[i + 1]:
                    masking[i].loc[n, node] = 1

    gene_in_network = list(masking[n_hidden].columns)
    train_x = train_x.loc[gene_in_network, :]
    test_x = test_x.loc[gene_in_network, :]
    for i in range(len(masking) - 1):
        masking[i] = np.array(masking[i + 1])
    del masking[len(masking) - 1]
    return masking, layers_node, train_x, test_x
