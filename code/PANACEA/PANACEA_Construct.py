import PANACEA_Distance as dist
import PANACEA_PEN_Dist as PEN_distance
from PANACEA_constant import *
from PANACEA_plot import plot_features
import pandas as pd


# find genes for cancer subtypes:
def find_subgene(network, tumour_type, cancer_name):  # tumour_type: 'prostate', cancer_name:'Prostate Cancer'
    onco_ = []
    target_ = []
    # find oncogenes for each cancer subtype from cancer gene census dataframe ()
    for index, row in cancer_df.iterrows():
        if tumour_type.lower() in str(row['Tumour Types(Somatic)']).lower() or tumour_type.lower() in str(
                row['Tumour Types(Germline)']).lower():
            if row['Gene Symbol'] in network.nodes():  # remove those not in signaling network
                onco_.append(row['Gene Symbol'])

    # find target genes for each cancer subtype
    for index, row in target_df.iterrows():
        if cancer_name.lower() in str(row['Indications']).lower():
            temp = str(row['Targets']).split('; ')
            for gene in temp:
                if gene in network.nodes():
                    target_.append(gene)

    target_ = list(set(target_))  # remove redundant elements from the list
    nontarget_ = []  # genes in signaling network which are not prostate cancer targets
    for u in network.nodes():
        if not u in target_:
            nontarget_.append(u)
    return onco_, target_, nontarget_


# construct a smallest connected network containing all the input cancer genes
def constructNetwork(cancer_genes, nodetype, targetset, g, cancer_name):  # construct cancer network
    for i in [f'{output_path}/{cancer_name}_{nodetype}/ppr',
              f'{output_path}/{cancer_name}_{nodetype}/dppr',
              f'{output_path}/{cancer_name}_{nodetype}/PEN_distance',
              f'{output_path}/{cancer_name}_{nodetype}/stats',
              f'{output_path}/{cancer_name}_{nodetype}/distance']:
        if not os.path.exists(i):
            os.makedirs(i)
    if not os.path.exists(f'{output_path}/{cancer_name}_{nodetype}/{cancer_name}_{nodetype}.gexf'):
        temp1 = cancer_genes.copy()  # cancer_genes is the set of cancer_genes of this kind of cancer type
        cancer_genes1 = []
        for u in temp1:
            if u in g.nodes():
                cancer_genes1.append(u)  # remove the cancer_genes which are not in the network
        temp2 = targetset.copy()  # target set is the set of drug targets for the input cancer type
        targetset1 = []
        for v in temp2:
            if v in g.nodes():
                targetset1.append(v)  # remove the targets which are not in the network
        print('There are {} cancer genes and {} targets'.format(len(cancer_genes1), len(targetset1)))
        new_nodes = set()
        count = 0
        # if there is at least one path from a target u to cancer_gene v
        # if the shortest path length <= 5, add the paths with the length <=5 to the network
        # else, add all the paths from u to v to the network.
        for u in targetset1:
            for v in cancer_genes1:
                count += 1
                print(count)
                if nx.has_path(g, u, v):
                    lp = nx.shortest_path_length(g, u, v)
                    # print('length of the shortest path between {} and {} is {}'.format(u,v,lp))
                    if lp <= 5:
                        temp = nx.all_simple_paths(g, u, v, 5)
                    else:
                        temp = nx.all_simple_paths(g, u, v, lp)
                    for p in temp:
                        new_nodes.update(p)
                # if there is no path from target u to cancer_gene v, add u and v to network only.
                else:
                    new_nodes.update({u, v})

                    # construct cancer network
        cancer_network = nx.subgraph(g, new_nodes)
        cancer_network_1 = cancer_network.copy()

        # compute the PEN_distance of the network
        PEN_distance.PEN_distance_alpha(cancer_network_1,
                                        f'{output_path}/{cancer_name}_{nodetype}/ppr',
                                        f'{output_path}/{cancer_name}_{nodetype}/stats',
                                        cancer_name, f'{output_path}/{cancer_name}_{nodetype}/dppr',
                                        f'{output_path}/{cancer_name}_{nodetype}/PEN_distance', iter=2)
        # compute the distance of the netowrk
        dist.compute_shortest_distance(cancer_network_1, f'{output_path}/{cancer_name}_{nodetype}/distance')

        # write the network to gexf file
        nx.write_gexf(cancer_network_1, f'{output_path}/{cancer_name}_{nodetype}/{cancer_name}_{nodetype}.gexf')
    else:
        cancer_network_1 = nx.read_gexf(f'{output_path}/{cancer_name}_{nodetype}/{cancer_name}_{nodetype}.gexf')

    # plot the histogram of ppr, dppr, and PEN-distance of the cancer signaling network
    pprdf = pd.read_csv(f'{output_path}//{cancer_name}_{nodetype}//ppr//alpha = 0.2//ppr.txt',
                        sep='\t', header=0, index_col=0)
    PEN_distancedf = pd.read_csv(
        f'{output_path}//{cancer_name}_{nodetype}//PEN_distance//alpha = 0.2//PEN_distance0.txt',
        sep='\t', header=0, index_col=0)
    plot_features(pprdf, cancer_name, nodetype, 'ppr')
    plot_features(PEN_distancedf, cancer_name, nodetype, 'PEN distance')

    return cancer_network_1
