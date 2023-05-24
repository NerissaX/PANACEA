from sklearn import preprocessing
import numpy as np
from PANACEA_constant import *


# the algorithm to compute ppr, dppr and PEN_distance
def compute_PEN_distance(G, ppr_path, stat_path, filename, dppr_path, PEN_distance_path, numIter=100, alpha=0.2):
    out_degree_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted out_degree.txt'),
        sep='\t', header=None)  # read out-degree of the network
    node = list(out_degree_df.iloc[:, 0])  # nodes in the network
    value = list(out_degree_df.iloc[:, 1])  # out-degree values of the nodes
    out_degree_dict = {}
    for i in range(len(node)):
        out_degree_dict[node[i]] = value[i]  # create an out-degree dictionary
    # create empty dataframes to store ppr, dppr and PEN_distance
    ppr_df = pd.DataFrame(index=pd.Index(sorted(G.nodes())), columns=pd.Index(sorted(G.nodes())))
    dppr_df = pd.DataFrame(index=pd.Index(sorted(G.nodes())), columns=pd.Index(sorted(G.nodes())))
    PEN_distance_df = pd.DataFrame(index=pd.Index(sorted(G.nodes())), columns=pd.Index(sorted(G.nodes())))

    # for each node in the network, compute ppr, dppr and PEN_distance
    for src in range(G.number_of_nodes()):  # src is the position of the current node in network
        nodename = sorted(G.nodes())[src]
        # the algorithm for computing ppr
        ppr_dict = {}
        A = nx.adjacency_matrix(G, nodelist=sorted(G.nodes()))  # adjacency matrix of the network
        P = preprocessing.normalize(A, norm='l1',
                                    axis=1).T  # Scale input vectors individually to unit norm (vector length)
        e = np.zeros(G.number_of_nodes())
        e[src] = 1.0
        ppr = e.copy()  # initialize the ppr list
        for i in range(numIter):  # compute ppr
            ppr = (1 - alpha) * P.dot(ppr) + e
        ppr_ls = alpha * ppr
        count = 0
        for j in sorted(G.nodes()):
            ppr_df.loc[nodename, j] = ppr_ls[count]  # store ppr values
            dppr_df.loc[nodename, j] = ppr_ls[count] * out_degree_dict[
                str(nodename)]  # store dppr values which is ppr * outdegree
            temp = round((1 - np.log(dppr_df.loc[nodename, j] + (1e-05))), 4)  # compute PEN_distance
            # if temp < 0:
            #     temp = 0
            PEN_distance_df.loc[nodename, j] = temp
            count += 1
        print(src)
        src += 1
    # write ppr, dppr, PEN_distance to txt files
    ppr_df.to_csv(f'{ppr_path}/ppr.txt', header=True, index=True, sep='\t')
    dppr_df.to_csv(f'{dppr_path}/dppr.txt', header=True, index=True, sep='\t')
    PEN_distance_df.to_csv(f'{PEN_distance_path}/PEN_distance.txt', header=True, index=True, sep='\t')
    PEN_distance_df0 = PEN_distance_df.copy()
    PEN_distance_df0[PEN_distance_df0 < 0] = 0
    PEN_distance_df0.to_csv(f'{PEN_distance_path}/PEN_distance0.txt', header=True, index=True, sep='\t')


# compute ppr, dppr and PEN_distance for a specific network
# iter is the number of different alphas
# for example, when iter = 3, the function will compute PEN_distance for alpha = 0.1, 0.2 and 0.3
# when iter = 9, the function will compute PEN_distance for alpha = 0.1 to 0.9
def PEN_distance_alpha(G, network_ppr_path, stat_network_path, network_name, network_dppr_path,
                       network_PEN_distance_path,
                       iter=2):
    alpha = 0
    for i in range(iter):
        alpha = (i + 1) / 10  # when iter = 1, alpha = 0.1
        dppr_path = os.path.join(network_dppr_path, f'alpha = {alpha}')
        PEN_distance_path = os.path.join(network_PEN_distance_path, f'alpha = {alpha}')
        ppr_path = os.path.join(network_ppr_path, f'alpha = {alpha}')
        for path in [dppr_path, ppr_path, PEN_distance_path]:
            if not os.path.exists(path):
                os.makedirs(path)
        compute_PEN_distance(G, network_name, stat_network_path, ppr_path,
                             dppr_path, PEN_distance_path, numIter=100, alpha=alpha)
