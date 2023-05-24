import networkx as nx
import pandas as pd
import numpy as np


# compute the shortest distance between any pairs of nodes in a network
def compute_shortest_distance(G, save_to):
    distance_df = pd.DataFrame(index=pd.Index(sorted(G.nodes())), columns=pd.Index(sorted(G.nodes())))
    # create an empty dataframe
    for u in sorted(G.nodes()):
        for v in sorted(G.nodes()):
            try:
                distance_df.loc[u, v] = nx.shortest_path_length(G, u, v)
                # store the shortest path length between each pair of nodes
            except nx.NetworkXNoPath:
                distance_df.loc[u, v] = np.NaN  # if no connection, store the value as NaN

    distance_df.to_csv(f'{save_to}\\distance.txt', index=True, header=True, sep='\t')  # write to file
