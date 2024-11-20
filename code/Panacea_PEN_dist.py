# Panacea_PEN_dist.py
# This script contains functions to compute Personalized PageRank (PPR), differential PPR (dPPR),
# and the PEN distance for a given network.

import os
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.preprocessing import normalize
from Panacea_constant import *


def compute_PEN_distance(G, filename, ppr_path, stat_path, dppr_path, PEN_distance_path, num_iter=100, alpha=0.2):
    """
    Computes the Personalized PageRank (PPR), differential PPR (dPPR), and PEN distance for a graph G.

    Parameters:
    - G: networkx.Graph, the input graph.
    - filename: str, base name for reading and writing files.
    - ppr_path: str, directory path to save PPR results.
    - stat_path: str, directory path where network statistics are stored.
    - dppr_path: str, directory path to save dPPR results.
    - PEN_distance_path: str, directory path to save PEN distance results.
    - num_iter: int, number of iterations for the power method (default 100).
    - alpha: float, teleportation parameter for PPR (default 0.2).
    """
    ppr_file = os.path.join(ppr_path, 'ppr.txt')
    if not os.path.exists(ppr_file):
        # Read out-degree of the network
        out_degree_file = os.path.join(stat_path, f'{filename}_sorted out_degree.txt')
        out_degree_df = pd.read_csv(out_degree_file, sep='\t', header=None)
        nodes = out_degree_df.iloc[:, 0].astype(str).tolist()  # Node names
        degrees = out_degree_df.iloc[:, 1].values  # Out-degree values
        out_degree_dict = dict(zip(nodes, degrees))

        # Prepare DataFrames to store results
        node_list = sorted(G.nodes())
        ppr_df = pd.DataFrame(index=node_list, columns=node_list)
        dppr_df = pd.DataFrame(index=node_list, columns=node_list)
        pen_distance_df = pd.DataFrame(index=node_list, columns=node_list)

        # Compute adjacency matrix and transition probability matrix
        A = nx.adjacency_matrix(G, nodelist=node_list)
        P = normalize(A, norm='l1', axis=1).T  # Column-stochastic matrix

        # Compute PPR, dPPR, and PEN distance for each node
        for idx, src_node in enumerate(node_list):
            # Initialize personalization vector
            e = np.zeros(len(node_list))
            e[idx] = 1.0
            ppr = e.copy()

            # Power iteration method for PPR computation
            for _ in range(num_iter):
                ppr = (1 - alpha) * P.dot(ppr) + e
            ppr_values = alpha * ppr

            # Store PPR, dPPR, and PEN distance values
            for count, tgt_node in enumerate(node_list):
                ppr_df.loc[src_node, tgt_node] = ppr_values[count]
                dppr_value = ppr_values[count] * out_degree_dict[src_node]
                dppr_df.loc[src_node, tgt_node] = dppr_value

                # Compute PEN distance
                pen_distance = round(1 - np.log(dppr_value + 1e-5), 4)
                pen_distance_df.loc[src_node, tgt_node] = pen_distance

            print(f'Processed node {idx + 1}/{len(node_list)}: {src_node}')

        # Save results to files
        ppr_df.to_csv(ppr_file, sep='\t')
        dppr_df.to_csv(os.path.join(dppr_path, 'dppr.txt'), sep='\t')
        pen_distance_df.to_csv(os.path.join(PEN_distance_path, 'PEN_distance.txt'), sep='\t')

        # Replace negative PEN distances with zero and save
        pen_distance_df_zeroed = pen_distance_df.clip(lower=0)
        pen_distance_df_zeroed.to_csv(os.path.join(PEN_distance_path, 'PEN_distance0.txt'), sep='\t')
        print('PPR, dPPR, and PEN distance computations completed and saved.')


def PEN_distance_alpha(G, network_name, network_ppr_path, stat_network_path, network_dppr_path,
                       network_PEN_distance_path, iterations=2):
    """
    Computes PPR, dPPR, and PEN distance for different alpha values.

    Parameters:
    - G: networkx.Graph, the input graph.
    - network_name: str, base name for reading and writing files.
    - network_ppr_path: str, base directory path for PPR results.
    - stat_network_path: str, directory path where network statistics are stored.
    - network_dppr_path: str, base directory path for dPPR results.
    - network_PEN_distance_path: str, base directory path for PEN distance results.
    - iterations: int, number of different alpha values to compute (default 2).

    Note:
    - Alpha values are computed as alpha = (i + 1) / 10 for i in range(iterations).
    """
    for i in range(iterations):
        alpha = (i + 1) / 10  # Alpha values: 0.1, 0.2, ..., depending on iterations
        alpha_str = f'alpha_{alpha}'

        # Define paths for this alpha value
        dppr_path = os.path.join(network_dppr_path, alpha_str)
        ppr_path = os.path.join(network_ppr_path, alpha_str)
        pen_distance_path = os.path.join(network_PEN_distance_path, alpha_str)

        # Create directories if they don't exist
        for path in [dppr_path, ppr_path, pen_distance_path]:
            os.makedirs(path, exist_ok=True)

        # Compute and save PPR, dPPR, and PEN distance
        compute_PEN_distance(
            G, network_name, ppr_path, stat_network_path,
            dppr_path, pen_distance_path, num_iter=100, alpha=alpha
        )
        print(f'Computed PEN distance for alpha = {alpha}')

    print('All computations for specified alpha values are completed.')