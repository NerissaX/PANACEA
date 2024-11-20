# Panacea_distance.py
# This script computes the shortest distances between all pairs of nodes in a given network and saves the results to a file.

import os
import networkx as nx
import pandas as pd
import numpy as np

def compute_shortest_distance(G, save_to):
    """
    Computes the shortest path lengths between all pairs of nodes in the graph G and saves the results to a file.

    Parameters:
    - G: networkx.Graph, the input network.
    - save_to: str, the directory path where the distance file will be saved.

    The function checks if the distance file already exists. If not, it computes the shortest path lengths,
    stores them in a DataFrame, and saves the DataFrame to a tab-separated text file named 'distance.txt'.
    """
    distance_file = os.path.join(save_to, 'distance.txt')
    if not os.path.exists(distance_file):
        nodes = sorted(G.nodes())
        # Initialize an empty DataFrame with nodes as both index and columns
        distance_df = pd.DataFrame(index=nodes, columns=nodes, dtype=float)

        # Compute shortest paths between all pairs of nodes
        for source in nodes:
            # Compute shortest paths from the source node to all other nodes
            lengths = nx.single_source_shortest_path_length(G, source)
            for target, length in lengths.items():
                distance_df.at[source, target] = length
            # For nodes not reachable from the source, the value remains NaN

        # Save the DataFrame to a file
        distance_df.to_csv(distance_file, sep='\t', index=True, header=True)
        print(f'Shortest distances computed and saved to {distance_file}')
    else:
        print(f'Distance file already exists at {distance_file}')
