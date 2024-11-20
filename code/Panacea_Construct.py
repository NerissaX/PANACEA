# Panacea_Construct.py
# This script contains functions to find genes associated with specific cancer subtypes
# and to construct the smallest connected network containing all input cancer genes.

import os
import pandas as pd
import networkx as nx
import Panacea_distance as dist
import Panacea_PEN_dist as PEN_distance
from Panacea_constant import *
from Panacea_plot import plot_features
import PanaceaStats as stats


def find_subgene(network, cancer_type, cancer_name):
    """
    Find oncogenes, target genes, and non-target genes for a specific cancer subtype.

    Parameters:
    - network: networkx.Graph, the signaling network
    - cancer_type: str, code for the tumor type (e.g., 'prostate')
    - cancer_name: str, full name of the cancer (e.g., 'Prostate Cancer')

    Returns:
    - oncogenes: list of oncogenes associated with the cancer subtype
    - targets: list of known drug targets associated with the cancer subtype
    - nontargets: list of genes in the network that are not targets
    """
    cancer_type = cancer_type.lower()
    cancer_name = cancer_name.lower()
    network_nodes = set(network.nodes())

    # Find oncogenes for the cancer subtype from cancer_df
    oncogenes = set()
    for _, row in cancer_df.iterrows():
        cancer_types_somatic = str(row['Tumour Types(Somatic)']).lower()
        cancer_types_germline = str(row['Tumour Types(Germline)']).lower()
        gene_symbol = row['Gene Symbol']

        if cancer_type in cancer_types_somatic or cancer_type in cancer_types_germline:
            if gene_symbol in network_nodes:
                oncogenes.add(gene_symbol)

    # Find target genes for the cancer subtype from target_df
    targets = set()
    for _, row in target_df.iterrows():
        indications = str(row['Indications']).lower()
        if cancer_name in indications:
            target_genes = str(row['Targets']).split('; ')
            for gene in target_genes:
                if gene in network_nodes:
                    targets.add(gene)

    # Find non-target genes in the network
    nontargets = network_nodes - targets

    return list(oncogenes), list(targets), list(nontargets)


def constructNetwork(cancer_genes, nodetype, target_set, network, cancer_name):
    """
    Construct the smallest connected network containing all input cancer genes and their targets.
    It also computes network statistics, PEN distance, and plots features.

    Parameters:
    - cancer_genes: list of cancer genes associated with the subtype
    - nodetype: str, type of nodes to use (e.g., 'oncogenes')
    - target_set: list of known drug targets for the cancer subtype
    - network: networkx.Graph, the whole signaling network
    - cancer_name: str, name of the cancer subtype

    Returns:
    - cancer_network: networkx.Graph, the constructed cancer network
    """
    # Define output directories
    base_path = f'{output_path}/{cancer_name}_{nodetype}'
    output_dirs = [
        f'{base_path}/ppr',
        f'{base_path}/dppr',
        f'{base_path}/PEN_distance',
        f'{base_path}/stats',
        f'{base_path}/distance'
    ]

    # Create directories if they do not exist
    for dir_path in output_dirs:
        os.makedirs(dir_path, exist_ok=True)

    # Path to the network file
    network_file = f'{base_path}/{cancer_name}_{nodetype}.gexf'

    if not os.path.exists(network_file):
        # Filter cancer genes and targets to those present in the network
        cancer_genes_in_network = [gene for gene in cancer_genes if gene in network.nodes()]
        targets_in_network = [gene for gene in target_set if gene in network.nodes()]

        print(f'There are {len(cancer_genes_in_network)} cancer genes and {len(targets_in_network)} targets')

        new_nodes = set()
        count = 0

        # For each pair of target and cancer gene, find paths and add nodes to new_nodes
        for target_gene in targets_in_network:
            for cancer_gene in cancer_genes_in_network:
                count += 1
                print(f'Processing pair {count}: {target_gene}, {cancer_gene}')
                if nx.has_path(network, target_gene, cancer_gene):
                    try:
                        # Find the shortest path length between target and cancer gene
                        shortest_length = nx.shortest_path_length(network, target_gene, cancer_gene)
                        if shortest_length <= 5:
                            # Collect all simple paths up to length 5
                            paths = nx.all_simple_paths(network, target_gene, cancer_gene, cutoff=5)
                            for path in paths:
                                new_nodes.update(path)
                        else:
                            # For longer paths, add nodes from the shortest path only
                            path = nx.shortest_path(network, target_gene, cancer_gene)
                            new_nodes.update(path)
                    except nx.NetworkXNoPath:
                        # If no path exists, add the genes themselves
                        new_nodes.update([target_gene, cancer_gene])
                else:
                    new_nodes.update([target_gene, cancer_gene])

        # Construct the cancer network
        cancer_network = network.subgraph(new_nodes).copy()

        # Compute network statistics
        stats.compute_network_stats(
            cancer_network,
            cancer_name,
            f'{base_path}/stats'
        )

        # Compute the PEN distance of the network
        PEN_distance.PEN_distance_alpha(
            G=cancer_network,
            network_name=cancer_name,
            network_ppr_path=f'{base_path}/ppr',
            stat_network_path=f'{base_path}/stats',
            network_dppr_path=f'{base_path}/dppr',
            network_PEN_distance_path=f'{base_path}/PEN_distance',
            iterations=2
        )

        # Compute the shortest distances in the network
        dist.compute_shortest_distance(
            cancer_network,
            f'{base_path}/distance'
        )

        # Write the network to a GEXF file
        nx.write_gexf(cancer_network, network_file)
    else:
        # Read the existing network from file
        cancer_network = nx.read_gexf(network_file)

    # Plot the histograms of PPR and PEN distance
    ppr_file = f'{base_path}/ppr/alpha = 0.2/ppr.txt'
    pen_distance_file = f'{base_path}/PEN_distance/alpha = 0.2/PEN_distance0.txt'

    ppr_df = pd.read_csv(ppr_file, sep='\t', header=0, index_col=0)
    pen_distance_df = pd.read_csv(pen_distance_file, sep='\t', header=0, index_col=0)
    plot_features(ppr_df, cancer_name, nodetype, 'ppr')
    plot_features(pen_distance_df, cancer_name, nodetype, 'PEN distance')

    return cancer_network
