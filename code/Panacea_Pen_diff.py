# Panacea_Pen_diff.py
# This script computes the differences in PEN-distance, Distance, and PPR between cancer genes and other nodes
# for known drug target combinations and candidate combinations.

import os
import itertools
import pandas as pd
from Panacea_constant import *


def remove_nodes(input_nodes, cancer_network):
    """
    Remove nodes from input_nodes that are not present in the cancer_network.

    Parameters:
    - input_nodes: list or set of node names
    - cancer_network: networkx.Graph, the cancer-specific signaling network

    Returns:
    - List of nodes present in both input_nodes and cancer_network
    """
    network_nodes = set(cancer_network.nodes())
    return [node for node in input_nodes if node in network_nodes]


def PEN_diff(candidates, k, cancer_genes, nodetype, target_nodes, cancer_network, cancer_name):
    """
    Compute differences in PEN-distance, Distance, and PPR between cancer genes and other nodes
    for known drug target combinations and candidate combinations.

    Parameters:
    - candidates: list or set of candidate nodes
    - k: int, number of nodes in a combination
    - cancer_genes: list or set of cancer genes (e.g., oncogenes, TSG, biomarkers)
    - nodetype: str, the type of cancer_genes (e.g., 'oncogenes')
    - target_nodes: list or set of target nodes related to a specific cancer
    - cancer_network: networkx.Graph, the cancer-specific signaling network
    - cancer_name: str, name of the cancer subtype

    Outputs:
    - Writes known target combinations data and candidate combinations data to files
    """

    # Clean the drug targets data
    drug_targets = target_df.copy()

    # Process 'Targets' column
    targets_list = []
    for targets in drug_targets['Targets']:
        if isinstance(targets, str):
            targets_list.append(targets.split('; '))
        else:
            targets_list.append(targets)

    # Process 'Indications' column
    indications_list = []
    for indication in drug_targets['Indications']:
        if isinstance(indication, str):
            indications_list.append(indication.split('; '))
        else:
            indications_list.append(indication)

    # Remove unwanted columns
    columns_to_drop = [
        'Targets', 'Indications', 'EMA', 'FDA', 'EN', 'Other', 'WHO', 'Year',
        'Generic', 'DrugBank ID', 'ATC', 'ChEMBL', 'Last Update'
    ]
    drug_targets.drop(columns=columns_to_drop, inplace=True)

    # Assign new 'Indications' and 'Targets' columns
    drug_targets['Indications'] = indications_list
    drug_targets['Targets'] = targets_list
    drug_targets.dropna(subset=['Targets'], inplace=True)  # Remove rows with no targets
    drug_targets.set_index('Product', inplace=True)

    drug_df = drug_targets.copy()

    # Select only drugs with the cancer type we want
    cancer_name_lower = cancer_name.lower()
    filtered_drug_indices = []
    for idx in drug_df.index:
        indications = drug_df.at[idx, 'Indications']
        if isinstance(indications, list):
            include_drug = False
            for indication in indications:
                indication_lower = indication.lower()
                if 'added' in cancer_name_lower or 'removed' in cancer_name_lower:
                    temp_str = ' '.join(cancer_name.split(' ')[:2]).lower()
                    if temp_str in indication_lower:
                        include_drug = True
                        break
                else:
                    if cancer_name_lower in indication_lower:
                        include_drug = True
                        break
            if include_drug:
                filtered_drug_indices.append(idx)
        else:
            continue  # Skip if 'Indications' is not a list

    # Filter drug_df to include only relevant drugs
    drug_df = drug_df.loc[filtered_drug_indices]

    # Save the filtered drug targets data
    output_dir = f'{output_path}/{cancer_name}_{nodetype}'
    os.makedirs(output_dir, exist_ok=True)
    drug_targets_file = os.path.join(output_dir, f'{cancer_name}_{nodetype}_drug_targets.csv')
    drug_df.to_csv(drug_targets_file)
    print(f'Filtered drug targets saved to {drug_targets_file}')

    # Remove nodes not in the cancer network
    cancer_genes = remove_nodes(cancer_genes, cancer_network)
    target_nodes = remove_nodes(target_nodes, cancer_network)

    # Import and clean the datasets
    base_path = f'{output_path}/{cancer_name}_{nodetype}'
    pen_distance_path = os.path.join(base_path, 'PEN_distance', 'alpha = 0.2')
    ppr_path = os.path.join(base_path, 'ppr', 'alpha = 0.2')
    distance_path = os.path.join(base_path, 'distance')

    PEN_distance_df = pd.read_csv(
        os.path.join(pen_distance_path, 'PEN_distance0.txt'), sep='\t', index_col=0
    )
    ppr_df = pd.read_csv(
        os.path.join(ppr_path, 'ppr.txt'), sep='\t', index_col=0
    )
    dist_df = pd.read_csv(
        os.path.join(distance_path, 'distance.txt'), sep='\t', index_col=0
    )

    # Identify other nodes (non-cancer genes)
    cancer_genes_set = set(cancer_genes)
    other_nodes = [node for node in cancer_network.nodes() if node not in cancer_genes_set]

    # Group the dataframes by cancer_genes and other_nodes
    PEN_cancer_df = PEN_distance_df.loc[:, cancer_genes]
    PEN_other_df = PEN_distance_df.loc[:, other_nodes]
    ppr_cancer_df = ppr_df.loc[:, cancer_genes]
    ppr_other_df = ppr_df.loc[:, other_nodes]
    dist_cancer_df = dist_df.loc[:, cancer_genes]
    dist_other_df = dist_df.loc[:, other_nodes]

    # Compute known target combinations
    known_data = []
    for drug in drug_df.index:
        targets = drug_df.at[drug, 'Targets']
        # Ensure targets are in the cancer network
        targets_in_network = remove_nodes(targets, cancer_network)
        # Generate combinations of size k
        for subset in itertools.combinations(targets_in_network, k):
            subset_list = list(subset)
            # Calculate mean values
            mean_PEN_cancer = PEN_cancer_df.loc[subset_list].mean(axis=1).mean()
            mean_PEN_other = PEN_other_df.loc[subset_list].mean(axis=1).mean()
            mean_dist_cancer = dist_cancer_df.loc[subset_list].mean(axis=1).mean()
            mean_dist_other = dist_other_df.loc[subset_list].mean(axis=1).mean()
            mean_ppr_cancer = ppr_cancer_df.loc[subset_list].mean(axis=1).mean()
            mean_ppr_other = ppr_other_df.loc[subset_list].mean(axis=1).mean()
            # Compute differences
            pen_diff = mean_PEN_other - mean_PEN_cancer
            dist_diff = mean_dist_other - mean_dist_cancer
            ppr_diff = mean_ppr_other - mean_ppr_cancer
            # Store results
            known_data.append({
                'Combination': subset,
                'PEN Distance(cancer genes)': mean_PEN_cancer,
                'PEN Distance(other genes)': mean_PEN_other,
                'PEN-diff': pen_diff,
                'Distance(cancer genes)': mean_dist_cancer,
                'Distance(other genes)': mean_dist_other,
                'Distance-diff': dist_diff,
                'ppr(cancer genes)': mean_ppr_cancer,
                'ppr(other genes)': mean_ppr_other,
                'ppr-diff': ppr_diff
            })

    # Create dataframe from known target combinations
    known_targets_df = pd.DataFrame(known_data)
    known_targets_df.set_index('Combination', inplace=True)

    # Save known target combinations to file
    known_targets_file = os.path.join(output_dir, f'{cancer_name}_{nodetype}_known_targets.txt')
    known_targets_df.to_csv(known_targets_file, sep='\t')
    print(f'Known targets data saved to {known_targets_file}')

    # Compute candidate combinations if not already done
    candidates_file = os.path.join(output_dir, f'{cancer_name}_{nodetype}_{k}set_combo.txt')
    if not os.path.exists(candidates_file):
        candidate_data = []
        total_combinations = itertools.combinations(candidates, k)
        total_count = 0
        for subset in total_combinations:
            subset_list = list(subset)
            # Calculate mean values
            mean_PEN_cancer = PEN_cancer_df.loc[subset_list].mean(axis=1).mean()
            mean_PEN_other = PEN_other_df.loc[subset_list].mean(axis=1).mean()
            mean_dist_cancer = dist_cancer_df.loc[subset_list].mean(axis=1).mean()
            mean_dist_other = dist_other_df.loc[subset_list].mean(axis=1).mean()
            mean_ppr_cancer = ppr_cancer_df.loc[subset_list].mean(axis=1).mean()
            mean_ppr_other = ppr_other_df.loc[subset_list].mean(axis=1).mean()
            # Compute differences
            pen_diff = mean_PEN_other - mean_PEN_cancer
            dist_diff = mean_dist_other - mean_dist_cancer
            ppr_diff = mean_ppr_other - mean_ppr_cancer
            # Store results
            candidate_data.append({
                'Combination': subset,
                'PEN Distance(cancer genes)': mean_PEN_cancer,
                'PEN Distance(other genes)': mean_PEN_other,
                'PEN-diff': pen_diff,
                'Distance(cancer genes)': mean_dist_cancer,
                'Distance(other genes)': mean_dist_other,
                'Distance-diff': dist_diff,
                'ppr(cancer genes)': mean_ppr_cancer,
                'ppr(other genes)': mean_ppr_other,
                'ppr-diff': ppr_diff
            })
            total_count += 1
            if total_count % 1000 == 0:
                print(f'Processed {total_count} candidate combinations')

        # Create dataframe from candidate combinations
        candidates_df = pd.DataFrame(candidate_data)
        candidates_df.set_index('Combination', inplace=True)

        # Save to file
        candidates_df.to_csv(candidates_file, sep='\t')
        print(f'Candidate combinations data saved to {candidates_file}')
    else:
        print(f'Candidate combinations data already exists at {candidates_file}')