# This file lists all constants, datasets, and paths used throughout the project.
# It centralizes configurations so that we don't have to hardcode strings multiple times.

import os
import pandas as pd
import networkx as nx

# list the path we will use in the project
project_path = os.path.dirname(os.getcwd())
output_path = os.path.join(project_path, 'Output')  # output path
input_path = os.path.join(project_path, 'Input')  # input path

# Load datasets and define constants
try:
    # Load drug targets dataset
    target_df = pd.read_csv(os.path.join(input_path, 'Drug targets.csv'), encoding="unicode_escape")

    # Extract and process the 'Targets' column
    target_raw = target_df['Targets'].dropna().tolist()
    target_ls = []
    for targets in target_raw:
        target_ls.extend(targets.split('; '))
    target_ls = list(set(target_ls))  # Remove duplicates

    # Load oncogenes dataset
    cancer_df = pd.read_csv(os.path.join(input_path, 'oncogene.csv'), encoding='unicode_escape')
    cancer_ls = cancer_df['Gene Symbol'].dropna().tolist()  # List of cancer genes

    # Load the whole signaling network
    whole_signaling = nx.read_gexf(f'{input_path}//whole_signaling.gexf')

except Exception as e:
    print(f"An error occurred while loading datasets: {e}")
