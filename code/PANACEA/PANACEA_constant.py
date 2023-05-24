# this file list all the constant, dataset, paths we will use in the whole project
# so that we will  not have to write the whole string every time when we want to use it

import os
import pandas as pd
import networkx as nx

# list the path we will use in the project
project_path = os.path.dirname(os.path.dirname(os.getcwd()))
output_path = os.path.join(project_path, 'Output')  # output path
input_path = os.path.join(project_path, 'Input')  # input path

for ps in [output_path, f'{output_path}//plot']:
    if not os.path.exists(ps):
        os.makedirs(ps, exist_ok=True)

# define the constant and dataset we will use in the project
target_df = pd.read_csv(os.path.join(input_path, 'Drug targets.csv'), encoding="unicode_escape")
target_raw = list(target_df['Targets'])
target_ls = []
for tg in target_raw:
    if type(tg) == str:
        temp = tg.split('; ')
        target_ls.extend(temp)
target_ls = list(set(target_ls))
cancer_df = pd.read_csv(os.path.join(input_path, 'oncogene.csv'), encoding='unicode_escape')
cancer_ls = list(cancer_df['Gene Symbol'])  # list of cancer genes
whole_signaling = nx.read_gexf(f'{input_path}//whole_signaling.gexf')