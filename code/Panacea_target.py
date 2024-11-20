# Panacea_target.py
# This script contains the function to write the target combinations within the constraint range to a file.

import os
import pandas as pd
from Panacea_constant import output_path


def target(target_combination, cancer_name, nodetype, fn):
    """
    Writes the target combinations within the constraint range to a file.

    Parameters:
    - target_combination: pandas.DataFrame, the target combinations to write.
    - cancer_name: str, name of the cancer subtype.
    - nodetype: str, type of cancer nodes used (e.g., 'oncogenes').
    - fn: str, feature name used in the constraint (e.g., 'PEN-diff').
    """
    # Construct the directory path
    dir_path = os.path.join(output_path, f'{cancer_name}_{nodetype}')
    os.makedirs(dir_path, exist_ok=True)

    # Construct the file path
    file_name = f'{cancer_name}_{nodetype}_{fn}_constraint_range.txt'
    file_path = os.path.join(dir_path, file_name)

    # Write the DataFrame to a tab-separated text file
    target_combination.to_csv(file_path, header=True, index=True, sep='\t')
    print(f'Target combinations saved to {file_path}')
