from Panacea_constant import *
import pandas as pd
import numpy as np
from math import floor
from Panacea_plot import plot_deltaHist
from Panacea_plot import plot_diff_distribution
from Panacea_target import *


# Panacea_profiling_comb.py
# This script contains the function to perform delta histogram profiling of a cancer network.
def delta_histogram(cancer_name, nodetype, k, bucket_no=5):
    """
    Performs delta histogram profiling of a cancer network.

    Parameters:
    - cancer_name: str, name of the cancer subtype.
    - nodetxype: str, type of cancer nodes used (e.g., 'oncogenes').
    - k: int, number of nodes in a combination.
    - bucket_no: int, number of buckets to divide the data into (default is 5).

    Returns:
    - deltamin: dict, minimum delta values for each feature.
    - deltamax: dict, maximum delta values for each feature.
    """

    print('**********************************************************************************')
    print(f'{cancer_name}:')

    # Load known target combinations and candidate combinations data
    known_target_file = os.path.join(
        output_path,
        f'{cancer_name}_{nodetype}',
        f'{cancer_name}_{nodetype}_known_targets.txt'
    )
    candidate_target_file = os.path.join(
        output_path,
        f'{cancer_name}_{nodetype}',
        f'{cancer_name}_{nodetype}_{k}set_combo.txt'
    )

    known_target = pd.read_csv(known_target_file, header=0, index_col=[0, 1], sep='\t')
    candidate_target = pd.read_csv(candidate_target_file, header=0, index_col=[0, 1], sep='\t')

    # Plot the distribution of differences for each feature
    for feature_name in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
        feature_values = candidate_target[feature_name].tolist()
        plot_diff_distribution(feature_values, feature_name, cancer_name, nodetype)

    # Rank the data and save to file by each feature
    pen_diff_ranked_file = os.path.join(
        output_path,
        f'{cancer_name}_{nodetype}',
        f'{cancer_name}_{nodetype}_{k}_PEN-diff_rankedresults.txt'
    )

    if not os.path.exists(pen_diff_ranked_file):
        # Initialize 'known targets' column in candidate_target DataFrame
        candidate_target['known targets'] = 0

        # Mark known target combinations in candidate_target DataFrame
        known_indices = known_target.index.tolist()
        candidate_indices = candidate_target.index.tolist()

        count = 0
        not_found_count = 0

        for ind in known_indices:
            count += 1
            print(count)
            if ind in candidate_indices:
                candidate_target.at[ind, 'known targets'] = 1
            elif ind[::-1] in candidate_indices:
                candidate_target.at[ind[::-1], 'known targets'] = 1
            else:
                not_found_count += 1

        # Prepare dictionaries to store dataframes
        known_target_dict = {}
        candidate_target_dict = {}
        ranked_results_dict = {}

        for feature_name in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
            # Extract the feature column for known targets
            known_target_dict[feature_name] = known_target[[feature_name]]

            # Extract the feature and 'known targets' columns for candidates
            candidate_target_dict[feature_name] = candidate_target[[feature_name, 'known targets']]

            # Sort candidate combinations by the feature in descending order
            ranked_results = candidate_target_dict[feature_name].sort_values(by=[feature_name], ascending=False)
            ranked_results_dict[feature_name] = ranked_results

            # Save the ranked results to file
            ranked_file = os.path.join(
                output_path,
                f'{cancer_name}_{nodetype}',
                f'{cancer_name}_{nodetype}_{k}_{feature_name}_rankedresults.txt'
            )
            ranked_results.to_csv(ranked_file, header=True, index=True, sep='\t')

    else:
        # Load existing ranked results
        known_target_dict = {}
        ranked_results_dict = {}

        for feature_name in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
            # Extract the feature column for known targets
            known_target_dict[feature_name] = known_target[[feature_name]]

            # Load ranked results from file
            ranked_file = os.path.join(
                output_path,
                f'{cancer_name}_{nodetype}',
                f'{cancer_name}_{nodetype}_{k}_{feature_name}_rankedresults.txt'
            )
            ranked_results = pd.read_csv(ranked_file, header=0, index_col=[0, 1], sep='\t')
            ranked_results_dict[feature_name] = ranked_results

    # Initialize dictionaries to store delta min and max values
    deltamin = {}  # Minimum delta values for each feature
    deltamax = {}  # Maximum delta values for each feature

    # For each feature, find delta histogram and calculate constraints
    for feature_name in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
        ranked_results = ranked_results_dict[feature_name]
        known_target_df = known_target_dict[feature_name]

        # Determine the range of feature values and calculate bucket boundaries
        max_value = ranked_results[feature_name].max()
        min_value = ranked_results[feature_name].min()
        step = (max_value - min_value) / bucket_no

        # Create bucket ranges
        bucket_ranges = [min_value + i * step for i in range(bucket_no + 1)]
        bucket_ranges_rev = bucket_ranges[::-1]

        # Assign group names based on bucket ranges
        group_names = []
        for i in range(bucket_no):
            lower_bound = round(bucket_ranges_rev[i + 1], 4)
            upper_bound = round(bucket_ranges_rev[i], 4)
            group_name = f'{lower_bound} - {upper_bound}'
            group_names.append(group_name)

            # Assign group labels to ranked_results and known_target_df
            ranked_results.loc[ranked_results[feature_name] <= bucket_ranges_rev[i], 'group'] = group_name
            known_target_df.loc[known_target_df[feature_name] <= bucket_ranges_rev[i], 'group'] = group_name

        group_names.reverse()  # Reverse to match original order

        # Collect data for all candidates and known targets in each bucket
        all_groups = []
        known_groups = []

        for group_name in group_names:
            group_candidates = ranked_results[ranked_results['group'] == group_name][feature_name].sort_values(
                ascending=False)
            all_groups.append(group_candidates)

            group_known = known_target_df[known_target_df['group'] == group_name][feature_name].sort_values(
                ascending=False)
            known_groups.append(group_known)

        # Calculate coverage of known combinations in top m% for each bucket
        delta_percent_groups = []
        all_percentage = []
        candidate_counts = []
        known_counts = []

        for i in range(bucket_no):
            total_known = len(known_target_df)
            known_in_bucket = len(known_groups[i])
            candidates_in_bucket = len(all_groups[i])

            # Percentage of known combinations in this bucket
            percentage_known = known_in_bucket / total_known if total_known > 0 else 0
            all_percentage.append(percentage_known)
            candidate_counts.append(candidates_in_bucket)
            known_counts.append(known_in_bucket)

            temp_dict = {}
            # For top m% in the bucket, compute the coverage of known combinations
            for m in [0.01, 0.10, 0.20, 0.50]:
                threshold_index = floor(candidates_in_bucket * m)
                if candidates_in_bucket > threshold_index:
                    threshold_value = all_groups[i].iloc[threshold_index]
                    count = (known_groups[i] >= threshold_value).sum()
                    coverage = count / known_in_bucket if known_in_bucket > 0 else 0
                else:
                    coverage = 0
                temp_dict[f'{int(m * 100)}%'] = coverage
            delta_percent_groups.append(temp_dict)

        # Create a DataFrame to store coverage information
        coverage_df = pd.DataFrame(delta_percent_groups, index=range(bucket_no))
        coverage_df['range'] = group_names.copy()

        # Find the bucket with the highest coverage at 50%
        coverage_at_50 = coverage_df['50%'].tolist()
        max_coverage_index = coverage_at_50.index(np.nanmax(coverage_at_50))

        delta_range_str = coverage_df.loc[max_coverage_index, 'range']
        constraint_known = known_groups[max_coverage_index]

        # Call target function with the known combinations in the constraint range
        target(constraint_known, cancer_name, nodetype, feature_name)

        # Extract delta_min and delta_max from the range string
        delta_min_str, delta_max_str = delta_range_str.split(' - ')
        delta_min = float(delta_min_str)
        delta_max = float(delta_max_str)

        # Print the results
        print('----------------------------------------------------------------------------------------')
        print(f'The delta_min of {feature_name} is {delta_min}')
        print(f'The delta_max of {feature_name} is {delta_max}')
        print(f'The coverage of {feature_name} is {list(delta_percent_groups[max_coverage_index].values())}')
        print(f'There are {candidate_counts[max_coverage_index]} candidate combinations in the constraint range.')

        # Plot the delta histogram
        plot_deltaHist(coverage_df, cancer_name, nodetype, feature_name)

        # Store delta_min and delta_max values
        deltamin[feature_name] = delta_min
        deltamax[feature_name] = delta_max

    return deltamin, deltamax
