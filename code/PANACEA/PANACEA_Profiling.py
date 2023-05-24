from PANACEA_constant import *
import pandas as pd
import numpy as np
from math import floor
from PANACEA_plot import plot_deltaHist
from PANACEA_plot import plot_diff_distribution


def delta_histogram(cancer_name, nodetype, k):
    print('**********************************************************************************')
    print(f'{cancer_name}:')
    # known_target is the dataframe of ppr-diff, pen-diff, distance-diff, average closeness centrality
    # for all known target combinations
    # candidate_target is for all candidate target combinations
    known_target = pd.read_csv(f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_known_targets.csv',
                               header=0, index_col=[0, 1], sep=',')
    candidate_target = pd.read_csv(f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}set_combo.txt',
                                   sep='\t', header=0, index_col=[0, 1])

    # for PEN-diff and each baseline, find their delta histogram
    for fn in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
        ls = list(candidate_target[fn])  # extract the column of the corresponding feature of candidates
        plot_diff_distribution(ls,fn,cancer_name,nodetype)

    # rank the data and save to file by each feature
    if not os.path.exists(
            f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}_PEN-diff_rankedresults.txt'):
        count = 0
        not_count = 0
        # a new column to show whether the combination is a known target combination
        # if yes, the value is 1, else 0
        if_known = [0] * len(candidate_target)
        candidate_target['known targets'] = if_known
        for ind in list(known_target.index):
            count += 1
            print(count)
            if ind in list(candidate_target.index):
                tar_ind = list(candidate_target.index).index(ind)
                candidate_target['known targets'][tar_ind] = 1
            elif ind[::-1] in list(candidate_target.index):
                tar_ind = list(candidate_target.index).index(ind[::-1])
                candidate_target['known targets'][tar_ind] = 1
            else:
                not_count += 1

        known_target_dict = {}
        candidate_target_dict = {}
        ranked_results_dict = {}
        for fn in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
            known_target_dict[fn] = known_target.loc[:, [fn]]
            candidate_target_dict[fn] = candidate_target.loc[:, [fn, 'known targets']]
            ranked_results_dict[fn].to_csv(
                f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}_{fn}_rankedresults.txt',
                header=True,
                index=True, sep='\t')

    else:  # if the ranked results already exist, skip the ranking step
        known_target_dict = {}
        ranked_results_dict = {}

        for fn in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
            known_target_dict[fn] = known_target.loc[:, [fn]]
            ranked_results_dict[fn] = pd.read_csv(
                f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}_{fn}_rankedresults.txt',
                header=0,
                index_col=[0, 1], sep='\t')

    deltamin = {}  # deltamin in descending order
    deltamax = {}

    bucket_no = 5  # fix the number of buckets
    # for each features in PEN-diff, ppr-diff, distance-diff, dppr-diff and avg closeness, find delta histogram
    for fn in ['PEN-diff', 'Distance-diff', 'ppr-diff']:
        ranked_results = ranked_results_dict[fn]
        known_target_df = known_target_dict[fn]
        max_p = ranked_results[fn].max()  # max value of the data
        min_p = ranked_results[fn].min()  # min value of the data
        step_p = (max_p - min_p) / bucket_no  # range of a bucket
        # group the data by buckets
        range_p = []
        for i in range(bucket_no + 1):
            range_p.append(min_p + i * step_p)
        reverse_range_p = range_p.copy()
        reverse_range_p.reverse()
        group_name = []
        for i in range(bucket_no):
            gname = f'{round(reverse_range_p[i + 1], 4)} - {round(reverse_range_p[i], 4)}'
            group_name.append(gname)
            # add a column to record the group name
            ranked_results.loc[ranked_results[fn] <= reverse_range_p[i], 'group'] = gname
            known_target_df.loc[known_target_df[fn] <= reverse_range_p[i], 'group'] = gname
        group_name.reverse()

        allgroups = []
        knowngroups = []
        # find delta histogram, deltamin, deltamax in descending order
        # all the steps are the same as above
        for i in range(bucket_no):
            allgroups.append(
                ranked_results_dict[fn].loc[ranked_results_dict[fn]['group'] == group_name[i], [fn]][
                    fn].sort_values(
                    ascending=False))
            knowngroups.append(
                known_target_dict[fn].loc[known_target_dict[fn]['group'] == group_name[i], [fn]][fn].sort_values(
                    ascending=False))
        delta_percent_groups = []
        all_percentage = []
        max_count = []
        known_count = []
        for i in range(bucket_no):
            all_percent = len(knowngroups[i]) / len(known_target_df)
            all_percentage.append(all_percent)
            max_count.append(len(allgroups[i]))
            known_count.append(len(knowngroups[i]))
            temp_dict = {}
            for m in [0.01, 0.10, 0.20, 0.50]:
                count = (knowngroups[i] >= allgroups[i][floor(len(allgroups[i]) * m)]).sum()
                temp = count / len(knowngroups[i])
                # temp = count / len(known_target_df)
                temp_dict[m] = temp
            delta_percent_groups.append(temp_dict)

        df = pd.DataFrame(columns=['1%', '10%', '20%', '50%'], index=[i for i in range(bucket_no)])
        for i in range(bucket_no):
            for j in range(4):
                df.iloc[i, j] = list(delta_percent_groups[i].values())[j]
        df['range'] = group_name.copy()
        s = [df.loc[i, '50%'] for i in range(bucket_no)]

        delta_range = df.loc[s.index(np.nanmax(s)), 'range']
        constraint_known = knowngroups[s.index(np.nanmax(s))]
        constraint_known.to_csv(
            f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{fn}_constraint_range.txt',
            header=True,
            index=True, sep='\t')
        delta_range = delta_range.split(' - ')
        delta_min = float(delta_range[0])
        delta_max = float(delta_range[1])

        # print the results on the
        print('----------------------------------------------------------------------------------------')
        print(f'The delta_min of {fn} is {delta_min}')
        print(f'The delta_max of {fn} is {delta_max}')
        print(f'The coverage of {fn} is {list((delta_percent_groups[s.index(np.nanmax(s))]).values())}')
        print(
            f'There are {max_count[s.index(np.nanmax(s))]} candidate combinations in the constraint range.')
        plot_deltaHist(df,cancer_name,nodetype,fn)
        deltamin[fn] = delta_min
        deltamax[fn] = delta_max

    return deltamin, deltamax
