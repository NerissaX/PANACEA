import itertools
from PANACEA_constant import *


# remove the nodes in inputnodes that are not in the given network
def remove_nodes(inputnodes, cancer_network):
    temp_c = inputnodes.copy()
    for u in temp_c:
        if u not in cancer_network.nodes():
            inputnodes.remove(u)
    return inputnodes


def PEN_diff(candidates, k, cancer_genes, nodetype, target_nodes, cancer_network,
             cancer_name):  # sort the candidate combination

    # clean the drug targets data
    drug_targets = target_df.copy()
    target_raw = list(target_df['Targets'])
    targets_in_list = []
    indications = []
    for tg in target_raw:
        if type(tg) == str:
            temp = tg.split('; ')
            targets_in_list.append(temp)
        else:
            targets_in_list.append(tg)
    for ind in list(target_df['Indications']):
        if type(ind) == str:
            temp = ind.split('; ')
            indications.append(temp)
        else:
            indications.append(ind)
    del drug_targets['Targets']  # delete the original Targets column
    del drug_targets['Indications']  # delete the original Indications column
    drug_targets['Indications'] = indications  # assign the new indications column
    drug_targets['Targets'] = targets_in_list  # assign the new targets column
    drug_targets.dropna(subset=['Targets'], inplace=True)  # delete the rows with no targets
    drug_targets.set_index(['Product'], inplace=True)

    # drug-indication summary
    if not os.path.exists(f'{output_path}/drug_targets.csv'):
        drugs = pd.DataFrame(0, index=list(drug_targets.index))
        ind = 0
        indi = []
        # drug's hallmark score = no. of targets with this hallmark/no. of targets of the drug
        for drug in drug_targets.index:
            indi.append(drug_targets.loc[drug, 'Indications'])
            ind += 1
        drugs['Indications'] = indi
        pd.DataFrame.to_csv(drugs, f'{output_path}/drug_targets.csv')  # save as a file

    else:
        # drugs is a dataframe with all the drugs for cancer and their hallmark score
        drugs = pd.read_csv(f'{output_path}/drug_targets.csv', sep=',',
                            index_col=0, header=0)

    # remove the nodes not in the cancer network
    cancer_genes = remove_nodes(cancer_genes, cancer_network)
    target_nodes = remove_nodes(target_nodes, cancer_network)

    # import and clean the dataset
    cancer_PEN_distance_path = f'{output_path}/{cancer_name}_{nodetype}/PEN_distance/alpha = 0.2'
    cancer_ppr_path = f'{output_path}/{cancer_name}_{nodetype}/ppr/alpha = 0.2'
    cancer_dist_path = f'{output_path}/{cancer_name}_{nodetype}/distance'
    PEN_distance_df = pd.read_csv(f'{cancer_PEN_distance_path}/PEN_distance0.txt', sep='\t', index_col=0, header=0)
    ppr_df = pd.read_csv(f'{cancer_ppr_path}/ppr.txt', sep='\t', index_col=0, header=0)
    dist_df = pd.read_csv(f'{cancer_dist_path}/distance.txt', sep='\t', index_col=0, header=0)
    othernodes = []  # other nodes means non-oncogenes
    for u in cancer_network.nodes():
        if u not in cancer_genes:
            othernodes.append(u)

    PEN_cancer_df = PEN_distance_df.loc[:, cancer_genes]
    PEN_other_df = PEN_distance_df.loc[:, othernodes]
    ppr_cancer_df = ppr_df.loc[:, cancer_genes]
    ppr_other_df = ppr_df.loc[:, othernodes]
    dist_cancer_df = dist_df.loc[:, cancer_genes]
    dist_other_df = dist_df.loc[:, othernodes]
    drug_df = drugs.copy()

    # select only rows with the cancer type we want
    for i in drug_df.index:
        if type(drug_df.loc[i, 'Indications']) == str:
            if cancer_name.lower() not in drug_df.loc[i, 'Indications'].lower():
                drug_df.drop(index=i, inplace=True)  # drugs with k targets for the input cancer type
        else:
            drug_df.drop(index=i, inplace=True)

    # summary of known target combinations
    # if k nodes are all from the same drug, then they can make up a known target combination
    PEN_distance_cancergenes_dict = {}
    PEN_distance_othergenes_dict = {}
    PEN_diff_dict = {}
    dist_cancergenes_dict = {}
    dist_othergenes_dict = {}
    dist_diff_dict = {}
    ppr_cancergenes_dict = {}
    ppr_othergenes_dict = {}
    ppr_diff_dict = {}
    for drug in drug_df.index:  # for each drug
        temp_target1 = drug_targets.loc[drug, 'Targets']
        temp_target = remove_nodes(temp_target1, cancer_network)
        for subset in itertools.combinations(temp_target, k):  # for each target combination in the drug
            # calculate the mean ppr,distance,PEN-distance
            mean_po = (PEN_cancer_df.loc[list(subset), :]).mean(axis=1).mean()  # avg PEN-distance to oncogenes
            mean_pr = (PEN_other_df.loc[list(subset), :]).mean(axis=1).mean()  # avg PEN-distance to other nodes
            mean_do = (dist_cancer_df.loc[list(subset), :]).mean(axis=1).mean()  # avg distance to oncogenes
            mean_dr = (dist_other_df.loc[list(subset), :]).mean(axis=1).mean()  # avg distance to other nodes
            mean_ppro = (ppr_cancer_df.loc[list(subset), :]).mean(axis=1).mean()  # avg ppr to oncogenes
            mean_pprr = (ppr_other_df.loc[list(subset), :]).mean(axis=1).mean()  # ppr to other
            pen_diff = mean_pr - mean_po  # PEN-diff
            distance_diff = mean_dr - mean_do  # distance-diff
            ppr_diff = mean_pprr - mean_ppro  # ppr-diff
            PEN_distance_cancergenes_dict[subset] = mean_po
            PEN_distance_othergenes_dict[subset] = mean_pr
            PEN_diff_dict[subset] = pen_diff
            dist_cancergenes_dict[subset] = mean_do
            dist_othergenes_dict[subset] = mean_dr
            dist_diff_dict[subset] = distance_diff
            ppr_cancergenes_dict[subset] = mean_ppro
            ppr_othergenes_dict[subset] = mean_pprr
            ppr_diff_dict[subset] = ppr_diff

    # combine the above dataset into one dataframe
    known_targets = pd.concat([pd.Series(d) for d in
                               [PEN_distance_cancergenes_dict,
                                PEN_distance_othergenes_dict,
                                PEN_diff_dict,
                                dist_cancergenes_dict,
                                dist_othergenes_dict,
                                dist_diff_dict,
                                ppr_cancergenes_dict,
                                ppr_othergenes_dict,
                                ppr_diff_dict
                                ]], axis=1)
    known_targets.columns = ['PEN Distance(cancer genes)',
                             'PEN Distance(other genes)',
                             'PEN-diff',
                             'Distance(cancer genes)',
                             'Distance(other genes)',
                             'Distance-diff',
                             'ppr(cancer genes)',
                             'ppr(other genes)',
                             'ppr-diff'
                             ]

    # write to file
    known_targets.to_csv(f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_known_targets.csv',
                         header=True, index=True, sep=',')
    print('known targets done')

    # use the same way as above to profile the candidate combinations
    if not os.path.exists(f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}set_combo.txt'):
        PEN_distance_cancergenes_dict = {}
        PEN_distance_othergenes_dict = {}
        PEN_diff_dict = {}
        dist_cancergenes_dict = {}
        dist_othergenes_dict = {}
        dist_diff_dict = {}
        ppr_cancergenes_dict = {}
        ppr_othergenes_dict = {}
        ppr_diff_dict = {}
        count = 0
        for subset in itertools.combinations(candidates, k):
            mean_po = (PEN_cancer_df.loc[list(subset), :]).mean(axis=1).mean()
            mean_pr = (PEN_other_df.loc[list(subset), :]).mean(axis=1).mean()
            mean_do = (dist_cancer_df.loc[list(subset), :]).mean(axis=1).mean()
            mean_dr = (dist_other_df.loc[list(subset), :]).mean(axis=1).mean()
            mean_ppro = (ppr_cancer_df.loc[list(subset),]).mean(axis=1).mean()
            mean_pprr = (ppr_other_df.loc[list(subset), :]).mean(axis=1).mean()
            pen_diff = mean_pr - mean_po  # difference between PEN_distance to cancergenes and other nodes
            distance_diff = mean_dr - mean_do
            ppr_diff = mean_pprr - mean_ppro
            PEN_distance_cancergenes_dict[subset] = mean_po
            PEN_distance_othergenes_dict[subset] = mean_pr
            PEN_diff_dict[subset] = pen_diff
            dist_cancergenes_dict[subset] = mean_do
            dist_othergenes_dict[subset] = mean_dr
            dist_diff_dict[subset] = distance_diff
            ppr_cancergenes_dict[subset] = mean_ppro
            ppr_othergenes_dict[subset] = mean_pprr
            ppr_diff_dict[subset] = ppr_diff
            count += 1
            print(f'{count},{subset}')

        cs_df = pd.concat([pd.Series(d) for d in
                           [PEN_distance_cancergenes_dict,
                            PEN_distance_othergenes_dict,
                            PEN_diff_dict,
                            dist_cancergenes_dict,
                            dist_othergenes_dict,
                            dist_diff_dict,
                            ppr_cancergenes_dict,
                            ppr_othergenes_dict,
                            ppr_diff_dict]], axis=1)

        cs_df.columns = ['PEN Distance(cancer genes)',
                         'PEN Distance(other genes)',
                         'PEN-diff',
                         'Distance(cancer genes)',
                         'Distance(other genes)',
                         'Distance-diff',
                         'ppr(cancer genes)',
                         'ppr(other genes)',
                         'ppr-diff']

        cs_df.to_csv(f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}set_combo.txt', header=True,
                     index=True, sep='\t')
