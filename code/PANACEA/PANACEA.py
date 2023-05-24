from PANACEA_Construct import *
from PANACEA_Profiling import delta_histogram
from PANACEA_PEN_diff import PEN_diff
from PANACEA_constant import *


def panacea(cancer_name, k, cancer_genes, target1, nodetype):
    # create the cancer signaling network using construcNetwork
    cancer_network = constructNetwork(cancer_genes, nodetype, target1, whole_signaling, cancer_name)
    # some targets may not in the signaling network because the human signaling network and drug targets are from
    # different data sources, we need to remove them first
    target_in_network = target1.copy()
    t_in_network = 0
    for u in target1:
        if u not in cancer_network.nodes():
            target_in_network.remove(u)
        else:
            t_in_network += 1

    # candidate3 is the list of all nodes in the cancer signaling network
    candidate3 = list(cancer_network.nodes())

    # compute PEN_diff, ppr_diff, distance_diff
    # this function will take several hours for each cancer network
    PEN_diff(candidate3, k, cancer_genes, nodetype, target1, cancer_network, cancer_name)

    t3_in_network = 0
    for u in candidate3:
        if u in target1:
            t3_in_network += 1

    # print basic information about the cancer
    with open(f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_basic_info.csv', 'w') as f:
        f.write('Cancer Subtype:{}'.format(cancer_name))
        f.write('\n')
        f.write('k-set:{}'.format(k))
        f.write('\n')
        f.write('There are {} known targets for this cancer subtype'.format(len(target1)))
        f.write('\n')
        f.write('There are {} known targets in the cancer network'.format(t_in_network))
        f.write('\n')
        f.write('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
        f.write('\n')
        f.write('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(cancer_network)))
        f.write('\n')
        f.write('\n')
    f.close()
    # find the delta histogram and range constraint of PEN-diff, ppr-diff, average closeness and distance-diff
    deltamin, deltamax = delta_histogram(cancer_name, nodetype, k)


if __name__ == "__main__":
    nodetype = 'oncogenes'  # the user can choose other cancer nodes like Biomarkers, disease nodes
    k = 2
    # prostate cancer
    # use find_subgene to find oncogenes, target genes and non-taragets of a specific cancer
    pc_onco1, pc_target1, pc_nontarget1 = find_subgene(whole_signaling, 'prostate', 'Prostate Cancer')
    # if a user wants to use other cancer genes instead of oncogenes, replace pc_onco1 in the input with other gene list
    panacea('Prostate Cancer', k, pc_onco1, pc_target1, nodetype)

    # breast cancer
    bc_onco1, bc_target1, bc_nontarget1 = find_subgene(whole_signaling, 'breast', 'Breast Cancer')
    panacea('Breast Cancer', k, bc_onco1, bc_target1, nodetype)

    # bladder cancer
    blc_onco1, blc_target1, blc_nontarget1 = find_subgene(whole_signaling, 'bladder', 'Bladder Cancer')
    panacea('Bladder Cancer', k, blc_onco1, blc_target1, nodetype)

    # colorectal cancer
    cc_onco1, cc_target1, cc_nontarget1 = find_subgene(whole_signaling, 'colorectal', 'Colorectal Cancer')
    panacea('Colorectal Cancer', k, cc_onco1, cc_target1, nodetype)
