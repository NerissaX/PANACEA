from Panacea_Construct import *
from Panacea_profiling_comb import delta_histogram
from Panacea_Pen_diff import PEN_diff
from Panacea_constant import *


def panacea(cancer_name, k, cancer_genes, targets, nodetype, cancer_network):
    """
    Perform analysis for a given cancer subtype.

    Parameters:
    - cancer_name: str, name of the cancer subtype
    - k: int, size of the k-set
    - cancer_genes: list, genes associated with the cancer
    - targets: list, known drug targets
    - nodetype: str, type of nodes to use (e.g., 'oncogenes')
    - cancer_network: networkx.Graph, the cancer signaling network
    """

    # Remove targets not present in the cancer network
    targets_in_network = [t for t in targets if t in cancer_network.nodes()]
    num_targets_in_network = len(targets_in_network)

    # Get list of all nodes in the cancer network
    candidate_nodes = list(cancer_network.nodes())

    # Compute PEN_diff, ppr_diff, and distance_diff
    # Note: This function may take several hours for each cancer network
    PEN_diff(candidate_nodes, k, cancer_genes, nodetype, targets, cancer_network, cancer_name)

    # Write basic information about the cancer to a CSV file
    output_file = f'{output_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_basic_info.csv'
    with open(output_file, 'w') as f:
        f.write(f'Cancer Subtype:{cancer_name}\n')
        f.write(f'k-set:{k}\n')
        f.write(f'There are {len(targets)} known targets for this cancer subtype\n')
        f.write(f'There are {num_targets_in_network} known targets in the cancer network\n')
        f.write(f'Number of nodes in signaling network:{nx.number_of_nodes(whole_signaling)}\n')
        f.write(f'Number of nodes in the cancer network:{nx.number_of_nodes(cancer_network)}\n\n')

    # Find the delta histogram and range constraints
    deltamin, deltamax = delta_histogram(cancer_name, nodetype, k, bucket_no=5)


if __name__ == "__main__":
    nodetype = 'oncogenes'  # You can choose other node types like 'Biomarkers' or 'disease nodes'
    k = 2

    # List of cancer types to analyze
    cancers = [
        ('prostate', 'Prostate Cancer'),
        ('breast', 'Breast Cancer'),
        ('bladder', 'Bladder Cancer'),
        ('colorectal', 'Colorectal Cancer')
    ]

    # Analyze each cancer type
    for cancer_code, cancer_name in cancers:
        # Find oncogenes, target genes, and non-targets for the cancer
        cancer_oncogenes, cancer_targets, cancer_nontargets = find_subgene(
            whole_signaling, cancer_code, cancer_name
        )

        # Create the cancer signaling network
        cancer_network = constructNetwork(
            cancer_oncogenes, nodetype, cancer_targets, whole_signaling, cancer_name
        )
        print("construct network done")

        # Analyze the cancer network
        panacea(cancer_name, k, cancer_oncogenes, cancer_targets, nodetype, cancer_network)
