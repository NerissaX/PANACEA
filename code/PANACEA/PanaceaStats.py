import matplotlib.pyplot as plt
from Panacea_constant import *
import networkx as nx


def rank(input_dict):
    """
    Ranks the values in a dictionary, handling ties appropriately.

    Parameters:
    - input_dict: dict
        A dictionary where keys are items (e.g., genes) and values are numerical scores to rank.

    Returns:
    - ranking_df: pandas.DataFrame
        A DataFrame containing 'Gene', 'Value', and 'Rank' columns, sorted by 'Value' in descending order.
        Tied values receive the same rank.
    """
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(list(input_dict.items()), columns=['Gene', 'Value'])

    # Sort the DataFrame by 'Value' in descending order
    df.sort_values(by='Value', ascending=False, inplace=True)

    # Compute ranks, assigning the same rank to equal values
    df['Rank'] = df['Value'].rank(method='dense', ascending=False).astype(int)

    # Reset index for neatness
    df.reset_index(drop=True, inplace=True)

    return df


def add_legend():
    """
    Adds a legend to the current plot, removing duplicate labels.
    """
    handles, labels = plt.gca().get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    plt.legend(unique.values(), unique.keys())


def compute_network_stats(graph, filename, save_to):
    """
    Computes various network statistics and writes them to files.

    Parameters:
    - graph: networkx.Graph, the network graph to analyze.
    - filename: str, the base name for output files.
    - save_to: str, the directory path to save output files.
    """
    stats_file = os.path.join(save_to, f'{filename}.txt')
    if not os.path.exists(stats_file):
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write(str(nx.info(graph)) + '\n')

            # Density
            density = nx.density(graph)
            f.write(f'Density of the network: {density:.4f}\n')
            print('Density computed')

            # Average clustering coefficient
            avg_clustering = nx.average_clustering(graph)
            f.write(f'Average clustering coefficient of the network: {avg_clustering:.4f}\n')
            print('Average clustering coefficient computed')

            # Strongly connected components
            scc = list(nx.strongly_connected_components(graph))
            num_scc = len(scc)
            largest_scc = max(scc, key=len)
            f.write(f'Number of strongly connected components: {num_scc}\n')
            f.write(f'Size of the largest strongly connected component: {len(largest_scc)}\n')

            # Weakly connected components
            wcc = list(nx.weakly_connected_components(graph))
            num_wcc = len(wcc)
            largest_wcc = max(wcc, key=len)
            f.write(f'Number of weakly connected components: {num_wcc}\n')
            f.write(f'Size of the largest weakly connected component: {len(largest_wcc)}\n')
            print('Components computed')

            # Centrality measures to compute
            centrality_measures = {
                'degree centrality': nx.degree_centrality(graph),
                'closeness centrality': nx.closeness_centrality(graph),
                'betweenness centrality': nx.betweenness_centrality(graph, weight='weight'),
                'eigenvector centrality': nx.eigenvector_centrality(graph, max_iter=600, weight='weight'),
                'degree': dict(graph.degree()),
                'in-degree': dict(graph.in_degree()),
                'out-degree': dict(graph.out_degree()),
                'pagerank': nx.pagerank(graph, weight='weight')
            }

            # Compute and save centrality measures
            for measure_name, measure_dict in centrality_measures.items():
                sorted_measure = rank(measure_dict)
                measure_filename = f'{filename}_sorted {measure_name.replace(" ", "_")}.txt'
                measure_filepath = os.path.join(save_to, measure_filename)

                # Save full centrality ranking
                sorted_measure.to_csv(measure_filepath, sep='\t', header=False, index=False)
                print(f'{measure_name.capitalize()} computed and saved to {measure_filename}')

                # Write top 10 nodes to the stats file
                f.write(f"\nTop 10 nodes by {measure_name}:\n")
                for _, row in sorted_measure.head(10).iterrows():
                    f.write(f"\t{row['Gene']} {row['Value']:.4f}\n")

            print('All centrality measures computed')


def plot_histogram(data, dataname, filename, save_to, mark_genes=True):
    """
    Plots a histogram of the provided data and saves the figure.

    Parameters:
    - data: pandas.DataFrame, the data to plot (assumes two columns: Gene and Value).
    - dataname: str, the name of the data (for labels and titles).
    - filename: str, the base name for the saved plot.
    - save_to: str, the directory path to save the plot.
    - mark_genes: bool, whether to mark specific genes (not implemented in this snippet).
    """
    values = data.iloc[:, 1].values

    plt.figure(figsize=(10, 6))
    plt.hist(values, bins=50, edgecolor='black')
    plt.yscale('log')
    plt.ylabel("Count")
    plt.xlabel(dataname)
    plt.title(f'{dataname} Distribution')
    plt.tight_layout()
    plot_filename = f'{filename} {dataname} Distribution.png'
    plt.savefig(os.path.join(save_to, plot_filename))
    plt.close()
    print(f'{dataname} histogram saved as {plot_filename}')


def read_stats(filename, stat_path):
    """
    Reads the centrality measure files and returns them as DataFrames.

    Parameters:
    - filename: str, the base name of the files to read.
    - stat_path: str, the directory path where the files are located.

    Returns:
    - A tuple of pandas.DataFrames containing the centrality measures.
    """
    measures = [
        'betweenness_centrality',
        'degree_centrality',
        'closeness_centrality',
        'eigenvector_centrality',
        'degree',
        'in_degree',
        'out_degree',
        'pagerank'
    ]
    dataframes = []
    for measure in measures:
        filepath = os.path.join(stat_path, f'{filename}_sorted {measure}.txt')
        df = pd.read_csv(filepath, sep='\t', header=None)
        dataframes.append(df)
    return tuple(dataframes)


def plot_network_stats(filename, stat_path, save_to):
    """
    Plots histograms for all centrality measures.

    Parameters:
    - filename: str, the base name of the files to read.
    - stat_path: str, the directory path where the centrality measure files are located.
    - save_to: str, the directory path to save the plots.
    """
    centrality_data = read_stats(filename, stat_path)
    measure_names = [
        'Betweenness Centrality',
        'Degree Centrality',
        'Closeness Centrality',
        'Eigenvector Centrality',
        'Degree',
        'In-Degree',
        'Out-Degree',
        'Pagerank'
    ]
    for data, name in zip(centrality_data, measure_names):
        plot_histogram(data, name, filename, save_to)
