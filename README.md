Overview

The Panacea Project is a computational framework designed to analyze cancer-specific signaling networks using network science methods.
It aims to identify potential drug targets and understand the complex interactions within
cancer networks by leveraging various network analysis techniques, including the computation of network statistics,
distances, and profiling of gene combinations.

Project Structure

The project consists of several Python modules, each responsible for specific tasks within the analysis pipeline.
Below is a summary of each module, their purposes, and how they interact with one another.

1. Panacea_constant.py

	•	Purpose: Defines constants, global variables, datasets, and file paths used throughout the project.
	•	Key Data Loaded:
        •	Drug Targets Dataset: Information on known drug targets.
        •	Oncogenes Dataset: List of known oncogenes.
        •	Signaling Network Dataset: The comprehensive signaling network used for analysis.
	•	Usage: Must be run first to set up the environment. It is imported by other modules that require access to the constants and datasets.

2. Panacea_plot.py

	•	Purpose: Provides functions for generating plots and visualizations.
	•	Key Functions:
        •	plot_features(...): Plots histograms and distributions of network features.
        •	plot_diff_distribution(...): Plots the distribution of differences in features.
        •	plot_deltaHist(...): Plots delta histograms for profiling.
	•	Usage: Used by multiple modules to generate plots for analysis and visualization.

3. PanaceaStats.py

	•	Purpose: Computes and saves network statistics such as degree centrality, closeness centrality, betweenness centrality, and eigenvector centrality.
	•	Key Functions:
        •	compute_network_stats(graph, filename, save_to): Calculates various network statistics and writes them to files.
        •	plot_network_stats(...): Generates histograms for the computed statistics.
	•	Usage: Called by Panacea_Construct.py to compute statistics of the constructed cancer networks.

4. Panacea_PEN_dist.py

	•	Purpose: Computes the Personalized PageRank (PPR), differential PPR (dPPR), and PEN distance for a given network.
	•	Key Functions:
        •	compute_PEN_distance(...): Computes PPR, dPPR, and PEN distance for the network.
        •	PEN_distance_alpha(...): Computes these metrics for different alpha values.
	•	Usage: Used by Panacea_Construct.py to compute distances after constructing or modifying networks.

5. Panacea_distance.py

	•	Purpose: Computes the shortest path lengths between all pairs of nodes in a given network.
	•	Key Function:
	    •	compute_shortest_distance(G, save_to): Calculates shortest distances and saves them.
	•	Usage: Called by Panacea_Construct.py to compute shortest distances in the network.

6. Panacea_Construct.py

	•	Purpose:
        •	Identifies genes associated with specific cancer subtypes.
        •	Constructs the smallest connected network containing all input cancer genes.
        •	Computes network statistics, distances, and plots features.
	•	Key Functions:
        •	find_subgene(...): Finds oncogenes, target genes, and non-target genes for a specific cancer subtype.
        •	constructNetwork(...): Builds the cancer network and computes required metrics.
	•	Usage: Called from the main script to construct cancer networks for each cancer type.

7. Panacea_Pen_diff.py

	•	Purpose: Computes differences in PEN-distance, Distance, and PPR between cancer genes and other nodes for known drug target combinations and candidate combinations.
	•	Key Function:
        •	PEN_diff(...): Performs the computation and saves results.
	•	Usage: Used in the main script after constructing the network to compute differences and prepare data for profiling.

8. Panacea_profiling_comb.py

	•	Purpose: Performs delta histogram profiling of the cancer network to identify the best constraints for selecting target combinations.
	•	Key Function:
        •	delta_histogram(...): Profiles combinations and generates delta histograms.
	•	Usage: Called from the main script to analyze and profile the combinations based on computed differences.

9. Panacea_target.py

	•	Purpose: Writes the target combinations within the constraint range to a file for further analysis or reporting.
	•	Key Function:
        •	target(...): Saves the constrained target combinations identified during profiling.
	•	Usage: Called by Panacea_profiling_comb.py.

10. Main Script (panacea.py)

	•	Purpose: Orchestrates the entire analysis process by coordinating the use of all the modules above.
	•	Key Steps:
        •	Imports necessary modules.
        •	Loops over cancer types to perform analysis.
        •	Calls functions to find genes, construct networks, compute differences, and profile combinations.
	•	Usage: Execute this script to run the entire analysis pipeline.

Prerequisites

	•	Programming Language: Python 3.x
	•	Required Python Libraries:
	•	pandas
	•	numpy
	•	networkx
	•	matplotlib
	•	seaborn
	•	scipy
	•	sklearn

Data Preparation

	•	Place the datasets (e.g., drug targets, oncogenes, signaling network) in the appropriate directories as expected by Panacea_constant.py.

Running the Analysis

	1.	Initialize Constants and Datasets:
        Ensure that Panacea_constant.py is correctly set up with paths to your datasets.
	2.	Execute the Main Script:
	    python panacea.py

	This script will perform the full analysis pipeline.

Analysis Workflow

For each cancer type:
	1.	Find Associated Genes:
	    •	Uses find_subgene to identify oncogenes and target genes.
	2.	Construct Cancer Network:
	    •	Uses constructNetwork to build the network and compute statistics.
	3.	Compute Differences:
	    •	Uses PEN_diff to calculate differences in distances between genes.
	4.	Profile Combinations:
	    •	Uses delta_histogram to perform profiling and identify optimal constraints.
	5.	Identify Targets:
	    •	Uses target to save target combinations within constraints.

Results

	•	Output Files:
	•	Network statistics and metrics are saved in specified directories.
	•	Plots and histograms are saved as image files.
	•	Target combinations and profiling results are saved as text files.

Project Dependencies

The following modules are interdependent:
	•	Panacea_constant.py is foundational and required by all other modules.
	•	Panacea_plot.py is used by multiple modules for visualization.
	•	PanaceaStats.py, Panacea_PEN_dist.py, and Panacea_distance.py are critical for network computations.
	•	Panacea_Construct.py relies on the above modules to build and analyze networks.
	•	Panacea_Pen_diff.py and Panacea_profiling_comb.py are used sequentially for profiling combinations.
	•	Panacea_target.py saves the results from profiling.

panacea-project//
├── Panacea_constant.py
├── PanaceaOtherFunc.py
├── Panacea_plot.py
├── PanaceaStats.py
├── Panacea_PEN_dist.py
├── Panacea_distance.py
├── Panacea_Construct.py
├── Panacea_Pen_diff.py
├── Panacea_profiling_comb.py
├── Panacea_target.py
├── panacea.py
├── datasets/
│   ├── drug_targets.csv
│   ├── oncogenes.csv
│   └── signaling_network.gexf
├── outputs/
│   ├── cancer_type_nodetype/
│   │   ├── stats/
│   │   ├── ppr/
│   │   ├── dppr/
│   │   ├── PEN_distance/
│   │   ├── distance/
│   │   └── plots/
│   └── ...
└── README.txt
