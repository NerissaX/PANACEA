import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from PANACEA_constant import *

plt.rc('font', size=12)  # controls default text size
plt.rc('axes', titlesize=12)  # fontsize of the title
plt.rc('axes', labelsize=12)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=12)  # fontsize of the y tick labels
plt.rc('legend', fontsize=12)  # fontsize of the legend
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['figure.dpi'] = 600


# plot the distribution of ppr, dppr, PEN-distance (featurename)
def plot_features(df, cancername, nodetype, featurename):
    df.stack().plot.hist(bins=50, edgecolor='black')
    plt.yscale('log')
    mn, mx = plt.xlim()
    plt.xlim(mn, mx)
    plt.ylabel("Count")
    plt.xlabel(featurename)
    plt.title(f'{cancername}',weight = 'bold')
    plt.savefig(f'{output_path}//{cancername}_{nodetype}//plot//{cancername}_{featurename} Distribution.png')
    plt.close('all')


def plot_deltaHist(df, cancer_name, nodetype, fn):
    plt.figure()
    df.plot(x='range',
            kind='bar',
            stacked=False)
    plt.title(f'{cancer_name}',weight = 'bold')
    plt.xlabel('buckets')
    plt.ylabel('percentage of known target combinations')
    plt.xticks(rotation=25)
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
    plt.tight_layout()
    plt.savefig(
        f'{output_path}//{cancer_name}_{nodetype}//plot//{cancer_name}_{nodetype}_{fn}_percentage_plot_des.png')
    plt.close()


def plot_diff_distribution(ls, fn, cancer_name, nodetype):
    plt.figure()  # print the distribution histogram
    # find the maximum size of bin in the distribution and print
    counts, edges, bars = plt.hist(ls, bins=100, edgecolor='black')
    counts = list(counts)
    print(f'There are {np.nanmax(counts)} combinations in the highest bar of {fn}')
    plt.xlabel(fn)
    plt.ylabel('count')
    # plt.yscale('log')
    plt.title(f'{cancer_name}',weight = 'bold')
    plt.tight_layout()
    plt.savefig(f'{output_path}//{cancer_name}_{nodetype}//plot//{cancer_name}_{nodetype}_{fn}_distribution.png')
