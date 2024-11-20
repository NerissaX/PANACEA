import os.path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from Panacea_constant import *


plt.rc('font', size=12)  # controls default text size
plt.rc('axes', titlesize=14)  # fontsize of the title
plt.rc('axes', labelsize=12)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=12)  # fontsize of the y tick labels
plt.rc('legend', fontsize=12)  # fontsize of the legend
plt.rcParams['figure.dpi'] = 600


# plot the distribution of ppr, dppr, PEN-distance (featurename)
# df: a dataframe of a specific feature of the network
# cancername: string: the name of cancer (e.g.:"Breast Cancer")
# featurename: string: "ppr", "PEN-distance", etc.
def plot_features(df, cancername, nodetype, featurename):
    print(f'plot_features,{cancername}')
    if not os.path.exists(f'{output_path}//{cancername}_{nodetype}//plot'):
        os.makedirs(f'{output_path}//{cancername}_{nodetype}//plot')
    cancername2 = cancername
    for per in [0.01,0.05,0.07,0.1,0.15]:
        change = per * 100
        if f'{change}' in cancername:
            print(True)
            cancername2 = cancername.replace(f'{change}',f'{int(change)}')
            break
    df.stack().plot.hist(bins=50, edgecolor='black')
    plt.yscale('log')
    mn, mx = plt.xlim()
    plt.xlim(mn, mx)
    plt.ylabel("Count")
    plt.xlabel(featurename)
    plt.title(f'{cancername2}',weight = 'bold')  # bold the text of the title to make it more visible
    plt.savefig(f'{output_path}//{cancername}_{nodetype}//plot//{cancername2}_{featurename} Distribution.png')
    plt.close('all')


# plot the delta histogram of a cancer specific signaling network
# df: dataframe of PEN-diff, ppr-diff, distance-diff
# fn: the name of the feature: "PEN-diff", "ppr-diff", etc.
def plot_deltaHist(df, cancer_name, nodetype, fn):
    print(f'plot_deltaHist,{cancer_name}')
    cancer_name2 = cancer_name
    for per in [0.01,0.05,0.07,0.1,0.15]:
        change = per * 100
        if f'{change}' in cancer_name:
            print(True)
            cancer_name2 = cancer_name.replace(f'{change}',f'{int(change)}')
            break
    plt.figure()
    df.plot(x='range',
            kind='bar',
            stacked=False)
    plt.title(f'{cancer_name2}',weight = 'bold')
    plt.xlabel('buckets')
    plt.ylabel('percentage of known target combinations')
    plt.xticks(rotation=25)
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
    plt.tight_layout()
    plt.savefig(
        f'{output_path}//{cancer_name}_{nodetype}//plot//{cancer_name2}_{nodetype}_{fn}_percentage_plot_des.png')
    plt.close()


# plot the distribution of PEN-diff, ppr-diff or distance-diff
# ls: a list of PEN-diff, ppr-diff or distance-diff for a cancer network
def plot_diff_distribution(ls, fn, cancer_name, nodetype):
    cancer_name2 = cancer_name
    print(f'plot diff distribution:{cancer_name}')
    for per in [0.01, 0.05, 0.07, 0.1, 0.15]:
        change = per * 100
        if f'{change}' in cancer_name:
            print(True)
            cancer_name2 = cancer_name.replace(f'{change}', f'{int(change)}')
            break
    plt.figure()  # print the distribution histogram
    # find the maximum size of bin in the distribution and print
    counts, edges, bars = plt.hist(ls, bins=100, edgecolor='black')
    counts = list(counts)
    print(f'There are {np.nanmax(counts)} combinations in the highest bar of {fn}')  # maximum count of candidates
    plt.xlabel(fn)
    plt.ylabel('count')
    # plt.yscale('log')
    plt.title(f'{cancer_name2}',weight = 'bold')
    plt.tight_layout()
    plt.savefig(f'{output_path}//{cancer_name}_{nodetype}//plot//{cancer_name2}_{nodetype}_{fn}_distribution.png')
