from Bio import Cluster
# Cluster.__version__
import os, sys, datetime, copy
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def cluster_and_assign_all_cuts(file, dist):

    with open(file) as handle:
        record = Cluster.read(handle)
    row_tree = record.treecluster(transpose=False, method='a', dist=dist)
    row_tree.scale()
    
    record.save(name.split('.')[0], geneclusters=row_tree)
    
    # cut tree into n clusters, for 1 <= n <= n_elements
    n_elements = len(record.geneid)
    mat = np.zeros((n_elements, n_elements))

    for i in range(n_elements):
        mat[:, i] = row_tree.cut(i+1) # cut into 1 <= i+1 <= n_elements

    return pd.DataFrame(data = mat,
                        dtype = int,
                        index = record.geneid,
                        columns = range(1, n_elements+1))
        
        
def rand_index(c1, c2):
    '''
    This function computes the rand index score for two
        clusterings of the same dataset (note that this
        assumes len(c1) = len(c2))
    
    The rand index is the fraction of agreements between two
        clusterings. This is equal to 1 minus the fraction of
        disagreements.
        
    Adapted from the fossil package in R, code walkthrough in
    https://davetang.org/muse/2017/09/21/the-rand-index/
    
    inputs:
        c1 - array-like object containing integers assigning
                element i to a cluster denoted by value at i
        c2 - array-like object containing integers assigning
                element i to a cluster denoted by value at i
    
    returns:
        rand_index - a score between 0 and 1
    ''' 
        
    
    # coerce the vectors to numpy arrays
    c1 = np.asarray(c1)
    c2 = np.asarray(c2)
    
    # for both clusterings, make a matrix showing pairwise
    # cluster membership. i.e. cell [i,j] = 0 if i and j are
    # in the same cluster
    x = abs(c1[:, np.newaxis] - c1)
    x[np.where(x > 1)] = 1
    y = abs(c2[:, np.newaxis] - c2)
    y[np.where(y > 1)] = 1    
    
    # sum to get the number of disagreements
    # divide by two bc the matrix is symmetric
    n_disagreements = np.sum(abs(x-y))/2
    
    # get the total number of pairs
    n_total_pairs = sp.special.comb(len(c1), 2)
    
    return(1 - n_disagreements/n_total_pairs)




    # make list of tuples of form (name, file), where name will be used in graphs
datasets = [('full_emap', 'gsp1_emap_for_clustering-GAPGEF.txt'),
            ('interface', 'core_deltarASA_by_mutant_for_clustering-GAPGEF.txt'), # TODO: this clustering is probably bad bc of the NAs
            ('GAP_GEF_effic', 'GAP_GEF_effic_for_clustering.txt')]

# add the emap subsets
fdir = 'clustered_emaps'
files = [fdir + '/' + file for file in os.listdir(fdir) if (
            file[-4:] == '.txt' and 'GAPGEF' in file)]
names = [file.split('_')[0].split('-')[0] for file in os.listdir(fdir) if (
            file[-4:] == '.txt' and 'GAPGEF' in file)]
datasets += [tup for tup in zip(names, files)]

# compute all possible cluster definitions for all datasets
clusterings = {}
for name, file in datasets:
    if name == 'GAP_GEF_effic':
        dist = 'e'
    else:
        dist = 'c'
    clusterings[name] = cluster_and_assign_all_cuts(file, dist = dist)



    ref_clusterings= ['interface', 'GAP_GEF_effic']
ri_vecs_dict = dict.fromkeys([name for name in clusterings.keys()
                         if name not in ref_clusterings])                   
ri_sums_dict = copy.deepcopy(ri_vecs_dict)

for dset in ri_vecs_dict.keys():

    n = len(clusterings[dset])

    d1 = clusterings[dset]
    d2 = clusterings['interface']
    d3 = clusterings['GAP_GEF_effic']

    # cut the tree into i+1 clusters, with 0<=i<=n
    ri1 = [rand_index(d1[i+1], d2[i+1]) for i in range(n)] # interface
    ri2 = [rand_index(d1[i+1], d3[i+1]) for i in range(n)] # GAPGEF ratio

    ri_vecs_dict[dset] = (ri1, ri2)
    ri_sums_dict[dset] = [np.sum(ri1)/n, np.sum(ri2)/n]




pp = PdfPages('RandIndex_vs_full_emap.pdf')

for dset in ri_vecs_dict:

    interface, GAPGEF = ri_vecs_dict[dset]
    full_emap_interface, full_emap_GAPGEF = ri_vecs_dict['full_emap']

    interface_sum, GAPGEF_sum = ri_sums_dict[dset]
    full_emap_interface_sum, full_emap_GAPGEF_sum = ri_sums_dict['full_emap']
    
    delta1 = round(interface_sum - full_emap_interface_sum, 3)
    delta2 = round(GAPGEF_sum - full_emap_GAPGEF_sum, 3)
        
    n = len(interface)

    fig, ax = plt.subplots(nrows=1, ncols=2,
                           sharex='all', sharey='all')  
    plt.axis([0, n, 0, 1.1])

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off',
                    bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel('number of clusters when tree was cut')
    plt.ylabel('Rand Index')
    plt.title('dataset = {}\n\n'.format(dset))
                  
    ax[0].plot(range(n), interface, 'b', label='subset vs. interface')
    ax[0].plot(range(n), full_emap_interface, 'r', label='full_emap vs. interface')
    ax[0].set_title('Interface\nΔ to full_emap: {}'.format(delta1))  
    ax[0].legend(fontsize='small')
              
    ax[1].plot(range(n), GAPGEF, 'b', label='subset vs. GAP_GEF_effic')
    ax[1].plot(range(n), full_emap_GAPGEF, 'r', label='full_emap vs. GAP_GEF_effic')
    ax[1].set_title('GAP/GEF effic\nΔ to full_emap: {}'.format(delta2))
    ax[1].legend(fontsize='small')
    fig.tight_layout()

    pp.savefig()
    plt.close()
    
pp.close()


pp = PdfPages('interface_vs_GAPGEF.pdf')

for dset in ri_vecs_dict:

    interface, GAPGEF = ri_vecs_dict[dset]
    interface_sum, GAPGEF_sum = ri_sums_dict[dset]        
    n = len(interface)

    plt.figure(figsize=(5,5))
    plt.plot(range(n), interface, 'b', label='subset vs. interface')
    plt.plot(range(n), GAPGEF, 'r', label='subset vs. GAP_GEF_effic')
    plt.axis([0, n, 0, 1.1])
    plt.xlabel('number of clusters when tree was cut')
    plt.ylabel('Rand Index')
    plt.title('dataset = {}'.format(dset))
    plt.tight_layout()
    plt.legend(loc=4)
    pp.savefig()
    plt.close()
    
pp.close()



ri_df = pd.DataFrame.from_dict(ri_sums_dict,
                               orient='index',
                               columns=('interface','GAP_GEF_effic'))

ri_df['Δfull_emap, interface'] = ri_df['interface'] - ri_df.loc['full_emap'].interface
ri_df['Δfull_emap, GAP_GEF_effic'] = ri_df['GAP_GEF_effic'] - ri_df.loc['full_emap'].GAP_GEF_effic
ri_df['Δinterface - ΔGAP_GEF_effic'] = ri_df['Δfull_emap, interface'] - ri_df['Δfull_emap, GAP_GEF_effic']

ax = ri_df[['Δfull_emap, GAP_GEF_effic',
            'Δfull_emap, interface',
            'Δinterface - ΔGAP_GEF_effic']].plot(kind="bar",
                                            title=title,
                                            fontsize=15,
                                            figsize=(10,5))

ax.set_xlabel('Complex / Process', size=15)
ax.set_ylabel('ΔAUC compared to full emap', size=15)
title='Does interface or GAP/GEF effic better match\nthe full emap or specific complexes?'
ax.set_title(title, size=15)
fig = ax.get_figure()
fig.savefig('delta_to_full_emap.png')