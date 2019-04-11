# Set up Python

from Bio import Cluster
# Cluster.__version__
import numpy as np
import pandas as pd
import os, sys, datetime

from biocluster_fns import cluster_dataset, cut_and_get_clusters

# SECTION 1: BASIC CLUSTERING OF GSP1 EMAP DATA
dataset = 'gsp1_emap_for_clustering.txt'
method = 'a' # average linkage
dist = 'c' # centered pearson correlation as distance metric

with open('gsp1_mutants_ordered_by_residue.txt') as f:
    row_order = f.read().splitlines()

# Cluster mutants (rows) and arrays (cols)
name = 'full_emap_mut_array_clustered'
res = cluster_dataset(dataset, name=name, method=method,
                      dist=dist, row_order=row_order)
record, row_tree, col_tree = res

# For another visualization, keep mutants ordered, cluster arrays
name = 'full_emap_array_clustered'
cluster_dataset(dataset, name=name, method=method, dist=dist, row_order=row_order, rows=False)



SGD_descriptions = pd.read_csv('SGD_descriptions.txt', sep='\t')

n_clusters = 1200
cluster_dict = cut_and_get_clusters(col_tree, n_clusters, record.expid)
cluster_tuples_by_gene_name = [(v,k) for k,v in cluster_dict.items()]

gene_cluster_dict = {}
for cl in cluster_tuples_by_gene_name:
    for gene in cl[0]:
        gene_cluster_dict[gene] = cl[1]

cluster_df = pd.DataFrame(gene_cluster_dict.items(),
                          columns=['name', 'cluster'])

df = pd.merge(SGD_descriptions, cluster_df, how = 'left', on = 'name')


clusterdir = 'clusters'
min_n_to_print = 3

if not os.path.exists(clusterdir):
    os.makedirs(clusterdir)

df.to_csv('unbiased_array_clusters_with_descriptions.txt',
          sep='\t', na_rep='', index=False)

for cl, cluster in df.groupby('cluster'):
    if len(cluster) >= min_n_to_print:
        cluster.to_csv(clusterdir + '/cluster_{}.csv'.format(cl))

with open(clusterdir + '/cluster_info.txt', 'w') as f:
    f.write('date and time of clustering: {}\n'.format(datetime.datetime.now()))
    f.write('number of clusters cut: {}\n'.format(n_clusters))
    f.write('min members per cluster for printing: {}\n'.format(min_n_to_print))




