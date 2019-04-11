
from Bio import Cluster
from collections import defaultdict


def cut_and_get_clusters(tree, n_clusters, list_of_ids):

    assigned_clusters = tree.cut(n_clusters) # cut tree to get cluster assignments

    # make a dict of cluster:list pairs, where the list contains member names
    clusters = defaultdict(list)
    for i, cluster in enumerate(assigned_clusters):
        clusters[cluster].append(list_of_ids[i])
        
    return dict(clusters)


# Python function cluster_dataset


# we read in the dataset (EMAP) as a handle, then read the handle into a record object
# our record contains:
# -- the emap scores in record.data
# -- the mutant names in record.geneid
# -- the library gene names in record.expid
# -- a "mask" matrix showing 0's for missing data in record.mask
# -- gorder and eorder are lists giving the preferred order for rows/cols
    
def cluster_dataset(filename, name=None, rows=True, cols=True,
                    method='a', dist='c', row_order=None, col_order=None):


    with open(filename) as handle:
        record = Cluster.read(handle)
    
    if row_order:
        record.gorder = [row_order.index(row) for row in record.geneid]
    if col_order:
        record.eorder = [col_order.index(col) for col in record.expid]
    
    # initialize trees as 'None' in case only one axis is being clustered
    row_tree = None
    col_tree = None
    
    # cluster rows (mutants)
    if rows:
        row_tree = record.treecluster(transpose=False, method=method, dist=dist)
        row_tree.scale() # scale to [0,1] for ease of viewing in Java TreeView
    
    # cluster columns (library genes)
    if cols:
        col_tree = record.treecluster(transpose=True, method=method, dist=dist)
        col_tree.scale()
    
    if not name:
        name = filename.split('.')[0]
    
    record.save(name, row_tree, col_tree)

    return record, row_tree, col_tree