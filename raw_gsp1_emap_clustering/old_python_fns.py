# Old python functions for Bio.Cluster


def cut_and_get_clusters(tree, nclusters, list_of_ids):
    
    # cut tree to get cluster assignments
    assigned_clusters = tree.cut(nclusters)
    
    # count the number of members in each cluster
    _, n_in_cluster = np.unique(assigned_clusters, return_counts=True)
    
    # count the number of clusters of a given size, return as a dict
    n_members, n_clusters_of_that_size = np.unique(n_in_cluster, return_counts=True)
    cluster_count_by_size = dict(zip(n_members, n_clusters_of_that_size))
    
    # make a dict of cluster:list pairs, where the list contains member names
    clusters = defaultdict(list)
    for i, cluster in enumerate(assigned_clusters):
        clusters[cluster].append(list_of_ids[i])
        
    return dict(clusters), cluster_count_by_size

