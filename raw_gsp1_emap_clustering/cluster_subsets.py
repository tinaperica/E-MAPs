


f_dir = r.cluster_emaps_dir

for file in os.listdir(f_dir):

    filename = os.fsdecode(file)
    f_str = filename.split('_')[0]
    handle = open(f_dir + '/' + filename)
    record = Cluster.read(handle)
    
    row_tree = record.treecluster(transpose=False, method='a', dist='c')
    row_tree.scale()
    col_tree = record.treecluster(transpose=True, method='a', dist='c')
    col_tree.scale()
    record.save(f_dir + '/' + f_str +'_clustered', row_tree, col_tree)