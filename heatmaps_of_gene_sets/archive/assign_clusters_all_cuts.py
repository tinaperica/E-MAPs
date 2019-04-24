

# Goal here: write a single function that takes in a clustering of mutants (row_tree), and spits out a matrix of size 56 x 56,
# where each column j is filled with the cluster number of mutant i when the tree is cut into j clusters. Just order the matrix to ascend alphabetically and numerically for simplicity.

# make output dir:

outdir = 'cluster_assigment_matrices'

if not os.path.exists(outdir):
    os.makedirs(outdir)

# with open(outdir + '/cluster_assigment_log.txt', 'w') as f:
#     f.write('# File format:\n')
#     f.write('''Rows are mutants, column i contain the cluster assignments
#                when cutting tree into i clusters\n\n''')


datasets = ['gsp1_emap_for_clustering.txt']
fdir = 'clustered_emaps'
datasets += [fdir + '/' + file for file in os.listdir(fdir) if file[-4:] == '.txt']

for dataset in datasets:
    
    if len(dataset.split('/')) > 1:
        name = dataset.split('.')[0].split('/')[1]
    else:
        name = dataset.split('.')[0]

    # Cluster mutants (rows) and arrays (cols)
    with open(dataset) as handle:
        record = Cluster.read(handle)
    row_tree = record.treecluster(transpose=False, method='a', dist='c')
    row_tree.scale() # scale to [0,1] for ease of viewing in Java TreeView
    col_tree = record.treecluster(transpose=True, method='a', dist='c')
    col_tree.scale()

    # cut tree to get cluster assignments fo
    n_elements = len(record.geneid)
    mat = np.zeros((n_elements, n_elements))

    for i in range(n_elements):
        mat[:, i] = row_tree.cut(i+1) # cut into 1 <= i+1 <= n_elements

    assignment_mat = pd.DataFrame(data = mat, dtype = int,
                                  index = record.geneid,
                                  columns = range(1, n_elements+1))
    assignment_mat.to_csv(outdir + '/{}_assigned_mat.txt'.format(name))

    with open(outdir + '/cluster_assigment_log.txt', 'a+') as f:
        f.write('# Dataset Name\n{}\n'.format(name))
        f.write('# Date and time of clustering:\n{}\n\n'.format(datetime.datetime.now()))

