import numpy as np
import pandas as pd


def add_diagonal_placeholders(table):
    """
    This function adds pairwise interactions of the same gene
    to the list, giving them a score of "NaN". This ensures the dataframe
    has the same columns and rows, allowing it to be truly symmetric after
    unstacking.

    input:
        table - a pandas dataframe in long format, with each row containing
                two gene names in alphabetical order (Pert1 and Pert2)
    output:
        table - table has entires containing pairwise interactions of the same gene appended (i.e. "YAL002W" "YAL002W" "NaN" is one appended row)
    """

    geneA_set = set(table.geneA.unique())
    geneB_set = set(table.geneB.unique())
    all_genes = geneA_set.union(geneB_set)

    diagonal = pd.DataFrame(data={'geneA' : list(all_genes),
                                  'geneB' : list(all_genes),
                                  'score' : np.nan})
    return table.append(diagonal)


def make_symmetric_genetic_interaction_matrix(table):
    """
    Takes a long format genetic interaction table and returns a symmetric matrix

    input: table - a long format pandas dataframe
    output: mat - a wide format symmetrical pandas dataframe
    """

    # unstack to get upper triangle and lower triangle matrices
    # also convert every NaN to 0.00, allowing for algebra of dataframes
    upp = table.set_index(['geneA','geneB'])['score'].unstack().fillna(0.0)
    low = table.set_index(['geneB','geneA'])['score'].unstack().fillna(0.0)
    
    # sum the two triangular matrices to get the symmetric matrix, then return the 0.0s to NaNs
    return (low + upp).replace(0.0, np.nan)


# Read in datasets
emap = pd.read_csv('data/EMAP_for_scaling_long.txt', sep='\t')
sga = pd.read_csv('data/SGA_for_scaling_long.txt', sep='\t')

# Add diagonal entries in long format
emap = add_diagonal_placeholders(emap)
sga = add_diagonal_placeholders(sga)

# convert to symmetric matrix
emap = make_symmetric_genetic_interaction_matrix(emap)
sga = make_symmetric_genetic_interaction_matrix(sga)

# write out matrices to txt file
emap.to_csv('data/EMAP_for_scaling_symm.txt', sep='\t')
sga.to_csv('data/SGA_for_scaling_symm.txt', sep='\t')

# save SGA as a matlab struct in the format used by EMAP toolbox
sga_dict = {'rowlabels': np.array(sga.columns),
            'collabels': np.array(sga.columns),
            'data': sga.values}

scipy.io.savemat('data/sga_for_scaling_unaveraged.mat',
                 {'SGA': sga_dict},
                 oned_as = 'column')

