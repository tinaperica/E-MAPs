To perform the analysis in this directory, run the scripts in the following order:

(1) Run 'make_SGD_descriptions.py'
(2) Run 'make_emap_for_clustering.R'
(3) Run 'find_unbiased_clusters.py' => NOW MADE IN R
    This cuts clustered emap into 1200 clusters, and considers
    those clusters with at least 3 members. After examining the SGD 
    descriptions, about half (21/44) strongly enrich one biological process or complex. These clusters were manually labeled and the
    file 'cluster_names.txt' was made
(4) Run 'make_emap_subsets.R'
(5) Run 'compute_rand_indexes.py'



Versions for Python:

Versions for R:




# reticulate for Python setup, including environment
library(reticulate)
use_python("/Users/cjmathy/miniconda3/bin/python")
use_condaenv('lab')
use_python("/Users/cjmathy/miniconda3/bin/python")