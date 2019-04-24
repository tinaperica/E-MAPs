#!/bin/bash
Rscript gene_set_heatmaps.R pearson score
Rscript gene_set_heatmaps.R pearson correlations
Rscript gene_set_heatmaps.R euclidean score
Rscript gene_set_heatmaps.R euclidean correlations

