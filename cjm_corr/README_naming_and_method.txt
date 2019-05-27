Correlations of correlations plots

Query genes and mutants were first filtered by number of shared significant S-scores, then standard pairwise pearson correlations were computed. This matrix of correlations was filtered for p < 0.05. Then, a correlation of correlations value was computed for each pair of queries.

Four sets of correlations of correlations were computed:
allQs, withGsp1
allQs, woutGsp1
APMSonly, withGsp1
APMSonly, woutGsp1

If a heatmap is labeled "allQs", it means correlations of correlations were computed with all queries (after the filtering). If a heatmap is labeled with APMSonly, the correlations of correlations were only computed using queries from the set of AP-MS SAINT hits.

If a heatmap is labeled withGsp1, the Gsp1 mutants were included in the correlations of correlations. If a heatmap is labeled woutGsp1, the mutants were not included in the correlations of correlations

To make the heatmaps, dendrograms were computed using euclidean distance as the metric, and complete linkage