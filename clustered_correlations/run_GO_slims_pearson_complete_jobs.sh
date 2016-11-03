# !/bin/bash
#$ -cwd
#$ -o Output/GO_slims_pearson_complete/task-$TASK_ID.stout
#$ -e Output/GO_slims_pearson_complete/task-$TASK_ID.sterr
#$ -S /bin/bash
#$ -l h_rt=300:00:00
#$ -l mem_free=3G
#$ -t 1-8605026:56242
echo Starting task $SGE_TASK_ID
R CMD BATCH correlations_on_the_cluster_GO_slims_pearson_complete.R Output/GO_slims_pearson_complete/R-${SGE_TASK_ID}.Rout
