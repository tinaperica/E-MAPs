input_corr_dir <- "Output/correlations"
#input_corr_dir <- "~/Box Sync/kortemmelab/home/tina/Gsp1/E-MAP_analysis_backup_data/Aug2018_analysis/correlations"
outfilename <- "Output/preprocessed_correlations/weighted_corr"
#outfilename <- "~/Desktop/weighted_corr"
clust_i <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) )
#clust_i <- 1
# mutants <- c("CTER3XFLAG WT","D79A","D79S","E115A","E115I","F58A","F58L","G80A","GSP1-NAT","H141E",
#              "H141I","H141R","H141V","K101R","K129E","K129F","K129I","K129T","K132H","K143H", "K143W",
#              "K143Y","K154M","K169I","N102I", "N102K","N102M","N105L","N105V","N84Y",
#              "NTER3XFLAG WT","Q147E","Q147L","R108A","R108D","R108G","R108I","R108L","R108Q","R108S",
#              "R108Y","R112A","R112S","R78K","T137G","T139A","T139R","T34A","T34D","T34E",
#              "T34G","T34L","T34N","T34Q","T34S","T34Y","Y148I","Y157A", "A180T")
clusters <- c("endomembrane system", "protein targeting", "protein targeting_GO_15_2_mut", "protein targeting_GO_30_1_all",
              "protein targeting_GO_30_1_mut", "whole_library", "whole_library_30_5_all", "chromatin", "chromatin organization",
              "chromatin organization_GO_15_1_mut", "chromatin organization_GO_30_1_mut", "chromatin_GO_30_3_mut", "chromosome",
              "chromosome_GO_30_1_mut", "sig_Gsp1_GI", "sig_Gsp1_GI_30_1_mut", "sig_Gsp1_GI_30_2_all", "chromatin_GO_30_2_all",
              "chromosome_GO_30_2_all", "histone modification", "histone modification_GO_30_1_all", "response to heat",
              "response to heat_GO_15_1_all", "response to heat_GO_30_1_all", "transcription", "transcription and mRNA processing",
              "transcription from RNA polymerase II promoter", "whole_library_30_1_all",
              "cellular response to DNA damage stimulus", "cellular response to DNA damage stimulus_GO_30_1_mut", "DNA repair",
              "DNA repair_GO_15_1_all", "DNA repair_GO_30_1_mut", "transcription and mRNA processing_GO_15_14_mut",
              "transcription and mRNA processing_GO_30_2_mut", "transcription from RNA polymerase II promoter_GO_15_3_mut",
              "transcription from RNA polymerase II promoter_GO_30_1_all", "transcription from RNA polymerase II promoter_GO_30_1_mut",
              "whole_library_15_22_mut", "whole_library_30_14_mut", "budding", "cytoskeleton and microtubules", "cytoskeleton organization",
              "ribosomes and translation", "ribosomes and translation_GO_30_1_all", "RNA modification_GO_15_1_all",
              "transcription and mRNA processing_GO_15_8_all", "cell cycle", "cell cycle_GO_30_2_all",
              "DNA-templated transcription, elongation", "mitotic cell cycle", "mitotic cell cycle_GO_15_1_all",
              "regulation of cell cycle", "RNA catabolic process", "transcription and mRNA processing_GO_30_2_all",
              "whole_library_30_3_all", "proteolysis involved in cellular protein catabolic process",
              "proteolysis involved in cellular protein catabolic process_GO_15_1_all", "sig_Gsp1_GI_30_4_all",
              "sig_Gsp1_GI_30_5_mut", "chromatin_GO_15_1_all", "chromosome segregation", "chromosome segregation_GO_15_1_mut",
              "chromosome segregation_GO_30_1_all", "chromosome segregation_GO_30_1_mut", "mitotic cell cycle_GO_30_1_mut",
              "organelle fission", "regulation of cell cycle_GO_30_1_all", "regulation of cell cycle_GO_30_1_mut",
              "endomembrane system_GO_30_3_all", "ribosomal small subunit biogenesis",
              "ribosomal small subunit biogenesis_GO_15_1_all", "sig_Gsp1_GI_15_1_mut", "sig_Gsp1_GI_30_3_mut",
              "cellular response to DNA damage stimulus_GO_15_1_all", "cellular response to DNA damage stimulus_GO_30_3_all",
              "cell cycle_GO_15_3_all", "endomembrane system_GO_15_3_mut", "regulation of organelle organization_GO_15_1_all",
              "chromatin organization_GO_30_2_all", "DNA binding_GO_15_3_mut", "transcription_GO_15_1_mut",
              "transcription_GO_15_2_all", "transcription_GO_30_4_all", "nucleobase-containing compound transport",
              "nucleobase-containing compound transport_GO_15_1_mut", "nucleobase-containing compound transport_GO_15_2_all",
              "protein modification by small protein conjugation or removal", "protein modification by small protein conjugation or removal_GO_15_1_mut"
              , "protein modification by small protein conjugation or removal_GO_30_1_mut","sig_Gsp1_GI_15_11_all",
              "sig_Gsp1_GI_30_2_mut", "sig_Gsp1_GI_30_3_all", "cellular response to DNA damage stimulus_GO_15_3_mut", "chromatin organization_GO_15_5_all",
              "chromatin_GO_15_5_all", "DNA repair_GO_15_3_mut", "endomembrane system_GO_15_8_all", "lipid binding_GO_15_2_all","lipids_GO_15_2_all",
              "nuclear transport","nuclear transport and organization","nuclear transport and organization_GO_15_1_mut","nuclear transport and organization_GO_15_2_all",
              "nuclear transport and organization_GO_30_1_all","nuclear transport_GO_15_1_mut","nuclear transport_GO_15_2_all","peptidyl-amino acid modification_GO_15_2_mut",
              "protein targeting_GO_15_3_all","response to heat_GO_15_1_mut","whole_library_15_23_all","histone modification_GO_15_1_mut",
              "telomere organization","telomere organization_GO_15_1_mut","budding_GO_30_1_mut",
              "cytoskeleton organization_GO_15_1_mut","cytoskeleton organization_GO_15_3_all",
              "cytoskeleton organization_GO_30_1_all","cytoskeleton organization_GO_30_1_mut",
              "regulation of cell cycle_GO_15_2_mut","sig_Gsp1_GI_15_5_all","sig_Gsp1_GI_15_6_mut",
              "cytoskeleton_GO_15_2_mut","sig_Gsp1_GI_15_3_mut","chromatin_GO_15_3_mut",
              "DNA-templated transcription, elongation_GO_15_1_all","transcription and mRNA processing_GO_15_3_all",
              "RNA modification_GO_15_1_mut","transcription and mRNA processing_GO_30_1_all","whole_library_15_24_mut",
              "cytoplasmic translation","ribosome","ribosome_GO_15_1_all","structural constituent of ribosome",
              "sig_Gsp1_GI_15_4_mut","whole_library_15_14_mut","whole_library_30_12_mut","histone modification_GO_15_1_all",
              "cell cycle_GO_15_9_mut","cell cycle_GO_30_3_mut","proteolysis involved in cellular protein catabolic process_GO_15_2_mut",
              "sig_Gsp1_GI_15_1_all","cytoplasmic translation_GO_15_1_all","ribosomes and translation_GO_15_2_mut",
              "transcription and mRNA processing_GO_15_1_all","chromatin_GO_15_7_mut","chromosome segregation_GO_15_1_all",
              "cytoskeleton and microtubules_GO_15_2_mut","cytoskeleton and microtubules_GO_15_3_all",
              "cytoskeleton and microtubules_GO_30_3_all","nuclear transport and organization_GO_15_3_all",
              "nuclear transport_GO_15_1_all","nucleus organization","nucleus organization_GO_15_1_all",
              "nucleus organization_GO_15_1_mut","sig_Gsp1_GI_15_9_all","sig_Gsp1_GI_30_1_all","whole_library_15_36_all",
              "transcription and mRNA processing_GO_15_7_all","cell cycle_GO_15_5_all","RNA processing","budding_GO_15_1_all",
              "organelle fission_GO_30_1_mut","mRNA processing","RNA processing_GO_15_1_all","DNA recombination",
              "DNA recombination_GO_15_3_all","DNA recombination_GO_30_1_all","transcription from RNA polymerase II promoter_GO_15_4_all",
              "sig_Gsp1_GI_15_7_mut","mitotic cell cycle_GO_15_2_all","whole_library_15_3_all","whole_library_30_9_all",
              "sig_Gsp1_GI_15_3_all","sig_Gsp1_GI_15_9_mut","DNA recombination_GO_15_1_mut","DNA recombination_GO_30_1_mut",
              "organelle fission_GO_15_3_mut","DNA-templated transcription, elongation_GO_15_1_mut","whole_library_15_13_all",
              "cytoskeleton and microtubules_GO_15_5_all","chromosome_GO_15_2_mut","regulation of cell cycle_GO_15_4_all"  )
cluster <- clusters[clust_i]
correlations <- data.frame() ### contains the correlations per cluster
all_genes <- vector()
files_list <- list.files(path = input_corr_dir)
files_list <- files_list[grepl(files_list, pattern = "RData")]
#files_list <- c("20180702_7110586_correlations.RData", "20180702_7114206_correlations.RData", "20180702_1_correlations.RData", "20180702_7107871_correlations.RData", "20180702_1088716_correlations.RData")
#files_list <- files_list[1:30]
for (i in seq_along( files_list) ) {
  filename = files_list[i]
  print (filename)
  file_path <- file.path(input_corr_dir, filename)
  load(file_path)
  if (length(correlations_df) > 0) {
    if ( length(correlations_df[,1]) > 0 ) {
      correlations_df <- correlations_df[correlations_df$cluster == cluster,]
      if ( length(correlations_df[,1]) > 0 ) {
        temp_all_genes <- as.character(unique(c(correlations_df$Gene_uniq1, correlations_df$Gene_uniq2)))
        all_genes <- unique(append(all_genes, temp_all_genes))
        correlations_df$cluster <- NULL
        correlations <- rbind( correlations, correlations_df)
      }
    } 
  }
}
print("Done with reading files.")
corr_lst <- list()  ## final

#corr_lst[["correlations"]] <- correlations[[cluster]][! duplicated(correlations[[cluster]]), ]
corr_lst[["correlations"]] <- correlations
corr_lst[["cluster"]] <- cluster
save(corr_lst, file = paste0(outfilename, "_", clust_i, ".RData"))


