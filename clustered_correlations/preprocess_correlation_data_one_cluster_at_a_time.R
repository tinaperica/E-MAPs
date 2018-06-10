input_corr_dir <- "Output/files_to_sync"
outfilename <- "Output/preprocessed_correlations/corr_w_random_"
clust_i <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) )
mutants <- c("T34G","D79A","H141E","D79S","T34Q","R112S","R112A","R78K",
             "H141R","K101R","T34E","R108Y","NTER3XFLAG,WT","CTER3XFLAG,WT","R108G","R108Q",
             "Q147E","R108L","H141I","R108A","T34A","Y148I","G80A","Y157A",
             "R108S","R108I","K143Y","T34N","N84Y","E115I","K154M","T137G",
             "K143W","T139A","N105L","GSP1-NAT","K143H","K132H","K169I","K129F",
             "A180T","E115A","N105V","H141V","T34S","K129E","K129I","F58L",
             "N102I","T34D","T139R","N102K","T34L","T34Y","Q147L","F58A",
             "N102M","R108D","K129T", "GSP1-NAT", "NTER3XFLAG WT", "CTER3XFLAG WT")
clusters <- c("rRNA processing","transcription from RNA polymerase II promoter","RNA binding",
              "mRNA processing","hydrolase activity","regulation of cell cycle","mitochondrion",
              "ribosome","structural constituent of ribosome","mitochondrion organization","mitochondrial envelope",
              "oxidoreductase activity","transmembrane transporter activity","transferase activity","DNA recombination",
              "ion transport","transmembrane transport","organelle assembly","tRNA processing","DNA binding",
              "cytoplasmic translation","protein folding","nucleolus","protein targeting","chromosome",
              "nucleic acid binding transcription factor activity","transcription factor activity, protein binding",
              "cytoplasmic vesicle","endomembrane system","endosomal transport","cytoskeleton organization",
              "regulation of organelle organization","plasma membrane","vacuole","nuclear transport",
              "protein modification by small protein conjugation or removal",
              "proteolysis involved in cellular protein catabolic process","endoplasmic reticulum",
              "Golgi vesicle transport","phosphatase activity","lipid metabolic process","chromatin organization",
              "DNA replication","carbohydrate metabolic process","histone modification","regulation of DNA metabolic process",
              "response to heat","DNA repair","cellular response to DNA damage stimulus","cellular bud","cytoskeleton",
              "site of polarized growth","chromosome segregation","mitotic cell cycle","organelle fission","kinase activity",
              "protein phosphorylation","cell wall organization or biogenesis","meiotic cell cycle","sporulation",
              "chromatin binding","RNA modification","DNA-templated transcription, elongation","RNA catabolic process",
              "nucleobase-containing compound transport","Golgi apparatus","conjugation","endocytosis",
              "response to osmotic stress","enzyme regulator activity","regulation of protein modification process",
              "peptidyl-amino acid modification","protein acylation","pseudohyphal growth","lipid binding","telomere organization",
              "mRNA binding","cell cortex","ribosomes and translation","transcription and mRNA processing","transcription",
              "RNA processing","Golgi and ER","Golgi","ER","vacuoles","mitochondria","chromatin","cytoskeleton and microtubules",
              "cell cycle","budding","lipids","nuclear transport and organization","metabolic",
              "transcription from RNA polymerase II promoter_GO_3","transcription from RNA polymerase II promoter_GO_2",
              "transcription from RNA polymerase II promoter_GO_1","RNA binding_GO_1","hydrolase activity_GO_4",
              "hydrolase activity_GO_3","hydrolase activity_GO_2","hydrolase activity_GO_1","regulation of cell cycle_GO_2",
              "regulation of cell cycle_GO_1","mitochondrion_GO_2","mitochondrion_GO_1","ribosome_GO_2","ribosome_GO_1",
              "mitochondrion organization_GO_1","transferase activity_GO_6","transferase activity_GO_5",
              "transferase activity_GO_4","transferase activity_GO_3","transferase activity_GO_2","transferase activity_GO_1",
              "DNA recombination_GO_1","ion transport_GO_1","DNA binding_GO_3","DNA binding_GO_2","DNA binding_GO_1",
              "nucleolus_GO_1","protein targeting_GO_2","protein targeting_GO_1","chromosome_GO_3","chromosome_GO_2",
              "chromosome_GO_1","nucleic acid binding transcription factor activity_GO_1","cytoplasmic vesicle_GO_1",
              "endomembrane system_GO_6","endomembrane system_GO_5","endomembrane system_GO_4","endomembrane system_GO_3",
              "endomembrane system_GO_2","endomembrane system_GO_1","endosomal transport_GO_1","cytoskeleton organization_GO_2",
              "cytoskeleton organization_GO_1","regulation of organelle organization_GO_2",
              "regulation of organelle organization_GO_1","plasma membrane_GO_1","vacuole_GO_2","vacuole_GO_1",
              "nuclear transport_GO_1","protein modification by small protein conjugation or removal_GO_2",
              "protein modification by small protein conjugation or removal_GO_1","proteolysis involved in cellular protein catabolic process_GO_1",
              "endoplasmic reticulum_GO_3","endoplasmic reticulum_GO_2","endoplasmic reticulum_GO_1","Golgi vesicle transport_GO_1",
              "lipid metabolic process_GO_1","chromatin organization_GO_3","chromatin organization_GO_2","chromatin organization_GO_1",
              "DNA replication_GO_1","carbohydrate metabolic process_GO_1","histone modification_GO_1","regulation of DNA metabolic process_GO_1",
              "DNA repair_GO_2","DNA repair_GO_1","cellular response to DNA damage stimulus_GO_2","cellular response to DNA damage stimulus_GO_1",
              "cellular bud_GO_1","cytoskeleton_GO_2","cytoskeleton_GO_1","site of polarized growth_GO_2","site of polarized growth_GO_1",
              "chromosome segregation_GO_1","mitotic cell cycle_GO_3","mitotic cell cycle_GO_2","mitotic cell cycle_GO_1",
              "organelle fission_GO_3","organelle fission_GO_2","organelle fission_GO_1","kinase activity_GO_2","kinase activity_GO_1",
              "protein phosphorylation_GO_2","protein phosphorylation_GO_1","cell wall organization or biogenesis_GO_1",
              "meiotic cell cycle_GO_2","meiotic cell cycle_GO_1","sporulation_GO_1","chromatin binding_GO_1","Golgi apparatus_GO_1",
              "response to osmotic stress_GO_1","enzyme regulator activity_GO_2","enzyme regulator activity_GO_1",
              "regulation of protein modification process_GO_1","peptidyl-amino acid modification_GO_2","peptidyl-amino acid modification_GO_1",
              "mRNA binding_GO_1","cell cortex_GO_1","ribosomes and translation_GO_4","ribosomes and translation_GO_3","ribosomes and translation_GO_2",
              "ribosomes and translation_GO_1","transcription and mRNA processing_GO_6","transcription and mRNA processing_GO_5",
              "transcription and mRNA processing_GO_4","transcription and mRNA processing_GO_3","transcription and mRNA processing_GO_2",
              "transcription and mRNA processing_GO_1","transcription_GO_4","transcription_GO_3","transcription_GO_2","transcription_GO_1",
              "Golgi and ER_GO_7","Golgi and ER_GO_6","Golgi and ER_GO_5","Golgi and ER_GO_4","Golgi and ER_GO_3","Golgi and ER_GO_2",
              "Golgi and ER_GO_1","Golgi_GO_3","Golgi_GO_2","Golgi_GO_1","ER_GO_3","ER_GO_2","ER_GO_1","vacuoles_GO_1","mitochondria_GO_5",
              "mitochondria_GO_4","mitochondria_GO_3","mitochondria_GO_2","mitochondria_GO_1","chromatin_GO_5","chromatin_GO_4","chromatin_GO_3",
              "chromatin_GO_2","chromatin_GO_1","cytoskeleton and microtubules_GO_3","cytoskeleton and microtubules_GO_2",
              "cytoskeleton and microtubules_GO_1","cell cycle_GO_5","cell cycle_GO_4","cell cycle_GO_3","cell cycle_GO_2","cell cycle_GO_1",
              "budding_GO_5","budding_GO_4","budding_GO_3","budding_GO_2","budding_GO_1","lipids_GO_3","lipids_GO_2","lipids_GO_1",
              "nuclear transport and organization_GO_1","metabolic_GO_6","metabolic_GO_5","metabolic_GO_4","metabolic_GO_3","metabolic_GO_2",
              "metabolic_GO_1","whole_library_9","whole_library_8","whole_library_7","whole_library_6","whole_library_5","whole_library_4",
              "whole_library_3","whole_library_26","whole_library_25","whole_library_24","whole_library_23","whole_library_22","whole_library_21",
              "whole_library_20","whole_library_2","whole_library_19","whole_library_18","whole_library_17","whole_library_16",
              "whole_library_15","whole_library_14","whole_library_13","whole_library_12","whole_library_11","whole_library_10",
              "whole_library_1")
cluster <- clusters[clust_i]
correlations <- list() ### contains the correlations per cluster
all_genes <- vector()
files_list <- list.files(path = input_corr_dir)
for (i in seq_along( files_list) ) {
  filename = files_list[i]
  file_path <- file.path(input_corr_dir, filename)
  load(file_path)
  correlations_df <- correlations_df[correlations_df$cluster == cluster,]
  if ( length(correlations_df[,1]) > 0 ) {
    temp_all_genes <- as.character(unique(c(correlations_df$Gene_uniq1, correlations_df$Gene_uniq2)))
    all_genes <- unique(append(all_genes, temp_all_genes))
    correlations_df <- correlations_df[,c("Gene_uniq1", "Gene_uniq2", "correlation", "random_correlation", 
                      "random_correlation_2", "random_correlation_3", "random_correlation_4")]
    correlations[[cluster]] <- rbind( correlations[[cluster]], correlations_df)
  }
}
corr_lst <- list()  ## final

corr_lst[["correlations"]] <- correlations[[cluster]][! duplicated(correlations[[cluster]]), ]
corr_lst[["cluster"]] <- cluster
save(corr_lst, file = paste0(outfilename, "_", clust_i, ".RData"))


