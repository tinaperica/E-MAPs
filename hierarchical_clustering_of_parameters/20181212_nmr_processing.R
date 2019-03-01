library(tidyverse)

nmr_raw_data <- as.tibble(read_delim('hierarchical_clustering_of_parameters/20190204_nmr_data.csv', delim = ','))

experiments.to.ignore <- c(1,4,19,8,3,7,13,2)
data.to.plot <- nmr_raw_data %>%
    filter(!experiment_id %in% experiments.to.ignore)

write_tsv(data.to.plot, path = "hierarchical_clustering_of_parameters/NMR_data_clean.txt")
ggplot(data.to.plot, aes(reorder(gsp1_mutation, gamma2_percent), gamma2_percent, 
                fill = factor(ifelse(gsp1_mutation == "WT","Wild Type", "Mutant")))) +
    geom_bar(color = "black", stat = "identity") +
    scale_fill_manual(name = "variant", values = c("grey50","red")) + 
    geom_text(aes(label = sprintf("%0.2f", round(gamma2_percent, digits = 2))), vjust=-1, size = 6) +
    xlab("Gsp1 mutant") +
    ylab("Percent in State 2 Conformation") +
    ylim(0,100) +
    theme(text = element_text(size=25),
    axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("hierarchical_clustering_of_parameters/NMR_gamma2_percents.pdf", height = 10, width = 15)
