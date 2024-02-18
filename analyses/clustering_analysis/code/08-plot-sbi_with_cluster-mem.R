################################################################################
# 05-plot-sbi_with_cluster-mem.R
# script that takes in "splicing_index.total.txt" data file from previous module
# and intersects with cluster membership info from current module output file 
# and plots number of samples per cluster stratified by splicing burden
# 
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 05-plot-sbi_with_cluster-mem.R
################################################################################
library("tidyverse")

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")

plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")
output_dir <- file.path(analysis_dir, "output")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# define number of clusters
n_clusters <- 13

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
sbi_cluster_plot <- file.path(analysis_dir, "plots", 
                              "cluster_by_sbi.pdf")

cc_members_file = "ccp_output/non_expr_pan_cancer_splice_subset_pam_canberra_0_ccp.rds"
cc_members_df  <- readRDS(paste0(output_dir,"/",cc_members_file) )
cc_members   <- cc_members_df[[n_clusters]]$consensusClass 
cc_members <- tibble::rownames_to_column(as.data.frame(cc_members), "Sample")
  
sbi_file <- file.path(input_dir, "splicing_index.total.txt")
sbi_df <- read_tsv(sbi_file)

SI_total_high     <- quantile(sbi_df$SI, probs=.75, names=FALSE)
SI_total_low      <- quantile(sbi_df$SI, probs=.25, names=FALSE)

splice_index_high <- filter(sbi_df, sbi_df$SI >SI_total_high)
splice_index_low  <- filter(sbi_df, sbi_df$SI <SI_total_low)

## add column with "High or "Low" for SBI info
sbi_outliers <- sbi_df %>%
  filter(sbi_df$SI <SI_total_low | sbi_df$SI >SI_total_high) %>% 
  mutate(level=case_when(SI < SI_total_low ~ "Low",SI >SI_total_high  ~ "High" )) %>% 
  inner_join(cc_members, by="Sample")

summ_sbi_cl <- sbi_outliers %>% group_by(cc_members, level) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>%
  ungroup()

# are there any clusters without high or low SBI members? YES
clusters <- as.list(1:n_clusters)

missing_clusters <- setdiff(clusters, unique(as.character(summ_sbi_cl$cc_members)))

# Ensure missing_cluster_df has the correct structure with the specified columns
missing_cluster_df <- data.frame(cc_members = integer(),
                                 level = character(),
                                 n = integer(),
                                 Freq = numeric())

# Loop through each missing cluster to add a new row for it
for(missing in missing_clusters) {
  missing_cluster_df <- rbind(missing_cluster_df, 
                              data.frame(cc_members = missing,
                                         level = "High",
                                         n = NA,
                                         Freq = NA))
}

summ_sbi_cl <- summ_sbi_cl %>%
  bind_rows(missing_cluster_df)

# plot
sbi_vs_cl_barplot <- ggplot(summ_sbi_cl, aes(x = as.factor(cc_members), y = n, fill = level)) + 
  geom_bar(stat = "identity") +
  labs(x="Cluster", y="Number of Samples") +
  theme_Publication() + 
  scale_fill_manual(values= c("Low" = "#FFC20A", "High"="#0C7BDC"), name = "SBI Level") +

pdf(sbi_cluster_plot, height = 3, width = 8)
print(sbi_vs_cl_barplot)
dev.off()

