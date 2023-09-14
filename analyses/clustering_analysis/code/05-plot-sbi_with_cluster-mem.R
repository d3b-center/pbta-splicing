################################################################################
# 05-plot-sbi_with_cluster-mem.R
# script that takes in "splicing_index.total.txt" data file from previous module
# and intersects with cluster membership info from current module output file 
# and plots number of samples per cluster stratified by splicing burden
# 
# written by Ammar Naqvi
#
# usage: Rscript 05-plot-sbi_with_cluster-mem.R
################################################################################

library("vroom")
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

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
sbi_cluster_plot <- file.path(analysis_dir, "plots", 
                              "cluster_by_sbi.tiff")

cc_members_file = "ccp_output/non_expr_pan_cancer_splice_subset_pam_pearson_0_ccp.rds"
cc_members_df  <- readRDS(paste0(output_dir,"/",cc_members_file) )
cc_members   <- cc_members_df[[3]]$consensusClass 
cc_members <- tibble::rownames_to_column(as.data.frame(cc_members), "Sample")
  
sbi_file <- "splicing_index.total.txt"
sbi_df <- vroom(paste0(input_dir,"/",sbi_file),delim = "\t")

SI_total_high     <- quantile(sbi_df$SI, probs=.75, names=FALSE)
SI_total_low      <- quantile(sbi_df$SI, probs=.25, names=FALSE)

splice_index_high <- filter(sbi_df, sbi_df$SI >SI_total_high )
splice_index_low  <- filter(sbi_df, sbi_df$SI <SI_total_low  )

## add column with "High or "Low" for SBI info
sbi_outliers <- sbi_df%>%filter(sbi_df$SI <SI_total_low | sbi_df$SI >SI_total_high) %>% 
  mutate(level=case_when(SI < SI_total_low ~ "Low",SI >SI_total_high  ~ "High" )) %>% inner_join(cc_members, by="Sample")

summ_sbi_cl <- sbi_outliers %>% group_by(cc_members, level) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n))

sbi_vs_cl_barplot <- ggplot(summ_sbi_cl, aes(x = level, y = n, fill = level)) + geom_bar(stat = "identity") +
                            facet_wrap( ~ cc_members, ncol = 3) + 
                            labs(x="SBI Level", y="Num samples") +  
                            coord_flip()  +
                            theme_Publication() + 
                            theme(legend.position = "none") +
                            scale_fill_manual(values= c("Low" = "#FFC20A", "High"="#0C7BDC" ))
                                        

print(sbi_vs_cl_barplot)

  
ggsave(
  sbi_cluster_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =12,
  height = 4,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)


