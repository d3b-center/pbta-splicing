################################################################################
# survival_analysis.R
# script that takes 
# written by Ammar Naqvi
#
# usage: Rscript survival_analysis.R
################################################################################
library(dplyr)


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

file = "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/cc_bs_ids.txt"
cluster_mem  = read.delim(file, sep = "\t", header=TRUE)



cluster_not6 <- filter(cluster_mem, cluster_id %in% c(1,2,3,4,5))
cluster_not5 <- filter(cluster_mem, cluster_id %in% c(1,2,3,4,6))
cluster_not4 <- filter(cluster_mem, cluster_id %in% c(1,2,3,6,5))
cluster_not3 <- filter(cluster_mem, cluster_id %in% c(1,2,4,6,5))
cluster_not2 <- filter(cluster_mem, cluster_id %in% c(1,3,4,6,5))
cluster_not1 <- filter(cluster_mem, cluster_id %in% c(3,2,4,6,5))

cluster_1 <- filter(cluster_mem, cluster_id %in% 1)
cluster_2 <- filter(cluster_mem, cluster_id %in% 2)
cluster_3 <- filter(cluster_mem, cluster_id %in% 3)
cluster_4 <- filter(cluster_mem, cluster_id %in% 4)
cluster_5 <- filter(cluster_mem, cluster_id %in% 5)
cluster_6 <- filter(cluster_mem, cluster_id %in% 6)

# source function to compute survival
source(file.path("util", "survival_models.R"))
`%>%` <- dplyr::`%>%`

# set up directories
setwd("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results")
input_dir <- "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer" # you can set this to the directory which contains this code


data_dir  <- file.path("/Users/naqvia/Desktop/AS-DMG/data") # contains pbta-histologies file and all input .dat files
plots_dir <- file.path("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/plots")  # output directory for plot

# create output directories
dir.create(plots_dir, recursive = T, showWarnings = F)

# read pbta-histology file
metadata <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"))
#metadata <- metadata %>%
#  dplyr::filter(broad_histology == "Diffuse astrocytic and oligodendroglial tumor",
#                composition == "Solid Tissue")

# read in given list of files
file_list <- list.files(path = data_dir, pattern = ".dat", recursive = T, full.names = T)

splicing_types <- splicing_types %>%
  dplyr::left_join(metadata, by = "Kids_First_Biospecimen_ID")


# loop through each file and add metadata, generate survival plot, add to splots which is a list of plots
splots <- list()
for(i in 1:length(file_list)){
  print(i)
  splicing_types <- readr::read_tsv(file.path(file_list[i]))
  splicing_types <- splicing_types %>%
    dplyr::left_join(metadata, by = "Kids_First_Biospecimen_ID")
  
  # Patients with multiple different molecular subtype might affect downstream analysis for survival so using random selection of 1-1 sample-patient matches
  unique_match <- splicing_types %>%
    dplyr::group_by(Kids_First_Participant_ID) %>%
    dplyr::summarize(Kids_First_Biospecimen_ID = sample(Kids_First_Biospecimen_ID, 1))
  
  # select only independent primary WGS bs_ids to get overall survival
  metadata_sub <- metadata %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% unique_match$Kids_First_Biospecimen_ID) %>%
    dplyr::left_join(splicing_types[,c("sample_id","type")], suffix = c("_DNA","_RNA")) %>%
    dplyr::filter(!is.na(type))
  
  # Kaplan-Meier for all HGG/DMG subtypes
  kap_fit <- survival_analysis(metadata_sub,
                               ind_var = "type",
                               test = "kap.meier",
                               metadata_sample_col = "Kids_First_Biospecimen_ID")
  
  surv_pvalue(kap_fit$model, data = kap_fit$original_data)
  spval<- surv_pvalue(kap_fit$model, data = kap_fit$original_data)
  if(spval$pval > 0.05){
    next
  }
  
  sign_findings=file.path(file_list[i])
  write(sign_findings,file="/Users/naqvia/Desktop/DIPG/surv_analysis/sign_findings.txt",append=TRUE)
  
  # make plots (I removed the bottom table and x-axis limits)
  surv_plot <- survminer::ggsurvplot(kap_fit$model,
                                     pval = TRUE,
                                     data = kap_fit$original_data,
                                     risk.table = TRUE,
                                     break.time.by = 500,
                                     ggtheme = theme_minimal(),
                                     risk.table.y.text.col = TRUE,
                                     risk.table.y.text = FALSE)
  
  surv_plot$plot <- surv_plot$plot +
    ggtitle(gsub('.*/|_[0-9].*','',file_list[i])) +
    theme(legend.position = "right")
  
  # add to splots list
  splots[[i]] <- surv_plot
}

splots_mod = splots[-which(sapply(splots, is.null))]
res <- arrange_ggsurvplots(splots_mod, print = FALSE, ncol = 4, nrow = 31)

# output file name
kap_meier_plot_file <- file.path(plots_dir, "survival_curve.pdf")
ggsave(filename = kap_meier_plot_file, plot = res, device = "pdf", height = 150, width = 40, limitsize = F)

##theme for all plots
theme_Publication <- function(base_size=15, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            # legend.margin = unit(0.5, "cm"),
            legend.margin = margin(5,5,5,5),
            legend.title = element_text(face="bold"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            plot.margin=unit(c(10,5,5,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}
                    
