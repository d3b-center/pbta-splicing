# 06-plot_tmb_by_sbi.R
# written by Ammar Naqvi
#
# This script plots TMB based on high vs low splicing burden 
#
# usage: Rscript 06-plot_tmb_by_sbi.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
  library("ggplot2")
  library("ggpubr")
} )

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses","splicing_index")
input_dir   <- file.path(analysis_dir, "input")
data_dir   <- file.path(root_dir, "data")
results_dir   <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## define output files
boxplot_sbi_vs_tmb_file  <- file.path(plots_dir,"boxplot_tmb_vs_sbi-level.tiff")


## get and setup input
tmb_coding_file  <- file.path(input_dir,"tmb-coding.tsv")
tmb_coding_df  <-  vroom(tmb_coding_file, comment = "#",delim="\t")  %>% mutate(Kids_First_Biospecimen_ID=Tumor_Sample_Barcode)

sbi_coding_file  <- file.path(results_dir,"splicing_index.total.txt")
sbi_coding_df  <-  vroom(sbi_coding_file, comment = "#",delim="\t") %>% mutate(Kids_First_Biospecimen_ID=Sample) %>% filter(Histology=="HGG")

clin_file = file.path(data_dir, "histologies.tsv")
clin_df  <-  vroom(clin_file, comment = "#",delim="\t") 
sbi_ids_clin <- clin_df %>% inner_join(sbi_coding_df, by="Kids_First_Biospecimen_ID") %>% filter(tumor_descriptor=="Initial CNS Tumor")
tmb_ids_clin <- clin_df %>% inner_join(tmb_coding_df, by="Kids_First_Biospecimen_ID") %>% filter(tumor_descriptor=="Initial CNS Tumor")

## intersect tmb values with SBI tumors
sbi_vs_tmb_innerjoin_df <- left_join(sbi_ids_clin,tmb_ids_clin, by="Kids_First_Participant_ID") 

## compute quantiles to define high vs low Exon 4 SBI tumors
quartiles_sbi <- quantile(sbi_vs_tmb_innerjoin_df$SI, probs=c(.25, .75), na.rm = FALSE)
IQR_sbi <- IQR(sbi_vs_tmb_innerjoin_df$SI)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

## subset tmb values and samples by high vs low SBI tumors
high_sbi_df <- dplyr::filter(sbi_vs_tmb_innerjoin_df, SI > upper_sbi) %>% dplyr::mutate(SBI_level="high")
low_sbi_df  <- dplyr::filter(sbi_vs_tmb_innerjoin_df, SI < lower_sbi) %>% dplyr::mutate(SBI_level="low")
high_vs_low_df <- rbind(low_sbi_df,high_sbi_df)

boxplot_grp_sbi_tmb <- ggplot(high_vs_low_df,aes(SBI_level,log10(tmb))) + 
  geom_boxplot(aes(fill=SBI_level)) + 
  stat_compare_means(method = "t.test") + 
  xlab("Splicing Burden Level") + 
  geom_jitter() + theme_Publication() 

# save plot tiff version
tiff(boxplot_sbi_vs_tmb_file, height =2000, width = 1500, res = 300)
print(boxplot_grp_sbi_tmb)
dev.off()


