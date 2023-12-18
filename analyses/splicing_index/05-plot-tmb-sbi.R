# 05-plot-tmb-sbi.R
# written by Ammar Naqvi
#
# This script plots TMB based on high vs low splicing burden 
# usage: Rscript 05-plot-tmb-sbi.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("tidyverse")
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

## Call plot publication theme script 
source(file.path(root_dir, "figures", "theme_for_plots.R"))

## define output files
boxplot_sbi_vs_tmb_file  <- file.path(plots_dir,"boxplot_sbi-tmb.pdf")
boxplot_tmb_vs_sbi_file  <- file.path(plots_dir,"boxplot_tmb-sbi.pdf")

## input file
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_wgs_file <- file.path(data_dir, "independent-specimens.wgswxspanel.primary.prefer.wgs.tsv")
tmb_coding_file  <- file.path(input_dir,"tmb-coding.tsv")
sbi_coding_file  <- file.path(results_dir,"splicing_index.SE.txt")

indep_rna_df <- vroom(indep_rna_file) %>% filter(cohort == 'PBTA')
indep_wgs_df <- vroom(indep_wgs_file) %>% filter(cohort == 'PBTA')

## filter for samples that have both RNA and WGS
indep_sample_wgs_rna <- inner_join(indep_wgs_df,indep_rna_df, by='Kids_First_Participant_ID')

## get tmb file (source: open-pedcan)
tmb_coding_df  <-  vroom(tmb_coding_file)  %>% 
  mutate(Kids_First_Biospecimen_ID=Tumor_Sample_Barcode)

sbi_coding_df  <-  vroom(sbi_coding_file) %>% 
  mutate(Kids_First_Biospecimen_ID=Sample) 

## subset histology file for indep/primary samples 
clin_file = file.path(data_dir, "histologies.tsv")
clin_df  <-  vroom(clin_file) %>% filter(cohort == 'PBTA') %>% 
  inner_join(indep_sample_wgs_rna, by='Kids_First_Participant_ID')

sbi_ids_clin <- clin_df %>% inner_join(sbi_coding_df, by="Kids_First_Biospecimen_ID") 
tmb_ids_clin <- clin_df %>% inner_join(tmb_coding_df, by="Kids_First_Biospecimen_ID")

## intersect tmb values with SBI tumors
sbi_vs_tmb_innerjoin_df <- inner_join(sbi_ids_clin,tmb_ids_clin, by="Kids_First_Participant_ID") 

## identify samples by high vs low SBI tumors
quartiles_sbi <- quantile(sbi_vs_tmb_innerjoin_df$SI, probs=c(.25, .75), na.rm = FALSE)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]
high_sbi_df <- dplyr::filter(sbi_vs_tmb_innerjoin_df, SI > upper_sbi) %>% dplyr::mutate(SBI_level="High")
low_sbi_df  <- dplyr::filter(sbi_vs_tmb_innerjoin_df, SI < lower_sbi) %>% dplyr::mutate(SBI_level="Low")
high_vs_low_df <- rbind(low_sbi_df,high_sbi_df)

## identify samples by high vs low TMB tumors
quartiles_tmb <- quantile(sbi_vs_tmb_innerjoin_df$tmb, probs=c(.25, .75), na.rm = FALSE)
lower_tmb <- quartiles_tmb[1]
upper_tmb <- quartiles_tmb[2]
high_tmb_df <- dplyr::filter(sbi_vs_tmb_innerjoin_df, tmb > upper_tmb) %>% dplyr::mutate(TMB_level="High")
low_tmb_df  <- dplyr::filter(sbi_vs_tmb_innerjoin_df, tmb < lower_tmb) %>% dplyr::mutate(TMB_level="Low")
high_vs_low_TMB_df <- rbind(low_tmb_df,high_tmb_df)


## remove hyper/ultra-mutant samples
high_vs_low_hyper_rem_df <- high_vs_low_df %>% dplyr::filter(tmb < 10)
ggplot(high_vs_low_hyper_rem_df,aes(SBI_level,log10(tmb)) ) +  
  geom_violin(aes(SBI_level)) +
  ggforce::geom_sina(aes(color = SBI_level), size = 2,method="density") +
  stat_compare_means() + 
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  theme(legend.position="none")

## look at samples with low TMB
high_vs_low_lowTMB_df <- high_vs_low_df %>% dplyr::filter(tmb <= lower_tmb)
sbi_tmb_plot <- ggplot(high_vs_low_lowTMB_df,aes(SBI_level,log10(tmb)) ) +  
  geom_violin(aes(SBI_level)) +
  ggforce::geom_sina(aes(color = SBI_level), size = 2,method="density") +
  stat_compare_means() + 
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  theme(legend.position="none")

## look at samples with high TMB
high_vs_low_highTMB_df <- high_vs_low_df %>% dplyr::filter(tmb >= upper_tmb)
sbi_tmb_plot <- ggplot(high_vs_low_highTMB_df,aes(SBI_level,log10(tmb)) ) +  
  geom_violin(aes(SBI_level)) +
  ggforce::geom_sina(aes(color = SBI_level), size = 2,method="density") +
  stat_compare_means() + 
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  theme(legend.position="none")

## generate two boxplots (SBI and TMB)
sbi_tmb_plot <- ggplot(high_vs_low_df,aes(SBI_level,log10(tmb)) ) +  
  geom_violin(aes(SBI_level)) +
  ggforce::geom_sina(aes(color = SBI_level), size = 2,method="density") +
  stat_compare_means() + 
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  theme(legend.position="none")

## save pdf
pdf(boxplot_sbi_vs_tmb_file, 
    width = 5, height = 5)
sbi_tmb_plot
dev.off()

tmb_sbi_plot <- ggplot(high_vs_low_TMB_df,aes(TMB_level,SI) ) +  
  geom_violin(aes(TMB_level)) +
  ggforce::geom_sina(aes(color = TMB_level), size = 2,method="density") +
  stat_compare_means() + 
  scale_color_manual(name = "TMB_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + labs(y="Splicing Burden Index", x="Tumor Mutation Burden Level") + 
  theme(legend.position="none")

## save pdf
pdf(boxplot_tmb_vs_sbi_file, 
    width = 5, height = 5)
tmb_sbi_plot
dev.off()

