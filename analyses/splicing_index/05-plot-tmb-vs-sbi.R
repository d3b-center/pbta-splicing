# 05-plot-tmb-sbi.R
# written by Ammar Naqvi & Jo Lynne Rokita
#
# This script plots TMB based on high vs low splicing burden 
# usage: Rscript 05-plot-tmb-sbi.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("ggplot2")
  library("ggpubr")
} )

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
map_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
input_dir   <- file.path(analysis_dir, "input")
data_dir   <- file.path(root_dir, "data")
results_dir   <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## Call plot publication theme script 
source(file.path(root_dir, "figures", "theme_for_plots.R"))

## define output files
boxplot_sbi_vs_tmb_by_mutation_status_file  <- file.path(plots_dir, "boxplot_sbi-tmb-by-mutation-status.pdf")
boxplot_sbi_vs_tmb_by_cg_file  <- file.path(plots_dir, "boxplot_sbi-tmb-by-cg.pdf")
boxplot_sbi_vs_tmb_by_mutation_status_dmg_file <- file.path(plots_dir, "boxplot_sbi-tmb-by-dmg.pdf")
boxplot_sbi_vs_tmb_by_mutation_status_dipg_file <- file.path(plots_dir, "boxplot_sbi-tmb-by-dipg.pdf")

## input files
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_wgs_file <- file.path(data_dir, "independent-specimens.wgswxspanel.primary.prefer.wgs.tsv")
tmb_coding_file  <- file.path(input_dir,"snv-mutation-tmb-coding.tsv") # OPC v13 TMB file
sbi_coding_file  <- file.path(results_dir,"splicing_index.SE.txt")

# read in files
hist_pal <- read_tsv(palette_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cancer_group)

indep_rna_df <- read_tsv(indep_rna_file) %>% 
  filter(cohort == "PBTA") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_RNA = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_RNA)

indep_wgs_df <- read_tsv(indep_wgs_file) %>% 
  filter(cohort == "PBTA") %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_DNA)

## get tmb file (source: OpenPedCan v13)
tmb_coding_df  <-  read_tsv(tmb_coding_file)  %>% 
  mutate(mutation_status = case_when(tmb <10 ~ "Normal",
                                     tmb >=10 & tmb < 100 ~ "Hypermutant",
                                     tmb >=100 ~ "Ultra-hypermutant")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Tumor_Sample_Barcode) %>%
  dplyr::select(-experimental_strategy) %>%
  right_join(indep_wgs_df)

sbi_coding_df  <-  read_tsv(sbi_coding_file) %>% 
  dplyr::rename(Kids_First_Biospecimen_ID_RNA = Sample) %>%
  #dplyr::select(-Histology) %>%
  right_join(indep_rna_df)

## intersect tmb values with SBI tumors
sbi_vs_tmb_innerjoin_df <- tmb_coding_df %>%
  inner_join(sbi_coding_df, by="Kids_First_Participant_ID") %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA) #%>%
  # add histology
  #left_join(hist_pal)

## identify samples by high vs low SBI tumors
quartiles_sbi <- quantile(sbi_vs_tmb_innerjoin_df$SI, probs=c(.25, .75), na.rm = TRUE)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

high_vs_low_df <- sbi_vs_tmb_innerjoin_df %>%
  mutate(SBI_level = case_when(SI > upper_sbi ~ "High",
                               SI < lower_sbi ~ "Low")) %>%
  # remove samples without either TMB or SI data
  filter(!is.na(SBI_level),
         !is.na(mutation_status))

# relevel for plotting
high_vs_low_df$mutation_status <- factor(high_vs_low_df$mutation_status, levels = c("Normal", "Hypermutant", "Ultra-hypermutant"))
high_vs_low_df$SBI_level <- factor(high_vs_low_df$SBI_level, levels = c("Low", "High"))

sbi_tmb_plot <- ggplot(high_vs_low_df, aes(SBI_level, log10(tmb))) +  
  ggforce::geom_sina(aes(color = SBI_level, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  facet_wrap("mutation_status") +
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + 
  labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  ylim(c(-3,3)) +
  theme(legend.position="none")

## save pdf
pdf(boxplot_sbi_vs_tmb_by_mutation_status_file, width = 7, height = 4)
sbi_tmb_plot
dev.off()

## asses DIMG only
sbi_tmb_dmg <- sbi_vs_tmb_innerjoin_df %>%  filter(Histology=='Diffuse midline glioma')

## identify samples by high vs low SBI tumors
quartiles_sbi <- quantile(sbi_tmb_dmg$SI, probs=c(.25, .75), na.rm = TRUE)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

high_vs_low_df <- sbi_tmb_dmg %>%
  mutate(SBI_level = case_when(SI > upper_sbi ~ "High",
                               SI < lower_sbi ~ "Low")) %>%
  # remove samples without either TMB or SI data
  filter(!is.na(SBI_level),
         !is.na(mutation_status))

# relevel for plotting
high_vs_low_df$mutation_status <- factor(high_vs_low_df$mutation_status, levels = c("Normal", "Hypermutant", "Ultra-hypermutant"))
high_vs_low_df$SBI_level <- factor(high_vs_low_df$SBI_level, levels = c("Low", "High"))

sbi_tmb_dmgplot <- ggplot(high_vs_low_df, aes(SBI_level, log10(tmb))) +  
  ggforce::geom_sina(aes(color = SBI_level, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  stat_compare_means(label.y = 2.9) + 
  facet_wrap("mutation_status") +
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + 
  labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  ylim(c(-3,3)) +
  theme(legend.position="none")

## save pdf
pdf(boxplot_sbi_vs_tmb_by_mutation_status_dmg_file, width = 7, height = 4)
sbi_tmb_dmgplot
dev.off()

## asses DIPG only
sbi_tmb_dipg <- sbi_vs_tmb_innerjoin_df %>%  filter(cancer_group=='Diffuse intrinsic pontine glioma')

## identify samples by high vs low SBI tumors
quartiles_sbi <- quantile(sbi_tmb_dipg$SI, probs=c(.25, .75), na.rm = TRUE)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

high_vs_low_df <- sbi_tmb_dipg %>%
  mutate(SBI_level = case_when(SI > upper_sbi ~ "High",
                               SI < lower_sbi ~ "Low")) %>%
  # remove samples without either TMB or SI data
  filter(!is.na(SBI_level),
         !is.na(mutation_status))

# relevel for plotting
high_vs_low_df$mutation_status <- factor(high_vs_low_df$mutation_status, levels = c("Normal", "Hypermutant", "Ultra-hypermutant"))
high_vs_low_df$SBI_level <- factor(high_vs_low_df$SBI_level, levels = c("Low", "High"))

sbi_tmb_dipg_plot <- ggplot(high_vs_low_df, aes(SBI_level, log10(tmb))) +  
  ggforce::geom_sina(aes(color = SBI_level, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  stat_compare_means(label.y = 2.9) + 
  facet_wrap("mutation_status") +
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + 
  labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  ylim(c(-3,3)) +
  theme(legend.position="none")

## save pdf
pdf(boxplot_sbi_vs_tmb_by_mutation_status_dipg_file, width = 7, height = 4)
sbi_tmb_dipg_plot
dev.off()


