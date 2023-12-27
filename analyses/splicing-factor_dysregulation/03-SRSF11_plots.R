################################################################################
# 03-SRSF11_plots.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# This script uses SRSF11 specific input files to compute and plot correlations 
# between PSI vs expression and also makes barplot of splicing changes across 
# samples. 
#
# usage: Rscript 03-SRSF11_plots.R
################################################################################

suppressPackageStartupMessages({
  library("reshape2")
  library("tidyverse")
  library("ggpubr")
  library("ggplot2")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

##theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_SRSF11_plot <- file.path(analysis_dir, "plots", "SRSF11_hgg_stacked.pdf")
file_SRSF11_corr_plot <- file.path(analysis_dir, "plots", "corr_rna_vs_psi_SRSF11.pdf")

## get SRSF11 psi values in tumors and ctrls
psi_tumors_file <- file.path(input_dir, "se.srsf11.incl.txt")
psi_ctrl_file <- file.path(input_dir, "se.srsf11.ctrl.incl.txt")
clin_file <- file.path(data_dir,"histologies.tsv")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary-plus.tsv")

psi_tab  <-  read_tsv(psi_tumors_file, col_names = c("Kids_First_Biospecimen_ID","psi","type"))
psi_ctrl_tab  <-  read_tsv(psi_ctrl_file, col_names = c("Kids_First_Biospecimen_ID","psi","type"))
indep_df <- read_tsv(indep_file)

## intersect just midline HGGs, independent specimens
clin_df  <-  read_tsv(clin_file, guess_max = 100000) %>% 
  # Only include RNA-Seq samples 
  filter(experimental_strategy == "RNA-Seq",
         cohort == "PBTA",
         # Include samples that are "stranded"
         RNA_library == "stranded",
         # Include those samples originated from the Midline
         CNS_region == "Midline",
                # Include only HGAT sample histology
         short_histology == "HGAT",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID) 

psi_midline_hgg_and_ctrl <- psi_tab %>%
  bind_rows(psi_ctrl_tab) %>%
  # only midline hgg
  filter(Kids_First_Biospecimen_ID %in% c(clin_df$Kids_First_Biospecimen_ID, "control-BS")) %>%
  left_join(clin_df[,c("Kids_First_Biospecimen_ID", "sample_id")]) %>%
  mutate(sample_id = case_when(Kids_First_Biospecimen_ID == "control-BS" ~ "brainstem",
                               TRUE ~ sample_id))

# Make a column for exposure of interest to order samples by
samples_in_order <- psi_midline_hgg_and_ctrl %>%
  filter(type == "inclusion") %>%
  select(psi, sample_id) %>%
  arrange(-psi) %>%
  pull(sample_id)

## stacked barplot 
plot_df <- psi_midline_hgg_and_ctrl %>%
  mutate(type = case_when(Kids_First_Biospecimen_ID == "control-BS" & type == "inclusion" ~ "control inclusion",
                          Kids_First_Biospecimen_ID == "control-BS" & type == "skip" ~ "control skipping",
                          Kids_First_Biospecimen_ID != "control-BS" & type == "inclusion" ~ "tumor inclusion",
                          Kids_First_Biospecimen_ID != "control-BS" & type == "skip" ~ "tumor skipping",
                                 TRUE ~ type))
# reorder
plot_df_fct <- plot_df %>%
  mutate(sample_id = factor(sample_id,
                            levels = samples_in_order))

        
stacked_barplot_SRSF11 <- ggplot(plot_df_fct,
                                 aes(x = sample_id, y = psi, fill= type)) +
  geom_bar(position="stack", stat="identity", colour="black")    + 
  scale_fill_manual(values=c("lightyellow","lightblue", "#FFC20A","#0C7BDC"),
                    name = "Splice event")      + 
  theme_Publication() +  
  ylab(expression(bold(bolditalic("SRSF11")~" Isoform Fraction"))) + 
  xlab("Sample") + 
    
  theme(axis.text.x = element_text(angle = 75, hjust = 1))

# Save plot as pdf
pdf(file_SRSF11_plot, height = 4, width = 8, useDingbats = FALSE)
print(stacked_barplot_SRSF11)
dev.off()
