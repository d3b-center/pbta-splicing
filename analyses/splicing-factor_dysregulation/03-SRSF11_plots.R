################################################################################
# 03-SRSF11_plots.R
# written by Ammar Naqvi
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


psi_tab  <-  read_tsv(psi_tumors_file, col_names = c("Kids_First_Biospecimen_ID","psi","type")) %>%
  mutate(sample = fct_reorder(Kids_First_Biospecimen_ID, psi))

psi_ctrl_tab  <-  read_tsv(psi_ctrl_file, col_names = c("Kids_First_Biospecimen_ID","psi","type")) %>%
  filter(grepl("control-BS", Kids_First_Biospecimen_ID))

## intersect just midline HGGs
clin_df  <-  read_tsv(clin_file, guess_max = 100000)

psi_with_midline_filter <- clin_df %>% 
  # Only include RNA-Seq samples 
  filter(experimental_strategy == "RNA-Seq",
         cohort == "PBTA",
         # Include samples that are "stranded"
         RNA_library == "stranded",
         # Include those samples originated from the Midline
         CNS_region == "Midline",
                # Include only HGAT sample histology
         short_histology == "HGAT") %>% 
  inner_join(psi_tab, by="Kids_First_Biospecimen_ID") %>% 
  dplyr::select("Kids_First_Biospecimen_ID","psi","type") %>% 
  distinct()

psi_with_midline_filter_and_ctrl <- rbind(psi_ctrl_tab, psi_with_midline_filter)
## stacked barplot 

plot_df <- psi_with_midline_filter_and_ctrl %>%
  mutate(type = case_when(Kids_First_Biospecimen_ID == "control-BS" & type == "inclusion" ~ "control inclusion",
                          Kids_First_Biospecimen_ID == "control-BS" & type == "skip" ~ "control skipping",
                          Kids_First_Biospecimen_ID != "control-BS" & type == "inclusion" ~ "tumor inclusion",
                          Kids_First_Biospecimen_ID != "control-BS" & type == "skip" ~ "tumor skipping",
                                 TRUE ~ type))
        
stacked_barplot_SRSF11 <- ggplot(plot_df,
                                 aes(x = Kids_First_Biospecimen_ID, y = psi, fill= type)) +
  geom_bar(position="stack", stat="identity")    + 
  scale_fill_manual(values=c("black","darkgray", "#FFC20A","#0C7BDC"))      + 
  theme_Publication() +  
  ylab("Isoform %") + 
  xlab("Sample") + 
  theme(axis.text.x=element_blank(), 
        #text=element_text(size=20),
        axis.ticks.x=element_blank())

stacked_barplot_SRSF11

# Save plot as pdf
pdf(file_SRSF11_plot, height = 3, width = 6, useDingbats = FALSE)
print(stacked_barplot_SRSF11)
dev.off()
