################################################################################
# 04-CLK1_PSI_plots.R
# written by Ammar Naqvi
#
# This script generates stacked barplot of splicing changes across 
# samples for CLK1 exon 4 
#
# usage: Rscript 04-CLK1_PSI_plots.R
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## create plots dir if it doesn't exist
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
CLK1_plot_path <- file.path(plots_dir, "CLK1_exon4_inclustion_fraction_hgg_stacked.pdf")

## get SRSF11 psi values in tumors and ctrls
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_df <- read_tsv(indep_file) %>% dplyr::filter(cohort=='PBTA')

rmats_file <- file.path(data_dir,"rMATS_merged.comparison.tsv.gz")
hist_file <- file.path(data_dir,"histologies.tsv")

## load histologies info
histologies_df <- vroom(hist_file) %>% 
  dplyr::filter(Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

## load rmats input for CLK1
clk1_rmats <- vroom(rmats_file, comment = "#", delim="\t") %>%
  # filter for CLK1 and exon 4
  dplyr::filter(geneSymbol=="CLK1",
         exonStart_0base=="200860124", 
         exonEnd=="200860215",
         FDR < 0.05, 
         PValue < 0.05) %>% 
  mutate("Inclusion"=IncLevel2,
         "Skipping"=1-IncLevel2) %>%  
  dplyr::select(sample_id, Inclusion, Skipping) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=sample_id) %>% 
  ## reformat and rename to plot
  pivot_longer(!Kids_First_Biospecimen_ID, 
               names_to = "Type",
               values_to = "PSI") %>% 
  dplyr::filter(Kids_First_Biospecimen_ID %in% histologies_df$Kids_First_Biospecimen_ID) %>% 
  right_join(histologies_df, by='Kids_First_Biospecimen_ID') %>% 
  dplyr::select(sample_id, Type, PSI) 
  
# Make a column for exposure of interest to order samples by
samples_in_order <- clk1_rmats %>%
  dplyr::filter(Type == "Inclusion") %>%
  dplyr::select(PSI, sample_id) %>%
  dplyr::arrange(-PSI) %>%
  dplyr::pull(sample_id)

# create df for plotting and reorder
plot_df <- clk1_rmats %>%
mutate(sample_id = factor(sample_id,
                          levels = samples_in_order)) %>% 
                          na.omit()

stacked_barplot <- ggplot(plot_df, aes(x = sample_id, y = PSI, fill= Type )) +
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=c("#FFC20A","#0C7BDC"), name="Exon 4") + 
  theme_Publication() + ylab(expression(bold(bolditalic("CLK1")~" Isoform Fraction"))) + xlab("Sample") + 
  theme(
      axis.text.x = element_text(angle = 75, hjust = 1)) 

# Save plot as pdf
pdf(CLK1_plot_path, height = 6, width = 14, useDingbats = FALSE)
print(stacked_barplot)
dev.off()
