################################################################################
# 04-CLK1_PSI_plots.R
# written by Ammar Naqvi, Jo Lynne Rokita
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
  library("data.table")
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
hist_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")

## create plots dir if it doesn't exist
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## input files
rmats_file <- file.path(data_dir,"splice-events-rmats.tsv.gz")
clin_file  <- file.path(hist_dir,"histologies-plot-group.tsv")

## output files for final plots
hgg_plot_file <- file.path(plots_dir,"all_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf")
dmg_plot_file <- file.path(plots_dir,"dmg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf")
other_hgg_plot_file <- file.path(plots_dir,"other_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf")

## get CLK1 psi values in tumors and ctrls
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_df <- vroom(indep_file)

## load histologies info for HGG subty  
histologies_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

hgg_bs_id <- histologies_df %>%
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
  pull(Kids_First_Biospecimen_ID)

dmg_bs_id <- histologies_df %>%
  # Select only "RNA-Seq" samples
  filter(plot_group == "DIPG or DMG") %>%
  pull(Kids_First_Biospecimen_ID)

other_hgg_bs_id <- histologies_df %>%
  filter(plot_group == "Other high-grade glioma") %>%
  pull(Kids_First_Biospecimen_ID)

## load rmats input for CLK1
clk1_rmats <- fread(rmats_file) %>%
  # filter for CLK1 and exon 4, HGGs
  dplyr::filter(
                geneSymbol=="CLK1",
         exonStart_0base=="200860124", 
         exonEnd=="200860215",
         #FDR < 0.05, 
         #PValue < 0.05
         ) %>% 
  mutate("Inclusion"=IncLevel1,
         "Skipping"=1-IncLevel1) %>%  
  dplyr::select(sample_id, Inclusion, Skipping) %>%
  dplyr::rename(Kids_First_Biospecimen_ID=sample_id) %>% 
  ## reformat and rename to plot
  pivot_longer(!Kids_First_Biospecimen_ID, 
               names_to = "Type",
               values_to = "PSI") %>% 
  dplyr::select(Kids_First_Biospecimen_ID, Type, PSI) %>%
  # Join rmats data with clinical data
  inner_join(histologies_df, by='Kids_First_Biospecimen_ID') 

bs_list <- list("all_hgg" = hgg_bs_id, "dmg" = dmg_bs_id, "other_hgg" = other_hgg_bs_id)
names <- names(bs_list)

for (each in names) {
  # Assign correct plot file
  if (each == "all_hgg") {
    plot_file <- hgg_plot_file
  } else if (each == "dmg") {
    plot_file <- dmg_plot_file
  } else if (each == "other_hgg") {
    plot_file <- other_hgg_plot_file
  }
  
  # Filter the DataFrame based on current group's IDs
  new_df <- clk1_rmats %>%
    filter(Kids_First_Biospecimen_ID %in% bs_list[[each]])
  
  # Make a column for exposure of interest to order samples by
  samples_in_order <- new_df %>%
    dplyr::filter(Type == "Inclusion") %>%
    dplyr::select(PSI, Kids_First_Biospecimen_ID) %>%
    dplyr::arrange(-PSI) %>%
    dplyr::pull(Kids_First_Biospecimen_ID)
  

  # create df for plotting and reorder
  plot_df <- new_df %>%
    mutate(sample_id = factor(Kids_First_Biospecimen_ID,
                              levels = samples_in_order))
  plot_df_incl <- plot_df %>%
    filter(Type == "Inclusion")
  print(new_df)
  
  
  stacked_barplot <- ggplot(plot_df, aes(x = sample_id, y = PSI, fill= Type)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values=c("#FFC20A","#0C7BDC"), name="Exon 4") +
    theme_Publication() + 
    ylab(expression(bold(bolditalic("CLK1")~" Isoform Fraction"))) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    scale_y_continuous(expand = c(0, 0))  + # Set the expand argument to ensure the bottom line starts from 0
    geom_hline(yintercept = mean(plot_df_incl$PSI), color="black",linetype='dotted')
  
  # save mean CLK1 PSI for manuscript
  write_lines(mean(plot_df_incl$PSI), 
              file.path(results_dir, paste0(each,"-mean_clk1_psi.txt")) )
  
  # Save plot pdf
  pdf(plot_file, height = 3, width = 6)
  print(stacked_barplot)
  dev.off()
 
}
