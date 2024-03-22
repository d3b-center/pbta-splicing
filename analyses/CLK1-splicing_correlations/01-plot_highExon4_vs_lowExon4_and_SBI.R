################################################################################
# 01-plot_highExon4_vs_lowExon4_and_SBI.R
# written by Ammar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 01-plot_highExon4_vs_lowExon4_and_SBI.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(vroom)
  library(data.table)
  })


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
hist_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
figures_dir <- file.path(root_dir, "figures")
input_dir   <- file.path(root_dir, "analyses", "splicing_index", "results")


# Specify file paths
sbi_file <-  file.path(input_dir,"splicing_index.SE.txt")
clin_file  <- file.path(hist_dir,"histologies-plot-group.tsv")
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
  
# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(figures_dir, "theme_for_plots.R"))

# Path to save file as
hgg_plot_file <- file.path(plots_dir,"all_hgg_boxplot_SBI_high_vs_low_CLK1.pdf")
dmg_plot_file <- file.path(plots_dir,"dmg_boxplot_SBI_high_vs_low_CLK1.pdf")
other_hgg_plot_file <- file.path(plots_dir,"other_hgg_boxplot_SBI_high_vs_low_CLK1.pdf")

## Load clinical file
indep_df <- read_tsv(indep_file)

hist_rna_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

hgg_bs_id <- hist_rna_df %>%
  # Select only "RNA-Seq" samples
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
  pull(Kids_First_Biospecimen_ID)

dmg_bs_id <- hist_rna_df %>%
  # Select only "RNA-Seq" samples
  filter(plot_group == "DIPG or DMG") %>%
  pull(Kids_First_Biospecimen_ID)

other_hgg_bs_id <- hist_rna_df %>%
  filter(plot_group == "Other high-grade glioma") %>%
  pull(Kids_First_Biospecimen_ID)

## Load rmats file
rmats_df <-  vroom(rmats_file) %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>%
  # Select exon 4
  filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  select(sample_id, geneSymbol, IncLevel1) %>%
  # Join rmats data with clinical data
  inner_join(hist_rna_df, by=c('sample_id'='Kids_First_Biospecimen_ID')) 

## Load SBI file from previous module
sbi_vs_inclEx4_df <-  vroom(sbi_file) %>%
  inner_join(rmats_df, by=c("Sample"="sample_id"))

set.seed(45)

bs_list <- list("hgg" = hgg_bs_id, "dmg" = dmg_bs_id, "other" = other_hgg_bs_id)
names <- names(bs_list)

# Loop through groups
for (each in names) {
  # Assign correct plot file
  if (each == "hgg") {
    plot_file <- hgg_plot_file
  } else if (each == "dmg") {
    plot_file <- dmg_plot_file
  } else if (each == "other") {
    plot_file <- other_hgg_plot_file
  }
  
  # Filter the DataFrame based on current group's IDs
  new_df <- sbi_vs_inclEx4_df %>%
    filter(Sample %in% bs_list[[each]],
           RNA_library == "stranded")
  ## Compute quantiles to define high vs low Exon 4 PSI groups
  quartiles_psi <- quantile(new_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)
  # Calculate IQR
  IQR_psi <- IQR(new_df$IncLevel1)

  # Get lower quantile (25%)
  lower_psi <- quartiles_psi[1] 
  # Get upper quantile (75%)
  upper_psi <- quartiles_psi[2]

  # Create df with high/low PSI
  sbi_vs_inclEx4_by_extremePSI_df <- new_df %>%
    mutate(PSI = case_when(new_df$IncLevel1 > upper_psi ~ "high",
                           new_df$IncLevel1 < lower_psi ~ "low",
                         TRUE ~ NA_character_)
           ) %>%
    filter(!is.na(PSI)) %>%
    select(Sample,SI,PSI)

  # add a little extra room for wilcoxon test
  ylim_max <- max(sbi_vs_inclEx4_by_extremePSI_df$SI)+0.10*(max(sbi_vs_inclEx4_by_extremePSI_df$SI))

  ## Make box plot with stats
  boxplot_sbi_vs_incl <- ggboxplot(sbi_vs_inclEx4_by_extremePSI_df,x = "PSI", y = "SI") +
    xlab(expression(bold(bolditalic("CLK1")~"Exon 4 PSI Level"))) +
    ylab(expression(bold("Splicing Burden Index"))) +
    lims(y = c(0,ylim_max)) +
    ggforce::geom_sina(aes(color = PSI), size = 2, shape = 21, fill = NA, stroke = 1) +
    scale_color_manual(name = "PSI Level", values = c(high = "#FFC20A", low = "#0C7BDC")) +
    stat_compare_means(position = "identity", label.x = 1) +
    theme_Publication()

  # Save plot pdf
  pdf(plot_file, height = 4, width = 4)
  print(boxplot_sbi_vs_incl)
  dev.off()

}
