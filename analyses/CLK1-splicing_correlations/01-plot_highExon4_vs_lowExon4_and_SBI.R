################################################################################
# 01-plot_highExon4_vs_lowExon4_and_SBI.R
# written by Ammar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 01-plot_highExon4_vs_lowExon4_and_SBI.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(vroom)
  })


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_splicing_correlations")
figures_dir <- file.path(root_dir, "figures")
input_dir   <- file.path(root_dir, "analyses", "splicing_index", "results")


# Specify file paths
sbi_file <-  file.path(input_dir,"splicing_index.SE.txt")
clin_file  <- file.path(data_dir,"histologies.tsv")
rmats_file <- file.path(data_dir, "rMATS_merged.comparison.tsv.gz")

# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(figures_dir, "theme_for_plots.R"))

# Path to save file as
plot_file <- file.path(plots_dir,"boxplot_high_vs_low_SBI.pdf")

## Load clinical file
clin_HGG_midline_str_df  <-  read_tsv(clin_file) %>%
  # Select only "RNA-Seq" samples
  filter(experimental_strategy=="RNA-Seq",
         # Select only "HGAT" samples
         short_histology=="HGAT",
         # Select only "Midline" HGATs
         CNS_region %in% c("Midline", "Spine"),
         # Select the "stranded" RNA library samples
         #RNA_library=="stranded" 
         ) 

## Load rmats file
rmats_df <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>% 
  # Select exon 4
  filter(exonStart_0base=="200860124", exonEnd=="200860215") %>% 
  # Join rmats data with clinical data
  inner_join(clin_HGG_midline_str_df, by=c('sample_id'='Kids_First_Biospecimen_ID')) %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  select(sample_id, geneSymbol, IncLevel2) 

## Load SBI file from previous module
sbi_vs_inclEx4_df <-  read_tsv(sbi_file) %>% 
  inner_join(rmats_df, by=c("Sample"="sample_id")) 

## Compute quantiles to define high vs low Exon 4 PSI groups
quartiles_psi <- quantile(sbi_vs_inclEx4_df$IncLevel2, probs=c(.25, .75), na.rm = FALSE)
# Calculate IQR
IQR_psi <- IQR(sbi_vs_inclEx4_df$IncLevel2)
# Get lower quantile (25%)
lower_psi <- quartiles_psi[1] 
# Get upper quantile (75%)
upper_psi <- quartiles_psi[2] 

# Create df with high/low PSI
sbi_vs_inclEx4_by_extremePSI_df <- sbi_vs_inclEx4_df %>% 
  mutate(PSI = case_when(sbi_vs_inclEx4_df$IncLevel2 > upper_psi ~ "high",
                         sbi_vs_inclEx4_df$IncLevel2 < lower_psi ~ "low",
                         TRUE ~ NA_character_)
         ) %>% 
  filter(!is.na(PSI)) %>%
  select(Sample,SI,PSI)

# add a little extra room for wilcoxon test
ylim_max <- max(sbi_vs_inclEx4_by_extremePSI_df$SI)+0.10*(max(sbi_vs_inclEx4_by_extremePSI_df$SI))

## Make box plot with stats
set.seed(45)
boxplot_sbi_vs_incl <- ggboxplot(sbi_vs_inclEx4_by_extremePSI_df,x = "PSI", y = "SI") +
  xlab(expression(bold(bolditalic("CLK1")~"Exon 4 PSI Level"))) +
  ylab(expression(bold("Splicing Burden Index"))) +
  lims(y = c(0,ylim_max)) + 
  ggforce::geom_sina(aes(color = PSI), size = 2) +
  scale_color_manual(name = "PSI Level", values = c(high = "#FFC20A", low = "#0C7BDC")) +
  stat_compare_means(position = "identity", label.x = 1) +
  theme_Publication()


# Save plot tiff version
pdf(plot_file, height = 4, width = 4)
print(boxplot_sbi_vs_incl)
dev.off()

