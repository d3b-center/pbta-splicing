################################################################################
# 09-plot_highExon4_vs_lowExon4_and_SBI.R
# written by Ammar Naqvi
#
# usage: Rscript 09-plot_highExon4_vs_lowExon4_and_SBI.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  # Library for fast data loading
  library(vroom)
  # Library for data manipulation
  library(dplyr)
  # Library for plotting
  library(ggplot2)
  # Library for extra plotting functionality
  library(ggpubr)
  
})


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")
figures_dir <- file.path(root_dir, "figures")
input_dir   <- file.path(analysis_dir, "input")
# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(figures_dir, "theme_for_plots.R"))

## Load clinical file
# Specify clinical file path
clin_file  <- file.path(data_dir,"histologies.tsv")
# Load clinical file
clin_HGG_midlin_str_df  <-  vroom(clin_file, comment = "#",delim="\t") %>%
  # Select only "RNA-Seq" samples
  filter(experimental_strategy=="RNA-Seq",
         # Select only "HGAT" samples
         short_histology=="HGAT",
         # Select only "Midline" HGATs
         CNS_region=="Midline",
         # Select the "stranded" RNA library samples
         RNA_library=="stranded" 
         ) 

## Load rmats file
# Specify rmats file path
rmats_file <- file.path(data_dir, "rMATS_merged.single.SE.tsv.gz")
# Load rmats file
rmats_df <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>% 
  # Select exon 4
  filter(exonStart_0base=="200860124", exonEnd=="200860215") %>% 
  # Join rmats data with clinical data
  inner_join(clin_HGG_midlin_str_df, by=c('sample'='Kids_First_Biospecimen_ID')) %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  select(sample, geneSymbol, IncLevel1) 


## Load SBI file from previous module
# Specify SBI file path
sbi_file <-  file.path(input_dir,"splicing_index.total.txt")
# Load SBI file
sbi_df <-  vroom(sbi_file, comment = "#",delim="\t")
# Merge with rmats_df
sbi_vs_inclEx4_df <- inner_join(rmats_df, sbi_df, by=c("sample"="Sample")) 

## Compute quantiles to define high vs low Exon 4 PSI groups
# Compute quantiles
quartiles_psi <- quantile(sbi_vs_inclEx4_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)
# Calculate IQR
IQR_psi <- IQR(sbi_vs_inclEx4_df$IncLevel1)
# Get lower quantile (25%)
lower_psi <- quartiles_psi[1] 
# Get upper quantile (75%)
upper_psi <- quartiles_psi[2] 

# Create dataframe with low PSI values
sbi_vs_inclEx4_lowPSI_df <- sbi_vs_inclEx4_df %>%
  dplyr::filter(sbi_vs_inclEx4_df$IncLevel1 < lower_psi) %>%
  mutate(PSI="low") %>% 
  dplyr::select(sample, SI, PSI)

# Create dataframe with high PSI values
sbi_vs_inclEx4_highPSI_df <- sbi_vs_inclEx4_df %>% 
  dplyr::filter(sbi_vs_inclEx4_df$IncLevel1 > upper_psi) %>% 
  mutate(PSI="high") %>% 
  dplyr::select(sample,SI,PSI)

# Combine both high and low PSI dataframes
sbi_vs_inclEx4_by_extremePSI_df <- rbind(sbi_vs_inclEx4_lowPSI_df,sbi_vs_inclEx4_highPSI_df)


## Make box plot with stats
boxplot_sbi_vs_incl <- ggboxplot(sbi_vs_inclEx4_by_extremePSI_df, 
                                 # Specify x values
                                 x = "PSI",
                                 # Specify y values
                                 y = "SI",
                                 # Color in the box plot
                                 fill = "PSI",
                                 # Specify color palette
                                 palette = "jco", 
                                 # Add x-axis label
                                 xlab="Exon 4 PSI Level",
                                 # Add y-axis label
                                 ylab="Splicing Burden Index",
                                 # Add points
                                 add = "jitter") +
                       # Add p-value
                       stat_compare_means(method = "t.test") +
                       # Change theme
                       theme_Publication()

# View boxplot
print(boxplot_sbi_vs_incl)

# Path to save file as
plot_file <- file.path(plots_dir,"boxplot_high_vs_low_SBI.tiff")

# Save plot as tiff
tiff(plot_file, 
     res = 600, width = 6, height = 8, units = "in")

# Close the plotting device
dev.off()




