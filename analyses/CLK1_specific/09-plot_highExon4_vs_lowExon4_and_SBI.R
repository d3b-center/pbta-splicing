################################################################################
# 09-plot_highExon4_vs_lowExon4_and_SBI.R
# written by Ammar Naqvi
#
# usage: Rscript 09-plot_highExon4_vs_lowExon4_and_SBI.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("vroom")
  library("tidyverse")
  library("ggpubr")
  
})

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")
figures_dir <- file.path(root_dir, "figures")


input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# source function for theme for plots 
source(file.path(figures_dir, "theme_for_plots.R"))

##clinical file
clin_file  <- file.path(data_dir,"histologies.tsv")
clin_HGG_midlin_str_df  <-  vroom(clin_file, comment = "#",delim="\t") %>% filter(experimental_strategy=="RNA-Seq", short_histology=="HGAT",CNS_region=="Midline", RNA_library=="stranded" ) 

##rmats input
rmats_file <-  "/Users/naqvia/d3b_coding/pbta-splicing/data/rMATS_merged.single.SE.tsv.gz"
rmats_df <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  # select specific samples and extract CLK1 exon 4  
  dplyr::filter(geneSymbol=="CLK1") %>% dplyr::filter(exonStart_0base=="200860124", exonEnd=="200860215") %>% dplyr::select(sample, geneSymbol, IncLevel1) %>% 
  inner_join(clin_HGG_midlin_str_df, by=c('sample'='Kids_First_Biospecimen_ID')) %>%  dplyr::select(sample, geneSymbol, IncLevel1) 

## SBI file from previous module
sbi_file <-  file.path(input_dir,"splicing_index.total.txt")
sbi_df <-  vroom(sbi_file, comment = "#",delim="\t") 
sbi_vs_inclEx4_df <- inner_join(rmats_CLK1_ex4_df,sbi_df, by=c("sample"="Sample")) 

## compute quantiles to define high vs low Exon 4 PSI groups
quartiles_psi <- quantile(sbi_vs_inclEx4_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)
IQR_psi <- IQR(sbi_vs_inclEx4_df$IncLevel1)
lower_psi <- quartiles[1] 
upper_psi <- quartiles[2] 

## create dataframe with groupings
sbi_vs_inclEx4_lowPSI_df  <-dplyr::filter(sbi_vs_inclEx4_df, sbi_vs_inclEx4_df$IncLevel1 < lower_psi) %>% mutate(PSI="low") %>% dplyr::select(sample,SI,PSI)
sbi_vs_inclEx4_highPSI_df <-dplyr::filter(sbi_vs_inclEx4_df, sbi_vs_inclEx4_df$IncLevel1 > upper_psi) %>% mutate(PSI="high") %>% dplyr::select(sample,SI,PSI)
sbi_vs_inclEx4_by_extremePSI_df <- rbind(sbi_vs_inclEx4_lowPSI_df,sbi_vs_inclEx4_highPSI_df)

## box plot with stats
boxplot_sbi_vs_incl <- ggboxplot(sbi_vs_inclEx4_by_extremePSI_df, x = "PSI", y = "SI",
                                 color = "PSI", palette = "jco", xlab="Exon 4 PSI Level",ylab="Splicing Burden Index",
                                 add = "jitter")
# add p-value
boxplot_sbi_vs_incl + stat_compare_means()

#change method
boxplot_sbi_vs_incl + stat_compare_means(method = "t.test") + theme_Publication()

plot_file = file.path(plots_dir,"boxplot_high_vs_low_SBI.tiff")

# Save plot as tiff
tiff(plot_file, 
     res = 600, width = 6, height = 8, units = "in")
boxplot_sbi_vs_incl
dev.off()




