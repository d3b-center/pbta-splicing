################################################################################
# 06-identify_and_plot_histologies_by_SBI.R
# written by Ammar Naqvi
#
# script identifies high vs low SBI tumors and plots their histology makeup
#
# usage: Rscript 06-identify_and_plot_histologies_by_SBI.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
  library("ggplot2")
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

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## get and setup input
sbi_coding_file  <- file.path(results_dir,"splicing_index.total.txt")
sbi_coding_df  <-  vroom(sbi_coding_file, comment = "#",delim="\t") %>% mutate(Kids_First_Biospecimen_ID=Sample) 

clin_file = file.path(data_dir,"histologies.tsv")
clin_df_w_SBI  <-  vroom(clin_file, comment = "#",delim="\t") %>% inner_join(sbi_coding_df, by="Kids_First_Biospecimen_ID")


## compute quantiles to define high vs low Exon 4 SBI tumors
quartiles_sbi <- quantile(sbi_coding_df$SI, probs=c(.25, .75), na.rm = FALSE)
IQR_sbi <- IQR(sbi_coding_df$SI)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

mutate(failreason = case_when(P1memory == "FALSE" ~ "memory", 
                              P1choice == "FALSE" ~ "choice",
                              P1future == "FALSE" ~ "justification"))

clin_df_w_SBI <- clin_df_w_SBI %>% 
                 mutate(SBI_level = case_when(SI > upper_sbi ~ "high", 
                                              SI < lower_sbi ~ "low")) %>% 
                 filter(!is.na(SBI_level))

clin_df_w_highSBI <- clin_df_w_SBI %>% filter(SBI_level=="high")
clin_df_w_lowSBI <- clin_df_w_SBI %>% filter(SBI_level=="low")

histology_counts_by_sbi <- clin_df_w_SBI %>% count(SBI_level, short_histology)
  
## define output files
piechart_hist_by_sbi_file  <- file.path(plots_dir,"piechart_hist_by_sbi-level.tiff")

# create color palette for short histology
palette_dir   <- file.path(root_dir, "palettes")
palette_df <- file.path(palette_dir, "short_histology_color_palette.tsv") %>% read_tsv()

short_histology_palettes <- palette_df %>%
  dplyr::select(short_histology, hex_code) %>%
  unique()

mycolors <- list()
mycolors[['short_histology']] <- short_histology_palettes$hex_code
names(mycolors[['short_histology']]) <- short_histology_palettes$short_histology

# make piechart
piechart_hist_by_sbi <- ggplot(data = histology_counts_by_sbi, aes(x = "", y = n, fill = short_histology )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  facet_wrap(~ SBI_level, ncol=1)  +
  scale_fill_manual(name = "short_histology",values = mycolors[['short_histology']]) +
  xlab("") + ylab("Numer of samples") + 
  theme_Publication()

# save plot tiff version
tiff(piechart_hist_by_sbi_file, height = 3000, width = 2000, res = 300)
print(piechart_hist_by_sbi)
dev.off()  



