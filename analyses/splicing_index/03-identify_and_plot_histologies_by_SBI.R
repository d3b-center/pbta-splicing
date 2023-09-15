################################################################################
# 03-identify_and_plot_histologies_by_SBI.R
# written by Ammar S Naqvi
#
# script identifies high vs low SBI tumors and creates piechart of their hist 
#
# usage: Rscript 03-identify_and_plot_histologies_by_SBI.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("tidyverse")
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

## define output files
piechart_hist_by_sbi_file  <- file.path(plots_dir,"piechart_hist_by_sbi-level.tiff")

## get and setup input
sbi_coding_file  <- file.path(results_dir,"splicing_index.SE.txt")
sbi_coding_df  <-  vroom(sbi_coding_file, comment = "#",delim="\t") %>% mutate(Kids_First_Biospecimen_ID=Sample) 

clin_file = file.path(data_dir,"histologies.tsv")
clin_df_w_SBI  <-  vroom(clin_file, comment = "#",delim="\t",show_col_types = FALSE) %>% inner_join(sbi_coding_df, by="Kids_First_Biospecimen_ID")

## compute quantiles to define high vs low Exon 4 SBI tumors
quartiles_sbi <- quantile(sbi_coding_df$SI, probs=c(.25, .75), na.rm = FALSE)
IQR_sbi <- IQR(sbi_coding_df$SI)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]


clin_df_w_SBI <- clin_df_w_SBI %>% 
  mutate(SBI_level = case_when(SI > upper_sbi ~ "High SBI", 
                               SI < lower_sbi ~ "Low SBI")) %>% 
  filter(!is.na(SBI_level))

clin_df_w_highSBI <- clin_df_w_SBI %>% filter(SBI_level=="High SBI")
clin_df_w_lowSBI <- clin_df_w_SBI %>% filter(SBI_level=="Low SBI")

histology_counts_by_sbi <- clin_df_w_SBI %>% dplyr::count(SBI_level, short_histology) %>%
  dplyr::mutate(SBI_level = fct_relevel(SBI_level,
                                        c("Low SBI", "High SBI")))


# create color palette for short histology
palette_dir   <- file.path(root_dir, "palettes")
palette_df <- file.path(palette_dir, "short_histology_color_palette.tsv") %>% read_tsv()

short_histology_palettes <- palette_df %>%
  dplyr::select(short_histology, hex_code) %>%
  unique()

short_histology_palettes <- palette_df %>%
  dplyr::select(short_histology,plot_group_display, hex_code) %>%
  unique()

mycolors <- list()
mycolors[['short_histology']] <- short_histology_palettes$hex_code
names(mycolors[['short_histology']]) <- short_histology_palettes$short_histology

plot_labels <- list()
plot_labels[['short_histology']] <- short_histology_palettes$plot_group_display
names(plot_labels[['short_histology']]) <- short_histology_palettes$short_histology

# make piechart
piechart_hist_by_sbi <- ggplot(data = histology_counts_by_sbi, aes(x = "", y = n, fill = short_histology )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5),size=3) +
  coord_polar(theta = "y",start=0) +
  facet_wrap(~ SBI_level, ncol=1)  +
  scale_fill_manual(name = "Histology",values = mycolors[['short_histology']], labels=plot_labels[['short_histology']]) +
  xlab("") + ylab("Number of Samples") +
  theme_Publication() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8), 
        axis.title=element_text(size=10))

# save plot tiff version
tiff(piechart_hist_by_sbi_file, height = 1500, width = 1850, res = 300)
print(piechart_hist_by_sbi)
dev.off()  