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
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses","splicing_index")
input_dir   <- file.path(analysis_dir, "input")
results_dir   <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## define output files
plot_path  <- file.path(plots_dir,"hist_by_sbi-level_lolliplot.pdf")

## get and setup input
sbi_coding_file  <- file.path(results_dir,"splicing_index.SE.txt")
sbi_coding_df  <-  vroom(sbi_coding_file, comment = "#",delim="\t") %>% 
  dplyr::rename(Kids_First_Biospecimen_ID=Sample) 

clin_file = file.path(results_dir,"histologies-plot-group.tsv") 
histologies_df <- vroom(clin_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID, plot_group) %>% 
  dplyr::rename(Histology = plot_group)

sbi_coding_df %>% inner_join(sbi_coding_df, by="Kids_First_Biospecimen_ID")

clin_df_w_SBI  <-  vroom(clin_file, comment = "#",delim="\t",show_col_types = FALSE) %>% 
  inner_join(sbi_coding_df, by="Kids_First_Biospecimen_ID") %>%
  dplyr::select(Kids_First_Biospecimen_ID, SI, plot_group) %>% 
  dplyr::rename(Histology = plot_group)

## compute quantiles to define high vs low SBI tumors
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

histology_counts_by_sbi <- clin_df_w_SBI %>% dplyr::count(SBI_level, Histology) %>%
  dplyr::mutate(SBI_level = fct_relevel(SBI_level,
                                        c("Low SBI", "High SBI")))


## create color palette 
palette_file <- file.path(results_dir, "histologies-plot-group.tsv")

## make plot colors list
palette_df <- read_tsv(palette_file) %>%
  dplyr::rename(Histology = plot_group) %>%
  select(Histology, plot_group_hex) %>%
  unique()

plot_colors <- list()
plot_colors[['Histology']] <- palette_df$plot_group_hex

plot_labels <- list()
plot_labels[['Histology']]  <- palette_df$Histology
#names(plot_labels) <- palette_df$Histology

lolliplot_plot <- ggplot(histology_counts_by_sbi, aes(x=Histology, y=n, color=Histology)) +
  geom_segment(aes(x=Histology, xend=Histology, y=0, yend=n)) +
  scale_color_manual(name = "Histology",values = plot_colors[['Histology']]) + 
  geom_point(size=4) +
  xlab("") + ylab("# of Samples") +
  facet_wrap(SBI_level ~ .) +
  theme_Publication() + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=12), 
        axis.title=element_text(size=14)) 


# save plot pdf version
pdf(plot_path, height = 4, width = 14)
print(lolliplot_plot)
dev.off()