################################################################################
# 03-identify_and_plot_histologies_by_SBI.R
# written by Ammar S Naqvi, Jo Lynne Rokita
#
# script identifies high vs low SBI tumors and creates piechart of their hist 
#
# usage: Rscript 03-identify_and_plot_histologies_by_SBI.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses","splicing_index")
map_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
input_dir   <- file.path(analysis_dir, "input")
results_dir   <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))
if(!interactive()) pdf(NULL)

## define output files
barplot_path <- file.path(plots_dir, "hist_by_sbi_level_barplot.pdf")

## get and setup input files
sbi_coding_file  <- file.path(results_dir,"splicing_index.SE.txt")
palette_file <- file.path(map_dir,"histologies-plot-group.tsv") 

# read in files, join palette with sbi file
sbi_coding_df  <-  read_tsv(sbi_coding_file, comment = "#") %>% 
  dplyr::rename(Kids_First_Biospecimen_ID = Sample) 



palette_df <- read_tsv(palette_file) %>%
  dplyr::rename(Histology = plot_group) %>%
  dplyr::select(Histology, plot_group_hex) %>%
  unique() %>% 
  add_row(Histology="Diffuse midline glioma",plot_group_hex="#ff40d9") %>% 
  add_row(Histology="Diffuse intrinsic pontine glioma",plot_group_hex="#ffccf5") 


plot_colors <- palette_df$plot_group_hex
names(plot_colors) <- palette_df$Histology
plot_colors <- list(plot_colors)


sbi_coding_df <- sbi_coding_df %>%
  left_join(palette_df)

## compute quantiles to define high vs low SBI tumors
quartiles_sbi <- quantile(sbi_coding_df$SI, probs=c(.25, .75), na.rm = FALSE)
IQR_sbi <- IQR(sbi_coding_df$SI)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

clin_df_w_SBI <- sbi_coding_df %>% 
  mutate(SBI_level = case_when(SI > upper_sbi ~ "High SBI", 
                               SI < lower_sbi ~ "Low SBI")) %>% 
  filter(!is.na(SBI_level))

clin_df_w_highSBI <- clin_df_w_SBI %>% 
  filter(SBI_level=="High SBI") %>%
  dplyr::count(SBI_level, Histology) %>%
  dplyr::rename(High_SBI = n) %>%
  dplyr::select(-SBI_level)

clin_df_w_lowSBI <- clin_df_w_SBI %>% 
  filter(SBI_level=="Low SBI") %>%
  dplyr::count(SBI_level, Histology) %>%
  dplyr::rename(Low_SBI = n) %>%
  dplyr::select(-SBI_level)

plot_df <- clin_df_w_highSBI %>%
  full_join(clin_df_w_lowSBI) %>%
  left_join(palette_df) %>%
  mutate(Histology = fct_reorder(Histology, High_SBI)) %>%
  dplyr::rename(`Samples with High SBI (N)` = High_SBI,
                `Samples with Low SBI (N)` = Low_SBI)

# make colors
plot_colors <- palette_df$plot_group_hex
names(plot_colors) <- palette_df$Histology

g.mid <- ggplot(plot_df,aes(x=1,y=Histology)) +
  geom_text(aes(label=Histology), size = 4.25) +
  ylab(NULL)+
  xlab(NULL) +
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.06))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(2,1,10,2), "mm")) 

g1 <- ggplot(data = plot_df, aes(x = Histology, y = `Samples with Low SBI (N)`)) +
  geom_bar(stat = "identity", aes(fill = Histology), show.legend = FALSE) + 
  scale_fill_manual(values = plot_colors) +
  theme_Publication() +
  theme(#axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    plot.margin = unit(c(2,1,1,1), "mm"),
    axis.line.y.left = element_blank(),
    panel.grid.major.y = element_blank(),  # Remove major grid lines
    panel.grid.minor.y = element_blank())  + # Remove minor grid lines
  scale_y_reverse(limits = c(200,0)) + 
  coord_flip() +
  #labs(x = "Number of samples with Low SBI") +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", linewidth = 0.75)

g2 <- ggplot(data = plot_df, aes(x = Histology, y = `Samples with High SBI (N)`)) +
  geom_bar(stat = "identity", aes(fill = Histology), show.legend = FALSE) +
  scale_fill_manual(values = plot_colors) +
  theme_Publication() +
  theme(#axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(1,1,1,2), "mm"),
    axis.line.y.left = element_blank(),
    panel.grid.major.y = element_blank(),  # Remove major grid lines
    panel.grid.minor.y = element_blank())  + # Remove minor grid lines) +
  # labs(x = "Number of samples with High SBI") +
  scale_y_continuous(limits = c(0,200)) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", linewidth = 0.75)

gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

# save barplot 
pdf(barplot_path, height = 5, width = 11)
grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(3/8,2/8,3/8))
dev.off()
