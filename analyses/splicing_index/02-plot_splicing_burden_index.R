################################################################################
# 02-plot_splicing_burden_index.R
# script that takes in SBI tsv files and  computes relative proportion 
# (splicing burden index of aberrant splicing changes in samples
#
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 02-plot_splicing_burden_index.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("viridis")
  library("RColorBrewer")
  })

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
map_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
# source functions
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

# read in SI files
splice_index_SE_file   <- file.path(results_dir, "splicing_index.SE.txt")
splice_index_RI_file   <- file.path(results_dir, "splicing_index.RI.txt")
splice_index_A5SS_file <- file.path(results_dir, "splicing_index.A5SS.txt")
splice_index_A3SS_file <- file.path(results_dir, "splicing_index.A3SS.txt")

splice_index_SE_df   <- readr::read_tsv(splice_index_SE_file)
splice_index_RI_df   <- readr::read_tsv(splice_index_RI_file)
splice_index_A5SS_df <- readr::read_tsv(splice_index_A5SS_file)
splice_index_A3SS_df <- readr::read_tsv(splice_index_A3SS_file)

# read in color palette
palette_file <- file.path(map_dir, "histologies-plot-group.tsv")

palette_df <- read_tsv(palette_file) %>%
  dplyr::rename(Histology = plot_group) %>%
  select(Histology, plot_group_hex) %>%
  unique()

# Function to calculate medians and ranks
prepare_data_for_plot <- function(df, grouping_variable = NULL, min_samples = 5) {
  df %>%
    # Group by specified column
    group_by({{grouping_variable}}) %>%
    # Only keep groups with the specified minimum number of samples
    filter(n() > min_samples) %>%
    # Calculate group median
    mutate(
      group_median = median(SI, na.rm = TRUE),
      group_rank = rank(SI, ties.method = "first") / n(),
      sample_size = paste0("n = ", n())
    ) %>%
    ungroup() 
}

# create filenames for plots
file_si_SE_plot = "sbi-plot-SE.pdf"
file_si_RI_plot = "sbi-plot-RI.pdf"
file_si_A5SS_plot = "sbi-plot-A5SS.pdf"
file_si_A3SS_plot = "sbi-plot-A3SS.pdf"

plot_sbi <- function(sbi_df, plot_file) {
  
  si_cdf_plot <- sbi_df %>%
    as_tibble() %>%
    select(Sample, SI, Histology) %>%
    left_join(palette_df) %>%
    drop_na(Histology) %>%
    # Perform calculations needed for plot
    prepare_data_for_plot(grouping_variable = Histology) %>%
    # remove "Other" cancer group
   # filter(Histology != "Other") %>%
    # Order cancer groups by median TMB
    mutate(Histology = str_wrap(Histology, 18),
           Histology = fct_reorder(Histology, SI, .fun = median)
    ) 
  
  plot_colors <- si_cdf_plot$plot_group_hex
  names(plot_colors) <- si_cdf_plot$Histology
  
  
  p <- ggplot(si_cdf_plot) + 
    aes(
      x = group_rank,
      y = SI 
    ) 
  
  p <- p + geom_point(aes(fill=Histology, alpha = 0.8), 
               shape = 21, colour = "black", size = 3) + 
  
    # Add summary line for median
    geom_segment(
      x = 0, xend = 1, color = "black",
      aes(y = group_median, yend = group_median)
    ) +
    xlim(-1, 1.2) +
    scale_y_continuous(
      #trans = "log10",
      limits = c(0,0.6),
     # limits = c(0, 400),
      #breaks = c(0, 3, 10, 30, 100, 300)
    ) +
    
    # add labels
    labs(x = expression(bold("Histology")),
         y = expression(bold("Splicing Burden Index (SBI)"))
         ) +

    ggpubr::theme_pubr() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 12, angle = 90, hjust = 1),
      strip.background = element_rect(fill = NA, color = NA),
      legend.position = "none"
    ) +
    
    facet_wrap(~Histology + sample_size, nrow = 1, strip.position = "bottom")  +
    scale_fill_manual(values = plot_colors) 

  # Save plot
  ggsave(filename = plot_file, path = plots_dir, plot = p,
         height = 6, width = 10, useDingbats = FALSE)
}

## plot SBI for each splicing case
plot_sbi(splice_index_SE_df,file_si_SE_plot)
plot_sbi(splice_index_RI_df,file_si_RI_plot)
plot_sbi(splice_index_A5SS_df,file_si_A5SS_plot)
plot_sbi(splice_index_A3SS_df,file_si_A3SS_plot)

