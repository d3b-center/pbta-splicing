################################################################################
# splicing_index_tumors.R
# script that takes in "splicing_index.total.txt" data file
# and computes relative proportion (splicing burden index of aberrant splicing 
# changes in samples
#
# written by Ammar Naqvi
#
# usage: Rscript splicing_index_tumors.R
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

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

file_si_plot = "SI_total.tiff"

# theme for all plots
# source functions
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

# read in SI file
splice_index_file <- file.path(results_dir, "splicing_index.total.txt")
splice_index_df <- readr::read_tsv(splice_index_file)


splice_index <- splice_index_df %>%
  as.data.frame(stringsAsFactors = FALSE)

# Set up the data.frame for plotting
si_cdf_plot <- splice_index %>%
  # Group by specified column
  dplyr::group_by(Histology) %>%
  # Only keep groups with the specified minimum number of samples
  dplyr::filter(dplyr::n() > 1) %>%
  # Calculate group median
  dplyr::mutate(
    group_median = median(SI, na.rm = TRUE),
    group_rank = rank(SI, ties.method = "first") / dplyr::n(),
    sample_size = paste0("n = ", dplyr::n())
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Histology = reorder(Histology, group_median)) 

si_plot <- si_cdf_plot %>%
  # Now we will plot these as cumulative distribution plots
  ggplot2::ggplot(ggplot2::aes(
    x = group_rank,
    y = log(SI+1, 10)
  )) +
  
  ggplot2::geom_point(color = "black", alpha =0.7, shape = 1) +
  
  # Add summary line for median
  ggplot2::geom_segment(
    x = 0, xend = 1, color = "red",linetype=2,
    ggplot2::aes(y = log(group_median+1,10), yend = log(group_median+1,10))
  ) +
  
  # Separate by histology
  ggplot2::facet_wrap(~ Histology + sample_size, nrow = 1, strip.position = "bottom", labeller = ggplot2::label_wrap_gen(multi_line = FALSE)) +
  ggplot2::xlab("Histology") +
  ggplot2::ylab("Splicing Burden Index - log10(SI+1)") +
  
  # Making it pretty
  ggplot2::theme(legend.position = "none") +
  theme_Publication() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.placement = "outside"
  )  

dev.set(dev.next())

# Save plot
tiff(file.path(plots_dir, file_si_plot), height = 1200, width = 3000, res = 300)
print(si_plot)
dev.off()
