################################################################################
# 02-plot_splicing_burden_index.R
# script that takes in SBI tsv files and  computes relative proportion 
# (splicing burden index of aberrant splicing changes in samples
#
# written by Ammar Naqvi
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

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

file_si_SE_plot = "SBI-plot.SE.pdf"
file_si_RI_plot = "SBI-plot.RI.pdf"
file_si_A5SS_plot = "SBI-plot.A5SS.pdf"
file_si_A3SS_plot = "SBI-plot.A3SS.pdf"

# theme for all plots
# source functions
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

# read in SI file
splice_index_SE_file   <- file.path(results_dir, "splicing_index.SE.txt")
splice_index_RI_file   <- file.path(results_dir, "splicing_index.RI.txt")
splice_index_A5SS_file <- file.path(results_dir, "splicing_index.A5SS.txt")
splice_index_A3SS_file <- file.path(results_dir, "splicing_index.A3SS.txt")

splice_index_SE_df   <- readr::read_tsv(splice_index_SE_file)
splice_index_RI_df   <- readr::read_tsv(splice_index_RI_file)
splice_index_A5SS_df <- readr::read_tsv(splice_index_A5SS_file)
splice_index_A3SS_df <- readr::read_tsv(splice_index_A3SS_file)


plot_sbi <- function(sbi_df, plot_file) {
  
  # Set up the data.frame for plotting
  si_cdf_plot <- sbi_df %>%
    # Group by specified column
    dplyr::group_by(Histology) %>%
    # Only keep groups with the specified minimum number of samples
    dplyr::filter(dplyr::n() > 1) %>%
    # Calculate group median
    dplyr::mutate(
      group_mean = mean(SI, na.rm = TRUE),
      group_rank = rank(SI, ties.method = "first") / dplyr::n(),
      sample_size = paste0("n = ", dplyr::n())
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Histology = reorder(Histology, group_mean))

  si_plot <- si_cdf_plot %>%
    # Now we will plot these as cumulative distribution plots
    ggplot2::ggplot(ggplot2::aes(
      x = group_rank,
      y = SI 
    )) +
    ggplot2::geom_point(color = "red", alpha = 0.7, size = 3) +

    # Add summary line for median
    ggplot2::geom_segment(
      x = .3, xend = .7, color = "black", linetype = 2, size = 1,
      ggplot2::aes(y = group_mean, yend = group_mean)
    ) +
    scale_y_continuous(
      #trans = "log1p",
      limits = c(0,0.6),
    ) +

    # Separate by histology
    ggplot2::facet_wrap(~ Histology + sample_size, nrow = 1, 
                        strip.position = "bottom", labeller = ggplot2::label_wrap_gen(multi_line = TRUE)) +
    ggplot2::xlab("Histology") +
    ggplot2::ylab("Splicing Burden Index") +

    # Making it pretty
    ggplot2::theme(legend.position = "none") +
    ggpubr::theme_pubr() +
    theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 10, angle = 90, hjust = 1),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      legend.position = "none"
    ) 

  # Save plot
  ggsave(filename = plot_file, path = plots_dir, plot = si_plot,
         height = 10, width = 40)
}

## plot SBI for each splicing case
splice_index_SE_df   <- readr::read_tsv(splice_index_SE_file) %>% as.data.frame(stringsAsFactors = FALSE)
splice_index_RI_df   <- readr::read_tsv(splice_index_RI_file) %>% as.data.frame(stringsAsFactors = FALSE)
splice_index_A5SS_df <- readr::read_tsv(splice_index_A5SS_file) %>% as.data.frame(stringsAsFactors = FALSE)
splice_index_A3SS_df <- readr::read_tsv(splice_index_A3SS_file) %>% as.data.frame(stringsAsFactors = FALSE)

plot_sbi(splice_index_SE_df,file_si_SE_plot)
plot_sbi(splice_index_RI_df,file_si_RI_plot)
plot_sbi(splice_index_A5SS_df,file_si_A5SS_plot)
plot_sbi(splice_index_A3SS_df,file_si_A3SS_plot)
