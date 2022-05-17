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

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}



##theme for all plots
theme_Publication <- function(base_size=15, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            # legend.margin = unit(0.5, "cm"),
            legend.margin = margin(5,5,5,5),
            legend.title = element_text(face="bold"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            plot.margin=unit(c(10,5,5,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}


## newer version with just tumors
file <- "splicing_index.total.txt"
splice_index <- readr::read_tsv(file.path(results_dir, file))

splice_index <- splice_index %>%
  as.data.frame(stringsAsFactors = FALSE)

# Set up the data.frame for plotting
si_cdf <- splice_index %>%

  # We only really need these two variables from data.frame
  dplyr::transmute(
    group = Histology,
    number = log10(as.numeric(SI*100))
  ) %>%

  # Group by specified column
  dplyr::group_by(group) %>%

  # Only keep groups with the specified minimum number of samples
  dplyr::filter(dplyr::n() > 1) %>%

  # Calculate group median
  dplyr::mutate(
    group_median = median(number, na.rm = TRUE),
    group_rank = rank(number, ties.method = "first") / dplyr::n(),
    sample_size = paste0("n = ", dplyr::n())
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = reorder(group, group_median))

si_cdf %>%
  # Now we will plot these as cumulative distribution plots
  ggplot2::ggplot(ggplot2::aes(
    x = group_rank,
    y = number
  )) +

  ggplot2::geom_point(color = "black") +

  # Add summary line for median
  ggplot2::geom_segment(
    x = 0, xend = 1, color = "blue",
    ggplot2::aes(y = group_median, yend = group_median)
  ) +

  # Separate by histology
  ggplot2::facet_wrap(~ group + sample_size, nrow = 1, strip.position = "bottom") +
  ggplot2::theme_classic() +
  ggplot2::xlab("Histology") +
  ggplot2::ylab("Splicing Burden Index") +

  # Making it pretty
  #ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.text = ggplot2::element_text(size = 14, angle = 90, hjust = .5),
    strip.background = ggplot2::element_rect(fill = NA, color = NA)
  ) + theme_Publication() 

file_si_plot = "/SI_total.png"
filename = paste0(plots_dir, file_si_plot)
ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =6000,
  height = 2000,
  units = "px",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

## SBI plot for HGATs only
file <- "splicing_index.total.hgg_clusters.txt"
splice_index <- readr::read_tsv(file.path(results_dir, file))

splice_index <- splice_index %>%
  as.data.frame(stringsAsFactors = FALSE)

# Set up the data.frame for plotting
si_cdf <- splice_index %>%
  
  # We only really need these two variables from data.frame
  dplyr::transmute(
    group = Histology,
    number = log10(as.numeric(SI*100))
  ) %>%
  
  # Group by specified column
  dplyr::group_by(group) %>%
  
  # Only keep groups with the specified minimum number of samples
  dplyr::filter(dplyr::n() > 1) %>%
  
  # Calculate group median
  dplyr::mutate(
    group_median = median(number, na.rm = TRUE),
    group_rank = rank(number, ties.method = "first") / dplyr::n(),
    sample_size = paste0("n = ", dplyr::n())
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = reorder(group, group_median))

si_cdf %>%
  # Now we will plot these as cumulative distribution plots
  ggplot2::ggplot(ggplot2::aes(
    x = group_rank,
    y = number
  )) +
  
  ggplot2::geom_point(color = "black") +
  
  # Add summary line for median
  ggplot2::geom_segment(
    x = 0, xend = 1, color = "blue",
    ggplot2::aes(y = group_median, yend = group_median)
  ) +
  
  # Separate by histology
  ggplot2::facet_wrap(~ group + sample_size, nrow = 1, strip.position = "bottom") +
  ggplot2::theme_classic() +
  ggplot2::xlab("Histology") +
  ggplot2::ylab("Splicing Burden Index") +
  
  # Making it pretty
  #ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.text = ggplot2::element_text(size = 14, angle = 90, hjust = .5),
    strip.background = ggplot2::element_rect(fill = NA, color = NA)
  ) + theme_Publication() 

file_si_plot = "/SI_total.hgat.png"
filename = paste0(plots_dir, file_si_plot)
ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =6000,
  height = 2000,
  units = "px",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)
