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

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

splice_index_file <- file.path(results_dir, "splicing_index.total.txt")
splice_index_df <- readr::read_tsv(splice_index_file)

splice_index <- splice_index_df %>%
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

## survival based on high vs low SBI
file <- "splicing_index.total.txt"
splice_index <- read.delim(paste0(results_dir,"/", file), sep = "\t", header=TRUE) %>% rename(Kids_First_Biospecimen_ID = Sample)

SI_total_high     <- quantile(splice_index$SI, probs=.75, names=FALSE)
SI_total_low      <- quantile(splice_index$SI, probs=.25, names=FALSE)

splice_index_high <- filter(splice_index, splice_index$SI >SI_total_high )
splice_index_low  <- filter(splice_index, splice_index$SI <SI_total_low | splice_index$SI >SI_total_high  )

## add column with "High or "Low" for SBI info
splicing_index_outliers <- splice_index%>%filter(splice_index$SI <SI_total_low | splice_index$SI >SI_total_high) %>% 
  mutate(level=case_when(SI < SI_total_low ~ "Low",SI >SI_total_high  ~ "High" ))

# source function to compute survival
util_dir <- file.path(root_dir, "util")

source(file.path(util_dir, "survival_models.R"))
`%>%` <- dplyr::`%>%`

# read in clinical file
clin_tab <- readr::read_tsv(file.path(data_dir, "histologies.tsv")) %>% 
  dplyr::left_join(splicing_index_outliers, by = "Kids_First_Biospecimen_ID")

# Kaplan-Meier for all clusters
kap_fit <- survival_analysis(clin_tab,
                             ind_var = "level",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)


file_surv_plot = "/surv_si.png"
filename = paste0(plots_dir, file_surv_plot)
ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  #width =2800,
  #height = 1800,
  units = "px",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)
