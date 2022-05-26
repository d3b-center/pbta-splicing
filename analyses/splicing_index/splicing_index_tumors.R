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

## SI with survival in HGGs
file <- "splicing_index.total.hgg_clusters.surv.txt"
splice_index_surv = read.delim(paste0(dataDir, file), sep = "\t", header=TRUE)

# source function to compute survival
source(file.path("/Users/naqvia/Desktop/pbta-splicing/analyses/pan_cancer/util", "survival_models.R"))
`%>%` <- dplyr::`%>%`

# read pbta-histology file
metadata <- readr::read_tsv("/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/pbta-histologies.tsv")

splice_index_surv <- splice_index_surv %>%
  dplyr::left_join(metadata, by = "Kids_First_Biospecimen_ID")

# Kaplan-Meier for all clusters
kap_fit <- survival_analysis(splice_index_surv,
                             ind_var = "Level",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

surv_plot <- survminer::ggsurvplot(kap_fit$model,
                                   pval = TRUE,
                                   data = kap_fit$original_data,
                                   risk.table = FALSE,
                                   break.time.by = 500,
                                   ggtheme = theme_Publication(),
                                   risk.table.y.text.col = TRUE,
                                   risk.table.y.text = FALSE)

file_surv_plot = "surv_si.hgat.png"
filename = paste0(plotDir, file_surv_plot)
ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =2800,
  height = 1800,
  units = "px",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

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

## survival based on high vs low SBI
file <- "splicing_index.total.txt"
splice_index  <-  read.delim(paste0(resultsDir, file), sep = "\t", header=TRUE) %>% rename(Kids_First_Biospecimen_ID = Sample)


SI_total_high     <- quantile(splice_index$SI, probs=.75, names=FALSE)
SI_total_low      <- quantile(splice_index$SI, probs=.25, names=FALSE)

splice_index_high <- filter(splice_index, splice_index$SI >SI_total_high )
splice_index_low  <- filter(splice_index, splice_index$SI <SI_total_low | splice_index$SI >SI_total_high  )

## add column with "High or "Low" for SBI info
splicing_index_outliers <- filter(splice_index, splice_index$SI <SI_total_low | splice_index$SI >SI_total_high  ) %>% 
  mutate(level=case_when(splicing_index_outliers$SI < SI_total_low ~ "Low",
                         splicing_index_outliers$SI >SI_total_high  ~ "High" ))


splicing_index_outliers_HGAT <- filter(splice_index, Histology=="HGAT") 
SI_total_high_HGAT     <- quantile(splicing_index_outliers_HGAT$SI, probs=.75, names=FALSE)
SI_total_low_HGAT      <- quantile(splicing_index_outliers_HGAT$SI, probs=.25, names=FALSE)

splicing_index_outliers_HGAT <-  filter(splicing_index_outliers_HGAT, SI <SI_total_low_HGAT | SI >SI_total_high_HGAT  ) %>% 
  mutate(level=case_when(SI < SI_total_low_HGAT ~ "Low",
                         SI >SI_total_high_HGAT  ~ "High" ))

# source function to compute survival
util_dir <- file.path(root_dir, "util")

source(file.path(util_dir, "survival_models.R"))
`%>%` <- dplyr::`%>%`

# read in clinical file
clin_tab <- readr::read_tsv(file.path(data_dir, "v19_plus_20210311_pnoc_rna.tsv")) %>% 
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
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

## focus on HGAT
splicing_index_outliers_HGAT <- filter(splice_index, Histology=="HGAT") 
SI_total_high_HGAT     <- quantile(splicing_index_outliers_HGAT$SI, probs=.75, names=FALSE)
SI_total_low_HGAT      <- quantile(splicing_index_outliers_HGAT$SI, probs=.25, names=FALSE)

splicing_index_outliers_HGAT <-  filter(splicing_index_outliers_HGAT, SI <SI_total_low_HGAT | SI >SI_total_high_HGAT  ) %>% 
  mutate(level=case_when(SI < SI_total_low_HGAT ~ "Low",
                         SI >SI_total_high_HGAT  ~ "High" ))

# source function to compute survival
source(file.path("/Users/naqvia/Desktop/pbta-splicing/analyses/pan_cancer/util", "survival_models.R"))
`%>%` <- dplyr::`%>%`

# read in clinical file
clin_tab <- readr::read_tsv(file.path(data_dir, "v19_plus_20210311_pnoc_rna.tsv")) %>% 
  dplyr::left_join(splicing_index_outliers_HGAT, by = "Kids_First_Biospecimen_ID")


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
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

