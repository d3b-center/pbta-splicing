################################################################################
# 05-plot-flip_mixed_events.R
# script that plots a table of flip-like splicing events across functional sites
# written by Ammar Naqvi
#
# usage: Rscript 05-plot-flip_mixed_events.R 
################################################################################

## libraries needed
suppressPackageStartupMessages({
  library("ggstatsplot")
  library("dplyr")
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_events_functional_sites")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

##theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output file names for plots
file_flip_events_plot <- file.path(analysis_dir, "plots", "flip_barplots.pdf")

## retrieve psi values from tables
file_psi_pos_total <- "splicing_events.total.pos.tsv"
file_psi_neg_total <- "splicing_events.total.neg.tsv"

file_psi_pos_func <- "splicing_events.total.pos.intersectUnip.ggplot.txt"
file_psi_neg_func <- "splicing_events.total.neg.intersectUnip.ggplot.txt"

psi_pos_func_tab <-  read.delim(file.path(results_dir, file_psi_pos_func), sep = "\t", row.names = NULL, header=TRUE)
psi_pos_tab      <-  read.delim(file.path(results_dir, file_psi_pos_total), sep = "\t", row.names = NULL, header=TRUE)  %>% filter(flip == 1)

psi_neg_func_tab <-  read.delim(file.path(results_dir, file_psi_neg_func), sep = "\t", row.names = NULL, header=TRUE)
psi_neg_tab      <-  read.delim(file.path(results_dir, file_psi_neg_total), sep = "\t", row.names = NULL, header=TRUE)  %>% filter(flip == 1)


## num of sites that are undergo flip at functional sites
flip_event_func_skipping  <- psi_pos_tab[psi_pos_func_tab$SpliceID %in% psi_pos_tab$gene, ]
flip_event_func_inclusion <- psi_neg_tab[psi_neg_func_tab$SpliceID %in% psi_neg_tab$gene, ]

##calculate number of totals and percentages 
## get skipping counts
num_flip_skip_func <- c(nrow(flip_event_func_skipping))
num_skip_func      <- c(nrow(psi_pos_func_tab))
num_flip_skip_func_perc <- (num_flip_skip_func/(num_skip_func)) * 100
num_non_flip_skip_func_perc      <- ( (num_skip_func-num_flip_skip_func)/(num_skip_func)) * 100

## getting inclusion counts
num_flip_incl_func <- c(nrow(flip_event_func_inclusion))
num_incl_func <- c(nrow(psi_neg_func_tab))
num_flip_incl_func_perc <- (num_flip_incl_func/(num_incl_func)) * 100
num_non_flip_incl_func_perc      <- ( (num_incl_func-num_flip_incl_func)/(num_incl_func)) * 100

##create dataframe for plotting
type <- c("Skipping", "Skipping", "Inclusion", "Inclusion")
counts <- c(num_flip_incl_func_perc,num_non_flip_incl_func_perc,
            num_flip_skip_func_perc,num_non_flip_skip_func_perc)
event <- c("Flip", "Non-flip","Flip", "Non-flip")
num_of_hits_perc <- data.frame(type,counts, event) %>% 
  filter(event=="Flip")

##plot
plot_flip <- ggplot(num_of_hits_perc, aes(x = type, y = counts, fill = type)) + 
         geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=c("#FFC20A","#0C7BDC")) +
  xlab("Type ") +
  ylab("Splice Variants (%)") + 
  coord_flip() +
  theme_Publication() + 
  ggtitle("Flip Events") +
  theme(legend.position="none", legend.title = element_text(size=19), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))

# Save plot as PDF
pdf(file_flip_events_plot, 
    width = 5, height = 3)
plot_flip
dev.off()


