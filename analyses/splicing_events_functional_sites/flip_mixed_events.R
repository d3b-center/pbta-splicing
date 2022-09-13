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

## output file names for plots
file_flip_events_plot <- file.path(analysis_dir, "plots", "flip_barplots.png")

## retrieve psi values from tables
file_psi_pos_total <- "/splicing_events.total.pos.tsv"
file_psi_neg_total <- "/splicing_events.total.neg.tsv"

file_psi_pos_func <- "/splicing_events.total.pos.intersectUnip.ggplot.txt"
file_psi_neg_func <- "/splicing_events.total.neg.intersectUnip.ggplot.txt"

psi_pos_func_tab <-  read.delim(paste0(results_dir, file_psi_pos_func), sep = "\t", row.names = NULL, header=TRUE)
psi_pos_tab      <-  read.delim(paste0(results_dir, file_psi_pos_total), sep = "\t", row.names = NULL, header=TRUE)  %>% filter(flip == 1)

psi_neg_func_tab <-  read.delim(paste0(results_dir, file_psi_neg_func), sep = "\t", row.names = NULL, header=TRUE)
psi_neg_tab      <-  read.delim(paste0(results_dir, file_psi_neg_total), sep = "\t", row.names = NULL, header=TRUE)  %>% filter(flip == 1)

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
            legend.position = "top",
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


##num of sites that are undergo flip at functional sites
flip_event_func_skipping  <- psi_pos_tab[psi_pos_func_tab$SpliceID %in% psi_pos_tab$splice_event, ]
flip_event_func_inclusion <- psi_neg_tab[psi_neg_func_tab$SpliceID %in% psi_neg_tab$splice_event, ]

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

##create daatafram for plotting
type <- c("Skipping", "Skipping", "Inclusion", "Inclusion")
counts <- c(num_flip_incl_func_perc,num_non_flip_incl_func_perc,
            num_flip_skip_func_perc,num_non_flip_skip_func_perc)
event <- c("Flip", "Non-flip","Flip", "Non-flip")
num_of_hits_perc <- data.frame(type,counts, site)

##plot
ggplot(num_of_hits_perc, aes(x = type, y = counts, fill = event)) + 
         geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=c("red","blue")) + 
  coord_flip() +
  theme_Publication()

## save plot
ggsave(
  file_flip_events_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =410,
  height = 250,
  units = "mm",
  dpi = 800,
  limitsize = TRUE,
  bg = NULL
)

