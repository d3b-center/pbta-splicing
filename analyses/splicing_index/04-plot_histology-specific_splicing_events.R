################################################################################
# upsetR_of_hist-specific_splicing.R
# written by Ammar Naqvi
#
# usage: Rscript upsetR_of_hist-specific_splicing.R
################################################################################

##libraries 
suppressPackageStartupMessages({
  library("ggplot2")
  library("UpSetR")
  library("grid")
  library("cowplot")
  library("viridis")
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## output files for final plots
upsetR_es_plot_file     <- file.path(analysis_dir, "plots", "upsetR_histology-specific.es.pdf")
upsetR_ei_plot_file <- file.path(analysis_dir, "plots", "upsetR_histology-specific.ei.pdf")
upsetR_tiff_es_plot_file     <- file.path(analysis_dir, "plots", "upsetR_histology-specific.es.tiff")
upsetR_tiff_ei_plot_file <- file.path(analysis_dir, "plots", "upsetR_histology-specific.ei.tiff")

splice_event_df = vroom(paste0(results_dir,"/","splicing_events.hist-labeled_list.thr10freq.txt"), delim="\t", trim_ws = TRUE, col_names = TRUE)

## select and create list for exon skipping only events
splice_event_df_ATRT <- splice_event_df %>% filter(histology=='ATRT', type=='skipping') %>% select(splicing_event)
splice_event_df_CPG <- splice_event_df %>% filter(histology=='CPG', type=='skipping') %>% select(splicing_event)
splice_event_df_GNG <- splice_event_df %>% filter(histology=='GNG', type=='skipping') %>% select(splicing_event)
splice_event_df_EPN <- splice_event_df %>% filter(histology=='EPN', type=='skipping') %>% select(splicing_event)
splice_event_df_HGG <- splice_event_df %>% filter(histology=='HGG', type=='skipping') %>% select(splicing_event)
splice_event_df_MB <- splice_event_df %>% filter(histology=='MB', type=='skipping') %>% select(splicing_event)
splice_event_df_LGG <- splice_event_df %>% filter(histology=='LGG', type=='skipping') %>% select(splicing_event)

list_for_skipping_upsetR <- list("ATRT"=splice_event_df_ATRT$splicing_event,
                                  "CPG"=splice_event_df_CPG$splicing_event,
                                 "GNG"=splice_event_df_GNG$splicing_event,
                                 "EPN"=splice_event_df_EPN$splicing_event,
                                 "HGG"=splice_event_df_HGG$splicing_event,
                                 "MB"=splice_event_df_MB$splicing_event,
                                 "LGG"=splice_event_df_LGG$splicing_event)
                                 
                                 

es_events <- upset(fromList(list_for_skipping_upsetR), 
                   mainbar.y.label = "", sets=c("ATRT","CPG","GNG","EPN","HGG","MB","LGG"), sets.x.label = "Histology", order.by = "freq",
                   mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 2, line.size = 1.5, nsets = 7)


# Save plot as PDF
pdf(upsetR_es_plot_file, width = 16, height = 8)
es_events
dev.off()

# Save plot tiff version
tiff(upsetR_tiff_es_plot_file, height = 1200, width = 2400, res = 300)
print(es_events)
dev.off()

## select and create list for exon inclusion only events
splice_event_df_ATRT <- splice_event_df %>% filter(histology=='ATRT', type=='inclusion') %>% select(splicing_event)
splice_event_df_CPG <- splice_event_df %>% filter(histology=='CPG', type=='inclusion') %>% select(splicing_event)
splice_event_df_GNG <- splice_event_df %>% filter(histology=='GNG', type=='inclusion') %>% select(splicing_event)
splice_event_df_EPN <- splice_event_df %>% filter(histology=='EPN', type=='inclusion') %>% select(splicing_event)
splice_event_df_HGG <- splice_event_df %>% filter(histology=='HGG', type=='inclusion') %>% select(splicing_event)
splice_event_df_MB <- splice_event_df %>% filter(histology=='MB', type=='inclusion') %>% select(splicing_event)
splice_event_df_LGG <- splice_event_df %>% filter(histology=='LGG', type=='inclusion') %>% select(splicing_event)

list_for_inclusion_upsetR <- list("ATRT"=splice_event_df_ATRT$splicing_event,
                                 "CPG"=splice_event_df_CPG$splicing_event,
                                 "GNG"=splice_event_df_GNG$splicing_event,
                                 "EPN"=splice_event_df_EPN$splicing_event,
                                 "HGG"=splice_event_df_HGG$splicing_event,
                                 "MB"=splice_event_df_MB$splicing_event,
                                 "LGG"=splice_event_df_LGG$splicing_event)

## make upsetR plots
ei_events <- upset(fromList(list_for_inclusion_upsetR), 
                   mainbar.y.label = "", sets.x.label = "Histology", order.by = "freq",
                   mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 3, line.size = 1.5, nsets = 7)

# Save plot as PDF
pdf(upsetR_ei_plot_file, width = 16, height = 8)
ei_events
dev.off()

# Save plot tiff version
tiff(upsetR_tiff_ei_plot_file, height = 1200, width = 2400, res = 300)
print(ei_events)
dev.off()


