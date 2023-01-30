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

## exon skipping
ATRT_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.pos.ATRT.txt"), sep = "\t", header=FALSE)
CPG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.pos.CPG.txt"), sep = "\t", header=FALSE)
GNG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.pos.GNG.txt"), sep = "\t", header=FALSE)
EPN_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.pos.EPN.txt"), sep = "\t", header=FALSE)
HGG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.pos.HGG.txt"), sep = "\t", header=FALSE)
MB_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.pos.MB.txt"), sep = "\t", header=FALSE)
LGG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.pos.LGG.txt"), sep = "\t", header=FALSE)

listInput <- list("ATRT" =ATRT_events$V1, 
                  "CPG" =CPG_events$V1,
                  "GNG" =GNG_events$V1,
                  "EPN" =EPN_events$V1,
                  "HGG" =HGG_events$V1, 
                  "MB" = MB_events$V1,
                  "LGG"=LGG_events$V1)

es_events <- upset(fromList(listInput), 
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

## exon inclusion
ATRT_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.neg.ATRT.txt"), sep = "\t", header=FALSE)
CPG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.neg.CPG.txt"), sep = "\t", header=FALSE)
GNG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.neg.GNG.txt"), sep = "\t", header=FALSE)
EPN_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.neg.EPN.txt"), sep = "\t", header=FALSE)
HGG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.neg.HGG.txt"), sep = "\t", header=FALSE)
MB_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.neg.MB.txt"), sep = "\t", header=FALSE)
LGG_events = read.delim(paste0(results_dir, "/splicing_events.hist-labeled_list.thr10freq.neg.LGG.txt"), sep = "\t", header=FALSE)

listInput_ei <- list("ATRT" =ATRT_events$V1, 
                  "CPG" =CPG_events$V1,
                  "GNG" =GNG_events$V1,
                  "EPN" =EPN_events$V1,
                  "HGG" =HGG_events$V1, 
                  "MB" =MB_events$V1,
                  "LGG"=LGG_events$V1)


ei_events <- upset(fromList(listInput_ei), 
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


