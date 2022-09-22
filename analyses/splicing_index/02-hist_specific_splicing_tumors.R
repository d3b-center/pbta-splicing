################################################################################
# upsetR_of_hist-specific_splicing.R
# written by Ammar Naqvi
#
# usage: Rscript upsetR_of_hist-specific_splicing.R
################################################################################
library(UpSetR)
library(ggplot2)
library(gplots)
library(cowplot)
library("grid")

##set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
plots_dir <- file.path(root_dir, "analyses/splicing_index", "plots")
results_dir <- file.path(root_dir, "analyses/splicing_index", "results")

## output files for final plots
upsetR_es_plot     <- file.path(analysis_dir, "plots", "upsetR_histology-specific.es.pdf")
upsetR_ei_plot <- file.path(analysis_dir, "plots", "upsetR_histology-specific.ei.pdf")

## exon skipping
ATRT_events = read.delim(paste0(results_dir, "splicing_events.hist-labeled_list.thr10freq.pos.ATRT.txt"), sep = "\t", header=FALSE)
CPG_events = read.delim(paste0(results_dir, "splicing_events.hist-labeled_list.thr10freq.pos.Craniopharyngioma.txt"), sep = "\t", header=FALSE)
GNG_events = read.delim(paste0(results_dir, "splicing_events.hist-labeled_list.thr10freq.pos.Ganglioglioma.txt"), sep = "\t", header=FALSE)
EPN_events = read.delim(paste0(results_dir, "splicing_events.hist-labeled_list.thr10freq.pos.Ependymoma.txt"), sep = "\t", header=FALSE)
HGAT_events = read.delim(paste0(results_dir, "splicing_events.hist-labeled_list.thr10freq.pos.HGAT.txt"), sep = "\t", header=FALSE)
MB_events = read.delim(paste0(results_dir, "splicing_events.hist-labeled_list.thr10freq.pos.Medulloblastoma.txt"), sep = "\t", header=FALSE)
LGAT_events = read.delim(paste0(results_dir, "splicing_events.hist-labeled_list.thr10freq.pos.LGAT.txt"), sep = "\t", header=FALSE)

listInput <- list("ATRT" =ATRT_events$V1, 
                  "CPG" =CPG_events$V1,
                  "GNG" =GNG_events$V1,
                  "EPN" =EPN_events$V1,
                  "HGAT" =HGAT_events$V1, 
                  "MB" =MB_events$V1,
                  "LGAT"=LGAT_events$V1)

es_events <- upset(fromList(listInput), 
      mainbar.y.label = "", sets=c("ATRT","CPG","GNG","EPN","HGAT","MB","LGAT"), sets.x.label = "Histology", order.by = "freq",
      mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 2, line.size = 1.5, nsets = 7)

ggsave(upsetR_es_plot, width = 16, height = 4)


## exon inclusion
ATRT_events = read.delim(paste0(analysis_dir, "splicing_events.hist-labeled_list.thr10freq.neg.ATRT.txt"), sep = "\t", header=FALSE)
CPG_events = read.delim(paste0(analysis_dir, "splicing_events.hist-labeled_list.thr10freq.neg.Craniopharyngioma.txt"), sep = "\t", header=FALSE)
GNG_events = read.delim(paste0(analysis_dir, "splicing_events.hist-labeled_list.thr10freq.neg.Ganglioglioma.txt"), sep = "\t", header=FALSE)
EPN_events = read.delim(paste0(analysis_dir, "splicing_events.hist-labeled_list.thr10freq.neg.Ependymoma.txt"), sep = "\t", header=FALSE)
HGAT_events = read.delim(paste0(analysis_dir, "splicing_events.hist-labeled_list.thr10freq.neg.HGAT.txt"), sep = "\t", header=FALSE)
MB_events = read.delim(paste0(analysis_dir, "splicing_events.hist-labeled_list.thr10freq.neg.Medulloblastoma.txt"), sep = "\t", header=FALSE)
LGAT_events = read.delim(paste0(analysis_dir, "splicing_events.hist-labeled_list.thr10freq.neg.LGAT.txt"), sep = "\t", header=FALSE)

listInput_ei <- list("ATRT" =ATRT_events$V1, 
                  "CPG" =CPG_events$V1,
                  "GNG" =GNG_events$V1,
                  "EPN" =EPN_events$V1,
                  "HGAT" =HGAT_events$V1, 
                  "MB" =MB_events$V1,
                  "LGAT"=LGAT_events$V1)


ei_events <- upset(fromList(listInput_ei), 
      mainbar.y.label = "", sets.x.label = "Histology", order.by = "freq",
      mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 2, line.size = 1.5, nsets = 7)

ggsave(upsetR_ei_plot, width = 16, height = 4)

