################################################################################
# splicing_functional_sites.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi
#
# usage: Rscript splicing_functional_sites.R <file>
################################################################################

library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(ggplot2)
library(dplyr)
library(ggstatsplot)

## hardcoded for now
## table with splicing events and associated uniprot-annotaed site 
file="/Users/naqvia/Desktop/AS-DMG/analyses/splicing_events_functional_sites/results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt" 
dpsi_unip <- read.table(file, header=TRUE,sep = "\t")

set.seed(123)

ggstatsplot::ggbetweenstats(
  data = dpsi_unip, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  #outlier.label = Splice_ID, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = FALSE,
  mean.ci = TRUE,
  outlier.tagging = FALSE,
  #title = "Tumor supressors",
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication()

## table with splicing events and associated uniprot-annotaed site dvided by negative and positive
file_pos = "/Users/naqvia/Desktop/AS-DMG/analyses/splicing_events_functional_sites/results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt" 
file_neg = "/Users/naqvia/Desktop/AS-DMG/analyses/splicing_events_functional_sites/results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt" 
dpsi_unip_pos <- read.table(file_pos, header=TRUE,sep = "\t")
dpsi_unip_neg <- read.table(file_neg, header=TRUE,sep = "\t")

set.seed(123)
ggstatsplot::ggbetweenstats(
  data = dpsi_unip_pos, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  #nboot = 15,
  outlier.label = Splice_ID, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  #outlier.tagging = TRUE,
  #title = "Tumor supressors",
  #type = "robust",
  #xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication()

set.seed(123)
ggstatsplot::ggbetweenstats(
  data = dpsi_unip_neg, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  #nboot = 15,
  outlier.label = Splice_ID, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  #outlier.tagging = TRUE,
  #title = "Tumor supressors",
  type = "robust",
  #xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication()