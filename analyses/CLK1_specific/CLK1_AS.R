library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(ggplot2)
library(dplyr)
library(ggstatsplot)


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

file="/Users/naqvia/Desktop/pbta-splicing/analyses/splicing_events_functional_sites/scr/CLK1_dpsi_samples.v2.txt" 
dpsi<- read.table(file, header=TRUE,sep = "\t")


#JitterPlot = qplot(Type, dPSI, data=dpsi, fill = Type, geom=c("boxplot", "jitter"), 
#                  main="CLK1 Aberrant Splicing",xlab="Type", ylab="dPSI") + theme_Publication()

##set working dir for module
setwd("/Users/naqvia/Desktop/pbta-splicing_git/pbta_splicing/analyses/CLK1_specific") 

set.seed(123)
ggstatsplot::ggbetweenstats(
  data = dpsi, 
  x = Type, 
  y = dPSI,
  k = 3,
  nboot = 15,
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  type = "robust",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication() 

ggsave("plots/CLK1_dPSI_sign.pdf", width = 8, height = 15)
