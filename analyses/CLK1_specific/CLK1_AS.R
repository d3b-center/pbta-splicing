################################################################################
# CLK1_AS.R
# script 
#
# written by Ammar Naqvi
#
# usage: Rscript CLK1_AS.R
################################################################################

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

dataDir = "/Users/naqvia/Desktop/AS-DMG/analyses/CLK1_specific/results/"

## more accurate version -- dPSI version
file <- "CLK1_exon_inclusion.txt"
inc_levels  = read.delim(paste0(dataDir, file), sep = "\t", header=TRUE,row.names=1)

# results/CLK1_exon_inclusion.above10perc.txt
file_fil <- "CLK1_exon_inclusion.above20perc.txt"
inc_levels  = read.delim(paste0(dataDir, file_fil), sep = "\t", header=TRUE,row.names=1)


## violin plot version
e <- ggplot(inc_levels, aes(x = Hist, y = Inclusion))


e + geom_violin(trim = FALSE) + 
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black"
  )

e + geom_violin(aes(fill = Hist), trim = FALSE) + 
  geom_boxplot(width = 0.2) +
  #scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  geom_point(color="black", size=1, position = position_jitter(w=0.02)) +
  theme(legend.position = "none") + theme_Publication()

## results/CLK1_exon_psi.hgats_subtyped.txt
file_fil <- "CLK1_exon_psi.hgats_subtyped.txt"
dpsi_levels  = read.delim(paste0(dataDir, file_fil), sep = "\t", header=TRUE,row.names=1)


## box plot version for subtypes
e <- ggplot(dpsi_levels, aes(x = Hist, y = Inclusion,fill = Hist))

e + geom_boxplot(width = 0.4) +
  geom_point(color="black", size=.5, position = position_jitter(w=0.02)) +
  theme(legend.position = "none") + theme_Publication()
  
  



