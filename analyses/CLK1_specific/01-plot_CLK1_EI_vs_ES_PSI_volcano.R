################################################################################
# 01-plot_CLK1_EI_vs_ES_PSI_volcano.R
# written by Ammar Naqvi
#
# usage: Rscript 01-plot_CLK1_EI_vs_ES_PSI_volcano.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("ggrepel")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## get data files from rMATS run and make table
input      = file.path(input_dir,"dca735c2-6e0e-4239-8a68-10c6d2aa9015.CLK1_EI_vs_CLK1_ES.non_denovo.SE.MATS.JC.txt")
tab =read.table(input,header=TRUE,sep = "\t")


# add a column of NAs
tab$dPSI <- "Not Significant"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
tab$dPSI[tab$IncLevelDifference >= .20 & tab$PValue < 0.05] <- "Skipping"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tab$dPSI[tab$IncLevelDifference <= -.20 & tab$PValue < 0.05] <- "Inclusion"


# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Inclusion", "Skipping", "Not Signficant")

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
tab$delabel <- NA
tab$delabel[tab$dPSI != "Not Significant"] <- tab$geneSymbol[tab$dPSI != "Not Significant"]

# plot adding up all layers we have seen so far
plot_volcano <- ggplot(data=tab, aes(x=IncLevelDifference, y=-log10(PValue), col=dPSI, label=delabel, label.size=.05)) +
  geom_point() + 
  theme_Publication() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-0.20, 0.20), col="black", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="black",linetype = "longdash") 

plot_file = file.path(plots_dir,"dPSI_volcano_CLK1.pdf") 

# Save plot as PDF
pdf(plot_file, width = 10, height = 15)
plot_volcano 
dev.off()

