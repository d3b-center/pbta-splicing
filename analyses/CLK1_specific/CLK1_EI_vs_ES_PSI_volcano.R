library("ggplot2")

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
root_dir <- "/Users/naqvia/Desktop/pbta-splicing_git/pbta-splicing"
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## get data files and make table
input      = file.path(input_dir,"dca735c2-6e0e-4239-8a68-10c6d2aa9015.CLK1_EI_vs_CLK1_ES.non_denovo.SE.MATS.JC.txt")
tab=read.table(input,header=TRUE,sep = "\t")

p <- ggplot(data=tab, aes(x=IncLevelDifference, y=-log10(PValue))) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-.20, +.20), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
tab$dPSI <- "Not Significant"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
tab$dPSI[tab$IncLevelDifference >= .20 & tab$PValue < 0.05] <- "Skipping"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tab$dPSI[tab$IncLevelDifference <= -.20 & tab$PValue < 0.05] <- "Inclusion"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=tab, aes(x=IncLevelDifference, y=-log10(PValue), col=dPSI)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.20, 0.20), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

## Change point color 
# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Inclusion", "Skipping", "Not Signficant")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
tab$delabel <- NA
tab$delabel[tab$dPSI != "Not Significant"] <- tab$geneSymbol[tab$dPSI != "Not Significant"]

ggplot(data=tab, aes(x=IncLevelDifference, y=-log10(PValue), col=dPSI, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

library(ggrepel)

# plot adding up all layers we have seen so far
ggplot(data=tab, aes(x=IncLevelDifference, y=-log10(PValue), col=dPSI, label=delabel, label.size=.05)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.20, 0.20), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") 

plot_file = file.path(plots_dir,"dPSI_volcano_CLK1.pdf") 
ggsave(plot_file, width = 10, height = 15)
  