################################################################################
# 01-plot_CLK1_EI_vs_ES_PSI_volcano.R
# written byAmmar Naqvi and Jo Lynne Rokita
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

## set directories and file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# set file path
plot_file = file.path(plots_dir,"dPSI_volcano_CLK1.pdf")

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))


diff_splice_file <- file.path(data_dir,"CLK1_EI_vs_CLK1_ES_nondenovo.tsv")
diff_splice_df <- read_tsv(diff_splice_file)

# Add a few columns required for plotting
diff_splice_df_label <- diff_splice_df %>%
  mutate(neg_log10_p = -log10(PValue),
  # Determine whether dPSI value indicates a significant event and annotate
  dPSI = case_when(IncLevelDifference >= .20 & PValue < 0.05 ~ "Skipping",
                   IncLevelDifference <= -.20 & PValue < 0.05 ~ "Inclusion",
                   TRUE ~ "Not Signficant"),
  # Create a new column "point_label" that will contain the name of genes only if differentially spliced
  point_label = case_when(dPSI %in% c("Skipping", "Inclusion") ~ as.character(geneSymbol),
                          TRUE ~ NA_character_))


plot_volcano <- ggplot(data=diff_splice_df_label,
                      aes(x=IncLevelDifference, y=neg_log10_p, col=dPSI,
                      label=point_label, label.size=.05)) +
                      geom_point() +
                      theme_Publication() +
                      theme(legend.position="none") +
                      geom_text_repel() +
                      scale_color_manual(values=c("blue", "grey", "red")) +
                      geom_vline(xintercept=c(-0.20, 0.20), col="black", linetype = "longdash") +
                      geom_hline(yintercept=-log10(0.05), col="black",linetype = "longdash") +
                      xlab(bquote("dPSI")) +
                      ylab(bquote(bold(~-log[10] ~ "p-value")))

plot_file = file.path(plots_dir,"dPSI_volcano_CLK1.pdf")

# Save plot as PDF
pdf(plot_file, width = 10, height = 10)
plot_volcano
dev.off()
