################################################################################
# corr_plots_psi_vs_expr.R
# written by Ammar Naqvi
#
# usage: Rscript 05-plot_psi_vs_expr.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  
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

## function to get stats from to print on box plot
get_box_stats <- function(y, upper_limit = max(corr_df$V3) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

## output files for final plots
file_corr_plot <- file.path(analysis_dir, "plots", "corr_CLK1_spl_vs_expr.pdf")
file_tiff_corr_plot <- file.path(analysis_dir, "plots", "corr_CLK1_spl_vs_expr.tiff")

corr_expr_spl_file = file.path(results_dir, "expr_vs_psi_corr_res.txt")
corr_df=read.delim(corr_expr_spl_file,header=FALSE,sep = " ", row.names = NULL)

## make corr plot
boxplot_corr <- ggplot(corr_df, aes(V3, V2))  + stat_boxplot(fill = "lightgreen", colour = "black", outlier.colour = "red") +
  xlab("correlation coeff") + geom_jitter(width = 0.3, shape=16) + ylab("") + xlim(-1,1) + 
  ggtitle("Exon splicing vs gene expression") +  stat_summary(fun.data = get_box_stats, geom = "text", hjust = .4, vjust = .2) +
  theme_Publication()  

## save plot
ggsave(file_corr_plot, width = 15, height = 5)

## save plot in tiff format
# Save plot tiff version
tiff(file_tiff_corr_plot, height = 1200, width = 3600, res = 300)
print(boxplot_corr)
dev.off()