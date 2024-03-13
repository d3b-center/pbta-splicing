################################################################################
# 05-clk1-gtex.R
# written by Jo Lynne Rokita
#
# This script generates CLK1 expression boxplot for gtex brain and HGGs
#
# usage: Rscript 05-clk1-gtex.R
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("ggplot2")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")

## create plots dir if it doesn't exist
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}


## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## input files
hist_file <- file.path(data_dir, "histologies.tsv")
v15_hist_file <- file.path(data_dir, "v15-histologies.tsv")
pbta_counts_file <- file.path(data_dir, "gene-counts-rsem-expected_count-collapsed.rds") 
gtex_counts_file <- file.path(data_dir, "gtex_gene-counts-rsem-expected_count-collapsed.rds")
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

## output file for plot
gtex_plot_path <- file.path(plots_dir, "CLK1-gtex-hgg.pdf")

## load files
indep_df <- read_tsv(indep_file) %>% 
  dplyr::filter(cohort=='PBTA')

all_hgg_hist <- read_tsv(hist_file, guess_max = 100000) %>% 
  dplyr::filter(short_histology == "HGAT",
                Kids_First_Biospecimen_ID %in% 
                  indep_df$Kids_First_Biospecimen_ID)

gtex_hist <- read_tsv(v15_hist_file, guess_max = 100000) %>% 
  dplyr::filter(cohort == "GTEx")

library(tidyverse)

hgg_counts <- readRDS(pbta_counts_file) %>%
  select(all_hgg_hist$Kids_First_Biospecimen_ID) %>%
  rownames_to_column("Gene") %>%
  filter(Gene == "CLK1") %>%
  select(-Gene) %>%
  # Transpose and then immediately convert to a data frame or tibble, capturing row names as a column
  {data.frame(t(.))} %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::rename(CLK1 = `t...`) %>%
  mutate(plot_group = "pHGG")


gtex_counts <- readRDS(gtex_counts_file) %>%
  rownames_to_column("Gene") %>%
  filter(Gene == "CLK1") %>%
  select(-Gene) %>%
  {data.frame(t(.))} %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::rename(CLK1 = `t...`) %>%
  left_join(gtex_hist[,c("Kids_First_Biospecimen_ID", "gtex_group", "gtex_subgroup")]) %>%
  dplyr::rename(plot_group = gtex_subgroup) %>%
  filter(gtex_group == "Brain")

combined_df <- hgg_counts %>%
  bind_rows(gtex_counts)

# Calculate medians and reorder plot_group levels based on these medians
combined_df_ordered <- combined_df %>%
  group_by(plot_group) %>%
  mutate(median_CLK1 = median(CLK1, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(plot_group = factor(plot_group, levels = unique(plot_group[order(median_CLK1)])))

# plot
gtex_plot <- ggplot(combined_df_ordered, aes(x = plot_group, y = log2(CLK1))) +
  geom_boxplot(outlier.shape = NA) + # Remove outliers for cleaner look
  geom_jitter(width = 0.2, color = "black", alpha = 0.5) + # Add actual data points
  labs(title = "",
       x = "Group",
       y = "log2[CLK1 RSEM expected counts]") +
  theme_Publication() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 75, hjust = 1))

pdf(gtex_plot_path, height = 8, width = 8, useDingbats = FALSE)
print(gtex_plot)
dev.off()

  



