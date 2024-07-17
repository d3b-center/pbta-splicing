################################################################################
# 10-crispr-screen-intersection.R
# Intersect CRISPR scores with CLK1 targets
#
# written by Ammar Naqvi, Jo Lynne Rokita
# Usage: Rscript 10-crispr-screen-intersection.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("vroom")
  library("tidyverse")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")

input_dir   <- file.path(analysis_dir, "input")
output_dir   <- file.path(analysis_dir, "output")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## input files
clin_df <- file.path(data_dir,"histologies.tsv")
crispr_score <- file.path(input_dir,"CCMA_crispr_genedependency_042024.csv")
clk1_target_file <- file.path(results_dir,"common_genes_de_ds_functional.txt")

## ouput files
clk1_targets_crispr_file <- file.path(results_dir,"clk1-targets-crispr-cbtn-lines.txt")
crispr_score_plot_file <- file.path(plots_dir,"clk1-targets-crispr-cbtn-lines.pdf")
crispr_score_sign_plot_file <- file.path(plots_dir,"clk1-targets-crispr_cbtn_lines-sign.pdf")

clin_df <- vroom(clin_df) %>%
  dplyr::filter(short_histology == "HGAT",
                cohort == "PBTA") 

clk1_targets <- read_lines(clk1_target_file)

crispr <- read_csv(crispr_score) %>%
  filter(wald_p_value < 0.05 & wald_fdr < 0.10) %>% 
  mutate(sample_id = str_replace(sample, "_[^_]*$", "")) %>%
  mutate(sample_id = str_replace_all(sample_id, "_", "-")) %>%
  inner_join(clin_df, by="sample_id") %>%
  dplyr::filter(composition=="Solid Tissue") %>%
  dplyr::select(gene,sample_id,z,beta) %>%
  distinct()

clk1_targets_crispr <- intersect(clk1_targets, crispr$gene)

unique(clk1_targets_crispr) %>%
  write_lines(clk1_targets_crispr_file)

# Create swoosh plot
# Group by gene and calculate the mean z value
mean_z <- crispr %>%
  group_by(gene) %>%
  summarise(mean_z = mean(z, na.rm = TRUE)) %>%
  arrange(desc(mean_z))

## identify CLK1 targets
clk1_targets_scores <- mean_z %>%
  filter(gene %in% clk1_targets_crispr,
         mean_z < 1.5)

# Create the plot
crispr_scores_plot <- ggplot(mean_z, aes(x = reorder(gene, mean_z), y = mean_z)) +
  #geom_point(color = "blue", size = 3) +
  geom_point(size=3, colour="blue",size=3) + 
  geom_point(data=clk1_targets_scores, colour="red", size = 3) +
  geom_point(data=clk1_targets_scores, colour="black", size = 3, pch = 21) +
  geom_hline(yintercept = -1.5, color = "red", linetype = "dashed") +
  labs(title = "CRISPR Depedency Scores",
       x = "Gene",
       y = "Mean Z") +
  geom_text_repel(data = mean_z %>% filter(gene %in% clk1_targets_crispr & mean_z < 1.5), aes(label = gene), color = "black") +
  theme_Publication() +
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank()) # Remove x-axis tick

# Save plot as pdf
pdf(crispr_score_plot_file, height = 4, width = 7.5)
crispr_scores_plot
dev.off()

# Create swoosh plot (z with < 1.5 only)
# Group by gene and calculate the mean z value
mean_z_15 <- crispr %>%
  group_by(gene) %>%
  summarise(mean_z = mean(z, na.rm = TRUE)) %>%
  arrange(desc(mean_z)) %>% 
  filter(mean_z < -1.5)

## identify CLK1 targets
clk1_targets_scores_z15 <- mean_z_15 %>%
  filter(gene %in% clk1_targets_crispr)

# Create the plot
crispr_scores_z_plot <- ggplot(mean_z_15, aes(x = reorder(gene, mean_z), y = mean_z)) +
  geom_point(size=3, colour="blue",size=3) + 
  geom_point(data=clk1_targets_scores_z15, colour="red", size = 3) +
  geom_point(data=clk1_targets_scores_z15, colour="black", size = 3, pch = 21) +
  labs(title = "CRISPR Depedency Scores",
       x = "Gene",
       y = "Mean Z") +
  geom_text_repel(data = mean_z_15 %>% filter(gene %in% clk1_targets_crispr), aes(label = gene), color = "black") +
  theme_Publication() +
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank()) # Remove x-axis tick

# Save plot as pdf
pdf(crispr_score_sign_plot_file, height = 4, width = 7.5)
crispr_scores_z_plot
dev.off()
