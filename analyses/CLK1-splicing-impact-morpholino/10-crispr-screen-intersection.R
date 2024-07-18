################################################################################
# 09-crispr-screen-intersection.R
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
  library("ggrepel")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")

input_dir   <- file.path(analysis_dir, "input")
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
gene_dependent_crispr_file <- file.path(results_dir,"clk1-de-ds-crispr-targets.txt")

# read in clinical
clin_df <- read_tsv(clin_df, guess_max = 100000) %>%
  dplyr::filter(short_histology == "HGAT",
                cohort == "PBTA",
                composition == "Solid Tissue")

clk1_targets <- read_lines(clk1_target_file)

crispr_df <- read_csv(crispr_score) 

# filter for HGG cell lines only 
crispr_dep <- crispr_df %>%
  filter(wald_p_value < 0.05 & wald_fdr < 0.05) %>% 
  mutate(sample_id = str_replace(sample, "_[^_]*$", "")) %>%
  mutate(sample_id = str_replace_all(sample_id, "_", "-")) %>%
  dplyr::filter(sample_id %in% clin_df$sample_id) %>%
  dplyr::select(gene,sample,sample_id,z,beta) %>%
  group_by(gene, sample_id) %>%
  summarise(z = mean(z),
            beta = mean(beta)) %>%
  ungroup()

clk1_targets_crispr <- intersect(clk1_targets, crispr_dep$gene)

# write list of targets in any cell line
unique(clk1_targets_crispr) %>%
  write_lines(gene_dependent_crispr_file)

# Create swoosh plot (z with < -1.5 only)
# Group by gene and calculate the mean z value across all cell lines, add this to the df
mean_data <- crispr_dep %>%
  group_by(gene) %>%
  summarise(z = mean(z, na.rm = TRUE)) %>%
  arrange(desc(z)) %>% 
  filter(z < -1.5) %>%
  mutate(sample_id = "pHGG")

# add this to the full df to create an all hgg plot
crispr_dep_full <- crispr_dep %>%
  bind_rows(mean_data)

# Create cell line specific plots

for (each in unique(crispr_dep_full$sample_id)) {
  crispr_dep_each <- crispr_dep_full %>%
    filter(sample_id == each)
  clk1_targets_crispr_intersect <- crispr_dep_each %>%
    filter(gene %in% clk1_targets_crispr) 
  
crispr_scores_z_plot <- ggplot(crispr_dep_each, aes(x = reorder(gene, z), y = z)) +
  geom_point(size=3, colour="gray89") + 
  geom_point(size=3, colour = "gray50", pch = 21) + 
  geom_point(data=clk1_targets_crispr_intersect, colour="red", size = 3) +
  geom_point(data=clk1_targets_crispr_intersect, colour="black", size = 3, pch = 21) +
  geom_hline(yintercept = -1.5, linetype = "dashed", color = "gray50") +
  labs(title = paste0(each),
       x = "Gene",
       y = "CRISPR Dependencey Z-score") +
  #geom_text_repel(data = mean_z_15 %>% filter(gene %in% clk1_targets_crispr), aes(label = gene), color = "black") +
  geom_text_repel(data = crispr_dep_each %>% filter(gene %in% clk1_targets_crispr_intersect$gene), 
                  aes(label = gene), color = "black", 
                  nudge_y = 0.7, # Adjust the nudging value to avoid overlap
                  box.padding = 0.5, # Add padding around the label
                  point.padding = 0.5, # Add padding around the point
                  segment.size = 0.2, # Set the size of the line segment connecting the label to the point
                  max.overlaps = Inf) + # Allow for infinite overlaps, which ggrepel will handle
 # facet_grid(~sample_id, scales = "free_x") +
  theme_Publication() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  ) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 0.6)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 0.5)

# Save plot as pdf
pdf(file.path(paste0(plots_dir, "/clk1-crispr-swoosh-", each, ".pdf")), height = 4, width = 6)
print(crispr_scores_z_plot)
dev.off()

}
