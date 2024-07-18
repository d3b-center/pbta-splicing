################################################################################
# 09-crispr-screen-intersection.R
# Intersect CRISPR scores with CLK1 targets
#
# written by Ammar Naqvi, Jo Lynne Rokita
# Usage: Rscript 09-crispr-screen-intersection.R
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
gene_dependent_crispr_file <- file.path(results_dir,"sign-crispr-cbtn-lines.txt")

crispr_score_plot_file <- file.path(plots_dir,"clk1-targets-crispr-cbtn-lines.pdf")
crispr_score_sign_plot_file <- file.path(plots_dir,"clk1-targets-crispr_cbtn_lines-sign.pdf")

clin_df <- vroom(clin_df) %>%
  dplyr::filter(short_histology == "HGAT",
                cohort == "PBTA") 

clk1_targets <- read_lines(clk1_target_file)

crispr_df <- read_csv(crispr_score) 

# filter for HGG cell lines only 
crispr_dep <- crispr_df %>%
  filter(wald_p_value < 0.05 & wald_fdr < 0.05) %>% 
  mutate(sample_id = str_replace(sample, "_[^_]*$", "")) %>%
  mutate(sample_id = str_replace_all(sample_id, "_", "-")) %>%
  inner_join(clin_df, by="sample_id") %>%
  dplyr::filter(composition=="Solid Tissue") %>%
  dplyr::select(gene,sample,sample_id,z,beta) %>%
  group_by(gene, sample_id) %>%
  summarise(z = mean(z),
            beta = mean(beta)) %>%
  ungroup()

clk1_targets_crispr <- intersect(clk1_targets, crispr$gene)

unique(clk1_targets_crispr) %>%
  write_lines(clk1_targets_crispr_file)


# Create swoosh plot (z with < -1.5 only)
# Group by gene and calculate the mean z value
mean_z_15 <- crispr_dep %>%
  group_by(gene) %>%
  summarise(mean_z = mean(z, na.rm = TRUE)) %>%
  arrange(desc(mean_z)) %>% 
  filter(mean_z < -1.5)

unique(mean_z_15$gene) %>%
  write_lines(gene_dependent_crispr_file)

## identify CLK1 targets
clk1_targets_scores_z15 <- mean_z_15 %>%
  filter(gene %in% clk1_targets_crispr)

clk1_targets_scores_z15 <- crispr_dep %>%
  filter(gene %in% clk1_targets_crispr)

# Create the plot
crispr_scores_z_plot <- ggplot(crispr_dep, aes(x = reorder(gene, z), y = z)) +
  geom_point(size=3, colour="gray89") + 
  geom_point(size=3, colour = "gray50", pch = 21) + 
  geom_point(data=clk1_targets_scores_z15, colour="red", size = 3) +
  geom_point(data=clk1_targets_scores_z15, colour="black", size = 3, pch = 21) +
  geom_hline(yintercept = -1.5, linetype = "dashed") +
  labs(x = "Gene",
       y = "CRISPR Dependencey Z-score") +
  #geom_text_repel(data = mean_z_15 %>% filter(gene %in% clk1_targets_crispr), aes(label = gene), color = "black") +
  geom_text_repel(data = crispr_dep %>% filter(gene %in% clk1_targets_crispr), 
                  aes(label = gene), color = "black", 
                  nudge_y = 0.7, # Adjust the nudging value to avoid overlap
                  box.padding = 0.5, # Add padding around the label
                  point.padding = 0.5, # Add padding around the point
                  segment.size = 0.2, # Set the size of the line segment connecting the label to the point
                  max.overlaps = Inf) + # Allow for infinite overlaps, which ggrepel will handle
  facet_wrap(~sample_id,nrow = 2, scales = "free_x") +
  theme_Publication() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  ) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 0.6)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 0.5)
crispr_scores_z_plot

# Save plot as pdf
pdf(crispr_score_sign_plot_file, height = 4, width = 7.5)
crispr_scores_z_plot
dev.off()
