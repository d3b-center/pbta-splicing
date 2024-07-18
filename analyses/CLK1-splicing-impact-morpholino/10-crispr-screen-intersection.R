################################################################################
# 10-crispr-screen-intersection.R
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
  library('ggVennDiagram')
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
ds_genes_file <- file.path(results_dir,"differential_splice_by_goi_category.tsv")
de_genes_file <- file.path(results_dir, "de_genes.tsv")
clk1_target_file <- file.path(results_dir,"common_genes_de_ds_functional.txt")
categories_file <- file.path(results_dir, "gene_categories.tsv")

## ouput files
clk1_crispr_file <- file.path(results_dir,"clk1-de-ds-crispr-targets.txt")
gene_dep_file <- file.path(results_dir,"crispr-dependencies.txt")
venn_output_file <- file.path(results_dir, "ds-de-crispr-venn.pdf")
  
# read in clinical
clin_df <- read_tsv(clin_df, guess_max = 100000) %>%
  dplyr::filter(short_histology == "HGAT",
                cohort == "PBTA",
                composition == "Solid Tissue")

clk1_targets <- read_lines(clk1_target_file)
categories <- read_tsv(categories_file)

crispr_df <- read_csv(crispr_score) 
ds_genes <- read_tsv(ds_genes_file) %>%
  dplyr::rename(geneSymbol = gene) %>%
  filter(geneSymbol %in% categories$gene) %>%
  dplyr::select(SpliceID, geneSymbol, dPSI, Preference)
de_genes <- read_tsv(de_genes_file) %>%
  filter(geneSymbol %in% categories$gene) %>%
  dplyr::select(geneSymbol, baseMean, Preference)


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

crispr_dep_genes <- crispr_dep %>%
  filter(z < -1.5) %>%
  pull(gene) %>%
  unique() %>%
  write_lines(gene_dep_file)

clk1_targets_crispr <- intersect(clk1_targets, crispr_dep$gene)

# write list of targets in any cell line
unique(clk1_targets_crispr) %>%
  write_lines(clk1_crispr_file)

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
  
  if (each == "pHGG") {
    y_text <- "Dependency Mean Z-score"
  }
  if (each != "pHGG") {
    y_text <- "Dependency Z-score"
  }
  
crispr_scores_z_plot <- ggplot(crispr_dep_each, aes(x = reorder(gene, z), y = z)) +
  geom_point(size=3, colour="gray89") + 
  geom_point(size=3, colour = "gray50", pch = 21) + 
  geom_point(data=clk1_targets_crispr_intersect, colour="red", size = 3) +
  geom_point(data=clk1_targets_crispr_intersect, colour="black", size = 3, pch = 21) +
  geom_hline(yintercept = -1.5, linetype = "dashed", color = "gray50") +
  labs(title = paste0(each),
       x = "Gene",
       y = paste0(y_text)) +
  geom_text_repel(data = crispr_dep_each %>% filter(gene %in% clk1_targets_crispr_intersect$gene), 
                  aes(label = gene), color = "black", 
                  nudge_y = 0.7, # Adjust the nudging value to avoid overlap
                  box.padding = 0.5, # Add padding around the label
                  point.padding = 0.5, # Add padding around the point
                  segment.size = 0.2, # Set the size of the line segment connecting the label to the point
                  max.overlaps = Inf) + # Allow for infinite overlaps, which ggrepel will handle
  theme_Publication() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  )

# Save plot as pdf
pdf(file.path(paste0(plots_dir, "/clk1-crispr-swoosh-", each, ".pdf")), height = 4, width = 6)
print(crispr_scores_z_plot)
dev.off()

}


## plot venn diagram
venn_diag<- ggVennDiagram(x=list(de_genes$geneSymbol, unique(ds_genes$geneSymbol),crispr_dep_genes), 
                          edge_lty = "dashed", 
                          edge_size = 1,
                          label_size = 5,
                          set_size = 4,
                          category.names = c("DE" , "DS", "DG"),
                          label_percent_digit = 1) +  
  scale_fill_distiller(palette = "Blues", direction = 1, name = expression(bold("Gene count"))) + 
  #labs(title = expression(bold("Differentially expressed and spliced genes and dependent genes"))) +
  coord_flip()

ggplot2::ggsave(venn_output_file,
                plot=venn_diag,
                width=5.5,
                height=4,
                device="pdf",
                dpi=300)

# gather 3 way list to inspect
mean_data_renamed <- mean_data %>%
  dplyr::rename(geneSymbol = gene)
  
all_events <- ds_genes %>%
  full_join(de_genes, by='geneSymbol', relationship = "many-to-many",
            suffix = c("_psi", "_de")) %>%
  dplyr::select(SpliceID, geneSymbol, dPSI, Preference_psi, Preference_de) %>% 
  unique() %>%
  left_join(mean_data_renamed) %>%
  write_tsv(file.path(results_dir, "ds-de-crispr-events-functional-full.tsv"))

all_intersected <- all_events %>%
  filter(!is.na(z),
         !is.na(Preference_de),
         !is.na(Preference_psi)) %>%
  write_tsv(file.path(results_dir, "ds-de-crispr-events-functional-direction.tsv"))

