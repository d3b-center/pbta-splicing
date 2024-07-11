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
  library("data.table")
  library("grid")
  library("gridExtra")
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
all_hist_file <- file.path(data_dir, "histologies.tsv")
pbta_hist_file <- file.path(data_dir, "histologies-plot-group.tsv")
pbta_tpm_file <- file.path(data_dir, "rna-isoform-expression-rsem-tpm.rds") 
gtex_trans_tpm_file <- file.path(data_dir, "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz")
<<<<<<< HEAD
control_tpm_file <- file.path(data_dir, "control-rna-isoform-expression-rsem-counts-tpm.rds") 
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

## output file for plot
gtex_plot_path <- file.path(plots_dir, "CLK1ex4-hgg-normals.pdf")
=======
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")

## output file for plot
gtex_plot_path <- file.path(plots_dir, "CLK1ex4-hgg-gtex.pdf")
>>>>>>> main

## load files
indep_df <- read_tsv(indep_file) %>% 
  dplyr::filter(cohort=='PBTA')

hist <- read_tsv(all_hist_file, guess_max = 100000)

all_hgg_hist <- read_tsv(pbta_hist_file) %>% 
  dplyr::filter(plot_group %in% c("Other high-grade glioma", "DIPG or DMG"),
                Kids_First_Biospecimen_ID %in% 
                  indep_df$Kids_First_Biospecimen_ID)

gtex_brain <- hist %>% 
  dplyr::filter(cohort == "GTEx",
                gtex_group == "Brain")

# transcript is ENST00000321356.9 in pbta
hgg_counts <- readRDS(pbta_tpm_file) %>%
  select(c(transcript_id, all_hgg_hist$Kids_First_Biospecimen_ID)) %>%
  filter(transcript_id == "ENST00000321356.9") %>%
  select(-transcript_id) %>%
  # Transpose and then immediately convert to a data frame or tibble, capturing row names as a column
  {data.frame(t(.))} %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::rename(ENST00000321356 = `t...`) %>%
  left_join(all_hgg_hist[,c("Kids_First_Biospecimen_ID", "plot_group")])

# transcript is ENST00000321356.8 in gtex
gtex_clk1ex4_counts <- fread(gtex_trans_tpm_file, skip = 2) %>%
  select(c(transcript_id, gtex_brain$Kids_First_Biospecimen_ID)) %>%
  filter(transcript_id == "ENST00000321356.8") %>%
  select(-transcript_id) %>%
  {data.frame(t(.))} %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::rename(ENST00000321356 = `t...`) %>%
  left_join(gtex_brain[,c("Kids_First_Biospecimen_ID", "gtex_subgroup")]) %>%
  dplyr::rename(plot_group = gtex_subgroup) %>%
  # reduce labels
  dplyr::mutate(plot_group = gsub("Brain - ", "", plot_group))

<<<<<<< HEAD
control_counts <- readRDS(control_tpm_file) %>%
  filter(transcript_id == "ENST00000321356.9_CLK1-201",
         sample_id != "BM") %>%
  dplyr::mutate(plot_group = gsub("_", " ", sample_id)) %>%
  dplyr::rename(ENST00000321356 = TPM) %>%
  select(ENST00000321356, plot_group) %>%
  # add prefix
  dplyr::mutate(plot_group = ifelse(plot_group != "Fetal brain", paste("Pediatric", plot_group, sep = " "),
                                    plot_group))
  
combined_df <- hgg_counts %>%
  bind_rows(gtex_clk1ex4_counts, control_counts)
=======
combined_df <- hgg_counts %>%
  bind_rows(gtex_clk1ex4_counts)
>>>>>>> main

# Calculate medians and reorder plot_group levels based on these medians
combined_df_ordered <- combined_df %>%
  group_by(plot_group) %>%
  mutate(median_CLK1 = median(ENST00000321356, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(plot_group = factor(plot_group, levels = unique(plot_group[order(median_CLK1)])))

# color palette
pal <- combined_df %>%
  select(plot_group) %>%
  unique() %>%
  left_join(unique(all_hgg_hist[,c("plot_group", "plot_group_hex")])) %>%
  mutate(plot_group_hex = case_when(is.na(plot_group_hex) & plot_group %in% gtex_clk1ex4_counts$plot_group ~ "gray93",
<<<<<<< HEAD
                                    is.na(plot_group_hex) & plot_group %in% control_counts$plot_group ~ "lightskyblue1",
=======
>>>>>>> main
                                    TRUE ~ plot_group_hex))

colors <- setNames(pal$plot_group_hex, pal$plot_group)

# plot
gtex_plot <- ggplot(combined_df_ordered, aes(x = plot_group, y = ENST00000321356)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = plot_group), width = 0.2, size = 2, shape = 21, color = "black") + # Add actual data points
  scale_fill_manual(values = colors) +
  labs(title = "",
       x = "Tissue",
       y = expression(bold(bolditalic("CLK1 ENST00000321356")*" TPM (log"[2]*")"))) +
<<<<<<< HEAD
       theme_Publication() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 75, hjust = 1)) +
  scale_x_discrete(labels = function(x) sapply(x, function(l) str_wrap(l, width = 20))) # Wrap x-axis labels
=======
  theme_Publication() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 75, hjust = 1)) +
  scale_x_discrete(labels = function(x) sapply(x, function(l) str_wrap(l, width = 30))) # Wrap x-axis labels
>>>>>>> main


# Define legend colors
legend_colors <- c("DIPG or DMG" = "#ff40d9", 
                   "Other high-grade glioma" = "#ffccf5", 
<<<<<<< HEAD
                   "GTEx Brain" = "gray93", 
                   "Pediatric Brain" = "lightskyblue1")

# Manual legend creation
legend_grob <- grid::grid.legend(draw = TRUE, labels = c("Tissue", "DIPG or DMG", "Other high-grade glioma", "GTEx Brain", "Pediatric Brain"),
                                 pch = 21,
                                 gp = gpar(col = c("white", "black","black","black","black","black"), 
                                           fill = c("white", "#ff40d9", "#ffccf5", "gray93", "lightskyblue1"), 
                                           fontsize = 10))

pdf(gtex_plot_path, height = 8, width = 11)
grid.arrange(gtex_plot, legend_grob, ncol= 2, widths=c(8,2), heights=c(4,1))
dev.off()

=======
                   "GTEx Brain" = "gray93")

# Manual legend creation
legend_grob <- grid::grid.legend(draw = TRUE, labels = c("Tissue", "DIPG or DMG", "Other high-grade glioma", "GTEx Brain"),
                                 pch = 21,
                                 gp = gpar(col = c("white", "black","black","black","black"), 
                                           fill = c("white", "#ff40d9", "#ffccf5", "gray93"), 
                                           fontsize = 10))

pdf(gtex_plot_path, height = 8, width = 9)
grid.arrange(gtex_plot, legend_grob, ncol= 2, widths=c(8,2), heights=c(4,1))
dev.off()



>>>>>>> main
