################################################################################
# 02-plot_cell-proliferation-assay-res.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 02-plot_cell-proliferation-assay-res.R 
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("vroom")
  library("tidyverse")
  library("ggpubr")
  library("rstatix")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "KNS42-cell-line")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output for plot
file_line_plot = file.path(plots_dir,"cell_prolif-line.pdf")

## input file
cell_prolif_res_file <- file.path(input_dir,"cell_prolif_res.tsv")

cell_prolif_df <- vroom(cell_prolif_res_file, delim = "\t", col_names = TRUE) %>% 
  pivot_longer(
    c(-Time), 
    names_pattern = "(.+)_([0-9])", 
    names_to = c("Treatment", "Rep"),
    values_to = "Absorbance"
  ) %>%
  mutate(Treatment = case_when(Treatment == "CLK1" ~ "CKL1 Exon 4 morpholino",
                               Treatment == "Ctrl" ~ "Non-targeting morpholino",
                               Treatment == "Untreated" ~ "No morpholino"),
         Time = as.factor(Time)) 

# Subset the data to include only the two treatments of interest for stats
filtered_df <- cell_prolif_df %>%
  filter(Treatment %in% c("CKL1 Exon 4 morpholino", "Non-targeting morpholino"))

# Perform statistical tests
stat_results <- filtered_df %>%
  group_by(Time) %>%
  t_test(Absorbance ~ Treatment, paired = TRUE) %>%  # Adjust as per your study design
  mutate(p_sig = case_when(p < 0.05 ~ paste0("*", round(p, 2)),
                          TRUE ~ ""),
         y_pos = 20000,
         x_pos = Time)
# plot
ribbon_plot <- ggplot(cell_prolif_df, aes(x = Time, y = Absorbance)) + 
  stat_summary(aes(colour = Treatment, group = Treatment), fun = mean, geom = "line") +
  stat_summary(aes(fill = Treatment, group = Treatment), fun.data = mean_sd, geom = "ribbon", alpha = 0.25) +
  geom_text(data = stat_results, aes(x = Time, y = y_pos, label = p_sig), size = 3.5) + # Add p-values
  xlab("Time") + 
  ylab("Absorbance") +
  ylim(c(10000,60000)) +
  theme_Publication()

ribbon_plot

pdf(file_line_plot, width = 6.5, height = 3)
print(ribbon_plot)
dev.off()

