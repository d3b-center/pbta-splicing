################################################################################
# 02-plot_cell-viability-assay-res.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 02-plot_cell-viability-res.R 
################################################################################

## libraries used 
suppressPackageStartupMessages({
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
file_line_plot = file.path(plots_dir,"cell_viability-barplot.pdf")

## input file
cell_prolif_res_file <- file.path(input_dir,"cell_prolif_res.tsv")

cell_prolif_df <- read_tsv(cell_prolif_res_file) %>% 
  dplyr::select(Time, Ctrl_1, Ctrl_2, Ctrl_3, CLK1_1, CLK1_2, CLK1_3) %>%
  pivot_longer(
    c(-Time), 
    names_pattern = "(.+)_([0-9])", 
    names_to = c("Treatment", "Rep"),
    values_to = "Absorbance"
  ) %>%
  mutate(Treatment = case_when(Treatment == "CLK1" ~ "CLK1 Exon 4 morpholino",
                               Treatment == "Ctrl" ~ "Non-targeting morpholino"),
                                Time = as.factor(Time)) 

# Subset the data to include only the two treatments of interest for stats
filtered_df <- cell_prolif_df %>%
  filter(Treatment %in% c("CLK1 Exon 4 morpholino", "Non-targeting morpholino")) %>%
  mutate(Treatment = as.factor(Treatment),
         Treatment = fct_relevel(Treatment, "Non-targeting morpholino", "CLK1 Exon 4 morpholino"))

# Define mean_se function if not already defined
mean_se <- function(x) {
  n <- length(x)
  mean <- mean(x)
  sd <- sd(x)
  ymin <- mean - sd
  ymax <- mean + sd
  
  return(c(y = mean, ymin = ymin, ymax = ymax))
}

# Perform statistical tests using rstatix
stat_results <- filtered_df %>%
  group_by(Time) %>%
  t_test(Absorbance ~ Treatment, paired = TRUE) %>%  # Adjust for your study design
  mutate(p_sig = case_when(p < 0.05 ~ paste0("*", "p=", round(p, 2)),
                           TRUE ~ ""),
         y_pos = c(40000,42000,44000, 46000, 56000),
         #y_pos = max(filtered_df$Absorbance, na.rm = TRUE) + 1000,
         Treatment = "Non-targeting morpholino")

# Generate the bar plot
barplot <- ggplot(filtered_df, aes(x = Time, y = Absorbance, fill = Treatment)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.8), width = 0.25) +
  geom_text(data = stat_results, aes(x = Time, y = y_pos, label = p_sig), vjust = -0.5, size = 4) +
  xlab("Time (hours)") + 
  ylab("Luminescence (RLU)") +
  scale_y_continuous(labels = scales::scientific) +
  scale_fill_manual(values = c("lightgrey", "#0C7BDC")) +
  ylim(c(0, 6e4))+
  theme_Publication() 

pdf(file_line_plot, width = 6.5, height = 3)
print(barplot)
dev.off()






