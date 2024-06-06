################################################################################
# 02-plot_ctg-assay.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi
#
# usage: Rscript 02-plot_ctg-assay.R
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-inhibition-assays")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output for plot
file_plot = file.path(plots_dir,"cell_viability-barplot.pdf")

## input file
ctg_res_file <- file.path(input_dir,"kns42-ctg-3day.txt")

# Read the input file
ctg_data <- read.delim(ctg_res_file, header = TRUE, stringsAsFactors = FALSE)

ctg_tidy_df <- ctg_data %>%
  pivot_longer(cols = -1, names_to = "Treatment", values_to = "Measurement") %>%
  rename(Replicates = X) %>%
  mutate(Treatment = str_replace(Treatment, "X", ""), 
         Treatment = sub("\\.u", "u", Treatment) 
  ) %>%
  filter(!grepl("Free", Treatment)) 


# Calculate standard error for each Treatment
df_summary <- ctg_tidy_df %>%
  group_by(Treatment) %>%
  summarise(mean_measurement = mean(Measurement),
            se_measurement = sd(Measurement) / sqrt(n()))

# Perform t-tests to compare each treatment to the control (DMSO)
control_group <- ctg_tidy_df %>% filter(Treatment == "DMSO")

p_values <- ctg_tidy_df %>%
  group_by(Treatment) %>%
  summarise(p_value = t.test(Measurement, control_group$Measurement)$p.value) %>%
  mutate(p_value = case_when(p_value < 0.05 ~ paste0("*", "p=", format(p_value, scientific = TRUE, digits = 2)),
                           TRUE ~ "")) 
  

barplot <- ggplot(df_summary, aes(x = Treatment, y = mean_measurement)) +
  geom_bar(stat = "identity", fill = "#0C7BDC") +
  geom_errorbar(aes(ymin = mean_measurement - se_measurement,
                    ymax = mean_measurement + se_measurement),
                position = position_dodge(width = 18),
                width = 1, color = "black") +
  geom_text(data = p_values, aes(x = Treatment, y = max(ctg_tidy_df$Measurement+1000), label = p_value), vjust = -0.5, size = 3, color = "black") +
  
 
  labs(x = "Concentration", y = "Luminescence (RLU)") +
  theme_Publication()


pdf(file_plot, width = 8, height = 6)
print(barplot)
dev.off()