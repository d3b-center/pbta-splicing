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
ctg_res_file <- file.path(input_dir,"2024.05.30_KNS-42_3D_cirtuvivint_nogroups.txt")

ctg_data <- read_tsv(ctg_res_file, skip=10, col_names = TRUE) 

# Convert Date Time column to datetime
ctg_data <- ctg_data %>%
  mutate(`Date Time` = as.POSIXct(`Date Time`, format = "%m/%d/%Y %I:%M:%S %p"))

# Reshape data into long format for plotting
ctg_data_long <- ctg_data %>%
  pivot_longer(cols = -c(`Date Time`, Elapsed), names_to = "Treatment", values_to = "Measurement") %>%
  mutate(`Date Time` = format(`Date Time`, "%H:%M:%S")) %>%
  dplyr::rename("Time"=`Date Time`) %>%
  filter(!grepl("Std", Treatment)) %>%
  
  mutate(Treatment = sub(".*(?=CIrtuvivint)", "", Treatment, perl = TRUE),
         Treatment = sub("\\s*\\([^\\)]+\\)", "", Treatment),
         Treatment = str_replace(Treatment, "KNS-42 12K / well ", ""),
         Treatment = sub("\\s*\\([^\\)]+\\)", "", Treatment) 
  )

# Subset the data to ctg_data_long only the two treatments of interest for stats
filtered_df <- ctg_data_long %>%
  filter(Treatment %in% c("DMSO Vehicle 2%", "CIrtuvivint 10 µM", "CIrtuvivint 1 µM")) %>%
  mutate(Treatment = as.factor(Treatment),
         Treatment = fct_relevel(Treatment, "DMSO Vehicle 2%", "CIrtuvivint 10 µM", "CIrtuvivint 1 µM"))

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
  group_by(Elapsed) %>%
  t_test(Measurement ~ Treatment, paired = TRUE) %>%  # Adjust for your study design
  mutate(p_sig = case_when(p < 0.05 ~ paste0("*", "p=", round(p, 2)),
                           TRUE ~ "")
  )

barplot <- ggplot(filtered_df, aes(x = Elapsed, y = Measurement, fill = Treatment)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.8), width = 0.25) +
  #geom_text(data = stat_results, aes(x = Elapsed, label = p_sig), vjust = -0.5, size = 4) +
  xlab("Time (hours)") + 
  #ylab("Luminescence (RLU)") +
  #scale_y_continuous(labels = scales::scientific) +
  #scale_fill_manual(values = c("lightgrey", "#0C7BDC")) +
  #ylim(c(0, 6e4))+
  theme_Publication() 


pdf(file_plot, width = 8, height = 4)
print(barplot)
dev.off()
