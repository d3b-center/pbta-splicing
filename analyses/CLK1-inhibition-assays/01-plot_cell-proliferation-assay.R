################################################################################
# 01-plot_cell-proliferation-assay.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi
#
# usage: Rscript 01-plot_cell-proliferation-assay.R
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
file_line_plot = file.path(plots_dir,"cell-prolif.pdf")

## input file
cell_prolif_res_file <- file.path(input_dir,"2024.05.30_KNS-42_3D_cirtuvivint.txt")
cell_prolif_res_file <- file.path(input_dir,"2024.05.30_KNS-42_3D_cirtuvivint_nogroups.txt")

incucyte_data <- read_tsv(cell_prolif_res_file, skip=10, col_names = TRUE) 

# Convert Date Time column to datetime
incucyte_data <- incucyte_data %>%
  mutate(`Date Time` = as.POSIXct(`Date Time`, format = "%m/%d/%Y %I:%M:%S %p"))

# Reshape data into long format for plotting
incucyte_data_long <- incucyte_data %>%
  pivot_longer(cols = -c(`Date Time`, Elapsed), names_to = "Treatment", values_to = "Measurement") %>%
  mutate(`Date Time` = format(`Date Time`, "%H:%M:%S")) %>%
  dplyr::rename("Time"=`Date Time`) %>%
  filter(!grepl("Std", Treatment)) %>%
  
  mutate(Treatment = sub(".*(?=CIrtuvivint)", "", Treatment, perl = TRUE),
         Treatment = sub("\\s*\\([^\\)]+\\)", "", Treatment),
         Treatment = str_replace(Treatment, "KNS-42 12K / well ", ""),
         Treatment = sub("\\s*\\([^\\)]+\\)", "", Treatment) 
         )


means <- incucyte_data_long %>%
  group_by(Elapsed, Treatment) %>%
  summarise(mean_measurement = mean(Measurement),
            sd_measurement = sd(Elapsed),
            n = n())


# Perform statistical tests


# plot
# Plot the means with error bars as a line graph
plot <- ggplot(means, aes(x = Elapsed, y = mean_measurement, group = Treatment, color = Treatment)) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_measurement - sd_measurement/sqrt(n),
                    ymax = mean_measurement + sd_measurement/sqrt(n)),
                width = 0.2) +
  labs(x = "Elapsed Time", y = "KNS42 Confluence %") +
  theme_Publication()


pdf(file_line_plot, width = 8, height = 4)
print(plot)
dev.off()
