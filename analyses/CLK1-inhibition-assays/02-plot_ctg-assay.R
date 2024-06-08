################################################################################
# 02-plot_ctg-assay.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 02-plot_ctg-assay.R
################################################################################
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-inhibition-assays")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output for plot
day3_plot = file.path(plots_dir,"cell_viability-barplot-3day.pdf")
day6_plot = file.path(plots_dir,"cell_viability-barplot-6day.pdf")

## input file
day3_file <- file.path(input_dir,"kns42-ctg-3day.tsv")
day6_file <- file.path(input_dir,"kns42-ctg-6day.tsv")

# Define mean_se function if not already defined
mean_se <- function(x) {
  n <- length(x)
  mean <- mean(x)
  sd <- sd(x)
  ymin <- mean - sd
  ymax <- mean + sd
  
  return(c(y = mean, ymin = ymin, ymax = ymax))
}


for (eachfile in c(day3_file, day6_file)) {
  
  if (eachfile == day3_file) {
    plot_file <- day3_plot
  }
  if (eachfile == day6_file) {
    plot_file <- day6_plot
  }
   
  via_df <- read_tsv(eachfile) %>%
    dplyr::rename(Rep = `...1`) %>%
    dplyr::mutate(Rep = row_number()) %>%
    pivot_longer(cols = -1, names_to = "Dose", values_to = "Absorbance") %>%
    # Subset the data to include treatments of interest for stats
    filter(!grepl("Free", Dose)) %>%
    #separate(Rep, into = c("group", "rep"), sep = "-") %>%
    mutate(Dose = as.factor(gsub(" ", "", Dose))) %>%
    mutate(Dose = fct_relevel(Dose, c("DMSO", "0.01uM", "0.05uM", "0.1uM", "0.5uM", "1uM", "5uM", "10uM")))
  
  # Perform statistical tests using rstatix
  stat_results <- via_df %>%
    filter(Dose != "DMSO") %>%
    group_by(Dose) %>%
    do(tidy(t.test(Absorbance ~ Dose, data = bind_rows(via_df %>% filter(Dose == "DMSO"), .)))) %>%
    adjust_pvalue(method = "bonferroni") %>%
    mutate(p_sig = case_when(p.value.adj < 0.05 ~ paste0("*", "p=",
                                                         format(p.value.adj, scientific = TRUE, digits = 2)),
                             TRUE ~ ""),
           star = case_when(p.value.adj < 0.05 ~ paste0("*"),
                             TRUE ~ ""),
           y_pos = seq(max(via_df$Absorbance, na.rm = TRUE) + 500000, by = -500000, length.out = n())) # Stagger the p-values
  
  # Generate the bar plot
  barplot <- ggplot(via_df, aes(x = Dose, y = Absorbance, fill = Dose)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.8), width = 0.25) +
    geom_text(data = stat_results, aes(x = Dose, y = y_pos, label = star), vjust = -0.5, size = 7) +
    geom_text(aes(x = "5uM", y = 9e6, label = "*p<0.001"), vjust = -0.5, size = 4) +
    xlab("Treatment and Dose") + 
    ylab("Luminescence (RLU)") +
    scale_y_continuous(labels = scales::scientific) +
    scale_fill_manual(values = c("lightgrey", "#0C7BDC", "#0C7BDC", "#0C7BDC", "#0C7BDC", "#0C7BDC", "#0C7BDC", "#0C7BDC")) +
    ylim(c(0, 10e6))+
    theme_Publication() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) # Remove the legend
  
  pdf(plot_file, width = 4.5, height = 4)
  print(barplot)
  dev.off()
}

