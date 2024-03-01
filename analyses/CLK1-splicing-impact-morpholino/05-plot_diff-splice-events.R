################################################################################
# 05-plot_diff-splice-events.R
# written byAmmar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 03-plot_diff-splice-events.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggrepel")
  library("vroom")
  library("ggpubr")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories and file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")

results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}


##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## define output files
file_dpsi_plot <- file.path(plots_dir,"dPSI-distr-func.pdf")

## get and setup input

## retrieve psi values from tables
file_psi_func <- file.path(results_dir,"splicing_events.morpho.intersectUnip.ggplot.txt")

## read table of recurrent functional splicing (skipping)
dpsi_unip_incl <- vroom(file_psi_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  filter(dPSI<0) %>% 
  mutate(Preference='Inclusion',
         dPSI=abs(dPSI))

dpsi_unip_skp <- vroom(file_psi_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  filter(dPSI>0) %>% 
  mutate(Preference='Skipping')


psi_comb <- rbind(dpsi_unip_incl,dpsi_unip_skp) %>% 
  mutate(Uniprot = case_when(Uniprot == 'DisulfBond' ~ "Disulfide Bond",
                             Uniprot == 'LocSignal' ~ "Localization Signal",
                             .default = Uniprot),
         Uniprot_wrapped = stringr::str_wrap(Uniprot, width = 10)
  )


## ggstatplot across functional sites
set.seed(123)
counts_psi_comb <- psi_comb %>% 
  count(Preference, Uniprot_wrapped)
plot_dsp <-  ggplot(psi_comb, aes(Uniprot_wrapped, dPSI*100) ) +  
  ylab(expression(bold("dPSI"))) +
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 5, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  facet_wrap("Preference") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Disulfide\nBond", "Localization\nSignal"),
                                                                c("Disulfide\nBond", "Modifications"),
                                                                c("Disulfide\nBond", "Other"),
                                                                c("Localization\nSignal", "Modifications"),
                                                                c("Localization\nSignal", "Other"),
                                                                c("Modifications", "Other"))) + 
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A"))  + 
  theme_Publication() + 
  
  labs(y="Percent Spliced In (PSI)", x= "Uniprot-defined Functional Site") + 
  geom_text(data = counts_psi_comb, aes(label = paste("n =",n), x = Uniprot_wrapped, y = 0), vjust = 3, size = 4, hjust=.5) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +  # Angles x-axis text
  ylim(c(-20,170))

# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 8, height = 5)
print (plot_dsp)
dev.off()
