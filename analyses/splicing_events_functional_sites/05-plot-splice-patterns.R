################################################################################
# 05-plot-splice-patterns.R
# Script that plots a patterns of splicing events across total sites. 
# written by Ammar Naqvi
#
# usage: Rscript 05-plot-splice-patterns.R
################################################################################

## libraries needed
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_events_functional_sites")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

##theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output file names for plots
file_splice_pattern_plot <- file.path(analysis_dir, "plots", "splicing_pattern_plot.pdf")

file_psi <- file.path(results_dir,"splice_events.diff.SE.HGG.txt")
file_psi_func_incl <- file.path(results_dir,"splicing_events.total.HGG.pos.intersectUnip.ggplot.txt") 
file_psi_func_skip <-file.path(results_dir,"splicing_events.total.HGG.neg.intersectUnip.ggplot.txt") 

psi_tab      <-  read_tsv(file.path(file_psi)) 
psi_skip <- psi_tab %>% dplyr::filter(Type=="Skipping")
psi_incl <- psi_tab %>% dplyr::filter(Type=="Inclusion")

psi_skip_func <- read_tsv(file.path(file_psi_func_incl)) 
psi_incl_func <- read_tsv(file.path(file_psi_func_skip)) 


## find functional sites not mixed events or unidirectional
mixed_events_func_df <- inner_join(psi_incl_func,psi_skip_func, by="SpliceID",relationship = "many-to-many") %>%
  mutate(type="Mixed") %>% 
  dplyr::select("SpliceID", type) %>% 
  unique() %>% 
  mutate(Impact='Functional')

unidirectional_skip_events_func_df <-  psi_skip_func %>% 
  unique() %>% 
  mutate(type="Skipping") %>% 
  dplyr::select("SpliceID", type) %>%
  mutate(Impact='Functional')

unidirectional_incl_events_func_df <-  psi_incl_func %>% 
  unique() %>% 
  mutate(type="Inclusion") %>% 
  dplyr::select("SpliceID", type) %>% 
  mutate(Impact='Functional') 

## total non-functional
mixed_events_df <- inner_join(psi_incl,psi_skip, by="Splice ID") %>%
  mutate(type="Mixed") %>% 
  dplyr::select("Splice ID", type) %>% 
  dplyr::rename("SpliceID"="Splice ID") %>% 
  unique()   %>% 
  anti_join(mixed_events_func_df, by='SpliceID') %>%
  mutate(Impact='Non-functional')

unidirectional_skip_events_df <- anti_join(psi_skip,psi_incl, by="Splice ID") %>%
  mutate(type="Skipping") %>% 
  dplyr::select("Splice ID", type) %>% 
  dplyr::rename("SpliceID"="Splice ID") %>% 
  unique() %>% 
  anti_join(unidirectional_skip_events_func_df, by='SpliceID') %>%
  mutate(Impact='Non-functional')

unidirectional_incl_events_df <- anti_join(psi_incl,psi_skip, by="Splice ID") %>%
  mutate(type="Inclusion") %>% 
  dplyr::select("Splice ID", type) %>% 
  dplyr::rename("SpliceID"="Splice ID") %>% 
  unique() %>% 
  anti_join(unidirectional_incl_events_func_df, by='SpliceID') %>%
  mutate(Impact='Non-functional')

## combine into master table
psi_events_total <- rbind(mixed_events_df,
                          unidirectional_skip_events_df,
                          unidirectional_incl_events_df,
                          mixed_events_func_df,
                          unidirectional_skip_events_func_df,
                          unidirectional_incl_events_func_df)


plot_pattern <- ggplot(psi_events_total,
                       aes(x = type, fill= Impact)) +
  geom_bar(position="fill",stat="count", color="black")    + 
  scale_fill_manual(values=c("#FFC20A","#0C7BDC"),
                    name = "Predicted Impact")      + 
  theme(legend.position = "none",
        legend.title = element_text(size=19), 
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14)) + 
  xlab("Splicing Pattern") + ylab("Fraction of Variants") + 
  geom_text(stat='count',aes(label=..count..), position = position_fill(vjust = 0.5)) +
  theme_Publication() 

  
# Save plot as PDF
pdf(file_splice_pattern_plot, 
    width = 6, height = 4)
plot_pattern
dev.off()