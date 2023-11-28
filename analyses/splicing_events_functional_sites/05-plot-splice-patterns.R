################################################################################
# 05-plot-splice-patterns.R
# Script that plots a patterns of splicing events across total sites
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

## retrieve psi values from tables
file_psi_pos_func  <- "splicing_events.total.pos.intersectUnip.ggplot.txt"
file_psi_neg_func  <- "splicing_events.total.neg.intersectUnip.ggplot.txt"
file_psi_pos_total <- "splicing_events.total.pos.tsv"
file_psi_neg_total <- "splicing_events.total.neg.tsv"

psi_pos_tab      <-  read.delim(file.path(results_dir, file_psi_pos_total), sep = "\t", row.names = NULL, header=TRUE)
psi_neg_tab      <-  read.delim(file.path(results_dir, file_psi_neg_total), sep = "\t", row.names = NULL, header=TRUE)

psi_pos_func_tab      <-  read.delim(file.path(results_dir, file_psi_pos_func), sep = "\t", row.names = NULL, header=TRUE)
psi_neg_func_tab      <-  read.delim(file.path(results_dir, file_psi_neg_func), sep = "\t", row.names = NULL, header=TRUE)

mixed_events_df <- inner_join(psi_pos_tab,psi_neg_tab, by='gene') %>%
  mutate(type="Mixed") %>% 
  dplyr::select(gene, type) %>% 
  dplyr::rename("SpliceID"=gene)


psi_pos_non_mixed_df <-  psi_pos_tab %>% 
  filter(!gene %in% mixed_events_df$SpliceID) %>% 
  mutate(type = case_when(flip == 1 ~ "Flip",
                          flip == 0 ~ "Non-flip")) %>%
  dplyr::select(gene,type) %>% 
  dplyr::rename("SpliceID"=gene)

psi_neg_non_mixed_df <- psi_neg_tab %>% 
  filter(!gene %in% mixed_events_df$SpliceID) %>% 
  mutate(type = case_when(flip == 1 ~ "Flip",
                          flip == 0 ~ "Non-flip")) %>%
  dplyr::select(gene,type) %>% 
  dplyr::rename("SpliceID"=gene)


psi_anno_df <- rbind(mixed_events_df,psi_pos_non_mixed_df,psi_neg_non_mixed_df) %>% 
  mutate(Impact="Non-functional")

psi_anno_func_df <- rbind(psi_pos_func_tab,psi_neg_func_tab) %>%
  dplyr::select(SpliceID)

psi_anno_functional_df <- inner_join(psi_anno_df,psi_anno_func_df, by="SpliceID") %>% distinct() %>% 
  mutate(Impact="Functional")

psi_anno_df <- psi_anno_df %>% anti_join(psi_anno_functional_df, by='SpliceID')

psi_anno_total <- rbind(psi_anno_df,psi_anno_functional_df) 

plot_pattern <- ggplot(psi_anno_total,
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


