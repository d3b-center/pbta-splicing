################################################################################
# 01-plot_ont_vs_short_CLK1-Ex4_psi.R
# Plotting script that takes in Exon 4 PSI from short vs ONT reads
# written by Ammar Naqvi
#
# usage: Rscript 01-plot_ont_vs_short_CLK1-Ex4_psi.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("viridis")
  library("RColorBrewer")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "long-read-CLK1-validation")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output plot path
file_barplot = file.path(plots_dir,"isoform-stacked-barplot.pdf")

## retrive stringtie2 results for each cell line
cl_7316_1763_file = file.path(input_dir,"7316_1763.CLK1.processed.txt")
cl_7316_1769_file = file.path(input_dir,"7316_1769.CLK1.processed.txt")
cl_KNS42_file     = file.path(input_dir,"KNS42.CLK1.processed.txt")

## format and process cell line information 
cl_7316_1763_df <- vroom(cl_7316_1763_file,comment = "#", 
                         delim = "\t", 
                         col_names = c("transcript","tpm"), show_col_types = FALSE) %>% 
                   dplyr::mutate(cell_line="7316_1763") %>%
                   dplyr::filter(transcript=='STRG.1.2'| transcript=='STRG.1.1' ) %>% 
                   dplyr::mutate(Isoform = case_when(
                                transcript == "STRG.1.2" ~ "Inclusion",
                                TRUE ~ "Skipping")) %>% 
                  dplyr::mutate(type="long") %>% 
                  dplyr::mutate(PSI = (tpm / sum(tpm)) * 100 ) %>% 
                  # manual add rMATs results for short reads
                  add_row( cell_line = "7316_1763", type = "short", Isoform="Skipping", PSI = 59.8) %>%
                  add_row( cell_line = "7316_1763", type = "short", Isoform="Inclusion", PSI = 40.2) 
  
cl_7316_1769_df <- vroom(cl_7316_1769_file,comment = "#", 
                         delim = "\t", 
                         col_names = c("transcript","tpm"), show_col_types = FALSE) %>% 
                    mutate(cell_line="7316_1769") %>% 
                    dplyr::filter(transcript=='STRG.2.2'| transcript=='STRG.2.1' ) %>% 
                    dplyr::mutate(Isoform = case_when(
                                  transcript == "STRG.2.1" ~ "Inclusion",
                                  TRUE ~ "Skipping")) %>% 
                    dplyr::mutate(type="long") %>% 
                    dplyr::mutate(PSI = (tpm / sum(tpm)) * 100 ) %>% 
                    # manual add rMATs results for short reads
                    add_row( cell_line = "7316_1769", type = "short", Isoform="Skipping", PSI = 14.9) %>%
                    add_row( cell_line = "7316_1769", type = "short", Isoform="Inclusion", PSI = 85.1) 

cl_KNS42_df <- vroom(cl_KNS42_file,comment = "#", 
                         delim = "\t", 
                         col_names = c("transcript","tpm"), show_col_types = FALSE) %>% 
                         mutate(cell_line="KNS42")  %>% 
                         dplyr::mutate(Isoform = case_when(
                                       transcript == "STRG.1.2" ~ "Inclusion",
                                       TRUE ~ "Skipping")) %>% 
                        dplyr::mutate(type="long") %>% 
                        dplyr::mutate(PSI = (tpm / sum(tpm)) * 100 ) %>% 
                        # manual add rMATs results for short reads
                        add_row( cell_line = "KNS42", type = "short", Isoform="Skipping", PSI = 12.9) %>%
                        add_row( cell_line = "KNS42", type = "short", Isoform="Inclusion", PSI = 87.1) 


cell_lines_df <- rbind(cl_KNS42_df,cl_7316_1769_df,cl_7316_1763_df)

## make plot
stacked_barplot <- ggplot(cell_lines_df, aes(fill=Isoform, y=PSI, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#FFC20A","#0C7BDC")) +
  facet_wrap(~cell_line) +
  xlab("RNA-Seq sequencing strategy") + 
  theme_Publication()

pdf(file_barplot, 
    width = 7, height = 7)
stacked_barplot
dev.off()




