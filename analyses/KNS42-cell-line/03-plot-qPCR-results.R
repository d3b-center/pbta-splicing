################################################################################
# 03-plot-qPCR-res.R
# script that generates barplot of fold change for qRT-PCR
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 03-plot-qPCR-res.R
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("tidyverse")
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
source(file.path(figures_dir,"theme_for_plots.R"))

## output for plot
file_plot = file.path(plots_dir,"qPCR-morp.pdf")

qpcr_res_file = file.path(input_dir,"qpcr-results-raw-ct.csv")
qpcr_fc_df <- read_csv(qpcr_res_file) 

# pull out HPRT and average, this is the housekeeping gene
hprt_avg_ct <- qpcr_fc_df %>%
  filter(Primers == "HPRT") %>%
  select(CT) %>%
  colMeans(.)
  
# format and compute dct values from raw values
qpcr_dct <- qpcr_fc_df %>%
  # remove HPRT
  filter(Primers != "HPRT") %>%
  ## group by junction and calculate means
  group_by(Primers, Treatment) %>%
  summarise_all(mean) %>%
  ## calculate dct values
  mutate(dCT = CT - hprt_avg_ct) %>%
  arrange(Primers)

# calculate ddCT and Fold Change
control_34 <- qpcr_dct %>%
  filter(Treatment == "Control" & Primers == "Exons 3-4") %>%
  pull(dCT)

control_35_1 <- qpcr_dct %>%
  filter(Treatment == "Control" & Primers == "Exons 3-5 set 1") %>%
  pull(dCT)

control_35_2 <- qpcr_dct %>%
  filter(Treatment == "Control" & Primers == "Exons 3-5 set 2") %>%
  pull(dCT)

qpcr_ddct <- qpcr_dct %>%
  mutate(ddCT = case_when(Primers == "Exons 3-4" ~ dCT - control_34,
                          Primers == "Exons 3-5 set 1" ~ dCT - control_35_1,
                          Primers == "Exons 3-5 set 2" ~ dCT - control_35_2),
         FC = (2)^-ddCT,
         Primers = case_when(Primers == "Exons 3-5 set 1" ~ "Exons 3-5 (A)",
                             Primers == "Exons 3-5 set 2" ~ "Exons 3-5 (B)",
                             TRUE ~ Primers)) 

qpcr_ddct$Treatment <- factor(qpcr_ddct$Treatment, levels = c('Control','1uM','5uM','10uM'))

qpcr_plot <- ggplot(qpcr_ddct, aes(y=`FC`, x=Primers, fill=Treatment)) + 
  geom_bar(position="dodge", stat="identity",color = "black") +
  ylab("Fold Change") + 
  xlab(expression(bold(bolditalic("CLK1")~"exon-exon junction"))) +
  theme_Publication() + 
  scale_fill_manual(values = c("lightgrey", "lightblue", "#0C7BDC", "blue3")) +
  guides(fill = guide_legend(title = "Morpholino\nTreatment"))

# Save plot as PDF
pdf(file_plot, width = 6, height = 4)
print(qpcr_plot)
dev.off()

