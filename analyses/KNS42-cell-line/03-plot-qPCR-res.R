################################################################################
# 03-plot-qPCR-res.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi
#
# usage: Rscript 03-plot-qPCR-res.R
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("vroom")
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

qpcr_res_file = file.path(input_dir,"qpcr-raw-triplicates.tsv")
qpcr_fc_df <- read_tsv(qpcr_res_file) 

# format and compute ddct values from raw values
qpcr_ddct_transform <- as.data.frame(qpcr_fc_df) %>%
  ## group by junction and calculate means
  group_by(SJ) %>%
  summarise_all(mean) %>%
  
  ## calculate dct and ddct values
  mutate(across(Control:`10uM`,  ~ last(.) - . )) %>% 
  mutate(across(`1uM`:`10uM`, ~ Control - .  )) %>% 
  mutate(across(Control, ~ 1  )) %>% 
  
  ## reformat and rename to plot
  pivot_longer(cols=c('Control','1uM','5uM','10uM'), 
               names_to = "Concentration",
               values_to = "ddct") %>% 
  dplyr::filter(SJ != 'HPRT') %>% 
  mutate("ddct" = ddct) %>%
  mutate(`Fold Change` = (2)^-ddct) %>% 
  mutate(`Fold Change` = case_when(Concentration=='Control' ~ 1,
         .default = `Fold Change`))

qpcr_ddct_transform$Concentration <- factor(qpcr_ddct_transform$Concentration, levels = c('Control','1uM','5uM','10uM'))

# Save plot as PDF
pdf(file_plot, width =8, height = 5)

ggplot(qpcr_ddct_transform, aes(y=`Fold Change`, x=Concentration, fill=SJ)) + 
  geom_bar(position="dodge", stat="identity",color = "black") +
  ylab("Fold Change") + 
  theme_Publication() + 
  scale_fill_manual(name = "Splice Junction",values=c("#FFC20A","#0C7BDC","lightblue")) 
dev.off()

