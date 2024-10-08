################################################################################
# 02-plot_diff-splice-events.R
# written byAmmar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 02-plot_diff-splice-events.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggrepel")
  library("vroom")
  library("ggpubr")
  library("annoFuseData")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories and file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")
input_dir <- file.path(analysis_dir, "input")

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
file_dpsi_plot <- file.path(plots_dir,"dPSI_distr.pdf")

## get and setup input
## rmats file
rmats_merged_file  <- file.path(data_dir,"morpholno.merged.rmats.tsv")

## extract strong splicing changes 
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>% 
                 filter(FDR < 0.05 & PValue < 0.05) 

## extract strong differential splicing cases (dPSI >= |.10|) and use CLK1 high exon 4 as a reference (eg. -dPSI means there more inclusion in CLK1 exon 4)
splicing_df_ES <- splicing_df %>% filter(IncLevelDifference  >= .10) %>% mutate(Preference="Skipping")
splicing_df_EI <- splicing_df %>% filter(IncLevelDifference <= -.10) %>% mutate(Preference="Inclusion",
                                                                                IncLevelDifference = abs(IncLevelDifference) )

psi_comb <- rbind(splicing_df_EI,splicing_df_ES)

## ggstatplot across functional sites
set.seed(123)
counts_psi_comb <- psi_comb %>% 
  dplyr::count(splicing_case, Preference)

## ggstatplot across functional sites
set.seed(123)
plot_dsp <-  ggplot(psi_comb,aes(splicing_case, IncLevelDifference*100) ) +  
  ylab(expression(Delta*"PSI"))+
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  facet_wrap("Preference") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("A3SS", "A5SS"),
                                                                c("A3SS", "MXE"),
                                                                c("A3SS", "RI"),
                                                                c("A3SS", "SE"),
                                                                c("A5SS", "MXE"),
                                                                c("A5SS", "RI"),
                                                                c("A5SS", "SE"),
                                                                c("MXE", "RI"),
                                                                c("MXE", "SE"),
                                                                c("SE", "A5SS"),
                                                                c("SE", "A3SS"),
                                                                c("SE", "RI"))) + 
  scale_color_manual(name = "Preference (CLK1 exon 4 high)", values = c(Inclusion = "#FFC20A", Skipping = "#0C7BDC"))  + 
  theme_Publication() + 
  geom_text(data = counts_psi_comb, aes(label = paste("n =",n), x = splicing_case, y = 0), vjust = 2, size = 4, hjust=.5) +
  
  labs(x= "Splicing Case") + 
  theme(legend.position="none") +
  ylim(c(0,170))


# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 15, height = 7)
plot_dsp
dev.off()

# annotate significant splice events for TSG/Oncogene
annots <- read_tsv(system.file("extdata", "genelistreference.txt", package = "annoFuseData")) %>%
  dplyr::rename(geneSymbol = Gene_Symbol) %>%
  mutate(annotation = case_when(grepl("Oncogene|TumorSuppressorGene", type) ~ "Onco_TSG",
                              type %in% c("Kinase", "Kinase, CosmicCensus", "Kinase, TranscriptionFactor", "CosmicCensus, Kinase CosmicCensus, TranscriptionFactor") ~ "Kinase",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(annotation)) %>% 
  select(geneSymbol, annotation) %>%
  unique()

# export significant splice events
psi_comb_select <- psi_comb %>% 
  dplyr::select(splicing_case, geneSymbol, PValue, FDR, IncLevelDifference, exonStart_0base, exonEnd, 
                           "1stExonStart_0base",'1stExonEnd',
                           '2ndExonStart_0base','2ndExonEnd','riExonStart_0base', 'riExonEnd',"upstreamES", 
                           "upstreamEE","downstreamES","downstreamEE",
                           "longExonStart_0base","longExonEnd",
                           "shortES","shortEE",
                           "flankingES","flankingEE","Preference") %>%
  left_join(annots)

write_tsv(psi_comb_select, 
          file=file.path(results_dir,"splice-events-significant.tsv"))

