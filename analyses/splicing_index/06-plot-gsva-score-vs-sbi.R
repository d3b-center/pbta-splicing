################################################################################
# 06-plot-gsea-score-vs-sbi.R
# written by Ammar Naqvi 
#
# usage: Rscript 06-plot-gsea-score-vs-sbi.R
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("ggpubr")
  library("ggplot2")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Source function for plots theme
source(file.path(figures_dir, "theme_for_plots.R"))

## create plots dir if it doesn't exist
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

clin_file  <- file.path(data_dir,"histologies-plot-group.tsv")
sbi_file <- file.path(root_dir,"analyses/splicing_index/results/splicing_index.SE.txt")

## output files for final plots

## get CLK1 psi values in tumors and ctrls
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_df <- vroom(indep_file)

## load histologies info for HGG subty  
histologies_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID)

## Load rmats file
sbi_df <-  vroom(sbi_file) %>%
  # Join rmats data with clinical data
  inner_join(histologies_df, by=c('Sample'='Kids_First_Biospecimen_ID')) %>%
  dplyr::rename('Kids_First_Biospecimen_ID'='Sample')

## Load gsea score file
gsva_scores_df <- vroom(file.path(root_dir,"analyses/clustering_analysis/output/diff_pathways/non_expr_pan_cancer_splice_subset_pam_canberra_0_gsva_output.tsv")) %>%
  dplyr::rename('Kids_First_Biospecimen_ID'='sample_id') %>%
  inner_join(sbi_df ,by='Kids_First_Biospecimen_ID') %>% filter(geneset == 'KEGG_SPLICEOSOME') %>%
  select(Kids_First_Biospecimen_ID,SI,score) 

## create plot
scatterplot_score_sbi <- ggscatter(gsva_scores_df, 
                         x="SI", 
                         y="score",
                         add = "reg.line", 
                         color = "blue",
                         conf.int = TRUE, 
                         cor.coef = TRUE, 
                         cor.method = "spearman",
                         add.params = list(color = "red",
                                           fill = "pink"),
                         ticks = TRUE,
                         size = 2.5, alpha = 0.6) + 
  xlab("Splicing Burden Index") +
  ylab("Splicosome GSVA Score") +
  theme_Publication() 
  

# save plot
pdf(file.path(plots_dir,"corplot-CLK1-high-low-vs-gsva-spliceosome.pdf"),width = 4.5, height = 6)
print(scatterplot_score_sbi)
dev.off()
