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
rmats_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")

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
rmats_df <-  vroom(rmats_file) %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>%
  # Select exon 4
  filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%
  # Select "sample", "geneSymbol", and "IncLevel1" columns
  select(sample_id, geneSymbol, IncLevel1) %>%
  # Join rmats data with clinical data
  inner_join(histologies_df, by=c('sample_id'='Kids_First_Biospecimen_ID')) 

## Load gsea score file
gsva_scores_df <- vroom(file.path(root_dir,"analyses/clustering_analysis/output/diff_pathways/non_expr_pan_cancer_splice_subset_pam_canberra_0_gsva_output.tsv")) %>%
  inner_join(rmats_df, by='sample_id') %>% filter(geneset == 'KEGG_SPLICEOSOME') %>%
  select(sample_id,IncLevel1,score) 

## create plot
scatterplot_score_sbi <- ggscatter(gsva_scores_df, 
                         x="IncLevel1", 
                         y="score", 
                         add = "reg.line", 
                         conf.int = TRUE, 
                         cor.coef = TRUE, 
                         cor.method = "spearman",
                         add.params = list(color = "red",
                                           fill = "pink"),
                         ticks = TRUE) + 
  xlab("Splicing Burden Index") +
  ylab("Splicosome GSEA Score") +
  theme_Publication()    

# save plot
pdf(file.path(plots_dir,"sbi-score.pdf"), width = 4, height = 4)
print(scatterplot_score_sbi)
dev.off()

## Compute quantiles to define high vs low Exon 4 PSI groups
quartiles_sbi <- quantile(gsva_scores_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)

# Get lower quantile (25%)
lower_sbi <- quartiles_sbi[1] 

# Get upper quantile (75%)
upper_sbi <- quartiles_sbi[2]

# Create df with high/low PSI
gsva_scores_sbi_highlow_df <- gsva_scores_df %>%
  mutate(level = case_when(gsva_scores_df$IncLevel1 > upper_sbi ~ "high",
                         gsva_scores_df$IncLevel1 < lower_sbi ~ "low",
                         TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(level)) %>%
  select(sample_id,score,level)

## Make box plot with stats
boxplot_sbi_vs_score <- ggboxplot(gsva_scores_sbi_highlow_df,x = "level", y = "score") +
  ylab(expression(bold("Splicosome GSEA Score"))) +
  xlab(expression(bold("Splicing Burden Index Level"))) +
   lims(y = c(0,.37)) +
  ggforce::geom_sina(aes(color = level), size = 2, shape = 21, fill = NA, stroke = 1) +
  scale_color_manual(name = "level", values = c(high = "#FFC20A", low = "#0C7BDC")) +
  stat_compare_means(position = "identity", label.x = 1) +
  theme_Publication()

# save plot
pdf(file.path(plots_dir,"highlow-sbi-score.pdf"), width = 4.5, height = 6)
print(boxplot_sbi_vs_score)
dev.off()
