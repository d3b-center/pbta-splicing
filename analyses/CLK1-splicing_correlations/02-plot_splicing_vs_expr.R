################################################################################
# 02-plot_splicing_vs_expr.R 
# Plot scatter and compute R for CLK1 expression and PSI values for HGGs
# written by Ammar S Naqvi
# Usage: Rscript 02-plot_splicing_vs_expr.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("ggpubr")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))
source(paste0(analysis_dir, "/util/function-create-scatter-plot.E"))


## define input files
clin_file = file.path(data_dir,"histologies.tsv")
file_gene_counts = file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")
rmats_file <- file.path(data_dir, "rMATS_merged.comparison.tsv.gz")

## define output file
plot_total_hgg_path <- file.path(plots_dir, "CLK1-expr_vs_psi-totalHGG.tiff")
plot_midline_hgg_path <- file.path(plots_dir, "CLK1-expr_vs_psi-midlineHGG.tiff")

## read in histology file and count data
count_df <- readRDS(file_gene_counts)

## filter histology file for all HGG
clin_tab_hgg <- read_tsv(clin_file) %>% 
  filter(short_histology == 'HGAT',
         RNA_library == 'stranded',
         cohort == 'PBTA')

## filter histology file for midline HGG
clin_tab_midline_hgg <- clin_tab_hgg %>%
  filter(CNS_region == 'Midline')

## load rmats input
# create df for each case (midline HGG and all HGGs)
rmats_df_midline_hgg <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  # filter for CLK1 and exon 4
  filter(geneSymbol=="CLK1",
         exonStart_0base=="200860124", exonEnd=="200860215") %>% 
  # join with clinical data
  inner_join(clin_tab_midline_hgg, by=c('sample_id'='Kids_First_Biospecimen_ID')) %>%
  select(sample_id, geneSymbol, IncLevel2) 

rmats_df_hgg <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  # filter for CLK1 and exon 4
  filter(geneSymbol=="CLK1",exonStart_0base=="200860124", exonEnd=="200860215") %>% 
  # join with clinical data
  inner_join(clin_tab_hgg, by=c('sample_id'='Kids_First_Biospecimen_ID')) %>%
  select(sample_id, geneSymbol, IncLevel2) 

## combine with psi values for scatter plot/correlation
count_psi_midline_hgg_df <- filter(count_df,rownames(count_df) == 'CLK1')  %>% 
  select(any_of(rmats_df_midline_hgg$sample_id)) %>% 
  pivot_longer(cols = tidyselect::everything(),names_to=c("sample_id"), values_to="Expr") %>% inner_join(rmats_df_midline_hgg, by='sample_id')

count_psi_hgg_df <- filter(count_df,rownames(count_df) == 'CLK1')  %>% 
  select(any_of(rmats_df_hgg$sample_id)) %>% 
  pivot_longer(cols = tidyselect::everything(),names_to=c("sample_id"), values_to="Expr") %>% inner_join(rmats_df_hgg, by='sample_id')

## generate scatter plots and save to file
scatterplot_midline_hgg <- ggscatter(count_psi_midline_hgg_df, x="IncLevel2", y="Expr", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          add.params = list(color = "red",
                            fill = "pink"),
          ticks = TRUE,
          #xticks.by = .1, yticks.by = .1,
          xlab = "Exon 4 Inclusion (PSI)", ylab = "CLK1 Expr (TPM)") + theme_Publication()

ggplot2::ggsave(plot_midline_hgg_path,
                width=5,
                height=5,
                device="tiff",
                dpi=300)

scatterplot_total_hgg <- ggscatter(count_psi_hgg_df, x="IncLevel2", y="Expr", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          add.params = list(color = "red",
                            fill = "pink"),
          ticks = TRUE,
          #xticks.by = .1, yticks.by = .1,
          xlab = "Exon 4 Inclusion (PSI)", ylab = "CLK1 Expr (TPM)") + theme_Publication()

ggplot2::ggsave(plot_total_hgg_path,
                width=5,
                height=5,
                device="tiff",
                dpi=300)