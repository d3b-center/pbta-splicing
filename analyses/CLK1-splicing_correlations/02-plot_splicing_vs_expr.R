################################################################################
# 02-plot_splicing_vs_expr.R 
# Plot scatter and compute R for CLK1 expression and PSI values for HGGs
# written by Ammar S Naqvi, Jo Lynne Rokita
# Usage: Rscript 02-plot_splicing_vs_expr.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("ggpubr")
  library("vroom")
  library('data.table')
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
plots_dir <- file.path(analysis_dir, "plots")
figures_dir <- file.path(root_dir, "figures")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
source(file.path(figures_dir, "theme_for_plots.R"))
source(file.path(analysis_dir, "util", "function-create-scatter-plot.R"))

## define input files
clin_file <- file.path(data_dir,"histologies.tsv")
file_gene_counts <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")

## read in histology file and count data
## filter histology file for all HGG, only stranded samples
all_hgg_bsids <- read_tsv(clin_file, guess_max = 10000) %>% 
  filter(short_histology == "HGAT",
         RNA_library == "stranded",
         cohort == "PBTA") %>%
  select(Kids_First_Biospecimen_ID, CNS_region)

# keep only hgg ids
count_df <- readRDS(file_gene_counts) %>%
  select(all_hgg_bsids$Kids_First_Biospecimen_ID)

## load rmats input for CLK1
clk1_rmats <- fread(rmats_file) %>%
  # filter for CLK1 and exon 4
  filter(geneSymbol=="CLK1",
         exonStart_0base=="200860124", 
         exonEnd=="200860215") %>% 
  # select minimal info
  select(sample_id, IncLevel1) 

set.seed(2023)

## combine RNA with psi values for scatter plot/correlation for goi 
goi_list <- c("CLK1", "SRSF1")
region_list <- c("midline", "all")

for (gene in goi_list) {
  for (brain_region in region_list) {
    
    if (brain_region == "all"){
      # take all hgg bs ids
      bs_id_list <- all_hgg_bsids$Kids_First_Biospecimen_ID
    }
    
    else if (brain_region == "midline"){
      # take only midline bs ids
      bs_id_list <- all_hgg_bsids %>%
        filter(CNS_region == "Midline") %>%
        pull(Kids_First_Biospecimen_ID)
    }
    
  # filter RMATs file for the bs ids of interest
    rmats_df <- clk1_rmats %>%
    # filter for bs ids of interest and CLK1 and exon 4
    filter(sample_id %in% bs_id_list)
    
  # filter count data
    rmats_exp_df <- count_df %>%
      filter(rownames(.) == gene)  %>% 
      select(all_of(bs_id_list)) %>% 
      pivot_longer(cols = tidyselect::everything(),
                   names_to=c("sample_id"), 
                   values_to="Expr") %>% 
      inner_join(rmats_df, by= "sample_id") %>%
      # need to add this so that we can use this in the plot axis
      mutate(geneSymbol = paste0(gene))
    
  # make scatterplot
    p <- create_scatterplot(rmats_exp_df) 
    # save plot 
    pdf(file.path(paste(plots_dir, "/", gene, "_exp_vs_CLK1_psi_", brain_region, "_hgg.pdf", sep = "")), width = 4.5, height = 4.5)
    print(p)
    dev.off()
  }
}
