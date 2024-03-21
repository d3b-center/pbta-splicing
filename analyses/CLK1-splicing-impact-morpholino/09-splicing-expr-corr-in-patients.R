################################################################################
# 09-splicing-expr-corr-in-patients.R 
# Plot scatter and compute R for identified CLK1 targets expression and PSI 
# values in patient samples
#
# written by Ammar S Naqvi
# Usage: Rscript 09-splicing-expr-corr-in-patients.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")
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
rsem_counts <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")

## read in histology file and count data
## filter histology file for all HGG, only stranded samples
all_hgg_bsids <- read_tsv(clin_file, guess_max = 10000) %>% 
  filter(short_histology == "HGAT",
         RNA_library == "stranded",
         cohort == "PBTA") %>%
  select(Kids_First_Biospecimen_ID, CNS_region)

# keep only hgg ids
rsem_df <- readRDS(rsem_counts) %>%
  select(all_hgg_bsids$Kids_First_Biospecimen_ID)

## load rmats input for CLK1
clk1_rmats <- fread(rmats_file) %>%
  # filter for CLK1 and exon 4
  filter(sample_id %in% all_hgg_bsids$Kids_First_Biospecimen_ID,
         geneSymbol=="CLK1",
         exonStart_0base=="200860124", 
         exonEnd=="200860215") %>% 
  # select minimal info
  select(sample_id, IncLevel1) %>%
  unique()

set.seed(2023)

## combine RNA with psi values for scatter plot/correlation for goi 
## genes up in high exon 4 cells, should be positively correlated with exon 4 PSI
pos_genes <- vroom(file.path(analysis_dir,"results","ctrl_vs_treated.de.formatted.tsv")) %>%
  filter(log2FoldChange <= -2,
         padj <= 0.05) %>% 
  select(gene) %>% 
  pull()

## genes down in high CLK1 exon 4 cells, should be negatively correlated with exon 4 PSI
neg_genes <-vroom(file.path(analysis_dir,"results","ctrl_vs_treated.de.formatted.tsv")) %>%
  filter(log2FoldChange >= 2,
         padj <= 0.05) %>% 
  select(gene) %>% 
  pull()

plot_list_neg <- list()
plot_list_pos <- list()

region_list <- c("all")

for (gene in pos_genes) {
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
    
    # filter count data. 
    rmats_exp_df <- rsem_df %>%
      filter(rownames(.) == gene)  %>% 
      select(all_of(bs_id_list)) %>% 
      pivot_longer(cols = tidyselect::everything(),
                   names_to=c("sample_id"), 
                   values_to="Expr") %>% 
      inner_join(rmats_df, by= "sample_id") %>%
      # need to add this so that we can use this in the plot axis
      mutate(geneSymbol = paste0(gene))
    
    # make scatterplot
    if(nrow(rmats_exp_df) < 3 )
    {
      next
    }
    
    pval = cor.test(rmats_exp_df$IncLevel1, rmats_exp_df$Expr)$p.value 
    corr = cor.test(rmats_exp_df$IncLevel1, rmats_exp_df$Expr)$estimate 
    
    if(pval <= 0.05 & corr > .2)
    {
      p <- create_scatterplot(rmats_exp_df) 
      plot_list_pos <- append(plot_list_pos, list(p))
    }
  }
}

region_list <- c("all")


for (gene in neg_genes) {
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
    
    # filter count data. 
    rmats_exp_df <- rsem_df %>%
      filter(rownames(.) == gene)  %>% 
      select(all_of(bs_id_list)) %>% 
      pivot_longer(cols = tidyselect::everything(),
                   names_to=c("sample_id"), 
                   values_to="Expr") %>% 
      inner_join(rmats_df, by= "sample_id") %>%
      # need to add this so that we can use this in the plot axis
      mutate(geneSymbol = paste0(gene))
    
    # make scatterplot
    if(nrow(rmats_exp_df) < 3)
    {
      next
    }
    
    pval = cor.test(rmats_exp_df$IncLevel1, rmats_exp_df$Expr)$p.value 
    corr = cor.test(rmats_exp_df$IncLevel1, rmats_exp_df$Expr)$estimate 
    
    if(pval <= 0.05 & corr <= -.2)
    {
      p <- create_scatterplot(rmats_exp_df) 
      plot_list_neg <- append(plot_list_neg, list(p))
      
    }
  }
}

# save plots
plot_list_all = list()
plot_list_all = c(plot_list_pos,plot_list_neg)

# save plot 
pdf(file.path(paste(plots_dir, "/corr-CLK1-targets-patient.pdf", sep = "")), width = 21.3, height = 10.6)
plot_grid(plotlist = plot_list_all, ncol = 4, nrow=2)
dev.off()
