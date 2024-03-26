################################################################################
# 05-CLK-SRSF-expr-correlations.R
# Plot scatter and compute R for CLK1 expression and PSI values for HGGs
# written by Ryan Corbett
# Usage: Rscript 05-CLK-SRSF-expr-correlations.R
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
figures_dir <- file.path(root_dir, "figures")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
source(file.path(figures_dir, "theme_for_plots.R"))

## define input files
clin_file <- file.path(data_dir,"histologies.tsv")
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
rsem_counts <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")

cohort_file <- file.path(root_dir, "analyses", "cohort_summary",
                         "results", "histologies-plot-group.tsv")

## read in histology file and count data
## filter histology file for all HGG, only stranded samples
hist <- read_tsv(clin_file, guess_max = 10000)

cohort_df <- read_tsv(cohort_file)

indep_rna_df <- vroom(indep_rna_file) %>% 
  dplyr::filter(cohort == 'PBTA')

# extract hgg samples in indepenent specimens file, add plot_group
all_hgg_bsids <- hist %>% 
  filter(short_histology == "HGAT",
         RNA_library != "exome_capture",
         cohort == "PBTA",
      Kids_First_Biospecimen_ID %in% indep_rna_df$Kids_First_Biospecimen_ID) %>%
  select(Kids_First_Biospecimen_ID, RNA_library, CNS_region, match_id, molecular_subtype) %>%
  left_join(cohort_df %>% dplyr::select(Kids_First_Biospecimen_ID, plot_group))

# filter expr df for hgg samples
rsem_df <- readRDS(rsem_counts) %>%
  select(all_hgg_bsids$Kids_First_Biospecimen_ID)

# load rmats input for CLK1
clk1_rmats <- fread(rmats_file) %>%
  # filter for CLK1 and exon 4
  filter(sample_id %in% all_hgg_bsids$Kids_First_Biospecimen_ID,
         geneSymbol=="CLK1",
         exonStart_0base=="200860124",
         exonEnd=="200860215") %>%
  # select minimal info
  select(sample_id, IncLevel1) %>%
  unique()

# save CLK1 PSI df
write_tsv(clk1_rmats,
          file.path(results_dir, "clk1-exon4-psi-hgg.tsv"))

set.seed(2023)

## create lists of CLK and SRSF genes
srsf_list <- c("CLK1", "SRSF1", "SRSF2", "SRSF3", "SRSF4",
              "SRSF5", "SRSF6", "SRSF7", "SRSF8", "SRSF9",
              "SRSF10", "SRSF11")
clk_list <- c("CLK1", "CLK2", "CLK3", "CLK4")

goi_list <- list("SRSF" = srsf_list, "CLK" = clk_list)

region_list <- c("midline", "other", "all")

# Loop through SRSF and CLK gene lists to assess correlation between CLK1 exon4 inclusion and gene expression
 for (goi in names(goi_list)) {
   # assess correlations in all HGG, DMG, and other HGG separately
  for (brain_region in region_list) {
    
    if (brain_region == "all"){
      # take all hgg bs ids
      bs_id_list <- all_hgg_bsids$Kids_First_Biospecimen_ID
    }
    
    else if (brain_region == "midline"){
      # take only midline bs ids
      bs_id_list <- all_hgg_bsids %>%
        filter(plot_group == "DIPG or DMG") %>%
        pull(Kids_First_Biospecimen_ID)
    }
    
    else if (brain_region == "other"){
      # take only other hgg bs ids
      bs_id_list <- all_hgg_bsids %>%
        filter(plot_group == "Other high-grade glioma") %>%
        pull(Kids_First_Biospecimen_ID)
    }
    
  # filter RMATs file for the bs ids of interest
    rmats_df <- clk1_rmats %>%
    # filter for bs ids of interest and CLK1 and exon 4
    filter(sample_id %in% bs_id_list)
    
  # filter count data for GOI
    rmats_exp_df <- rsem_df %>%
      filter(rownames(.) %in% goi_list[[goi]])  %>%
      rownames_to_column("geneSymbol") %>%
      select(all_of(c("geneSymbol", bs_id_list))) %>% 
      gather(key = "sample_id", value = "Expr", -geneSymbol) %>%
      inner_join(rmats_df, by= "sample_id") %>%
      dplyr::mutate(logExp = log(Expr, 2))
    
  # make scatterplot
   if (goi == "SRSF"){
     
     # create scatterplots between ex4 inc and SRSF RNA expr
     p <-  ggplot(rmats_exp_df, aes(x = IncLevel1, y = logExp)) +
       geom_point() +
       stat_smooth(method = "lm", 
                   formula = y ~ x, 
                   geom = "smooth", 
                   colour = "red",
                   fill = "pink",
                   linetype="dashed") +
       labs(x = "CLK1 Exon 4 Inclusion (PSI)",
            y = "RSEM expected counts (log2)",
            colour = "Library type") + 
       stat_cor(method = "pearson", label.x = .1,
                label.y = 8.25, size = 3) +
       ylim(c(7.5, 15)) +
       facet_wrap(~geneSymbol, nrow = 4) + 
       theme_Publication()
     
     # save plot
     pdf(file.path(paste(plots_dir, "/", goi, "_exp_vs_CLK1_psi_", brain_region, "_hgg.pdf", sep = "")), width = 8, height = 10)
     print(p)
     dev.off()
     
   } else {
     
     # create scatterplots between ex4 inclusion and CLK RNA expr
     p <-  ggplot(rmats_exp_df, aes(x = IncLevel1, y = logExp)) +
       geom_point(colour = "black") +
       stat_smooth(method = "lm", 
                   formula = y ~ x, 
                   geom = "smooth", 
                   colour = "red",
                   fill = "pink",
                   linetype="dashed") +
       labs(x = "CLK1 Exon 4 Inclusion (PSI)",
            y = "RSEM expected counts (log2)",
            fill = "Library type") + 
       stat_cor(method = "pearson", label.x = .1,
                label.y = 14, size = 3) +
       ylim(c(7.5, 15)) +
       facet_wrap(~geneSymbol, nrow = 4) + 
       theme_Publication()
     
     # save plot
     pdf(file.path(paste(plots_dir, "/", goi, "_exp_vs_CLK1_psi_", brain_region, "_hgg.pdf", sep = "")), width = 4, height = 8)
     print(p)
     dev.off()
     
   }
   
  }
   
 }

# redefine srsf_list to exclude CLK1
srsf_list <- c("SRSF1", "SRSF2", "SRSF3", "SRSF4",
               "SRSF5", "SRSF6", "SRSF7", "SRSF8", "SRSF9",
               "SRSF10", "SRSF11")

# plot CLK expr vs SRSF expr
for (clk in clk_list) {
  for (brain_region in region_list) {
    
    if (brain_region == "all"){
      # take all hgg bs ids
      bs_id_list <- all_hgg_bsids$Kids_First_Biospecimen_ID
    }
    
    else if (brain_region == "midline"){
      # take only midline bs ids
      bs_id_list <- all_hgg_bsids %>%
        filter(plot_group == "DIPG or DMG") %>%
        pull(Kids_First_Biospecimen_ID)
    }
    
    else if (brain_region == "other"){
      # take only midline bs ids
      bs_id_list <- all_hgg_bsids %>%
        filter(plot_group == "Other high-grade glioma") %>%
        pull(Kids_First_Biospecimen_ID)
    }
    
    # Get CLK expr
    clk_expr <- rsem_df %>%
      filter(rownames(.) == clk) %>% 
      rownames_to_column("geneSymbol") %>%
      select(all_of(c("geneSymbol", bs_id_list))) %>% 
      gather(key = "sample_id", value = "Expr", -geneSymbol) %>%
      dplyr::mutate(logExp = log(Expr, 2)) %>%
      dplyr::select(sample_id, logExp) %>%
      dplyr::rename("clk_logExp" = logExp)
    
    # Get SRSF expr
    exp_df <- rsem_df %>%
      # filter(rownames(.) == gene)  %>%
      filter(rownames(.) %in% srsf_list)  %>%
      rownames_to_column("geneSymbol") %>%
      select(all_of(c("geneSymbol", bs_id_list))) %>% 
      gather(key = "sample_id", value = "Expr", -geneSymbol) %>%
      inner_join(clk_expr, by= "sample_id") %>%
      dplyr::mutate(logExp = log(Expr, 2))
    
    # Plot SRSF expr vs CLK expr
    if (clk == "CLK3") {
      
      p <-  ggplot(exp_df, aes(x = clk_logExp, y = logExp)) +
        geom_point(colour = "black") +
        stat_smooth(method = "lm", 
                    formula = y ~ x, 
                    geom = "smooth", 
                    colour = "red",
                    fill = "pink",
                    linetype="dashed") +
        labs(x = glue::glue("{clk} RSEM expected counts (log2)"),
             y = "RSEM expected counts (log2)") + 
        stat_cor(method = "pearson", label.x = 8,
                 label.y = 15, size = 3) +
        xlim(c(8,14)) +
        ylim(c(7.5, 16)) +
        facet_wrap(~geneSymbol, nrow = 4) + 
        theme_Publication()
      
    } else if (clk == "CLK4"){
      
      p <-  ggplot(exp_df, aes(x = clk_logExp, y = logExp)) +
        geom_point(colour = "black") +
        stat_smooth(method = "lm", 
                    formula = y ~ x, 
                    geom = "smooth", 
                    colour = "red",
                    fill = "pink",
                    linetype="dashed") +
        labs(x = glue::glue("{clk} RSEM expected counts (log2)"),
             y = "RSEM expected counts (log2)") + 
        stat_cor(method = "pearson", label.x = 8,
                 label.y = 15, size = 3) +
        xlim(c(8,12)) +
       ylim(c(7.5, 16)) +
        facet_wrap(~geneSymbol, nrow = 4) + 
        theme_Publication()
      
    } else {
      
      p <-  ggplot(exp_df, aes(x = clk_logExp, y = logExp)) +
        geom_point(colour = "black") +
        stat_smooth(method = "lm", 
                    formula = y ~ x, 
                    geom = "smooth", 
                    colour = "red",
                    fill = "pink",
                    linetype="dashed") +
        labs(x = glue::glue("{clk} RSEM expected counts (log2)"),
             y = "RSEM expected counts (log2)") + 
        stat_cor(method = "pearson", label.x = 8,
                 label.y = 15, size = 3) +
        ylim(c(7.5, 16)) +
        facet_wrap(~geneSymbol, nrow = 4) + 
        theme_Publication()
  
    }
      
      # save plot
      pdf(file.path(paste(plots_dir, "/", clk, "_exp_vs_SRSF_exp_", brain_region, "_hgg.pdf", sep = "")), width = 8, height = 10)
      print(p)
      dev.off()
    
  }
  
}
