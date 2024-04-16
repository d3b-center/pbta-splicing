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
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
rsem_tpm <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")
isoform_file <- file.path(data_dir, "rna-isoform-expression-rsem-tpm.rds")
rmats_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")

cohort_file <- file.path(root_dir, "analyses", "cohort_summary",
                         "results", "histologies-plot-group.tsv")

## read in histology file and count data
## filter histology file for all HGG, only stranded samples
cohort_df <- read_tsv(cohort_file)
indep_rna_df <- read_tsv(indep_rna_file)

# extract hgg samples in indepenent specimens file, add plot_group
all_hgg_bsids <- cohort_df %>% 
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma"),
         RNA_library == "stranded",
         cohort == "PBTA",
         Kids_First_Biospecimen_ID %in% indep_rna_df$Kids_First_Biospecimen_ID) %>%
  select(Kids_First_Biospecimen_ID, RNA_library, CNS_region, match_id, molecular_subtype, plot_group)

CLK1_ex4_rsem <- readRDS(isoform_file) %>%
  filter(transcript_id == "ENST00000321356.9")  %>%
  dplyr::select(-transcript_id) %>%
  gather(key = "sample_id", value = "Expr", -gene_symbol) %>%
  dplyr::mutate(logExp = log(Expr, 2)) %>%
  dplyr::rename("geneSymbol" = gene_symbol) %>%
  dplyr::filter(sample_id %in% all_hgg_bsids$Kids_First_Biospecimen_ID,
                !is.infinite(logExp)) 

srsf_list <- c("SRSF1", "SRSF2", "SRSF3", "SRSF4",
               "SRSF5", "SRSF6", "SRSF7", "SRSF8", "SRSF9",
               "SRSF10", "SRSF11")
clk_list <- c("CLK1-201", "CLK1", "CLK2", "CLK3", "CLK4")
srpk_list <- c("SRPK1", "SRPK2", "SRPK3")

# load rmats input for CLK1
clk1_rmats <- read_tsv(rmats_file) %>%
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

# load gene expression counts and merge CLK1 transcript counts, exon4 PSI values
rmats_exp_df <- readRDS(rsem_tpm) %>%
  select(all_hgg_bsids$Kids_First_Biospecimen_ID) %>%
  filter(rownames(.) %in% c(srsf_list, clk_list, srpk_list))  %>%
  rownames_to_column("geneSymbol") %>%
  #  select(all_of(c("geneSymbol", bs_id_list))) %>% 
  gather(key = "sample_id", value = "Expr", -geneSymbol) %>%
  dplyr::mutate(logExp = log(Expr, 2)) %>%
  bind_rows(CLK1_ex4_rsem) %>%
  inner_join(clk1_rmats, by= "sample_id")

set.seed(2023)

## create lists of CLK and SRSF genes
goi_list <- list("SRSF" = srsf_list, "CLK" = clk_list, "SRPK" = srpk_list)

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
    
    # make scatterplot
    if (goi == "SRSF"){
      
      # filter count data for GOI
      plot_df <- rmats_exp_df %>%
        dplyr::filter(sample_id %in% bs_id_list,
                      geneSymbol %in% srsf_list)
      
      # create scatterplots between ex4 inc and SRSF/SRPK RNA expr
      p <-  ggplot(plot_df, aes(x = IncLevel1, y = logExp)) +
        geom_point() +
        stat_smooth(method = "lm", 
                    formula = y ~ x, 
                    geom = "smooth", 
                    colour = "red",
                    fill = "pink",
                    linetype="dashed") +
        labs(x = "CLK1 Exon 4 Inclusion (PSI)",
             y = "RSEM TPM (log2)") + 
        stat_cor(method = "pearson",
                 label.x = 0, label.y = 8, size = 3) +
        # ylim(c(7.5, 15)) +
        facet_wrap(~geneSymbol, nrow = 4) + 
        theme_Publication()
      
      # save plot
      pdf(file.path(paste(plots_dir, "/", goi, "_exp_vs_CLK1_psi_", brain_region, "_hgg.pdf", sep = "")), width = 8, height = 10)
      print(p)
      dev.off()
      
    } else {
      
      for(i in 2:3){
        plot_df <- rmats_exp_df %>%
          dplyr::filter(sample_id %in% bs_id_list,
                        geneSymbol %in% goi_list[[i]])
        
        # create scatterplots between ex4 inclusion and CLK RNA expr
        p <-  ggplot(plot_df, aes(x = IncLevel1, y = logExp)) +
          geom_point(colour = "black") +
          stat_smooth(method = "lm", 
                      formula = y ~ x, 
                      geom = "smooth", 
                      colour = "red",
                      fill = "pink",
                      linetype="dashed") +
          labs(x = "CLK1 Exon 4 Inclusion (PSI)",
               y = "RSEM TPM (log2)") + 
          stat_cor(method = "pearson",
                   label.x = 0, label.y = 6, size = 3) +
          # ylim(c(NA, 15)) +
          facet_wrap(~geneSymbol, nrow = 5, scales = "free_y") + 
          theme_Publication()
        
        # save plot
        pdf(file.path(paste(plots_dir, "/", goi, "_exp_vs_CLK1_psi_", brain_region, "_hgg.pdf", sep = "")), width = 4, height = 8)
        print(p)
        dev.off()
      } 
    }
    
  }
  
}

# plot CLK expr vs SRSF/SRPK expr
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
    clk_expr <- rmats_exp_df %>%
      filter(geneSymbol == clk,
             sample_id %in% bs_id_list) %>% 
      dplyr::select(sample_id, logExp) %>%
      dplyr::rename("clk_logExp" = logExp)
    
    # Get SRSF/SRPK expr
    exp_df <- rmats_exp_df %>%
      filter(geneSymbol %in% c(srsf_list, srpk_list),
             sample_id %in% bs_id_list)  %>%
      inner_join(clk_expr, by= "sample_id") 
    
    # Plot SRSF expr vs CLK expr
    
    if (clk == "CLK1-201"){
      
      p <-  ggplot(exp_df, aes(x = clk_logExp, y = logExp)) +
        geom_point(colour = "black") +
        stat_smooth(method = "lm", 
                    formula = y ~ x, 
                    geom = "smooth", 
                    colour = "red",
                    fill = "pink",
                    linetype="dashed") +
        labs(x = glue::glue("{clk} RSEM TPM (log2)"),
             y = "RSEM TPM (log2)") + 
        stat_cor(method = "pearson",
                 label.x = -4, label.y = 8, size = 3) +
        # xlim(c(-5,6)) +
        # ylim(c(NA, 15.5)) +
        facet_wrap(~geneSymbol, nrow = 3, scales = "free_y") + 
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
        labs(x = glue::glue("{clk} RSEM TPM (log2)"),
             y = "RSEM TPM (log2)") + 
        stat_cor(method = "pearson",
                 label.x = 0, label.y = 8, size = 3) +
        # xlim(c(NA,13)) +
        # ylim(c(NA, 15.5)) +
        facet_wrap(~geneSymbol, nrow = 3, scales = "free_y") + 
        theme_Publication()
      
    }
    
    # save plot
    pdf(file.path(paste(plots_dir, "/", clk, "_exp_vs_SRSF_SRPK_exp_", brain_region, "_hgg.pdf", sep = "")), width = 12, height = 8)
    print(p)
    dev.off()
    
  }
  
}
