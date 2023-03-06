################################################################################
# 07-conduct_gsva_analysis.R
# This script conducts gene set enrichment analysis, specifically using the GSVA 
# method [1] for scoring hallmark human pathway enrichment from RNA-Seq results.
#
# The GSVA scores (i.e., enrichment scores) are calculated to produce a 
# **Gaussian distribution of scores** "under the null hypothesis of no change 
# in the pathway activity throughout the sample population."
#
# The authors claim a benefit to this approach:
#   + "Penalizes deviations that are large in *both* tails"
#   + "Provides a 'normalization' of the enrichment score by 
#      subtracting potential noise
#   + "Emphasizes genes in pathways that are concordantly activated in one 
#      direction only"
#   + "For pathways containing genes strongly acting in both directions, the 
#      deviations with cancel each other out and show little or no enrichment."
#
# Written by Stephanie J. Spielman for CCDL ALSF, 2020 and modified by Ammar S 
#             Naqvi
#
# usage: Rscript 07-conduct_gsva_analysis.R
#
# Reference:
# 1. Sonja Hänzelmann, Robert Castelo, and Justin Guinney. 2013. “GSVA: Gene Set 
# Variation Analysis for Microarray and RNA-Seq Data.” BMC Bioinformatics 14 
# (1): 7. https://doi.org/10.1186/1471-2105-14-7.
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("GSVA")
  library("dplyr")
  library("tidyverse")
  library("readr")
  library("msigdbr")
  library("tibble")
  library("vroom")
  library("ggpubr")
  library("ggplot2")
})

## Magrittr pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")
figures_dir <- file.path(root_dir, "figures")

# source function for theme for plots 
source(file.path(figures_dir, "theme_for_plots.R"))
  
## load input files
human_hallmark  <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") ## human hallmark genes from `migsdbr` package. The loaded data is a tibble.
#human_hallmark  <- msigdbr::msigdbr(species = "Homo sapiens", category = "C6") ## human hallmark genes from `migsdbr` package. The loaded data is a tibble.


## make histologies dataframe
histology_df <- readr::read_tsv( file.path(data_dir, "histologies.tsv")  , guess_max = 100000)

#### Prepare hallmark genes: Create a list of hallmarks, each of which is a list of genes -----------------------------------------------
human_hallmark_twocols <- human_hallmark %>% dplyr::select(gs_name, human_gene_symbol)
human_hallmark_list    <- base::split(human_hallmark_twocols$human_gene_symbol, list(human_hallmark_twocols$gs_name))

#### Perform gene set enrichment analysis --------------------------------------------------------------------
##rmats input

# Prepare expression data: log2 transform re-cast as matrix
histology_rna_df <- histology_df %>% 
  # Only include RNA-Seq, PBTA cohorts, stranded, Midline HGGs samples
  filter(experimental_strategy == "RNA-Seq", cohort == "PBTA", cohort == "PBTA", 
         CNS_region == 'Midline',short_histology == 'HGAT') 

rmats_file <-  file.path(data_dir,"rMATS_merged.single.SE.tsv.gz")
rmats_df <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  
  # select specific samples and extract CLK1 exon 4  
  dplyr::filter(geneSymbol=="CLK1") %>% dplyr::filter(exonStart_0base=="200860124", exonEnd=="200860215") %>% dplyr::select(sample, geneSymbol, IncLevel1) %>% 
  inner_join(histology_rna_df, by=c('sample'='Kids_First_Biospecimen_ID')) %>%  dplyr::select(sample, geneSymbol, IncLevel1) 


## compute quantiles to define high vs low Exon 4 PSI groups
quartiles_psi <- quantile(rmats_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)
IQR_psi <- IQR(rmats_df$IncLevel1)
lower_psi <- quartiles_psi[1] 
upper_psi <- quartiles_psi[2] 

rmats_subset_samples_lowEx4  <- rmats_df %>% dplyr::filter(rmats_df$IncLevel1 < lower_psi) %>% dplyr::select(sample)
rmats_subset_samples_highEx4 <- rmats_df %>% dplyr::filter(rmats_df$IncLevel1 > upper_psi) %>% dplyr::select(sample)
rmats_subset_samples_Ex4 <- rbind(rmats_subset_samples_lowEx4,rmats_subset_samples_highEx4)

histology_rna_df <- inner_join(histology_rna_df, rmats_subset_samples_Ex4,by=c('Kids_First_Biospecimen_ID'='sample')) #  %>% dplyr::select(Kids_First_Biospecimen_ID)

# filter expression data to exclude GTEx and TCGA
expression_data <- as.data.frame(readRDS(file.path(data_dir, "gene-counts-rsem-expected_count-collapsed.rds")  )) %>% 
  dplyr::select(matches(histology_rna_df$Kids_First_Biospecimen_ID) )
 
# for each type of the RNA library, we subset the expression matrix accordingly and run gsea scores for each RNA library 
rna_library_list <- histology_rna_df %>% pull(RNA_library) %>% unique()

# Further subset to each cohort to deal with size issues
cohort_list <- histology_rna_df %>% pull(cohort) %>% unique()

gsea_scores_df_tidy <- data.frame()

# iterate through each cohort and RNA library type 
for(i in 1:length(rna_library_list)){
  rna_library = rna_library_list[i]
  # get bs id for one particular rna library type
  rna_library_type_bs_id <- histology_rna_df %>% 
    dplyr::filter(RNA_library == rna_library) %>% 
    pull(Kids_First_Biospecimen_ID) %>%
    unique()
  
  # Filter the expression data to this RNA library type
  # Subset to the remaining samples 
  expression_data_each <- expression_data %>% 
    dplyr::select(rna_library_type_bs_id)
  
  ### Rownames are genes and column names are samples
  expression_data_each_log2_matrix <- as.matrix( log2(expression_data_each + 1) )
  
  #We then calculate the Gaussian-distributed scores
  gsea_scores_each <- GSVA::gsva(expression_data_each_log2_matrix,
                                 human_hallmark_list,
                                 method = "gsva",
                                 min.sz=1, max.sz=1500,## Arguments from K. Rathi
                                 parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                                 mx.diff = TRUE)        ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
  
  ### Clean scoring into tidy format
  gsea_scores_each_df <- as.data.frame(gsea_scores_each) %>%
    rownames_to_column(var = "hallmark_name")
  
  #first/last_bs needed for use in gather (we are not on tidyr1.0)
  first_bs <- head(colnames(gsea_scores_each), n=1)
  last_bs  <- tail(colnames(gsea_scores_each), n=1)
  
  rna_library<-gsub(" ", "_", rna_library)
  rna_library<-stringr::str_to_lower(gsub("-", "", rna_library))
  
  gsea_scores_each_df_tidy <- gsea_scores_each_df %>%
    tidyr::gather(Kids_First_Biospecimen_ID, gsea_score, !!first_bs : !!last_bs) %>%
    dplyr::select(Kids_First_Biospecimen_ID, hallmark_name, gsea_score) %>%
    dplyr::mutate(data_type = rna_library)
  
  gsea_scores_df_tidy <-  bind_rows(gsea_scores_df_tidy , gsea_scores_each_df_tidy)
}

#### Export GSEA scores to TSV --------------------------------------------------------------------
scores_output_file <- (file.path(results_dir, "gsea_out.tsv"))
write_tsv(gsea_scores_df_tidy, scores_output_file)

## compare high vs low SBI and their associated CLK1 expression 
expression_tidy_df <- as.data.frame(expression_data) %>%
  rownames_to_column(var = "gene")

#first/last_bs needed for use in gather (we are not on tidyr1.0)
first_bs <- head(colnames(expression_data), n=1)
last_bs  <- tail(colnames(expression_data), n=1)

expression_each_tidy_df <- expression_tidy_df %>%
  tidyr::gather(Kids_First_Biospecimen_ID, expr, !!first_bs : !!last_bs) 

expression_each_tidy_CLK1_df <- expression_each_tidy_df %>% filter(gene=='CLK1' | gene=="GALNT15" | gene=='DIAPH3')

ggplot(expression_each_tidy_CLK1_df, aes(y=gene, x=Kids_First_Biospecimen_ID, fill=expr)) + 
  geom_raster() + geom_tile() +  scale_fill_gradient(low="white", high="#DC3220") + xlab("Sample") +
  ylab("Gene") + theme(axis.text.x = element_text(size = 10,angle = 90,hjust = .7, vjust = .8),axis.text.y = element_text(size = 10)) 

## add column for CLK1 groupings and add in expression
rmats_subset_with_CLK1expr_df <- rmats_df %>% dplyr::rename(gene=geneSymbol) %>% 
                                              dplyr::rename(Kids_First_Biospecimen_ID=sample) %>% 
                                              dplyr::inner_join(expression_each_tidy_CLK1_df, by=c('gene','Kids_First_Biospecimen_ID')) %>%
                                              dplyr::mutate(CLK1_group = if_else(IncLevel1 < lower_psi, "Low", "High")) %>%  
                                              dplyr::select(Kids_First_Biospecimen_ID, gene,IncLevel1, expr, CLK1_group) 

##plot expression categorized by CLK1 expression
boxplot_grp_expr <- ggplot(rmats_subset_with_CLK1expr_df,aes(CLK1_group,expr)) + 
                      geom_boxplot(aes(fill=CLK1_group)) + 
                      stat_compare_means(method = "t.test") + 
                      geom_jitter() + theme_Publication()

file_tiff_plot = file.path(plots_dir,"boxplot_CLK1-expr.tiff")

# save plot tiff version
tiff(file_tiff_plot, height = 2000, width = 2000, res = 300)
print(boxplot_grp_expr)
dev.off()


