################################################################################
# 01-co-occurence-interactions.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 01-co-occurence-interactions.R
################################################################################

## libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(maftools)
  library(vroom)
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "mutation_vs_splicing")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")

## check and create plots dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## input files
maf_file <- file.path(data_dir,"snv-consensus-plus-hotspots.maf.tsv.gz")
clin_file <- file.path(data_dir, "histologies.tsv")
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_wgs_file <- file.path(data_dir, "independent-specimens.wgswxspanel.primary.prefer.wgs.tsv")
rmats_file <- file.path(data_dir, "rMATS_merged.comparison.tsv.gz")

indep_rna_df <- vroom(indep_rna_file)  %>% 
  dplyr::filter(cancer_group=='High-grade glioma',
         cohort == 'PBTA') 

indep_wgs_df <- vroom(indep_wgs_file)  %>% 
  filter(cancer_group=='High-grade glioma',
         cohort == 'PBTA') 

## filter for samples that have both RNA and WGS
indep_sample_rna_wgs <- inner_join(indep_wgs_df,indep_rna_df, by='Kids_First_Participant_ID') %>% 
  dplyr::rename('DNA_id'=Kids_First_Biospecimen_ID.x, 'RNA_id'=Kids_First_Biospecimen_ID.y)

histologies_df <- vroom(clin_file)  %>% 
  filter(short_histology == 'HGAT',
         cohort == 'PBTA',
         CNS_region == 'Midline'
        ) %>% 
  inner_join(indep_sample_rna_wgs, by='Kids_First_Participant_ID') %>% 
  filter(short_histology == 'HGAT',
         cohort == 'PBTA',
         CNS_region == 'Midline'
  ) %>%
  filter(experimental_strategy == 'RNA-Seq' |
         experimental_strategy == 'WGS' |
          experimental_strategy == 'WXS')

  

## get splicing events and subset based on samples with both wgs/rna
splice_df <-  vroom(rmats_file) %>%
  dplyr::filter(splicing_case == 'SE', 
               IncLevelDifference >= abs(.20) ) %>%
  dplyr::rename('RNA_id'=sample_id) %>% 
  inner_join(histologies_df, by='RNA_id', relationship = "many-to-many") %>%
  dplyr::mutate(geneSymbol=paste0(geneSymbol,"_spl")) %>%
  dplyr::select(Kids_First_Participant_ID,geneSymbol) %>%
  dplyr::rename('Tumor_Sample_Barcode'=Kids_First_Participant_ID, 'Hugo_Symbol'=geneSymbol) %>%
  mutate(Variant_Classification='Splicing', Variant_Type='Other')

## filter maf table for samples with RNA splicing + HGGs + midline + req'd cols
maf_df <- vroom(maf_file
                #, 
                #data.table = FALSE
                ) %>% 
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))

maf_fil_df <- maf_df %>%
  dplyr::filter(Tumor_Sample_Barcode %in% histologies_df$Tumor_Sample_Barcode.x) %>%
  dplyr::rename('Kids_First_Biospecimen_ID'=Tumor_Sample_Barcode) %>% 
  inner_join(histologies_df, by='Kids_First_Biospecimen_ID') %>% 
  dplyr::mutate('Tumor_Sample_Barcode'=Kids_First_Participant_ID) %>% 
  dplyr::filter(vaf > 0.05) %>% 
  dplyr::select(
    Hugo_Symbol,
    Entrez_Gene_Id,
    Chromosome,
    Start_Position,
    End_Position,
    Strand,
    Variant_Classification,
    Variant_Type,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    Tumor_Sample_Barcode,
    t_depth,
    t_ref_count,
    t_alt_count,
    Consequence
  ) %>%
  dplyr::mutate(Hugo_Symbol=paste0(Hugo_Symbol,"_mut"))


## rename histologies_df Kids_First_Participant_ID to Tumor_Sample_Barcode
histologies_df <- histologies_df %>%  
  dplyr::rename('Tumor_Sample_Barcode'=Kids_First_Participant_ID)
  
## combine mutation and splicing info
combined_maf <- bind_rows(maf_fil_df, splice_df)


maf = read.maf(maf = combined_maf, 
               clinicalData = histologies_df,
               vc_nonSyn = c("Frame_Shift_Del", 
                             "Frame_Shift_Ins", 
                             "Splice_Site", 
                             "Translation_Start_Site",
                             "Nonsense_Mutation",
                             "Nonstop_Mutation", 
                             "In_Frame_Del",
                             "In_Frame_Ins", 
                             "Missense_Mutation",
                             "Del",
                             "Amp",
                             "Splicing") )


## define colors
colors = c("Missense_Mutation" = "#35978f",
           "Nonsense_Mutation" = "#000000",
           "Frame_Shift_Del" = "#56B4E9",
           "Frame_Shift_Ins" = "#FFBBFF",
           "Splice_Site" = "#F0E442",
           "Translation_Start_Site" = "#191970",
           "Nonstop_Mutation" = "#545454",
           "In_Frame_Del" = "#CAE1FF",
           "In_Frame_Ins" = "#FFE4E1",
           "Multi_Hit" = "#f46d43",
           "Del" = "#0072B2",
           "Amp" = "#c51b7d",
           "Splicing" = "green")

variant_labels = c("Missense Mutation",
                   "Nonsense Mutation",
                   "Frameshift Deletion",
                   "Frameshift Insertion",
                   "Splice Site Mutation",
                   "Translation Start Site",
                   "Nonstop Mutation",
                   "Inframe Deletion",
                   "Inframe Insertion",
                   "Multi-Hit",
                   "Del",
                   "Amp",
                   "Splicing")

oncoplot(maf,annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, legendFontSize = 1,
         showTitle = F, logColBar = T, colors = colors, removeNonMutated= T,clinicalFeatures = c("molecular_subtype"))

somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
co_occurrence_df <- somaticInteractions(maf = maf) %>% 
  dplyr::filter(pValue < 0.05)
