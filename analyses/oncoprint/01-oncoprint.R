################################################################################
# 01-oncoprint.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 01-oncoprint.R 
################################################################################

## libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(maftools)
  library(vroom)
  library('data.table')
  library(ComplexHeatmap)
  library(circlize)
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "oncoprint")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")
results_dir   <- file.path(analysis_dir, "results")

# Load functions
#source(file.path(analysis_dir, "util/cooccur_functions.R"))

## check and create plots/results dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## input files
maf_df <- data.table::fread(file.path(data_dir,"snv-consensus-plus-hotspots.maf.tsv.gz"), data.table = FALSE)

clin_file <- file.path(data_dir, "histologies.tsv")
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_wgs_file <- file.path(data_dir, "independent-specimens.wgswxspanel.primary.prefer.wgs.tsv")
diff_psi_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")
goi_file <- file.path(input_dir,"oncoprint-goi-lists-OpenPedCan-gencode-v39.csv")

##output file
plot_out <- file.path(plots_dir,"oncoprint.pdf")

goi_df <- vroom(goi_file) %>% dplyr::select(HGAT) %>% dplyr::rename(gene=HGAT)

indep_rna_df <- vroom(indep_rna_file)  %>% 
  dplyr::filter(cohort == 'PBTA') 

indep_wgs_df <- vroom(indep_wgs_file)  %>% 
  dplyr::filter(cohort == 'PBTA') 

## filter for samples that have both RNA and WGS
indep_sample_rna_wgs <- inner_join(indep_wgs_df,indep_rna_df, by='Kids_First_Participant_ID') %>% 
  dplyr::rename('DNA_id'=Kids_First_Biospecimen_ID.x, 'RNA_id'=Kids_First_Biospecimen_ID.y)

histologies_df <- vroom(clin_file)  %>% 
  dplyr::filter(short_histology == 'HGAT',
         cohort == 'PBTA',
        # CNS_region == 'Midline'
        ) %>% 
  inner_join(indep_sample_rna_wgs, by='Kids_First_Participant_ID') %>% 
  dplyr::filter(short_histology == 'HGAT',
         cohort == 'PBTA',
         #CNS_region == 'Midline'
  ) %>%
  dplyr::filter(experimental_strategy == 'RNA-Seq' |
         experimental_strategy == 'WGS' |
          experimental_strategy == 'WXS')


splice_CLK1_df <-  fread(diff_psi_file) %>%
   dplyr::filter(geneSymbol=="CLK1",
                                 exonStart_0base=="200860124", 
                                 exonEnd=="200860215") %>% 
  dplyr::rename('RNA_id'=sample_id,
                'gene'=geneSymbol) %>% 
  inner_join(histologies_df, by='RNA_id', relationship = "many-to-many") 





## filter maf table for samples with RNA splicing + HGGs + midline + req'd cols
maf_df <- maf_df  %>% 
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))

## select required columns only
maf_df_fil <- maf_df %>%
  inner_join(histologies_df, by= c("Tumor_Sample_Barcode" ="DNA_id"),relationship = "many-to-many") %>% 
  dplyr::mutate('Tumor_Sample_Barcode'=Kids_First_Participant_ID) %>% 
  dplyr::filter(vaf > 0.05) %>% 
  dplyr::filter(Hugo_Symbol %in% goi_df$gene) %>% 
  dplyr::filter(Kids_First_Participant_ID  %in% splice_CLK1_df$Kids_First_Participant_ID) %>%
  dplyr::select(
    "Hugo_Symbol", 
    "Chromosome", 
    "Start_Position", 
    "End_Position", 
    "Reference_Allele", 
    "Tumor_Seq_Allele2", 
    "Variant_Classification", 
    "Variant_Type",
    "Tumor_Sample_Barcode"
  ) 

## color for barplot
source(file.path(input_dir, "mutation-colors.R"))


collapse_snv_dat <- dplyr::select(maf_df_fil,c(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification)) %>%
  #dplyr::unique() %>%
  dplyr::filter(Variant_Classification %in% names(colors)) %>%
  dplyr::group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  dplyr::summarise(count = as.double(length(Variant_Classification[!is.na(Variant_Classification)])),
            Variant_Classification=paste(Variant_Classification,collapse = ",")) 


# complex heatmap
gene_matrix<-reshape2::acast(collapse_snv_dat,
                             Hugo_Symbol~Tumor_Sample_Barcode,value.var = "Variant_Classification") %>%
  as.data.frame() %>%
  dplyr::mutate_if(is.character, ~replace_na(.,""))

# read in PSI
splice_CLK1_df <- splice_CLK1_df %>%
  dplyr::rename(PSI = IncLevel1) %>%
  dplyr::select(Kids_First_Participant_ID, PSI)

# mutate the hgg dataframe for plotting
histologies_df_sorted <- histologies_df %>%
  inner_join(splice_CLK1_df, by = "Kids_First_Participant_ID") %>%
  unique() %>%
  column_to_rownames("Kids_First_Biospecimen_ID") %>%
  distinct(RNA_id, .keep_all = TRUE) %>%
  arrange(PSI) %>%
  dplyr::rename("Gender"=reported_gender,
                "Molecular Subtype"=molecular_subtype,
                "CNS Region"=CNS_region, 
                "CLK1 Ex4 PSI"=PSI)

ha = HeatmapAnnotation(name = "annotation", df = histologies_df_sorted[histologies_df_sorted$Kids_First_Participant_ID %in% colnames(gene_matrix),c("Gender","Molecular Subtype","CLK1 Ex4 PSI")],
                       col=list(
                         "Gender" = c("Male" = "#56B4E9","Female" = "#CC79A7","Not Reported" = "white"),
                         "PSI" = colorRamp2(c(0, .5, 1.0), c("whitesmoke", "#CAE1FF","#0072B2")),
                         annotation_name_side = "right", 
                         annotation_name_gp = gpar(fontsize = 9),
                         na_col = "whitesmoke") )


col = colors
df = histologies_df_sorted[,c("Kids_First_Participant_ID", "Gender","CLK1 Ex4 PSI")]

colorder <- df$Kids_First_Participant_ID

plot_oncoprint <- oncoPrint(gene_matrix, get_type = function(x) strsplit(x, ",")[[1]],
          column_names_gp = gpar(fontsize = 9), show_column_names = F, #show_row_barplot = F,
          alter_fun = list(
            background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "whitesmoke",col="whitesmoke")),
            Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Missense_Mutation"]),col = NA)),
            Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonsense_Mutation"]),col = NA)),
            Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Del"]), col = NA)),
            Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Ins"]), col = NA)),
            Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Splice_Site"]), col = NA)),
            Translation_Start_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Translation_Start_Site"]), col = NA)),
            Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonstop_Mutation"]),col = NA)),
            In_Frame_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Del"]),col = NA)),
            In_Frame_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Ins"]), col = NA)),
            Stop_Codon_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Stop_Codon_Ins"]), col = NA)),
            Start_Codon_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Start_Codon_Del"]), col = NA)),
            Fusion = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Fusion"]),col = NA)),
            Multi_Hit_Fusion  = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Multi_Hit_Fusion"]),col = NA)),
            Multi_Hit = function(x, y, w, h) grid.rect(x, y, w*0.75, h*0.85, gp = gpar(fill = unname(col["Multi_Hit"]), col = NA)),
            Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Del"]), col = NA)),
            Amp = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Amp"]), col = NA)),
            `5'Flank` = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["5'Flank"]), col = NA))),
          col = col,
          top_annotation = ha,
          alter_fun_is_vectorized = FALSE,
          #bottom_annotation = ha1,
          column_order =  colnames(gene_matrix) )
                  
# Save plot as PDF
pdf(plot_out, width = 20, height = 10)
plot_oncoprint
dev.off()