################################################################################
# 01-plot_ont_vs_short_CLK1-Ex4_psi.R
# Plotting script that takes in Exon 4 PSI from short vs ONT reads
# written by Ammar Naqvi
#
# usage: Rscript 01-plot_ont_vs_short_CLK1-Ex4_psi.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("viridis")
  library("RColorBrewer")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "long-read-CLK1-validation")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output plot path
file_barplot = file.path(plots_dir,"isoform-stacked-barplot.pdf")

## retrive stringtie2 results for each cell line
cl_7316_1763_file = file.path(input_dir,"7316_1763.CLK1.processed.txt")
cl_7316_1769_file = file.path(input_dir,"7316_1769.CLK1.processed.txt")
cl_KNS42_file     = file.path(input_dir,"KNS42.CLK1.processed.txt")

## load KNS42 dataset
depmap_file = file.path(data_dir, "CLK1-CRISPR-DepMap-score.csv")
depmap_data <- vroom(depmap_file, show_col_types = FALSE) %>% 
  dplyr::rename("ModelID"=`Depmap ID`) %>% 
  filter( Lineage == 'CNS/Brain' ) 

depmap_data_KNS42 <- depmap_data %>% dplyr::filter(`Cell Line Name` =="KNS42")

omics_mappings_file <- file.path(data_dir,"OmicsDefaultModelProfiles.csv")
tpm_file <- file.path(data_dir,"OmicsExpressionTranscriptsTPMLogp1Profile.csv")

omics_id_mapping_df <- vroom(omics_mappings_file,show_col_types = FALSE) %>% 
  inner_join(depmap_data,by="ModelID")

depMap_transcr_expr <- vroom(tpm_file,show_col_types = FALSE) %>% 
  dplyr::rename("ProfileID"=`...1`) %>%
  inner_join(omics_id_mapping_df,by="ProfileID")

## load and store patient tumor cell line data
histology_file <- file.path(data_dir,'histologies.tsv')
rmats_file     <- file.path(data_dir, 'rMATS_merged.comparison.tsv.gz')

histology_df <- vroom(histology_file) %>%
  dplyr::filter(sample_id=='7316-1763' | sample_id=='7316-1769',
                experimental_strategy == 'RNA-Seq',
                composition == 'Derived Cell Line',
                RNA_library == 'stranded')  

rmats_df <- vroom(rmats_file) %>% 
            filter(splicing_case == 'SE',
                   sample_id %in% histology_df$Kids_First_Biospecimen_ID,
                   geneSymbol=="CLK1",
                   exonStart_0base=="200860124", 
                   exonEnd=="200860215",
                   row_number() <= n()-1
            ) 
psi_1763 <- rmats_df %>% dplyr::filter(sample_id =='BS_DRY58DTF') %>% dplyr::select(IncLevel2)
psi_1769 <- rmats_df %>% dplyr::filter(sample_id =='BS_40MP5BWR') %>% dplyr::select(IncLevel2)


## get psi for KNS42, 7316-1763, and 7316-1769  
CLK1_expr_KNS42 <- depMap_transcr_expr %>% 
  filter(ModelID=='ACH-000622') %>%
  dplyr::select(`CLK1 (ENST00000321356)`,`CLK1 (ENST00000409769)`,`CLK1 (ENST00000409769)` ,`CLK1 (ENST00000434813)`) %>% 
  mutate(Ex4_expr     = ( ( `CLK1 (ENST00000321356)`) / ( `CLK1 (ENST00000409769)` + `CLK1 (ENST00000434813)` + `CLK1 (ENST00000321356)`) ), 
         Non_Ex4_expr = ( ( `CLK1 (ENST00000409769)` + `CLK1 (ENST00000434813)` )  / (`CLK1 (ENST00000409769)` +`CLK1 (ENST00000434813)` + `CLK1 (ENST00000321356)`) )
  )

## format and process cell line information 
# STRG.1.2 is Ex.4 transcript by manual inspection
cl_7316_1763_df <- vroom(cl_7316_1763_file,comment = "#", 
                         delim = "\t", 
                         col_names = c("transcript","tpm"), show_col_types = FALSE) %>% 
                   dplyr::mutate(cell_line="7316-1763") %>%
                   dplyr::filter(transcript=='STRG.1.2'| transcript=='STRG.1.1' ) %>% 
                   dplyr::mutate(Isoform = case_when(
                                transcript == "STRG.1.2" ~ "Inclusion",
                                TRUE ~ "Skipping")) %>% 
                  dplyr::mutate(type="long") %>% 
                  dplyr::mutate(PSI = (tpm / sum(tpm)) * 100 ) %>% 
                  # manual add rMATs results for short reads
                  add_row( cell_line = "7316-1763", type = "short", Isoform="Inclusion", PSI = psi_1763$IncLevel2 * 100 ) %>%
                  add_row( cell_line = "7316-1763", type = "short", Isoform="Skipping", PSI = (1-psi_1763$IncLevel2) * 100  ) 

# STRG.2.2 is Ex.4 transcript by manual inspection
cl_7316_1769_df <- vroom(cl_7316_1769_file,comment = "#", 
                         delim = "\t", 
                         col_names = c("transcript","tpm"), show_col_types = FALSE) %>% 
                    mutate(cell_line="7316-1769") %>% 
                    dplyr::filter(transcript=='STRG.2.2'| transcript=='STRG.2.1' ) %>% 
                    dplyr::mutate(Isoform = case_when(
                                  transcript == "STRG.2.2" ~ "Inclusion",
                                  TRUE ~ "Skipping")) %>% 
                    dplyr::mutate(type="long") %>% 
                    dplyr::mutate(PSI = (tpm / sum(tpm)) * 100 ) %>% 
                    # manual add rMATs results for short reads
                    add_row( cell_line = "7316-1769", type = "short", Isoform="Inclusion", PSI = psi_1769$IncLevel2 * 100 ) %>%
                    add_row( cell_line = "7316-1769", type = "short", Isoform="Skipping", PSI = (1-psi_1769$IncLevel2) *100  )

## STRG 1.1 is Ex 4 transcript
cl_KNS42_df <- vroom(cl_KNS42_file,comment = "#", 
                         delim = "\t", 
                         col_names = c("transcript","tpm"), show_col_types = FALSE) %>% 
                         mutate(cell_line="KNS42")  %>% 
                         dplyr::mutate(Isoform = case_when(
                                       transcript == "STRG.1.1" ~ "Inclusion",
                                       TRUE ~ "Skipping")) %>% 
                        dplyr::mutate(type="long") %>% 
                        dplyr::mutate(PSI = (tpm / sum(tpm)) * 100 ) %>% 
                        # manual add rMATs results for short reads
                        add_row( cell_line = "KNS42", type = "short", Isoform="Skipping", PSI = CLK1_expr_KNS42$Non_Ex4_expr * 100)  %>%
                        add_row( cell_line = "KNS42", type = "short", Isoform="Inclusion", PSI = CLK1_expr_KNS42$Ex4_expr * 100 ) 


cell_lines_df <- rbind(cl_KNS42_df,cl_7316_1769_df,cl_7316_1763_df)

## make plot
stacked_barplot <- ggplot(cell_lines_df, aes(fill=Isoform, y=PSI, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#FFC20A","#0C7BDC")) +
  facet_wrap(~cell_line) +
  xlab("RNA-Seq sequencing strategy") + 
  theme_Publication()

pdf(file_barplot, 
    width = 7, height = 7)
stacked_barplot
dev.off()




