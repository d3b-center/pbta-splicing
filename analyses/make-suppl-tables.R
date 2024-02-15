library(tidyverse)
library(openxlsx)
library(survival)
library(vroom)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses")
data_dir <- file.path(root_dir, "data")

# output directory for supplementary tables
supp_tables_dir <- file.path(analysis_dir, "supp")
if(!dir.exists(supp_tables_dir)){
  dir.create(supp_tables_dir, recursive=TRUE)
}

# define input files
histology_file <- file.path(data_dir, "histologies.tsv")
histology_ei_splice_events <- file.path(analysis_dir, "histology-specific-splicing", "results", "unique_events-ei.tsv")
histology_es_splice_events <- file.path(analysis_dir, "histology-specific-splicing", "results", "unique_events-es.tsv")
CNS_match_json <- file.path(analysis_dir, "input", "CNS_primary_site_match.json")
deseq2_sf_file <- file.path(analysis_dir, "splicing-factor_dysregulation", "results", "diffSFs_sig_genes.txt")
func_sites_es_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.total.HGG.neg.intersectUnip.ggplot.txt") 
func_sites_ei_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.total.HGG.pos.intersectUnip.ggplot.txt") 
kinase_func_sites_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "kinases-functional_sites.txt")
deseq2_morph_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","ctrl_vs_treated.de.formatted.tsv")
rmats_tsv_file <- file.path(data_dir,"ctrl-vs-morpholino-merged-rmats.tsv")

# define suppl output files and sheet names, when appropriate
table_s1_file <- file.path(supp_tables_dir, "TableS1-histologies.xlsx")
table_s2_file <- file.path(supp_tables_dir, "TableS2-histology-specific-splice-events.xlsx")
table_s3_file <- file.path(supp_tables_dir, "TableS3-DeSeq2-sbi-SFs.xlsx")
table_s4_file <- file.path(supp_tables_dir, "TableS4-functional-sites.xlsx")
table_s5_file <- file.path(supp_tables_dir, "TableS5-CLK1-ex4-splicing-impact-morpholino.xlsx")

## write table for histologies
# Sheet 1: README tab
histology_df <- vroom(histology_file)
histology_df_tumor <- histology_df %>%
  filter(!is.na(pathology_diagnosis)) %>%
  mutate(broad_histology = case_when(broad_histology == "Ependymal tumor" ~ "Ependymoma",
                                     TRUE ~ broad_histology))

readme <- tribble(
  ~`Histology column`,~Definition,~`Possible values`,
  "age_at_diagnosis_days","Patient age at diagnosis in days","numeric",
  "age_last_update_days","Patient age at the last clinical event/update in days","numeric",
  "aliquot_id","External aliquot identifier","alphanumeric",
  "broad_histology","Broad WHO 2016 classification of cancer type",paste(unique(histology_df_tumor$broad_histology), collapse = "; "),
  "cancer_group","Harmonized cancer groupings for plots",paste(unique(histology_df$cancer_group), collapse = "; "),
  "cancer_predispositions","Reported cancer predisposition syndromes",paste(unique(histology_df$cancer_predispositions), collapse = "; "),
  "CNS_region","Harmonized brain region based on `primary_site`",paste(unique(histology_df$CNS_region), collapse = "; "),
  "cohort","Scientific cohort",paste(unique(histology_df$cohort), collapse = "; "),
  "cohort_participant_id","Scientific cohort participant ID","C#####-C######",
  "composition","Sample composition",paste(unique(histology_df$composition), collapse = "; "),
  "ethnicity","Patient reported ethnicity",paste(unique(histology_df$ethnicity), collapse = "; "),
  "experimental_strategy","Sequencing strategy",paste(unique(histology_df$experimental_strategy), collapse = "; "),
  # leaving this non-programmatic because of the duplicates that would come up (eg two selections in one patient, needing data cleanup)
  "extent_of_tumor_resection","Amount of tumor resected at time of surgical event","Biopsy only;Partial resection;Gross/Near total resection;Not Reported;Unavailable",
  "germline_sex_estimate","Predicted sex of patient based on germline X and Y ratio calculation (described in methods)",paste(unique(histology_df$germline_sex_estimate), collapse = "; "),
  "harmonized_diagnosis","`integrated_diagnosis` if exists or updated and harmonized diagnosis using pathology_free_text_diagnosis information","text",
  "integrated_diagnosis","2016 WHO diagnosis integrated from pathology diagnosis and molecular subtyping","text",
  "Kids_First_Biospecimen_ID","KidsFirst biospecimen identifier","BS_########",
  "Kids_First_Participant_ID","KidsFirst patient identifier","PT_########",
  "molecular_subtype","Molecular subtype defined by WHO 2016 guidelines","text",
  "normal_fraction","Theta2 normal DNA fraction estimate","numeric",
  "Notes","Free text field describing changes from `pathology_diagnosis` to `integrated_diagnosis` or manner in which molecular_subtype was determined","text",
  "OS_days","Overall survival in days","numeric",
  "OS_status","Overall survival status",paste(unique(histology_df$OS_status), collapse = "; "),
  "pathology_diagnosis","Reported and/or harmonized patient diagnosis from pathology reports","text",
  "pathology_free_text_diagnosis","Free text patient diagnosis from pathology reports","text",
  "PFS_days","Progression-free survival in days","numeric",
  "primary_site","Bodily site(s) from which specimen was derived","text",
  "race","Patient reported race",paste(unique(histology_df$race), collapse = "; "),
  "reported_gender","Patient reported gender",paste(unique(histology_df$reported_gender), collapse = "; "),
  "RNA_library","Type of RNA-Sequencing library preparation",paste(unique(histology_df$RNA_library), collapse = "; "),
  "sample_id","External biospecimen identifier","alphanumeric",
  "sample_type","Broad sample type",paste(unique(histology_df$sample_type), collapse = "; "),
  "seq_center","Sequencing center",paste(unique(histology_df$seq_center), collapse = "; "),
  "short_histology","Abbreviated `integrated_diagnosis` or `broad_histology` for plotting purposes",paste(unique(histology_df$short_histology), collapse = "; "),
  "tumor_descriptor","Phase of therapy from which tumor was derived",paste(unique(histology_df$tumor_descriptor), collapse = "; "),
  "tumor_fraction","Theta2 tumor DNA fraction estimate","numeric",
  "tumor_ploidy","Control-FREEC ploidy","numeric"
)

# Sheet 2: Histologies file (histology_df)

# Sheet 3: CNS region definition based on definitions from [Cassie Kline](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/tables/input/CNS_primary_site_match.json) with additional manual review of [HGG primary_site](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1025).
# This is integrated upstream of molecular subtyping

# CNS_region ~ primary_site matches
cns_regions_df <- purrr::imap_dfr(jsonlite::fromJSON(CNS_match_json),
                                  function(x, name) { tibble::tibble(CNS_region = name, primary_site = paste0(x, collapse = ';')) })


# Combine and output
list_s1_table <- list(README = readme,
                      histologies_file = histology_df,
                      CNS_region_definitions = cns_regions_df
)
write.xlsx(list_s1_table,
           table_s1_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 2 histology specific splicing events
## sheet 1, exon inclusion splicing
ei_events_df <- vroom(histology_ei_splice_events)

## sheet 2, exon inclusion splicing
es_events_df <- vroom(histology_es_splice_events)

# Combine and output
list_s2_table <- list(exon_inclusion = ei_events_df,
                      exon_skipping = es_events_df)

write.xlsx(list_s2_table,
           table_s2_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 3 DeSeq 2 results comparing high vs low SBI HGG tumors
## sheet 1, exon inclusion splicing
deseq_df <- vroom(deseq2_sf_file)


# Combine and output
list_s3_table <- list(deseq2=deseq_df)

write.xlsx(list_s3_table,
           table_s3_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 4 Differential HGG splicing events impacting functional sites
## sheet 1, exon inclusion events
ei_events_df <- vroom(func_sites_ei_file)

## sheet 2, exon skipping events
es_events_df <- vroom(func_sites_es_file)

##sheet 3 kinases
kinase_events_df <- vroom(kinase_func_sites_file)


# Combine and output
list_s4_table <- list(exon_inclusion = ei_events_df,
                      exon_skipping = es_events_df
                      kinases= kinase_events_df)

write.xlsx(list_s4_table,
           table_s4_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 5 morpholino vs ctrl DESeq2 and rMATs
deseq2_morpholino_df <- vroom(deseq2_morph_file)
rmats_df <-  vroom(rmats_tsv_file)

list_s5_table <- list(deseq2_morp = deseq2_morpholino_df,
                      rmats = rmats_df)

write.xlsx(list_s5_table,
           table_s5_file,
           overwrite=TRUE,
           keepNA=TRUE)

