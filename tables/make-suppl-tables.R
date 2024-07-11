library(tidyverse)
library(openxlsx)
<<<<<<< HEAD
library(survival)
=======
>>>>>>> main
library(vroom)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses")
table_dir <- file.path(root_dir, "tables")
data_dir <- file.path(root_dir, "data")
input_dir <- file.path(table_dir, "input")


# output directory for supplementary tables
supp_tables_dir <- file.path(table_dir, "output")
if(!dir.exists(supp_tables_dir)){
  dir.create(supp_tables_dir, recursive=TRUE)
}

# define input files
histology_file <- file.path(data_dir, "histologies.tsv")
histology_ei_splice_events <- file.path(analysis_dir, "histology-specific-splicing", "results", "unique_events-ei.tsv")
histology_es_splice_events <- file.path(analysis_dir, "histology-specific-splicing", "results", "unique_events-es.tsv")
optimal_cluster_tsv <- file.path(analysis_dir, "clustering_analysis", "output", "optimal_clustering", "lspline_output.tsv")
cluster_membership <- file.path(analysis_dir, "clustering_analysis", "output", "cluster_members_by_cancer_group_subtype.tsv")
CNS_match_json <- file.path(table_dir, "input", "CNS_primary_site_match.json")
<<<<<<< HEAD
deseq2_sf_file <- file.path(analysis_dir, "splicing-factor_dysregulation", "results", "diffSFs_sig_genes.txt")
func_sites_es_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.total.HGG.neg.intersectUnip.ggplot.txt") 
func_sites_ei_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.total.HGG.pos.intersectUnip.ggplot.txt") 
kinase_func_sites_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "kinases-functional_sites.tsv")
deseq2_morph_file <- file.path(root_dir,"analyses/CLK1-splicing-impact-morpholino","results","ctrl_vs_treated.de.tsv")
rmats_tsv_file <- file.path(data_dir,"ctrl-vs-morpholino-merged-rmats.tsv")
func_sites_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.intersectUnip.ggplot.txt")
=======
sf_list_file <- file.path(root_dir, "analyses","splicing-factor_dysregulation", "input","splicing_factors.txt")
hugo_file <- file.path(root_dir, "analyses", "oncoprint", "input", "hgnc-symbol-check.csv")
deseq2_sf_file <- file.path(analysis_dir, "splicing-factor_dysregulation", "results", "all_hgg-diffSFs_sig_genes.txt")
func_sites_es_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.SE.total.neg.intersectunip.ggplot.txt") 
func_sites_ei_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "splicing_events.SE.total.pos.intersectunip.ggplot.txt") 
kinase_func_sites_file <- file.path(analysis_dir, "splicing_events_functional_sites", "results", "kinases-functional_sites.tsv")
deseq2_morph_file <- file.path(root_dir,"analyses/CLK1-splicing-impact-morpholino","results","ctrl_vs_treated.de.tsv")
rmats_tsv_file <- file.path(data_dir,"ctrl-vs-morpholino-merged-rmats.tsv")
func_sites_SE_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.SE.intersectUnip.ggplot.txt")
func_sites_A5SS_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.A5SS.intersectUnip.ggplot.txt")
func_sites_A3SS_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.A3SS.intersectUnip.ggplot.txt")
func_sites_RI_morpho_tsv_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results","splicing_events.morpho.RI.intersectUnip.ggplot.txt")


func_sites_goi_file <- file.path(analysis_dir,"CLK1-splicing-impact-morpholino","results", "differential_splice_by_goi_category.tsv")
>>>>>>> main
primers_file <-  file.path(input_dir,"primers.tsv")

# define suppl output files and sheet names, when appropriate
table_s1_file <- file.path(supp_tables_dir, "TableS1-histologies.xlsx")
table_s2_file <- file.path(supp_tables_dir, "TableS2-histology-specific-splice-events.xlsx")
<<<<<<< HEAD
table_s3_file <- file.path(supp_tables_dir, "TableS3-DeSeq2-sbi-SFs.xlsx")
=======
table_s3_file <- file.path(supp_tables_dir, "TableS3-SF-dysreg.xlsx")
>>>>>>> main
table_s4_file <- file.path(supp_tables_dir, "TableS4-functional-sites.xlsx")
table_s5_file <- file.path(supp_tables_dir, "TableS5-CLK1-ex4-splicing-impact-morpholino.xlsx")

## write table for histologies
# Sheet 1: README tab
histology_df <- read_tsv(histology_file, guess_max = 100000)
histology_df_tumor <- histology_df %>%
  filter(!is.na(pathology_diagnosis))

readme <- tribble(
  ~`Histology column`,~Definition,~`Possible values`,
  "age_at_chemo_start","Patient age at chemotherapy start in days","numeric",
  "age_at_diagnosis_days","Patient age at diagnosis in days","numeric",
  "age_at_event_days","Patient age at sample collection event in days","numeric",
  "age_at_radiation_start","Patient age at radiation start in days","numeric",
  "age_last_update_days","Patient age at the last clinical event/update in days","numeric",
  "aliquot_id","External aliquot identifier","alphanumeric",
  "broad_histology","Broad WHO classification of cancer type",paste(unique(histology_df$broad_histology), collapse = "; "),
  "cancer_group","Harmonized cancer groupings for plots",paste(unique(histology_df$cancer_group), collapse = "; "),
  "cancer_predispositions","Reported cancer predisposition syndromes",paste(unique(histology_df$cancer_predispositions), collapse = "; "),
  "cell_line_composition","Cell line media",paste(unique(histology_df$cell_line_composition), collapse = "; "),
  "cell_line_passage","Cell line passage at collection","numeric",
  "clinical_status_at_event","Patient status at the time of sample collection", paste(unique(histology_df$clinical_status_at_event), collapse = "; "),
  "CNS_region","Harmonized brain region based on `primary_site`",paste(unique(histology_df$CNS_region), collapse = "; "),
  "cohort","Scientific cohort",paste(unique(histology_df$cohort), collapse = "; "),
  "cohort_participant_id","Scientific cohort participant ID","C#####-C######",
  "composition","Sample composition",paste(unique(histology_df$composition), collapse = "; "),
  "dkfz_v11_methylation_subclass","v11b6 DKFZ methylation-based CNS tumor subclass","text",
  "dkfz_v11_methylation_subclass_score","v11b6 DKFZ methylation-based CNS tumor subclass score","numeric",
  "dkfz_v12_methylation_subclass","v12b6 DKFZ methylation-based CNS tumor subclass score","text",
  "dkfz_v12_methylation_subclass_score","v12b6 DKFZ methylation-based CNS tumor subclass","numeric",
  "dkfz_v12_methylation_mgmt_status","v12b6 DKFZ MGMT promoter methylation status",paste(unique(histology_df$dkfz_v11_methylation_subclass), collapse = "; "),
  "dkfz_v12_methylation_mgmt_estimated","v12b6 DKFZ MGMT promoter methylation fraction","numeric",
  "EFS_days","Event-free survival in days","numeric",
  "EFS_event_type", "Event considered when calculating EFS", paste(unique(histology_df$EFS_event_type), collapse = "; "),
  "ethnicity","Patient reported ethnicity",paste(unique(histology_df$ethnicity), collapse = "; "),
  "experimental_strategy","Sequencing strategy",paste(unique(histology_df$experimental_strategy), collapse = "; "),
  # leaving this non-programmatic because of the duplicates that would come up (eg two selections in one patient, needing data cleanup)
  "extent_of_tumor_resection","Amount of tumor resected at time of surgical event","Biopsy only;Partial resection;Gross/Near total resection;Not Reported;Unavailable",
  "germline_sex_estimate","Predicted sex of patient based on germline X and Y ratio calculation (described in methods)",paste(unique(histology_df$germline_sex_estimate), collapse = "; "),
  "gtex_group","Tissue Type",paste(unique(histology_df$gtex_group), collapse = "; "),
  "gtex_subgroup","Tissue Subtype",paste(unique(histology_df$gtex_subgroup), collapse = "; "),
  "harmonized_diagnosis","`integrated_diagnosis` if exists or updated and harmonized diagnosis using pathology_free_text_diagnosis information","text",
  "integrated_diagnosis","WHO 2021 diagnosis integrated from pathology diagnosis and molecular subtyping","text",
  "Kids_First_Biospecimen_ID","Biospecimen identifier, Kids First or other cohort","BS_########",
  "Kids_First_Participant_ID","Patient identifier, Kids First or other cohort","PT_########",
  "match_id", "ID used to match experimental strategies within an event per sample composition", "Concatenation of sample_id, tumor descriptor, composition, and cell line composition and passage if applicable",
  "molecular_subtype","Molecular subtype defined by WHO 2021 guidelines","text",
  "molecular_subtype_methyl","DKFZ v12b6 or NIH v2 methylation class aligned to WHO 2021 subtypes","text",
  "normal_fraction","Theta2 normal DNA fraction estimate","numeric",
  "Notes","Free text field describing changes from `pathology_diagnosis` to `integrated_diagnosis` or manner in which molecular_subtype was determined","text",
  "OS_days","Overall survival in days","numeric",
  "OS_status","Overall survival status",paste(unique(histology_df$OS_status), collapse = "; "),
  "pathology_diagnosis","Reported and/or harmonized patient diagnosis from pathology reports","text",
  "pathology_free_text_diagnosis","Free text patient diagnosis from pathology reports","text",
  "primary_site","Bodily site(s) from which specimen was derived","text",
  "race","Patient reported race",paste(unique(histology_df$race), collapse = "; "),
  "reported_gender","Patient reported gender",paste(unique(histology_df$reported_gender), collapse = "; "),
  "RNA_library","Type of RNA-Sequencing library preparation",paste(unique(histology_df$RNA_library), collapse = "; "),
  "sample_id","Event id","alphanumeric",
  "sample_type","Broad sample type",paste(unique(histology_df$sample_type), collapse = "; "),
  "seq_center","Sequencing center",paste(unique(histology_df$seq_center), collapse = "; "),
  "short_histology","Abbreviated `cancer_group` or `broad_histology` for plotting purposes",paste(unique(histology_df$short_histology), collapse = "; "),
  "sub_cohort", "sub-cohort", paste(unique(histology_df$sub_cohort), collapse = "; "),
  "tumor_descriptor","Phase of therapy from which tumor was derived",paste(unique(histology_df$tumor_descriptor), collapse = "; "),
  "tumor_fraction","Theta2 tumor DNA fraction estimate","numeric",
  "tumor_fraction_LUMP","LUMP tumor DNA fraction estimate from methylation","numeric",
  "tumor_fraction_RFpurify_ABSOLUTE","RFpurify ABSOLUTE tumor DNA fraction estimate from methylation","numeric",
  "tumor_fraction_RFpurify_ESTIMATE","RFpurify ESTIMATE tumor DNA fraction estimate from methylation","numeric",
  "tumor_ploidy","Control-FREEC ploidy estimate","numeric"
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

## sheet 3, optimal clustering output/lspline
opt_cluster_df <- read_tsv(optimal_cluster_tsv)

# subset columns used for scoring methods
opt_cluster_df <- opt_cluster_df %>%
  dplyr::select(algorithm,
                distance,
                feature_selection,
                k,
                stretch,
                delta_auc,
                p_val,
                average.between,
                average.within,
                within.cluster.ss,
                dunn,
                entropy,
                avg_sil,
                cluster_qual,
                rank)
<<<<<<< HEAD
    
=======

>>>>>>> main
## sheet 4, cluster membership
cluster_membership_df <- read_tsv(cluster_membership)

# Combine and output
list_s2_table <- list(exon_inclusion = ei_events_df,
                      exon_skipping = es_events_df,
                      opt_cluster = opt_cluster_df,
                      clust_memb = cluster_membership_df)

write.xlsx(list_s2_table,
           table_s2_file,
           overwrite=TRUE,
           keepNA=TRUE)

<<<<<<< HEAD
## Table 3 DeSeq 2 results comparing high vs low SBI HGG tumors
## sheet 1, exon inclusion splicing
=======
## Table 3 Splicing factor dysregulation
## sheet 1, exon inclusion splicing
sf_list <- readLines(sf_list_file) 
hugo_list <- read_csv(hugo_file, skip = 1) %>%
           pull(`Approved symbol`) 
goi_list <- c(sf_list, hugo_list) %>%
           unique()

>>>>>>> main
deseq_df <- vroom(deseq2_sf_file) %>%
  dplyr::select(-Significant)

# Combine and output
<<<<<<< HEAD
list_s3_table <- list(deseq2=deseq_df)
=======
list_s3_table <- list(splicing_factors_and_splicosome=goi_list, deseq2=deseq_df)
>>>>>>> main

write.xlsx(list_s3_table,
           table_s3_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 4 Differential HGG splicing events impacting functional sites
## sheet 1, exon inclusion events
<<<<<<< HEAD
ei_events_df <- vroom(func_sites_ei_file)
=======
ds_events_es_df <- vroom(func_sites_es_file)
ds_events_ei_df <- vroom(func_sites_ei_file)
ds_events_A3SS_df <- vroom(func_sites_A3SS_morpho_tsv_file)
ds_events_RI_df <- vroom(func_sites_RI_morpho_tsv_file)

>>>>>>> main

## sheet 2, exon skipping events
es_events_df <- vroom(func_sites_es_file)

##sheet 3 kinases
kinase_events_df <- vroom(kinase_func_sites_file)


<<<<<<< HEAD
# Combine and output
list_s4_table <- list(exon_inclusion = ei_events_df,
                      exon_skipping = es_events_df,
=======

# Combine and output
list_s4_table <- list(ds_skipping = ds_events_es_df,
                      ds_inclusion = ds_events_ei_df,
>>>>>>> main
                      kinases= kinase_events_df)

write.xlsx(list_s4_table,
           table_s4_file,
           overwrite=TRUE,
           keepNA=TRUE)

## Table 5 morpholino vs ctrl DESeq2 and rMATs results
deseq2_morpholino_df <- vroom(deseq2_morph_file) %>%
  filter(padj < 0.05) %>%
  dplyr::select(Gene_Symbol,	
                baseMean,	
                log2FoldChange,	
                lfcSE,	
                stat,	
                pvalue,	
                padj)
<<<<<<< HEAD
  
rmats_df <-  vroom(rmats_tsv_file)
ds_func_df <- vroom(func_sites_morpho_tsv_file)
=======

rmats_df <-  vroom(rmats_tsv_file)
ds_events_SE_df <- vroom(func_sites_SE_morpho_tsv_file)
ds_events_A5SS_df <- vroom(func_sites_A5SS_morpho_tsv_file)
ds_events_A3SS_df <- vroom(func_sites_A3SS_morpho_tsv_file)
ds_events_RI_df <- vroom(func_sites_RI_morpho_tsv_file)
cancer_genes_func_df <- vroom(func_sites_goi_file)

>>>>>>> main
primers_df <- vroom(primers_file, delim = "\t")

list_s5_table <- list(deseq2_morp = deseq2_morpholino_df,
                      rmats = rmats_df,
<<<<<<< HEAD
                      func_sites = ds_func_df,
=======
                      ds_SE = ds_events_SE_df,
                      ds_A5SS = ds_events_A5SS_df,
                      ds_A3SS = ds_events_A3SS_df,
                      ds_RI = ds_events_RI_df,
                      ds_cancer_genes = cancer_genes_func_df,
>>>>>>> main
                      primers = primers_df)

write.xlsx(list_s5_table,
           table_s5_file,
           overwrite=TRUE,
           keepNA=TRUE)
<<<<<<< HEAD

=======
>>>>>>> main
