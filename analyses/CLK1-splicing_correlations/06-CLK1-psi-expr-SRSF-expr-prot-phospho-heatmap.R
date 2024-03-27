################################################################################
# 06-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R
# Generate heatmap of correlation coefficients between CLK1 features (PSI, Expr) and CLK1 and SRSF RNA, protein, and phosphoprotein expression
# written by Ryan Corbett
# Usage: Rscript 06-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggpubr")
  library("vroom")
  library("data.table")
  library("ComplexHeatmap")
  library("circlize")
  library("stringr")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
input_dir <- file.path(analysis_dir, "input")
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
cohort_file <- file.path(root_dir, "analyses", "cohort_summary",
                         "results", "histologies-plot-group.tsv")
rsem_counts <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")
clk1_psi_file <- file.path(results_dir, "clk1-exon4-psi-hgg.tsv")

cptac_proteo_file <- file.path(input_dir, "cptac-protein-imputed-prot-expression-abundance.tsv.gz")
hope_proteo_file <- file.path(input_dir, "hope-protein-imputed-prot-expression-abundance.tsv.gz")
hope_phospho_file <- file.path(input_dir, "hope-protein-imputed-phospho-expression-abundance.tsv.gz")

## read in histology, cohort, and independent specimens file
hist <- read_tsv(clin_file, guess_max = 10000)

cohort_df <- read_tsv(cohort_file)

indep_rna_df <- vroom(indep_rna_file) %>% 
  dplyr::filter(cohort == 'PBTA')

# extract HGG samples, retain stranded and polyA libraries, and filter for independent specimens
hist_rna <- hist %>% 
  filter(short_histology == "HGAT",
#         RNA_library == "stranded",
         cohort == "PBTA",
         Kids_First_Biospecimen_ID %in% indep_rna_df$Kids_First_Biospecimen_ID) %>%
  select(Kids_First_Biospecimen_ID, CNS_region, match_id, molecular_subtype) %>%
  left_join(cohort_df %>% dplyr::select(Kids_First_Biospecimen_ID, plot_group))

# filter expr df for hgg samples
rsem_df <- readRDS(rsem_counts) %>%
  select(hist_rna$Kids_First_Biospecimen_ID)

## load CLK1 psi file
clk1_rmats <- read_tsv(clk1_psi_file) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample_id) %>%
  left_join(hist_rna %>% dplyr::select(Kids_First_Biospecimen_ID, match_id))

# define proteomics and phosphoproteomics-specific histologies files
hist_proteo <- hist %>%
  dplyr::filter(short_histology == "HGAT",
                experimental_strategy %in% c("Whole Cell Proteomics")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, match_id) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_proteo = Kids_First_Biospecimen_ID)

hist_phospho <- hist %>%
  dplyr::filter(short_histology == "HGAT",
                experimental_strategy %in% c("Phospho-Proteomics")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, match_id) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_phospho = Kids_First_Biospecimen_ID)

# load proteomics and phosphoproteomics data
cptac_proteo <- read_tsv(cptac_proteo_file)
hope_proteo <- read_tsv(hope_proteo_file)
hope_phospho <- read_tsv(hope_phospho_file)

# prioritize HOPE samples over CPTAC
hist_proteo <- hist_proteo %>% 
  dplyr::mutate(proteo_cohort = case_when(
    Kids_First_Biospecimen_ID_proteo %in% colnames(hope_proteo) ~ "Hope",
    Kids_First_Biospecimen_ID_proteo %in% colnames(cptac_proteo) ~ "CPTAC"
  )) %>%
  dplyr::arrange(proteo_cohort) %>%
  distinct(match_id, .keep_all = T)

# filter phosphoproteomics samples
hist_phospho <- hist_phospho %>%
  dplyr::filter(Kids_First_Biospecimen_ID_phospho %in% colnames(hope_phospho))

# define lists of srsf and clk genes
srsf_list <- c("SRSF1", "SRSF2", "SRSF3", "SRSF4",
               "SRSF5", "SRSF6", "SRSF7", "SRSF8", "SRSF9",
               "SRSF10", "SRSF11")
clk_list <- c("CLK1", "CLK2", "CLK3", "CLK4")

# filter cptac proteo and phosphoproteo dfs for clk and srsf genes
cptac_proteo_df <- cptac_proteo %>%
  dplyr::select(-NP_id) %>%
  dplyr::filter(GeneSymbol %in% c(clk_list, srsf_list)) %>%
  gather(key = "Kids_First_Biospecimen_ID_proteo",
         value = "abundance",
         -GeneSymbol) %>%
  dplyr::filter(Kids_First_Biospecimen_ID_proteo %in% hist_proteo$Kids_First_Biospecimen_ID_proteo)

# Merge proteomics data frames and calculate z-scores
proteo_df <- hope_proteo %>%
  dplyr::select(-NP_id) %>%
  dplyr::filter(GeneSymbol %in% c(clk_list, srsf_list)) %>%
  gather(key = "Kids_First_Biospecimen_ID_proteo",
         value = "abundance",
         -GeneSymbol) %>%
  dplyr::filter(Kids_First_Biospecimen_ID_proteo %in% hist_proteo$Kids_First_Biospecimen_ID_proteo) %>%
  bind_rows(cptac_proteo_df) %>%
  group_by(GeneSymbol) %>%
  dplyr::mutate(proteo_zscore = scale(abundance)) %>%
  left_join(hist_proteo)

# calculate phosphoproteomics z-scores
phospho_df <- hope_phospho %>%
  dplyr::select(-NP_id, -Peptide_res_num,
                -Peptide_sequence) %>%
  dplyr::filter(GeneSymbol %in% c(clk_list, srsf_list)) %>%
  gather(key = "Kids_First_Biospecimen_ID_phospho",
         value = "abundance",
         -GeneSymbol, -Site) %>%
  dplyr::filter(Kids_First_Biospecimen_ID_phospho %in% hist_phospho$Kids_First_Biospecimen_ID_phospho)  %>%
  distinct(Kids_First_Biospecimen_ID_phospho, GeneSymbol, Site, .keep_all = TRUE) %>%
  group_by(GeneSymbol, Site) %>%
  dplyr::mutate(phospho_zscore = scale(abundance)) %>%
  left_join(hist_phospho) %>%
  dplyr::mutate(gene_site = glue::glue("{GeneSymbol}-{Site}"))

# create matrix of phosphorylation site z-scores with rownames as match ID
clk_srsf_phospho_df <- phospho_df %>%
  ungroup() %>%
  dplyr::filter(!grepl("NA", gene_site)) %>%
  dplyr::select(-abundance,
                -GeneSymbol,
                -Site, -Kids_First_Biospecimen_ID_phospho) %>%
  spread(gene_site, phospho_zscore)

# repeat for proteomics 
clk_srsf_proteo_df <- proteo_df %>%
  ungroup() %>%
  dplyr::mutate(GeneSymbol = glue::glue("{GeneSymbol} ")) %>%
  dplyr::select(GeneSymbol, match_id, proteo_zscore) %>%
  spread(GeneSymbol, proteo_zscore)
  
# merge CLK1 PSI and CLK and SRSF RNA, protein, and phosphoprotein expression
mol_df <- rsem_df %>%
  filter(rownames(.) %in% c(clk_list, srsf_list))  %>%
  rownames_to_column("GeneSymbol") %>%
  gather(key = "Kids_First_Biospecimen_ID", value = "Expr", -GeneSymbol) %>%
  dplyr::mutate(logExp = log(Expr, 2)) %>%
  left_join(hist_rna %>% dplyr::select(Kids_First_Biospecimen_ID, match_id)) %>%
  dplyr::select(-Expr, -Kids_First_Biospecimen_ID) %>%
  spread(GeneSymbol, logExp) %>%
  left_join(clk_srsf_proteo_df) %>%
  left_join(clk_srsf_phospho_df) %>%
  left_join(clk1_rmats)
  
# save df
write_tsv(mol_df, 
          file.path(results_dir, "hgg-dmg-clk-srsf-expression-phosphorylation.tsv"))

# Define DIPG/DMG and other HGG ids
dmg_ids <- hist_rna %>%
  dplyr::filter(plot_group == "DIPG or DMG") %>%
  pull(match_id)

hgg_ids <- hist_rna %>%
  dplyr::filter(plot_group == "Other high-grade glioma") %>%
  pull(match_id)

# extract molecular features from mol_df
features <- names(mol_df)[!names(mol_df) %in% c("match_id", "Kids_First_Biospecimen_ID",
                                                    "IncLevel1")]

# create empty df to store correlation coefficients between PSI and feature values
incl_cor_mat <- data.frame(row.names = features,
                           DMG = rep(0, length(features)),
                           HGG = rep(0, length(features)),
                           DMG_p = rep(0, length(features)),
                           HGG_p = rep(0, length(features)))  %>%
  dplyr::mutate(Feature = case_when(
    grepl("-", rownames(.)) ~ "Phospho-Protein z-score",
    grepl(" ", rownames(.)) ~ "Total Protein z-score",
    TRUE ~ "RNA log2Exp"
  ))

# create same empty df for CLK1 expression cor values
expr_cor_mat <- incl_cor_mat
phospho_cor_mat <- incl_cor_mat

# loop through subtypes and features to calculate pearson correlation coefficients
id_list <- list("DMG" = dmg_ids, "HGG" = hgg_ids)

for (subtype in c("DMG", "HGG")) {
  
  ids <- id_list[[subtype]]
  
  for (feature in rownames(incl_cor_mat)) {
    
    incl_cor_mat[feature, subtype] <- cor.test(mol_df$IncLevel1[mol_df$match_id %in% ids],
                                               unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                               method = "pearson")$estimate
    
    incl_cor_mat[feature, glue::glue("{subtype}_p")] <- cor.test(mol_df$IncLevel1[mol_df$match_id %in% ids],
                                                                 unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                                                 method = "pearson")$p.value
    
    expr_cor_mat[feature, subtype] <- cor.test(mol_df$CLK1[mol_df$match_id %in% ids],
                                               unlist(mol_df[mol_df$match_id %in% ids, feature]))$estimate
    
    expr_cor_mat[feature, glue::glue("{subtype}_p")] <- cor.test(mol_df$CLK1[mol_df$match_id %in% ids],
                                                                 unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                                                 method = "pearson")$p.value
    
    phospho_cor_mat[feature, subtype] <- cor.test(unlist(mol_df[mol_df$match_id %in% ids, colnames(mol_df) == "CLK1-S182"]),
                                               unlist(mol_df[mol_df$match_id %in% ids, feature]))$estimate
    
    phospho_cor_mat[feature, glue::glue("{subtype}_p")] <- cor.test(unlist(mol_df[mol_df$match_id %in% ids, colnames(mol_df) == "CLK1-S182"]),
                                                                 unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                                                 method = "pearson")$p.value  }
  
}

# define row_annot df to store feature type
row_annot <- expr_cor_mat %>%
  dplyr::select(Feature) %>%
  dplyr::mutate(Feature = factor(Feature,
                               levels = c("RNA log2Exp", "Total Protein z-score", "Phospho-Protein z-score")))

# add rownames
rownames(row_annot) <- rownames(expr_cor_mat)

# rename cor mat colnames for plotting
colnames(incl_cor_mat)[1:2] <- c("DIPG\nor DMG", "Other\nHGG")
colnames(expr_cor_mat)[1:2]  <- c("DIPG\nor DMG", "Other\nHGG")
colnames(phospho_cor_mat)[1:2]  <- c("DIPG\nor DMG", "Other\nHGG")

# create Heatmap row annotation for feature type
anno_col <- list(Feature = c("RNA log2Exp" = "#DC3220", "Total Protein z-score" = "#005AB5", "Phospho-Protein z-score" = "#40B0A6"))

row_anno = rowAnnotation(df = row_annot,
                         col = anno_col, show_legend = TRUE)

# define matrix indicating if correlation p < 0.05
incl_pval_mat <- ifelse(incl_cor_mat[,3:4] < 0.05, "*", "")
expr_pval_mat <- ifelse(expr_cor_mat[,3:4] < 0.05, "*", "")
phospho_pval_mat <- ifelse(phospho_cor_mat[,3:4] < 0.05, "*", "")

# create CLK1 ex4 psi correlation heatmap
incl_ht <- Heatmap(incl_cor_mat[,1:2],
                     name = "Pearson Correlation Coefficient",
                     col = colorRamp2(c(-1, 0, 1), c("#E66100", "white", "#5D3A9B")),
                     cluster_rows = FALSE,
                     row_split = row_annot$Feature,
                     column_gap = 0.5,
                     show_row_names = TRUE,
                   #  show_column_names = TRUE,
                     show_heatmap_legend=TRUE,
                     cluster_columns = FALSE,
                     right_annotation = row_anno,
                     row_title = NULL,
                     column_title = "CLK1 Exon4 PSI",
                     column_title_side = "top",
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(sprintf("%s", incl_pval_mat[i, j]), x, y, gp = gpar(fontsize = 14))
                     })

# create CLK1 expression correlation heatmap
expr_ht <- Heatmap(expr_cor_mat[,1:2],
                   name = "Pearson Correlation Coefficient",
                   col = colorRamp2(c(-1, 0, 1), c("#E66100", "white", "#5D3A9B")),
                   cluster_rows = FALSE,
                   row_split = row_annot$Feature,
                   column_gap = 0.5,
                   show_row_names = TRUE,
                   #  show_column_names = TRUE,
                   show_heatmap_legend=TRUE,
                   cluster_columns = FALSE,
                   right_annotation = row_anno,
                   row_title = NULL,
                   column_title = "CLK1 log2Expr",
                   column_title_side = "top",
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%s", expr_pval_mat[i, j]), x, y, gp = gpar(fontsize = 14))
                   })

# create CLK1 expression correlation heatmap
phospho_ht <- Heatmap(phospho_cor_mat[,1:2],
                   name = "Pearson Correlation Coefficient",
                   col = colorRamp2(c(-1, 0, 1), c("#E66100", "white", "#5D3A9B")),
                   cluster_rows = FALSE,
                   row_split = row_annot$Feature,
                   column_gap = 0.5,
                   show_row_names = TRUE,
                   #  show_column_names = TRUE,
                   show_heatmap_legend=TRUE,
                   cluster_columns = FALSE,
                   right_annotation = row_anno,
                   row_title = NULL,
                   column_title = "CLK1-S182 z-score",
                   column_title_side = "top",
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%s", phospho_pval_mat[i, j]), x, y, gp = gpar(fontsize = 14))
                   })

# save merged plot
pdf(file.path(plots_dir, "CLK1-psi-expr-correlation-heatmap.pdf"), width = 10, height = 12)
print(incl_ht + expr_ht + phospho_ht)
dev.off()
