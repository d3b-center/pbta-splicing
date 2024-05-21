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
cohort_file <- file.path(root_dir, "analyses", "cohort_summary",
                         "results", "histologies-plot-group.tsv")

cptac_proteo_file <- file.path(data_dir, "cptac-protein-imputed-prot-expression-abundance.tsv.gz")
hope_proteo_file <- file.path(data_dir, "hope-protein-imputed-prot-expression-abundance.tsv.gz")
hope_phospho_file <- file.path(data_dir, "hope-protein-imputed-phospho-expression-abundance.tsv.gz")

data_file <- file.path(results_dir, "clk1-nf1-psi-exp-df.rds")

## read in histology, cohort, and independent specimens file
cohort_df <- read_tsv(cohort_file)
indep_rna_df <- read_tsv(indep_rna_file)

data_df <- readRDS(data_file)


# define lists of nf1 and clk genes
goi <- c("CLK1", "NF1")
clk1_trans_list <- c("CLK1-201", "Total CLK1")
clk1_splice_list <- "CLK1-Exon4_PSI"
nf1_trans_list <- c("Total NF1", "NF1-202_PC", "NF1-215_RI",  "NF1-208_NMD")
nf1_splice_list <- c("NF1-Exon23a_PSI", "NF1-215_PSI")

all_lists <- c(clk1_trans_list, clk1_splice_list, nf1_trans_list, nf1_splice_list)

# extract HGG samples, retain stranded libraries, and filter for independent specimens
hist_rna <- cohort_df %>% 
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma"),
         RNA_library == "stranded",
         cohort == "PBTA",
         Kids_First_Biospecimen_ID %in% indep_rna_df$Kids_First_Biospecimen_ID) %>%
  select(Kids_First_Biospecimen_ID, RNA_library, CNS_region, match_id, molecular_subtype, plot_group)

# select only HGG/DIPG/DMG
data_df <- readRDS(data_file) %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  filter(Kids_First_Biospecimen_ID %in% hist_rna$Kids_First_Biospecimen_ID) %>%
  left_join(hist_rna[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  select(match_id, any_of(all_lists))

# define proteomics and phosphoproteomics-specific histologies files
hist_proteo <- cohort_df %>%
  dplyr::filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma"),
                experimental_strategy %in% c("Whole Cell Proteomics")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, match_id) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_proteo = Kids_First_Biospecimen_ID)

hist_phospho <- cohort_df %>%
  dplyr::filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma"),
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

# filter cptac proteo and phosphoproteo dfs for clk and NF1 genes
cptac_proteo_df <- cptac_proteo %>%
  dplyr::select(-NP_id) %>%
  dplyr::filter(GeneSymbol %in% goi) %>%
  gather(key = "Kids_First_Biospecimen_ID_proteo",
         value = "abundance",
         -GeneSymbol) %>%
  dplyr::filter(Kids_First_Biospecimen_ID_proteo %in% hist_proteo$Kids_First_Biospecimen_ID_proteo)

# Merge proteomics data frames and calculate z-scores
proteo_df <- hope_proteo %>%
  dplyr::select(-NP_id) %>%
  dplyr::filter(GeneSymbol %in% goi) %>%
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
  dplyr::filter(GeneSymbol %in% c(clk_list, nf1_list)) %>%
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
clk_nf1_phospho_df <- phospho_df %>%
  ungroup() %>%
  dplyr::filter(!grepl("NA", gene_site)) %>%
  dplyr::select(-abundance,
                -GeneSymbol,
                -Site, -Kids_First_Biospecimen_ID_phospho)

clk_nf1_phospho_df <- clk_nf1_phospho_df %>%
  spread(gene_site, phospho_zscore)

# repeat for proteomics 
clk_nf1_proteo_df <- proteo_df %>%
  ungroup() %>%
  dplyr::mutate(GeneSymbol = glue::glue("{GeneSymbol}")) %>%
  dplyr::select(GeneSymbol, match_id, proteo_zscore)

clk_nf1_proteo_df <- clk_nf1_proteo_df %>%
  spread(GeneSymbol, proteo_zscore)

# merge CLK1 PSI and CLK and NF1 RNA, protein, and phosphoprotein expression
mol_df <- data_df %>%
  left_join(clk_nf1_proteo_df) %>%
  left_join(clk_nf1_phospho_df) %>%
  select(match_id, everything(.))

# save df
write_tsv(mol_df, 
          file.path(results_dir, "hgg-dmg-clk-nf1-expression-phosphorylation.tsv"))

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
# put CLK-201 expr first
#features <- features[c(length(features), 1:(length(features)-1))]

# create empty df to store correlation coefficients between PSI and feature values
incl_cor_mat <- data.frame(row.names = features,
                           DMG = rep(0, length(features)),
                           HGG = rep(0, length(features)),
                           DMG_p = rep(0, length(features)),
                           HGG_p = rep(0, length(features)))  %>%
  dplyr::mutate(Feature = case_when(rownames(.) %in% names(clk_nf1_phospho_df) ~ "Phospho-Protein z-score",
                                    rownames(.) %in% names(clk_nf1_proteo_df) ~ "Total Protein z-score",
                                    rownames(.) %in% c(clk1_splice_list, nf1_splice_list) ~ "PSI",
                                    TRUE ~ "RNA log2Exp"
  ))


# create same empty df for CLK1 expression cor values
expr_cor_mat <- incl_cor_mat
ex4_expr_cor_mat <- incl_cor_mat

# loop through subtypes and features to calculate pearson correlation coefficients
id_list <- list("DMG" = dmg_ids, "HGG" = hgg_ids)

for (subtype in c("DMG", "HGG")) {
  
  ids <- id_list[[subtype]]
  
  for (i in 1:nrow(incl_cor_mat)) {
    
    feature <- rownames(incl_cor_mat)[i]
    
    incl_cor_mat[feature, subtype] <- cor.test(mol_df$`CLK1-Exon4_PSI`[mol_df$match_id %in% ids],
                                               unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                               method = ifelse(incl_cor_mat$Feature[i] == "RNA log2Exp", 
                                                               "pearson", "spearman"))$estimate
    
    incl_cor_mat[feature, glue::glue("{subtype}_p")] <- cor.test(mol_df$`CLK1-Exon4_PSI`[mol_df$match_id %in% ids],
                                                                 unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                                                 method = ifelse(incl_cor_mat$Feature[i] == "RNA log2Exp", 
                                                                                 "pearson", "spearman"))$p.value
    
    expr_cor_mat[feature, subtype] <- cor.test(mol_df$`Total CLK1`[mol_df$match_id %in% ids],
                                               unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                               method = ifelse(incl_cor_mat$Feature[i] == "RNA log2Exp", 
                                                               "pearson", "spearman"))$estimate
    
    expr_cor_mat[feature, glue::glue("{subtype}_p")] <- cor.test(mol_df$`Total CLK1`[mol_df$match_id %in% ids],
                                                                 unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                                                 method = ifelse(incl_cor_mat$Feature[i] == "RNA log2Exp", 
                                                                                 "pearson", "spearman"))$p.value
    
    ex4_expr_cor_mat[feature, subtype] <- cor.test(unlist(mol_df$`CLK1-201`[mol_df$match_id %in% ids]),
                                                   unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                                   method = ifelse(incl_cor_mat$Feature[i] == "RNA log2Exp", 
                                                                   "pearson", "spearman"))$estimate
    
    ex4_expr_cor_mat[feature, glue::glue("{subtype}_p")] <- cor.test(mol_df$`CLK1-201`[mol_df$match_id %in% ids],
                                                                     unlist(mol_df[mol_df$match_id %in% ids, feature]),
                                                                     method = ifelse(incl_cor_mat$Feature[i] == "RNA log2Exp", 
                                                                                     "pearson", "spearman"))$p.value
    
  }
  
}

# define row_annot df to store feature type
row_annot <- expr_cor_mat %>%
  dplyr::select(Feature) %>%
  dplyr::mutate(Feature = factor(Feature,
                                 levels = c("RNA log2Exp", "PSI", "Total Protein z-score", "Phospho-Protein z-score")))


# add rownames
rownames(row_annot) <- rownames(expr_cor_mat)

# rename cor mat colnames for plotting
colnames(incl_cor_mat)[1:2] <- c("DIPG or DMG", "Other HGG")
colnames(expr_cor_mat)[1:2]  <- c("DIPG or DMG", "Other HGG")
colnames(ex4_expr_cor_mat)[1:2]  <- c("DIPG or DMG", "Other HGG")

# create Heatmap row annotation for feature type
anno_col <- list(Feature = c("RNA log2Exp" = "#DC3220", "PSI" = "purple", "Total Protein z-score" = "#005AB5", "Phospho-Protein z-score" = "#40B0A6"))

row_anno = rowAnnotation(df = row_annot,
                         col = anno_col, show_legend = TRUE)

# define matrix indicating if correlation p < 0.05
incl_pval_mat <- ifelse(incl_cor_mat[,3:4] < 0.05, "*", "") 
expr_pval_mat <- ifelse(expr_cor_mat[,3:4] < 0.05, "*", "")
ex4_expr_pval_mat <- ifelse(ex4_expr_cor_mat[,3:4] < 0.05, "*", "")

# create CLK1 ex4 psi correlation heatmap
incl_ht <- Heatmap(as.matrix(incl_cor_mat[,1:2]),
                   width = unit(13, "mm"),
                   name = "Correlation Coefficient",
                   col = colorRamp2(c(-1, 0, 1), c("#E66100", "white", "#5D3A9B")),
                   cluster_rows = FALSE,
                   row_split = row_annot$Feature,
                   column_gap = 0.5,
                   show_row_names = TRUE,
                   show_heatmap_legend=TRUE,
                   cluster_columns = FALSE,
                   # right_annotation = row_anno,
                   row_title = NULL,
                   column_title = "CLK1 Exon4 PSI",
                   column_title_side = "top",
                   column_title_rot = 30,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%s", incl_pval_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                   })

# create CLK1 expression correlation heatmap
expr_ht <- Heatmap(as.matrix(expr_cor_mat[,1:2]),
                   width = unit(13, "mm"),
                   name = "Correlation Coefficient",
                   col = colorRamp2(c(-1, 0, 1), c("#E66100", "white", "#5D3A9B")),
                   cluster_rows = FALSE,
                   row_split = row_annot$Feature,
                   column_gap = 0.5,
                   show_row_names = TRUE,
                   show_heatmap_legend=TRUE,
                   cluster_columns = FALSE,
                   right_annotation = row_anno,
                   row_title = NULL,
                   column_title = "CLK1 Exp",
                   column_title_side = "top",
                   column_title_rot = 30,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%s", expr_pval_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                   })

# create CLK1 expression correlation heatmap
ex4_expr_ht <- Heatmap(as.matrix(ex4_expr_cor_mat[,1:2]),
                       width = unit(13, "mm"),
                       name = "Correlation Coefficient",
                       col = colorRamp2(c(-1, 0, 1), c("#E66100", "white", "#5D3A9B")),
                       cluster_rows = FALSE,
                       row_split = row_annot$Feature,
                       column_gap = 0.5,
                       show_row_names = TRUE,
                       #  show_column_names = TRUE,
                       show_heatmap_legend=TRUE, 
                       cluster_columns = FALSE,
                       # right_annotation = row_anno,
                       row_title = NULL,
                       column_title = "CLK1-201 Exp",
                       column_title_side = "top",
                       column_title_rot = 30,
                       na_col = "gray",
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%s", ex4_expr_pval_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                       })

# save merged plot
pdf(file.path(plots_dir, "CLK1-NF1-psi-expr-correlation-heatmap.pdf"), width = 7, height = 13)
print(incl_ht + ex4_expr_ht + expr_ht)
dev.off()

# create dfs for generating expr-prot and expr-prot scatter plots
proteo_scatter_df <- clk_nf1_proteo_df %>%
  left_join(mol_df %>% dplyr::select(match_id, `CLK1-Exon4_PSI`, `Total CLK1`, `CLK1-201`, `Total NF1`, `NF1-202_PC`, `NF1-Exon23a_PSI`, `NF1-215_RI`, `NF1-208_NMD`))


# generate CLK1 expr-CLK/NF1 prot and expr-protein scatterplots
for (each in all_lists) {
  
  for (subtype in c("DMG", "HGG")) {

    ids <- id_list[[subtype]]
    
    p_prot <- proteo_scatter_df %>%
      dplyr::filter(match_id %in% ids) %>%
      ggplot(aes(x = `NF1-215_RI`, y = each)) +
      geom_point(colour = "black") +
      stat_smooth(method = "lm", 
                  formula = y ~ x, 
                  geom = "smooth", 
                  colour = "red",
                  fill = "pink",
                  linetype="dashed") +
      labs(x = glue::glue(each),
           y = "NF1 Total protein abundance z-score") + 
      stat_cor(method = "spearman", cor.coef.name = "rho",  na.rm = F,
              size = 3) +
      theme_Publication()
  # label.x = 0, label.y = 3, 
    
  pdf(file.path(paste(plots_dir, "/", "NF1-215-cor-", each, "-", subtype, ".pdf", sep = "")), width = 8, height = 10)
  print(p_prot)
  dev.off()
  
  p_201_prot <- proteo_scatter_df %>%
    dplyr::filter(match_id %in% ids) %>%
    ggplot(aes(x = `NF1-208_NMD`, y = each)) +
    geom_point(colour = "black") +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth", 
                colour = "red",
                fill = "pink",
                linetype="dashed") +
    labs(x = glue::glue(each),
         y = "NF1 Total protein abundance z-score") + 
    stat_cor(method = "spearman", cor.coef.name = "rho", na.rm = F,
            size = 3) +
    theme_Publication()
  # label.x = 0, label.y = 4, 
  pdf(file.path(paste(plots_dir, "/", "NF1-208-cor-", each, "-", subtype, ".pdf", sep = "")), width = 8, height = 10)
  print(p_201_prot)
  dev.off()
  }
}

