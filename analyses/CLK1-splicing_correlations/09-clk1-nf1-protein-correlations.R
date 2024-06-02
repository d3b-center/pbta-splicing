################################################################################
# 09-clk1-nf1-protein-correlations.R
# Generate heatmap of correlation coefficients between CLK1 features (PSI, Expr) and CLK1 and NF1 RNA, protein, and phosphoprotein expression
# written by Ryan Corbett, Jo Lynne Rokita
# Usage: Rscript 09-clk1-nf1-protein-correlations.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggpubr")
  library("ComplexHeatmap")
  library("circlize")
  library("stringr")
  library("corrplot")
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
nf1_trans_list <- c("Total NF1", "NF1-202_PC", "NF1-215_RI")
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
  dplyr::filter(GeneSymbol %in% goi) %>%
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
  spread(gene_site, phospho_zscore) %>%
  # select only sites of interest
  select(match_id, `NF1-S864`, `NF1-S2796`)

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

# create matrix for cor
full_mat <- mol_df %>%
  column_to_rownames("match_id") %>%
  as.matrix()

# loop through subtypes and features to calculate pearson correlation coefficients
id_list <- list("DMG" = dmg_ids, "HGG" = hgg_ids)

# Initialize lists to store matrices for each subtype
cor_matrices <- list()
pval_matrices <- list()

meas_of_int <- c("CLK1-Exon4_PSI", "CLK1-201", "Total CLK1")

# Loop over each subtype and compute correlation and p-value matrices
for (subtype in c("DMG", "HGG")) {
  
  ids <- id_list[[subtype]]
  mat <- full_mat[ids,]
  
  # Compute the correlation matrix
  cor_mat <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  
  # Perform significance testing for the correlations
  p_mat <- cor.mtest(cor_mat, use = "pairwise.complete.obs", method = "spearman")$p

  # Define matrix indicating if correlation p < 0.05
  pval_mat <- ifelse(p_mat < 0.05, "*", "")
  
  # subset for only CLK1 cols
  cor_mat <- cor_mat[,meas_of_int]
  pval_mat <- pval_mat[,meas_of_int]
  
  # Append suffix to column names to avoid conflicts
  colnames(cor_mat) <- paste(colnames(cor_mat), subtype, sep = "_")
  colnames(pval_mat) <- paste(colnames(pval_mat), subtype, sep = "_")
  
  # Append matrices to the list
  cor_matrices[[subtype]] <- cor_mat
  pval_matrices[[subtype]] <- pval_mat
}

# Combine the matrices for all subtypes
combined_cor_mat <- do.call(cbind, cor_matrices)
combined_pval_mat <- do.call(cbind, pval_matrices)

# Print the combined matrices to check
print(combined_cor_mat)
print(combined_pval_mat)

# add rownames
rownames(combined_cor_mat) <- case_when(rownames(combined_cor_mat) == "CLK1-Exon4_PSI" ~ "CLK1-201 Exon4 PSI",
                                    rownames(combined_cor_mat) == "NF1-202_PC" ~ "NF1-202",
                                    rownames(combined_cor_mat) == "NF1-215_RI" ~ "NF1-215",
                                    rownames(combined_cor_mat) == "NF1-Exon23a_PSI" ~ "NF1-202 Exon23a PSI",
                                    rownames(combined_cor_mat) == "NF1-215_PSI" ~ "NF1-215 PSI",
                                    rownames(combined_cor_mat) == "NF1" ~ "Total NF1 protein",
                                    rownames(combined_cor_mat) == "NF1-S864" ~ "NF1 pS864",
                                    rownames(combined_cor_mat) == "NF1-S2796" ~ "NF1 pS2796",
                                    TRUE ~ rownames(combined_cor_mat))

# add rownames
rownames(combined_pval_mat) <- case_when(rownames(combined_pval_mat) == "CLK1-Exon4_PSI" ~ "CLK1-201 Exon4 PSI",
                                        rownames(combined_pval_mat) == "NF1-202_PC" ~ "NF1-202",
                                        rownames(combined_pval_mat) == "NF1-215_RI" ~ "NF1-215",
                                        rownames(combined_pval_mat) == "NF1-Exon23a_PSI" ~ "NF1-202 Exon23a PSI",
                                        rownames(combined_pval_mat) == "NF1-215_PSI" ~ "NF1-215 PSI",
                                        rownames(combined_pval_mat) == "NF1" ~ "Total NF1 protein",
                                        rownames(combined_pval_mat) == "NF1-S864" ~ "NF1 pS864",
                                        rownames(combined_pval_mat) == "NF1-S2796" ~ "NF1 pS2796",
                                        TRUE ~ rownames(combined_pval_mat))

# reorder rownames
ordered_rows <- c("CLK1-201 Exon4 PSI", "NF1-202 Exon23a PSI", "NF1-215 PSI",
                  "CLK1-201", "Total CLK1", "NF1-202", "NF1-215", "Total NF1",
                  "NF1 pS864", "NF1 pS2796", "Total NF1 protein")
combined_cor_mat <- combined_cor_mat[ordered_rows,]
combined_pval_mat <- combined_pval_mat[ordered_rows,]

# Heatmap annotation
row_annot <- as.data.frame(rownames(combined_cor_mat)) %>%
  dplyr::rename(ID = `rownames(combined_cor_mat)`) %>%
  mutate(Feature = case_when(ID == "Total NF1 protein" ~ "Whole Cell Protein",
                               ID %in% c("NF1 pS864", "NF1 pS2796") ~ "Phosphoprotein",
                               ID %in% c("NF1-202", "CLK1-201", "NF1-215",
                                         "Total NF1", "Total CLK1") ~ "RNA Expression",
                               ID %in% c("CLK1-201 Exon4 PSI", "NF1-202 Exon23a PSI", "NF1-215 PSI") ~ "mRNA splicing"))

row_anno_col <- list(Feature = c("RNA Expression" = "#DC3220",
                                   "mRNA splicing" = "orange",
                                   "Phosphoprotein" = "#005AB5", 
                                   "Whole Cell Protein" = "#40B0A6"))
  

row_anno <- rowAnnotation(df = row_annot["Feature"], col = row_anno_col, show_legend = TRUE,
                          show_annotation_name = FALSE)

# column annotation
# create single sample heatmap
col_annot <- cohort_df %>%
  dplyr::filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) %>%
  select(plot_group, plot_group_hex) %>%
  unique()

col_annot_df <- as.data.frame(colnames(combined_cor_mat)) %>%
  dplyr::rename(ID = `colnames(combined_cor_mat)`) %>%
  dplyr::mutate(Histology = case_when(grepl("DMG", ID) ~ "DIPG or DMG",
                                      grepl("HGG", ID) ~ "Other high-grade glioma")
                )

plot_colors <- setNames(col_annot$plot_group_hex, col_annot$plot_group)
plot_colors <- list(Histology = plot_colors)

ha <- HeatmapAnnotation(
  Histology = col_annot_df$Histology,
  col = plot_colors,
  annotation_name_side = "right")

column_labels_manual = rep(c("CLK1-201 Exon4 PSI", "CLK1 201", "Total CLK1"),2)

# create CLK1 ex4 psi correlation heatmap
heat_plot <- Heatmap(combined_cor_mat,
                   name = "Spearman's Rho",
                   col = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red")),
                   cluster_rows = FALSE,
                   row_split = row_annot$Feature,
                   column_split = col_annot_df$Histology,
                   show_parent_dend_line = FALSE,
                   column_labels = column_labels_manual,
                   column_gap = unit(2, "mm"), # Adjust the gap between columns
                   show_row_names = TRUE,
                   column_names_rot = 70, 
                   column_names_side = "top",
                   show_column_names = TRUE,
                   show_heatmap_legend = TRUE,
                   cluster_columns = FALSE,
                   right_annotation = row_anno,
                   top_annotation = ha,
                   row_title = NULL,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%s", combined_pval_mat[i, j]), x, y, gp = gpar(fontsize = 14))
                   })

# save merged plot
pdf(file.path(plots_dir, "CLK1-NF1-psi-expr-correlation-heatmap.pdf"), width = 8, height = 5)
print(heat_plot)
dev.off()




# create dfs for generating expr-prot and expr-prot scatter plots
proteo_scatter_df <- clk_nf1_proteo_df %>%
  left_join(mol_df %>% dplyr::select(match_id, `CLK1-Exon4_PSI`, `Total CLK1`, `CLK1-201`, 
                                     `Total NF1`, `NF1-202_PC`, `NF1-Exon23a_PSI`,
                                     `NF1-215_PSI`, `NF1-215_RI`)) %>%
  dplyr::rename(`CLK1-201 Exon4 PSI` = `CLK1-Exon4_PSI`,
                `NF1-202 Log2 TPM` = `NF1-202_PC`,
                `NF1-215 PSI` = `NF1-215_PSI`,
                `NF1-202 Exon23a PSI` = `NF1-Exon23a_PSI`,
                `NF1-215 Log2 TPM` = `NF1-215_RI`,
                `CLK1-201 Log2 TPM` = `CLK1-201`,
                `Total CLK1 Log2 TPM` = `Total CLK1`,
                `Total NF1 Log2 TPM` = `Total NF1`,
                `NF1 protein abundance z-score` = NF1)

# create dfs for generating expr-prot and expr-prot scatter plots
phos_scatter_df <- clk_nf1_phospho_df %>%
  # select only significant
  select(match_id, `NF1-S2796`, `NF1-S864`) %>%
  left_join(mol_df %>% dplyr::select(match_id, `CLK1-Exon4_PSI`, `Total CLK1`, `CLK1-201`, 
                                     `Total NF1`, `NF1-202_PC`, `NF1-Exon23a_PSI`,
                                     `NF1-215_PSI`, `NF1-215_RI`, NF1)) %>%
  dplyr::rename(`NF1 pS2796` = `NF1-S2796`,
                `NF1 pS864` = `NF1-S864`,
                `CLK1-201 Exon4 PSI` = `CLK1-Exon4_PSI`,
                `NF1-202 Log2 TPM` = `NF1-202_PC`,
                `NF1-215 PSI` = `NF1-215_PSI`,
                `NF1-202 Exon23a PSI` = `NF1-Exon23a_PSI`,
                `NF1-215 Log2 TPM` = `NF1-215_RI`,
                `CLK1-201 Log2 TPM` = `CLK1-201`,
                `Total CLK1 Log2 TPM` = `Total CLK1`,
                `Total NF1 Log2 TPM` = `Total NF1`,
                `NF1 protein abundance z-score` = NF1)


new_list <- names(proteo_scatter_df[,2:ncol(proteo_scatter_df)])

# generate CLK1 expr-CLK/NF1 prot and expr-protein scatterplots
for (each in new_list) {
  
  for (subtype in c("DMG", "HGG")) {

    ids <- id_list[[subtype]]
    
    p_prot <- proteo_scatter_df %>%
      dplyr::filter(match_id %in% ids,
                    !is.na(each)) %>%
      ggplot(aes(x = .data[[each]], y = `NF1 protein abundance z-score`)) +
      geom_point(colour = "black") +
      stat_smooth(method = "lm", 
                  formula = y ~ x, 
                  geom = "smooth", 
                  colour = "red",
                  fill = "pink",
                  linetype="dashed",
                  na.rm = TRUE) +
      labs(x = each,
           y = "NF1 protein abundance z-score") + 
      stat_cor(method = "spearman", cor.coef.name = "rho",  na.rm = F,
               label.y = 2.1, size = 3) +
      theme_Publication()

  pdf(file.path(paste(plots_dir, "/", "NF1-protein-cor-", each, "-", subtype, ".pdf", sep = "")), width = 4, height = 4)
  print(p_prot)
  dev.off()

  }
}


# for phos:
phos_list <- names(proteo_scatter_df[,2:ncol(proteo_scatter_df)])


# generate CLK1 expr-CLK/NF1 prot and expr-protein scatterplots
for (phos_sites in c("NF1 pS2796", "NF1 pS864")) {
    
  for (each in phos_list) {
    
    for (subtype in c("DMG", "HGG")) {
      
      ids <- id_list[[subtype]]
      
      p_prot <- phos_scatter_df %>%
        dplyr::filter(match_id %in% ids,
                      !is.na(each)) %>%
        ggplot(aes(x = .data[[each]], y = .data[[phos_sites]])) +
        geom_point(colour = "black") +
        stat_smooth(method = "lm", 
                    formula = y ~ x, 
                    geom = "smooth", 
                    colour = "red",
                    fill = "pink",
                    linetype="dashed",
                    na.rm = TRUE) +
        labs(x = each,
             y = paste0(phos_sites, " abundance z-score")) + 
        stat_cor(method = "spearman", cor.coef.name = "rho",  na.rm = F,
                 label.y = 2.1, size = 3) +
        theme_Publication()

      pdf(file.path(paste(plots_dir, "/", phos_sites, "-phos-cor-", each, "-", subtype, ".pdf", sep = "")), width = 4, height = 4)
      print(p_prot)
      dev.off()
      
    }
  }
}

