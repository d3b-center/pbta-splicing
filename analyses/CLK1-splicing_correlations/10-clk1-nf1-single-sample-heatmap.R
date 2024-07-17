################################################################################
# 10-clk1-nf1-single-sample-heatmap.R
# Generate heatmap of CLK1 features (PSI, Expr) and CLK1 and NF1 RNA, protein, and phosphoprotein expression
# written by Jo Lynne Rokita
# Usage: Rscript 10-clk1-nf1-single-sample-heatmap.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(data.table)
  library(Hmisc)
  library(corrplot)
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
gtf_file <- file.path(data_dir, "gencode.v39.primary_assembly.annotation.gtf.gz")
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
cohort_file <- file.path(root_dir, "analyses", "cohort_summary",
                         "results", "histologies-plot-group.tsv")
rsem_tpm <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")
rsem_transc_counts <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")
rmats_clk1_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")
rmats_nf1_file <- file.path(data_dir, "nf1-splice-events-rmats.tsv")

hope_proteo_file <- file.path(data_dir, "hope-protein-imputed-prot-expression-abundance.tsv.gz")
hope_phospho_file <- file.path(data_dir, "hope-protein-imputed-phospho-expression-abundance.tsv.gz")

# output file
heatmap_output_file <- file.path(plots_dir, "clk1-nf1-single-sample-exp-protein-heatmap-dmg.pdf")



# define lists of nf1 and clk genes
goi <- c("CLK1", "NF1")
clk1_trans_list <- c("CLK1-201", "Total CLK1")
clk1_splice_list <- "CLK1-201 (Exon4) PSI"
nf1_trans_list <- c("Total NF1", "NF1-201", "NF1-202", "NF1-215")
nf1_splice_list <- c("NF1-202 (Exon23a) PSI", "NF1-215 PSI")
phos_list <- "NF1 pS864"

trans_lists <- c(clk1_trans_list, nf1_trans_list)


## read in histology, cohort, and independent specimens file
cohort_df <- read_tsv(cohort_file)
indep_rna_df <- read_tsv(indep_rna_file)

# extract HGG samples, retain stranded libraries, and filter for independent specimens
rna_ids <- cohort_df %>% 
  filter(plot_group %in% c("DIPG or DMG"),
         RNA_library == "stranded",
         cohort == "PBTA",
         Kids_First_Biospecimen_ID %in% indep_rna_df$Kids_First_Biospecimen_ID) %>%
select(Kids_First_Biospecimen_ID, match_id)

# match ids for rna, now get those with proteomics
match_ids <-  rna_ids$match_id
prot_match_ids <- cohort_df %>%
  filter(match_id %in% match_ids,
         grepl("Proteo", experimental_strategy),
         sub_cohort == "HOPE") %>%
  select(Kids_First_Biospecimen_ID, match_id)

match_ids_to_use <- cohort_df %>%
  filter(match_id %in% prot_match_ids$match_id) %>%
  select(Kids_First_Biospecimen_ID, match_id, plot_group, plot_group_hex)

# splice tables

# read in rmats file once with genes of interest
rmats_nf1  <- fread(rmats_nf1_file)

rmats_goi  <- fread(rmats_clk1_file) %>%
  bind_rows(rmats_nf1) %>%
  # annotated the NF1 and CLK1 transcripts of interest
  mutate(ENST_id = case_when(geneSymbol == "NF1" & splicing_case == "SE" & exonStart_0base == "31252937" & exonEnd == "31253000" & downstreamES == "31258343" & downstreamEE == "31258502" ~ "ENST00000358273",
                             geneSymbol == "NF1" & splicing_case == "RI" & riExonStart_0base == "31229024" & riExonEnd == "31229974" ~ "ENST00000493220",
                             geneSymbol == "CLK1" & splicing_case == "SE" & exonStart_0base == "200860124" & exonEnd == "200860215" ~ "ENST00000321356",
                             TRUE ~ NA_character_),
         event_type = case_when(ENST_id == "ENST00000358273" ~ "NF1-202 Exon23a PSI", # NF1-202, also canonical
                                ENST_id == "ENST00000493220" ~ "NF1-215 PSI", # NF1-215
                                ENST_id == "ENST00000321356" ~ "CLK1-201 Exon4 PSI", # CLK1-201
                                TRUE ~ NA_character_)) %>%
  filter(!is.na(ENST_id))

rmats_df <- rmats_goi %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample_id) %>%
  filter(Kids_First_Biospecimen_ID %in% match_ids_to_use$Kids_First_Biospecimen_ID) %>%
  left_join(match_ids_to_use[,c("Kids_First_Biospecimen_ID", "match_id")], by = 'Kids_First_Biospecimen_ID') %>%
  select(event_type, IncLevel1, match_id) %>%
  #zscore
  group_by(event_type) %>%
  dplyr::mutate(psi_zscore = as.numeric(scale(IncLevel1))) %>%
  ungroup() %>%
  select(-IncLevel1) %>%
  spread(key = event_type, value = psi_zscore)


## gene tpm table 
gene_tpm <- readRDS(rsem_tpm) %>%
  dplyr::select(any_of(match_ids_to_use$Kids_First_Biospecimen_ID)) %>%
  rownames_to_column(var = "gene_symbol") %>%
  filter(gene_symbol %in% c("NF1", "CLK1")) %>%
  mutate(gene_symbol = paste0("Total ", gene_symbol, " mRNA")) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.numeric, ~ as.numeric(scale(.))) %>%
  rownames_to_column(var = "Kids_First_Biospecimen_ID") %>%
  left_join(match_ids_to_use[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  select(-Kids_First_Biospecimen_ID) 


# transcript table
transcript_df <- readRDS(rsem_transc_counts) %>%
  filter(gene_symbol %in% trans_lists) %>%
  select(gene_symbol, any_of(match_ids_to_use$Kids_First_Biospecimen_ID)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.numeric, ~ as.numeric(scale(.))) %>%
  rownames_to_column(var = "Kids_First_Biospecimen_ID") %>%
  left_join(match_ids_to_use[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  select(-Kids_First_Biospecimen_ID)


# prot data - zscore for only DMGs
prot_df <- read_tsv(hope_proteo_file) %>%
  dplyr::select(-NP_id) %>%
  dplyr::filter(GeneSymbol %in% goi) %>%
  gather(key = "Kids_First_Biospecimen_ID",
         value = "abundance",
         -GeneSymbol) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% all_ids$Kids_First_Biospecimen_ID)  %>%
  distinct(Kids_First_Biospecimen_ID, GeneSymbol, .keep_all = TRUE) %>%
  group_by(GeneSymbol) %>%
  dplyr::mutate(prot_zscore = as.numeric(scale(abundance))) %>%
  ungroup() %>%
  left_join(match_ids_to_use[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  dplyr::rename(`Total NF1 protein` = prot_zscore) %>%
  select(match_id, `Total NF1 protein`)


# phos data - zscore for only DMGs
phos_df <- read_tsv(hope_phospho_file) %>%
  dplyr::select(-NP_id, -Peptide_res_num,
                -Peptide_sequence) %>%
  dplyr::filter(GeneSymbol %in% goi) %>%
  gather(key = "Kids_First_Biospecimen_ID",
         value = "abundance",
         -GeneSymbol, -Site) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% all_ids$Kids_First_Biospecimen_ID)  %>%
  distinct(Kids_First_Biospecimen_ID, GeneSymbol, Site, .keep_all = TRUE) %>%
  group_by(GeneSymbol, Site) %>%
  dplyr::mutate(phospho_zscore = as.numeric(scale(abundance))) %>%
  ungroup() %>%
  left_join(match_ids_to_use[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  dplyr::mutate(gene_site = glue::glue("{GeneSymbol} p{Site}"))  %>%
  select(-Kids_First_Biospecimen_ID, -GeneSymbol, -Site, -abundance) %>%
  spread(key = gene_site, value = phospho_zscore) %>%
  select(match_id, any_of(phos_list))


full_data <- prot_df %>%
  left_join(phos_df) %>%
  left_join(rmats_df) %>%
  left_join(gene_tpm) %>%
  left_join(transcript_df) %>%
  # order by CLK1 exon 4 psi
  arrange(`CLK1-201 Exon4 PSI`) %>%
  column_to_rownames("match_id") %>%
  t() %>%
  as.data.frame()

mat <- full_data %>%
  as.matrix()

storage.mode(mat) <- "numeric"
class(mat)

anno_names <- names(full_data) %>%
  as.data.frame()
names(anno_names) <- "match_id"

# create single sample heatmap
col_annot <- anno_names %>%
  left_join(unique(match_ids_to_use[,c("match_id", "plot_group", "plot_group_hex")])) %>%
  column_to_rownames("match_id") %>%
  select(plot_group) %>%
  dplyr::rename(Histology = plot_group) %>%
  as.data.frame()

# create anno colors
col_anno_col <- list(Histology = c("DIPG or DMG" = "#ff40d9"))

ha = HeatmapAnnotation(
  Histology = col_annot$Histology, col = col_anno_col,
  annotation_name_side = "right")


# Heatmap annotation
row_annot <- as.data.frame(rownames(full_data)) %>%
  dplyr::rename(ID = `rownames(full_data)`) %>%
  mutate(Abundance = case_when(ID == "Total NF1 protein" ~ "Whole Cell Protein",
                           ID %in% c("NF1 pS864", "NF1 pS2796") ~ "Phosphoprotein",
                           ID %in% c("NF1-201", "NF1-202", "CLK1-201", "NF1-215") ~ "Transcript",
                           ID %in% c("Total NF1 mRNA", "Total CLK1 mRNA") ~ "Gene",
                           ID %in% c("CLK1-201 Exon4 PSI", "NF1-202 Exon23a PSI", "NF1-215 PSI") ~ "mRNA splicing"))

row_anno_col <- list(Abundance = c("Transcript" = "#DC3220",
                                   "Gene" = "purple",
                                   "mRNA splicing" = "orange",
                                   "Phosphoprotein" = "#005AB5", 
                                   "Whole Cell Protein" = "#40B0A6"))

row_anno = rowAnnotation(Abundance = row_annot$Abundance,
                         col = row_anno_col, 
                         show_legend = TRUE,
                         show_annotation_name = FALSE)

# Make heatmap
heat_plot <- Heatmap(mat,
                     name = "Z-score",
                     col = colorRamp2(c(-2, 0, 2), c("darkblue", "white", "red")),
                     #col = colorRamp2(c(-2, 0, 2), c("#E66100", "white", "#5D3A9B")),
                     cluster_rows = TRUE,
                     column_gap = 0.5,
                     show_row_names = TRUE,
                     show_column_names = FALSE,
                     show_heatmap_legend=TRUE,
                     cluster_columns = FALSE,
                    # clustering_distance_rows = "spearman",
                     clustering_distance_rows = "euclidean",
                     right_annotation = row_anno,
                     top_annotation = ha,
                     na_col = "lightgrey",
                     rect_gp = gpar(col = "white"),
                     row_title = NULL,
                     column_title = NULL,
                     column_title_side = "top")

pdf(heatmap_output_file, width = 8, height = 3)
print(heat_plot)
dev.off()


# corplot
rev_mat <- full_data %>%
  t()
gp_mat <- cor(rev_mat, use = "pairwise.complete.obs", method = "spearman")

# Perform significance testing for the correlations
p_gp <- cor.mtest(gp_mat, method = "spearman")$p

# Create correlogram RdYlBu
col <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988", "#BB4444"))

# Create the correlation plot with blank for non-significant correlations
pdf(file.path(plots_dir, "single-sample-cor-dmg.pdf"), width = 8, height = 8)
par(mar = c(5, 4, 4, 6) + 0.1)  # Adjust the right margin to bring the text closer
corrplot(gp_mat, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p_gp, sig.level = 0.05, insig = "blank", number.cex = 0.8, number.font = 1,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         outline = TRUE,
         addgrid.col = "black"
)
# Create a new plot overlay to add the y-axis title
mtext("Spearman's Rho", side = 4, line = 2, cex = 1.2, font = 2, las = 0, padj = -3)

dev.off()









