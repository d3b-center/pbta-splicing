################################################################################
# 01-plotting-depMap-data.R
#
# Plot and visualize dependency scores across brain tumors and CLK1 expr
#
# author: Ammar S Naqvi
# usage: Rscript 01-plotting-depMap-data.R
################################################################################# 

suppressPackageStartupMessages({
  library("tidyverse")
  library("dplyr")
  library("vroom")
  library("ggpubr")
  library("ggthemes")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

analysis_dir <- file.path(root_dir, "analyses", "KNS42_cell-line")
input_dir   <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots", "KNS42_cell-line")


if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}


## output for plots
file_depmap_score_plot <- file.path(plots_dir, "depmap_score_cell-lines.tiff")
file_expr_vs_score_plot <- file.path(plots_dir, "depmap_score_CLK1_vs_score_KNS42.tiff")

## call plot publication theme script 
source(file.path(root_dir, "figures", "theme_for_plots.R"))

## load dataset
depmap_file = file.path(data_dir, "CLK1-CRISPR-DepMap-score.csv")
depmap_data <- vroom(depmap_file, show_col_types = FALSE) %>% 
  dplyr::rename("ModelID"=`Depmap ID`) %>% 
  filter( Lineage == 'CNS/Brain' ) 

## general Depmap CLK1 Pertubation
depmap_data_KNS42 <- depmap_data %>% dplyr::filter(`Cell Line Name` =="KNS42")

gene_score_plot <- ggplot(data=depmap_data, aes(reorder(`Cell Line Name`,`CRISPR (DepMap Public 23Q2+Score, Chronos)`),`CRISPR (DepMap Public 23Q2+Score, Chronos)`,  group=1)) +
  #geom_line() + 
  geom_point(size=3, colour="grey") + 
  theme_Publication() + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  geom_point(data=depmap_data_KNS42, colour="red", size = 3) +
  xlab("Cell Line") + ylab("Dependency Score") + ggtitle("CLK1 Perturbation") 


# Save plot as tiff
tiff(file_depmap_score_plot, 
     res = 300, width = 1800, height = 1000)
gene_score_plot
dev.off()

## compare distributions based on CLK1 expression 
omics_mappings_file <- file.path(data_dir,"OmicsDefaultModelProfiles.csv")
tpm_file <- file.path(data_dir,"OmicsExpressionTranscriptsTPMLogp1Profile.csv")

omics_id_mapping_df <- vroom(omics_mappings_file,show_col_types = FALSE) %>% 
inner_join(depmap_data,by="ModelID")
depMap_transcr_expr <- vroom(tpm_file,show_col_types = FALSE) %>% 
  dplyr::rename("ProfileID"=`...1`) %>%
  inner_join(omics_id_mapping_df,by="ProfileID")

depMap_transcr_CLK1_pct_expr <- depMap_transcr_expr %>% rename("CRISPR_score"=`CRISPR (DepMap Public 23Q2+Score, Chronos)`) %>%
  dplyr::select(ModelID, `CLK1 (ENST00000321356)`,`CLK1 (ENST00000409769)`,`CLK1 (ENST00000409769)` ,`CLK1 (ENST00000434813)`, CRISPR_score) %>% 
  mutate(PCT_Expr =`CLK1 (ENST00000321356)`) 

quartiles_expr <- quantile(depMap_transcr_CLK1_pct_expr$PCT_Expr, probs=c(.25, .75), na.rm = TRUE)
IQR_expr <- IQR(depMap_transcr_CLK1_pct_expr$PCT_Expr, na.rm = TRUE)

lower_expr <- quartiles_expr[1] # low expression cut-off
upper_expr <- quartiles_expr[2] # high expression cut-off

depmap_high_expr_df <- dplyr::filter(depMap_transcr_CLK1_pct_expr, PCT_Expr > upper_expr) %>% dplyr::mutate(Expression_level="high")
depmap_low_expr_df  <- dplyr::filter(depMap_transcr_CLK1_pct_expr, PCT_Expr < lower_expr) %>% dplyr::mutate(Expression_level="low")

depmap_high_low <- rbind(depmap_low_expr_df,depmap_high_expr_df)

set.seed(42)
boxplot_expr_vs_score <- ggplot(depmap_high_low, 
                                aes(x = Expression_level, 
                                    y = CRISPR_score, 
                                    color = Expression_level)) + 
  geom_boxplot(outlier.size = 0,
               size = 0.5,
               alpha = 0,
               color = "black",
               coef = 0) +
  ggforce::geom_sina(aes(color = Expression_level) , size = 3) +
  scale_color_manual(name = "Expression level", values = c(high = "#FFC20A", low = "#0C7BDC")) +
  stat_compare_means(position = "identity", label.x = 1.5) + 
  
  ylab("Dependency Score") + 
  xlab(expression(bold(bolditalic("CLK1")~"Exon 4 Transcript Expression"))) +
  theme_Publication()

boxplot_expr_vs_score

# save plot tiff version
tiff(file_expr_vs_score_plot, height =1500, width = 1700, res = 300)
print(boxplot_expr_vs_score)
dev.off()
