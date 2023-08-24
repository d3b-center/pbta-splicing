# 01-plotting-depMap-data.R
#
#

suppressPackageStartupMessages({
  library("tidyverse")
  library("dplyr")
  library("vroom")
  library("ggpubr")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data/")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_splicing_impact")
input_dir   <- file.path(analysis_dir, "input")

## load dataset
depmap_file = "CLK1_CRISPR_depmap_score.csv"
depmap_data <- vroom(paste0(file.path(input_dir),"/", depmap_file), show_col_types = FALSE) %>% 
  dplyr::rename("ModelID"=`Depmap ID`) %>% 
  filter( Lineage == 'CNS/Brain' ) 

## general Depmap CLK1 Pertubation
depmap_data_KNS42 <- depmap_data %>% dplyr::filter(`Cell Line Name` =="KNS42")

ggplot(data=depmap_data, aes(reorder(`Cell Line Name`,`CRISPR (DepMap Public 23Q2+Score, Chronos)`),`CRISPR (DepMap Public 23Q2+Score, Chronos)`,  group=1)) +
  #geom_line() + 
  geom_point(size=4, colour="grey") + 
  theme_Publication() + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  geom_point(data=depmap_data_KNS42, colour="red", size = 4) +
  xlab("Cell Line") + ylab("Dependency Score") + ggtitle("DepMap CLK1 Pertubation") 


## compare distrubutions based on CLK1 expression 
omics_id_mapping_df <- vroom("~/Downloads/OmicsDefaultModelProfiles.csv") %>% inner_join(depmap_data,by="ModelID")

depMap_transcr_expr <- vroom("~/Downloads/OmicsExpressionTranscriptsTPMLogp1Profile.csv") %>% 
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

ggplot(depmap_high_low,aes(Expression_level,CRISPR_score)) + 
  geom_boxplot(aes(fill=Expression_level)) + 
  stat_compare_means() + 
  geom_jitter() 
