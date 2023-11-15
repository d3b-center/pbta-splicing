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

analysis_dir <- file.path(root_dir, "analyses", "KNS42-cell-line")
input_dir   <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")


if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}


## output for plots
file_depmap_all_score_plot <- file.path(plots_dir, "depmap_score_all_cell_lines.pdf")
file_expr_vs_score_plot <- file.path(plots_dir, "depmap_score_CLK1_vs_score_KNS42.pdf")
file_depmap_cns_score_plot <- file.path(plots_dir, "depmap_score_cns_cell_lines.pdf")

## call plot publication theme script 
source(file.path(root_dir, "figures", "theme_for_plots.R"))

## load dataset
depmap_file = file.path(data_dir, "CLK1-CRISPR-DepMap-score.csv")
depmap_data <- vroom(depmap_file, show_col_types = FALSE) %>% 
  dplyr::rename("ModelID" = `Depmap ID`,
                "Cell line"=`Cell Line Name`,
                Histology = `Lineage Subtype`) 
depmap_glioma <- depmap_data %>% 
  filter(Histology == "Diffuse Glioma")

## general Depmap CLK1 Pertubation
depmap_data_KNS42 <- depmap_glioma %>% 
  filter(`Cell line` =="KNS42")

gene_score_plot <- ggplot(data=depmap_glioma, 
                          aes(reorder(`Cell line`,`CRISPR (DepMap Public 23Q2+Score, Chronos)`),`CRISPR (DepMap Public 23Q2+Score, Chronos)`,  
                              group=1)) +
  #geom_line() + 
  geom_point(size=3, colour="gray89") + 
  geom_point(size=3, colour = "gray50", pch = 21) + 
  geom_point(data=depmap_data_KNS42, colour="red", size = 3) +
  geom_point(data=depmap_data_KNS42, colour="black", size = 3, pch = 21) +
  xlab("Cell Line") + 
  ylab("CRISPR Dependency Score") + 
  ggtitle("CLK1 Perturbation") +
  theme_Publication() + 
  theme(axis.text.x=element_text(angle = 75, hjust = 1, size = 8)) 

# Save plot as pdf
pdf(file_expr_vs_score_plot, height = 4, width = 7.5)
gene_score_plot
dev.off()

## compare distributions based on CLK1 expression 
# remove histologies with < 6 samples, need at least 3 per group later in boxplots
rm_low_n <- as.data.frame(table(depmap_data$Lineage)) %>%
  filter(Freq >= 6)

depmap_rm_low_n <- depmap_data %>%
  filter(Lineage %in% rm_low_n$Var1)

omics_mappings_file <- file.path(data_dir,"OmicsDefaultModelProfiles.csv")
tpm_file <- file.path(data_dir,"OmicsExpressionTranscriptsTPMLogp1Profile.csv")

omics_id_mapping_df <- vroom(omics_mappings_file, show_col_types = FALSE) %>% 
inner_join(depmap_rm_low_n, by="ModelID")

depMap_transcr_expr <- vroom(tpm_file, show_col_types = FALSE) %>% 
  dplyr::rename("ProfileID"=`...1`) %>%
  inner_join(omics_id_mapping_df, by="ProfileID")

depMap_transcr_CLK1_pct_expr <- depMap_transcr_expr %>% 
  dplyr::rename("CRISPR_score" = `CRISPR (DepMap Public 23Q2+Score, Chronos)`) %>%
  dplyr::select(ModelID, 
                Lineage,
                `CLK1 (ENST00000321356)`,
                `CLK1 (ENST00000409769)`,
                `CLK1 (ENST00000409769)`,
                `CLK1 (ENST00000434813)`, 
                CRISPR_score) %>% 
  mutate(CLK1_ENST00000321356_Expr =`CLK1 (ENST00000321356)`) 

# generate quantiles by lineage
quartiles_by_histology <- depMap_transcr_CLK1_pct_expr %>%
  group_by(Lineage) %>%
  summarise(
    lower_quantile = quantile(CLK1_ENST00000321356_Expr, probs = 0.25, na.rm = TRUE),
    upper_quantile = quantile(CLK1_ENST00000321356_Expr, probs = 0.75, na.rm = TRUE),
    IQR = IQR(CLK1_ENST00000321356_Expr, na.rm = TRUE)
  ) 

# create high/low expression groups and remove those not in highest/lowest quantiles
clk1_by_quartiles <- depMap_transcr_CLK1_pct_expr %>%
  left_join(quartiles_by_histology) %>%
  mutate(Expression_level = case_when(CLK1_ENST00000321356_Expr > upper_quantile ~ "high",
                                      CLK1_ENST00000321356_Expr < lower_quantile ~ "low",
                                      TRUE ~ NA_character_)) %>%
  # remove NAs
  filter(!is.na(Expression_level))

# remove groups with < 5 
exp_low_n_rm <- table(clk1_by_quartiles$Lineage, clk1_by_quartiles$Expression_level) %>%
  as.data.frame() %>%
  filter(Freq < 5)

clk1_rm_low_n <- clk1_by_quartiles %>%
  filter(!Lineage %in% exp_low_n_rm$Var1)

set.seed(42)
boxplot_expr_vs_score <- ggplot(clk1_rm_low_n, 
                                aes(x = Expression_level, 
                                    y = CRISPR_score, 
                                    color = Expression_level)) + 
  geom_boxplot(outlier.size = 0,
               size = 0.5,
               alpha = 0,
               color = "black",
               coef = 0) +
  ggforce::geom_sina(aes(color = Expression_level) , size = 3) +
  facet_wrap(~Lineage, nrow = 3, labeller = label_wrap_gen(15)) +
  scale_color_manual(values = c(high = "#FFC20A", low = "#0C7BDC")) +
  stat_compare_means(position = "identity", label.y = 0.2, label.sep = " ") + 
  ylab("CRISPR Dependency Score") + 
  xlab(expression(bold(bolditalic("CLK1")~"ENST00000321356 Expression"))) +
  theme_Publication() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10))


# save plot pdf version
pdf(file_depmap_all_score_plot, height = 10, width = 12)
print(boxplot_expr_vs_score)
dev.off()

# save brain plot separately for main figure
clk1_by_quartiles_brain <- clk1_by_quartiles %>%
  filter(Lineage == "CNS/Brain")

clk1_box <- ggplot(clk1_by_quartiles_brain, 
                                aes(x = Expression_level, 
                                    y = CRISPR_score, 
                                    color = Expression_level)) + 
  geom_boxplot(outlier.size = 0,
               size = 0.5,
               alpha = 0,
               color = "black",
               coef = 0) +
  ggforce::geom_sina(aes(color = Expression_level) , size = 3) +
  scale_color_manual(values = c(high = "#FFC20A", low = "#0C7BDC")) +
  stat_compare_means(position = "identity", label.y = 0.2, label.sep = " ") + 
  ylab("CRISPR Dependency Score") + 
  xlab(expression(bold(atop(bolditalic("CLK1"), 
                            "ENST00000321356 Expression")))) +
  theme_Publication() +
  theme(legend.position = "none")

# save plot pdf version
pdf(file_depmap_cns_score_plot, height = 4, width = 4)
print(clk1_box)
dev.off()

