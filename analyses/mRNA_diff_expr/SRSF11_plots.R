################################################################################
# volcano_plot_mRNA.R
# written by Ammar Naqvi
#
# usage: Rscript SRSF11_plot.R
################################################################################


suppressPackageStartupMessages({
  library("sva")
  library("EnhancedVolcano")
  library("DESeq2")
})


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "mRNA_diff_expr")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## output files for final plots
file_SRSF11_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_ctrl_hgat_SFs.png")
file_SRSF11_corr_plot <- file.path(analysis_dir, "plots", "corr_rna_vs_psi_SRSF11.png")

## get SRSF11 psi table
file <- "/SRSF11_psi.txt"
psi_tab  <-  read.delim(paste0(input_dir, file), sep = "\t", header=TRUE)

## transform and melt data for ggplot
melt_psi_tab <- melt(psi_tab, id=c("sample","id","histology"), variable.name =c("Type"))

## stacked barplot 
melt_psi_tab %>% 
  #arrange(desc(value)) %>%
  mutate(sample=fct_reorder(sample,value)) %>% 
  ggplot(aes(x = sample,y = value, fill= Type ))+ 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c("red",
                                                                           "blue")) 
ggsave(
  file_SRSF11_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =2047,
  height = 800,
  units = "px",
  dpi = 100,
  limitsize = TRUE,
  bg = NULL
)

## scatter plot of psi vs pct expression

file <- "/SRSF11_splicing_vs_expr.txt"
tab  <-  read.delim(paste0(input_dir, file), sep = "\t", header=TRUE)
tab = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/SRSF11_splicing_vs_expr.txt", header=TRUE, sep="\t")
ggscatter(tab, x="Expr", y="dPSI", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          add.params = list(color = "blue",
                            fill = "lightgray"),
          ticks = TRUE,
          xticks.by = .1, yticks.by = .1,
          xlab = "Expr", ylab = "dPSI") + theme_Publication()

"blue")) 

ggsave(
  file_SRSF11_corr_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =400,
  height = 400,
  units = "px",
  dpi = 100,
  limitsize = TRUE,
  bg = NULL
)

