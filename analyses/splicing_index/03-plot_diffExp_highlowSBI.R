################################################################################
# diffExp_highlowSBI.R
#
# generates volcanto plot of differentail expression between high vs low 
# splicing burden tumors using output from generate_splicing_index_tab_using_tumors.pl
#
# written by Ammar Naqvi
#
# usage: Rscript diffExp_highlowSBI.R
################################################################################

library("dplyr")
library("EnhancedVolcano")
library("DESeq2")

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`




## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_volc_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_hgat_sbi.png")
file_volc_non_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_nonhgat_sbi.png")

## retrieve and store input
# splicing index table
file <- "/splicing_index.total.txt"
splice_index  <-  read.delim(paste0(results_dir, file), sep = "\t", header=TRUE) %>% rename(Kids_First_Biospecimen_ID = Sample)

#count table for HGAT 
input      = file.path(input_dir,"tab_rsem.str.sbi.hgat.txt")
tab_rsem <- read.delim(input, header=TRUE, row.names=1)

#count table for non-HGAT  
input      = file.path(input_dir,"tab_rsem.str.sbi.non-hgat.txt")
tab_rsem <- read.delim(input, header=TRUE, row.names=1)

## grab HGAT samples and compute high vs low splicing burden index (SBI) values from splicing index table
# fitler for HGAT
splicing_index_outliers_HGAT <- filter(splice_index, Histology!="HGAT") 

# compute high vs low SBI values
SI_total_high_HGAT     <- quantile(splicing_index_outliers_HGAT$SI, probs=.75, names=FALSE)
SI_total_low_HGAT      <- quantile(splicing_index_outliers_HGAT$SI, probs=.25, names=FALSE)

## filter and get only HGAT outliers based on high and low SBI
splicing_index_outliers_HGAT <-  filter(splicing_index_outliers_HGAT, SI <SI_total_low_HGAT | SI >SI_total_high_HGAT  ) %>% 
  mutate(level=case_when(SI < SI_total_low_HGAT ~ "Low",
                         SI >SI_total_high_HGAT  ~ "High" ))


head(tab_rsem)

## HGAT differential gene expression analysis
# remove low expression genes
filtered.counts <- tab_rsem[rowSums(tab_rsem>=10) >= 38, ]
countTable <- filtered.counts


## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("High",36), rep("Low",40)))

condition = c(rep("High",36), rep("Low",40) )

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)

res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

file_volc_hgat_plot <- EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,8),
                xlim = c(-2,2),
                title = 'High vs Low SBI in HGATs',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                #drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black')

ggsave(
  file_volc_hgat_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =2.74,
  height = 2.38,
  units = "in",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

##non HGAT differential gene expression analysis
#filter low expressed genes
filtered.counts <- tab_rsem[rowSums(tab_rsem>=10) >= 145, ]
countTable <- filtered.counts

## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("High",145), rep("Low",145)))

condition = c(rep("High",145), rep("Low",145) )

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)

res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,30),
                xlim = c(-2,2),
                title = 'High vs Low SBI in non-HGATs',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                #drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black')

ggsave(
  file_volc_non_hgat_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =2.74,
  height = 2.38,
  units = "in",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)
