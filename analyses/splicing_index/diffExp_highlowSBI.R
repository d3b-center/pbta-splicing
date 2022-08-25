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

##theme for all plots
theme_Publication <- function(base_size=15, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            # legend.margin = unit(0.5, "cm"),
            legend.margin = margin(5,5,5,5),
            legend.title = element_text(face="bold"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            plot.margin=unit(c(10,5,5,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

## set directories

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
#root_dir <- "/Users/naqvia/Desktop/pbta-splicing_git/pbta-splicing"
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## get splicing index table
file <- "/splicing_index.total.txt"
splice_index  <-  read.delim(paste0(results_dir, file), sep = "\t", header=TRUE) %>% rename(Kids_First_Biospecimen_ID = Sample)


splicing_index_outliers_HGAT <- filter(splice_index, Histology!="HGAT") 
SI_total_high_HGAT     <- quantile(splicing_index_outliers_HGAT$SI, probs=.75, names=FALSE)
SI_total_low_HGAT      <- quantile(splicing_index_outliers_HGAT$SI, probs=.25, names=FALSE)

splicing_index_outliers_HGAT <-  filter(splicing_index_outliers_HGAT, SI <SI_total_low_HGAT | SI >SI_total_high_HGAT  ) %>% 
  mutate(level=case_when(SI < SI_total_low_HGAT ~ "Low",
                         SI >SI_total_high_HGAT  ~ "High" ))


## get data files and make table
input      = file.path(input_dir,"tab_rsem.str.sbi.hgat.txt")
tab_rsem <- read.delim(input, header=TRUE, row.names=1)

head(tab_rsem)

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

file_volc_hgat_plot = "/enhancedVolcano_hgat_sbi.png"
filename = paste0(plots_dir, file_volc_hgat_plot)
ggsave(
  filename,
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

##non HGAT
## get data files and make table
input      = file.path(input_dir,"tab_rsem.str.sbi.non-hgat.txt")
tab_rsem <- read.delim(input, header=TRUE, row.names=1)

head(tab_rsem)

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

file_volc_non_hgat_plot = "/enhancedVolcano_nonhgat_sbi.png"
filename = paste0(plots_dir, file_volc_non_hgat_plot)
ggsave(
  filename,
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
