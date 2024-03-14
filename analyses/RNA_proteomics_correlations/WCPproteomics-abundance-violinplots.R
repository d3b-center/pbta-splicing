################################################################################
# WCPproteomics-abundance-violinplots.R
# Author: Joseph M. Dybas
#
# Plot violin plots of proteins abundance distribution for each sample
# and highlight proteins of interest in each distribution.
# 
# Requires a text file of proteins to plot.
#
# usage: Rscript WCPproteomics-abundance-violinplots.R name_of_inputfile
################################################################################

##load libraries
library(tidyverse)
library(openxlsx)

##Magrittr pipe
`%>%` <- dplyr::`%>%`

##get comand line argument for input file
args = commandArgs(trailingOnly=TRUE)
if ( length(args) != 1 ){
  stop("must provide proteins to plot file name")
}

##set directory variables
#root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
#general data directory
data_dir <- file.path(root_dir, "data")
#directory to store plots from analysis
plot_dir <- "./plots"

##read whole cell proteomics data to data frame
#HOPE whole cell proteome
HOPEtotalprotfile <- "hope-protein-imputed-prot-expression-abundance.tsv.gz"
WCPdata_HOPE <- readr::read_tsv(file.path(data_dir, HOPEtotalprotfile))
#head(WCPdata_HOPE)
#
#GBM whole cell proteome
GBMtotalprotfile <- "gbm-protein-imputed-prot-expression-abundance.tsv.gz"
WCPdata_GBM <- readr::read_tsv(file.path(data_dir, GBMtotalprotfile))
#head(WCPdata_GBM)
#
#join HOPE and GBM tibbles on GeneSymbol and NP_id columns
WCPdata<-full_join(WCPdata_HOPE, WCPdata_GBM, join_by("GeneSymbol","NP_id"))
#head(WCPdata)

##get all WCP samples ids
WCP_KFid <- colnames(WCPdata)[4:ncol(WCPdata)]
#

##read file of proteins to be plotted. file should be formatted as a list of protein IDs with no header
protsfile <- args[1]
prots_data <- readr::read_table(protsfile, col_names=FALSE)
#head(prots_data)

##function to draw violin plots of abundance distributions with selected gene abundance highlighted
MakeAbundanceViolinPlots <- function(select_prot, WCPdata, WCP_KFid){
  
  ##create violin plot to show relative abundances of selected gene across samples in which it is identified
  #tibble to store abundances and genes for violin plot
  AllSamples_AbdTotDataframe <- tibble("Gene"=character(),"Abd"=numeric(),"Condition"=character())
  gene_symbols<- "GeneSymbol"
  #parse all samples to add to violin plot data frame
  for ( i in 1:length(WCP_KFid) ){
    
    #current sample ID
    sampleID<-WCP_KFid[i]
    
    #skip this sample if the selected gene is not identified
    if ( is.na(dplyr::filter(WCPdata, GeneSymbol %in% select_prot) %>% select(sampleID)) ){
      next
    }
    
    #select the genes and abundances for the current sample ID
    #change the column names to correspond to the plotting tibble
    #add the "Condition" column to store the sample ID
    SampleData <- WCPdata %>%
      dplyr::select({{gene_symbols}},{{sampleID}}) %>%
      dplyr::rename('Gene' = {{gene_symbols}},
                    'Abd' = {{sampleID}}) %>%
      mutate(Condition = {{sampleID}})
    #combine the sample data to the plotting tibble
    AllSamples_AbdTotDataframe<-bind_rows(AllSamples_AbdTotDataframe,SampleData)
    
  }
  #
  ##draw violin plot 
  #file name
  file_name<-paste("SampleAbundanceDist", select_prot, "violin_NEW.pdf", sep="_")
  pdf(file.path(".", plot_dir, file_name), height=6, width=12)
  #ggplot
  AllSamples_AbdTotDataframe_violin<-ggplot(AllSamples_AbdTotDataframe, aes(x=Condition, y=Abd, fill=Condition)) +
    geom_violin(alpha=.3, width=0.5) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=1, col="grey") +
    stat_summary(fun.y=median, geom="point", shape=23, size=1, col="grey") +
    theme_classic() + theme(legend.position = 'none', axis.line = element_line(colour = 'black', size = 1), axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=rep("grey60",length(unique(AllSamples_AbdTotDataframe$Condition)))) +
    geom_point(data = dplyr::filter(AllSamples_AbdTotDataframe, Gene %in% select_prot), color = "blue")
  print(AllSamples_AbdTotDataframe_violin)
  dev.off()
  #
  
  return()
  
} #end MakeAbundanceViolinPlots function


##draw violin plots of abundance distributions with selected gene abundance highlighted
#get proteins of interest from input list
prots_to_plot <- prots_data %>% dplyr::pull(X1)
#print(prots_to_plot)
#parse prots to plot
for ( i in 1:length(prots_to_plot) ){

  #selected gene to highlight in violin plot
  select_prot <- prots_to_plot[i]
  #skip if gene is not in WCP data
  if ( dim(dplyr::filter(WCPdata, GeneSymbol %in% select_prot))[1]==0 ){
    next
  }
  #call function to draw violin plots of abundance distributions with selected gene abundance highlighted
  MakeAbundanceViolinPlots(select_prot, WCPdata, WCP_KFid)

} #end parse prots_to_plot
