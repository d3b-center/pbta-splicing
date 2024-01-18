

##load libraries
library(tidyverse)
library(openxlsx)

##Magrittr pipe
`%>%` <- dplyr::`%>%`

##set directory variables
#root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
#general data directory
data_dir <- file.path(root_dir, "data")
#directory containing input files for splicing files
input_dir <- "./input"
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

##read splicing events data to data frame
#splicing events file prepared by naqvia with proteins of interest
splicingfile <- "splicing_events-psi.tsv"
SPLICEdata <- readr::read_tsv(file.path(input_dir, splicingfile))
#

##get sample IDs from splicing file. these are the IDs used in the HGG splicing analysis
SPLICE_KFid <- colnames(SPLICEdata)[9:ncol(SPLICEdata)]
#

##function to draw violin plots of abundance distributions with selected gene abundance highlighted
MakeAbundanceViolinPlots <- function(select_gene, WCPdata, WCP_KFid){
  
  ##create violin plot to show relative abundances of selected gene across samples in which it is identified
  #tibble to store abundances and genes for violin plot
  AllSamples_AbdTotDataframe <- tibble("Gene"=character(),"Abd"=numeric(),"Condition"=character())
  gene_symbols<- "GeneSymbol"
  #parse all samples to add to violin plot data frame
  for ( i in 1:length(WCP_KFid) ){
    
    #current sample ID
    sampleID<-WCP_KFid[i]
    
    #skip this sample if the selected gene is not identified
    if ( is.na(dplyr::filter(WCPdata, GeneSymbol %in% select_gene) %>% select(sampleID)) ){
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
  file_name<-paste("SampleAbundanceDist", select_gene, "violin.pdf", sep="_")
  pdf(file.path(".", plot_dir, file_name), height=6, width=12)
  #ggplot
  AllSamples_AbdTotDataframe_violin<-ggplot(AllSamples_AbdTotDataframe, aes(x=Condition, y=Abd, fill=Condition)) +
    geom_violin(alpha=.3, width=0.5) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=1, col="grey") +
    stat_summary(fun.y=median, geom="point", shape=23, size=1, col="grey") +
    theme_classic() + theme(legend.position = 'none', axis.line = element_line(colour = 'black', size = 1), axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=rep("grey60",length(unique(AllSamples_AbdTotDataframe$Condition)))) +
    geom_point(data = dplyr::filter(AllSamples_AbdTotDataframe, Gene %in% select_gene), color = "blue")
  print(AllSamples_AbdTotDataframe_violin)
  dev.off()
  #
  
  return()
  
} #end MakeAbundanceViolinPlots function


##draw violin plots of abundance distributions with selected gene abundance highlighted
#get genes of interest from splice data
genes_to_plot<- SPLICEdata %>% dplyr::pull(gene)
#parse genes to plot
for ( i in 1:length(genes_to_plot) ){
  #selected gene to highlight in violin plot
  select_gene <- genes_to_plot[i]
  #skip if gene is not in WCP data
  if ( dim(dplyr::filter(WCPdata, GeneSymbol %in% select_gene))[1]==0 ){
    next
  }
  #call function to draw violin plots of abundance distributions with selected gene abundance highlighted
  MakeAbundanceViolinPlots(select_gene, WCPdata, WCP_KFid)
}
