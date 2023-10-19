

#load libraries
library(tidyverse)
library(openxlsx)

#Magrittr pipe
`%>%` <- dplyr::`%>%`

#directory containing input files
input_dir <- "./input"
#directory to store plots
plot_dir <- "./plots"

#read histology file to data frame
histologyfile <- "Hope-GBM-histologies-base.tsv"
HISTdata <- readr::read_tsv(file.path(input_dir, histologyfile))
#

#read whole cell proteomics data to data frame
totalprotfile <- "Hope_proteome_imputed_data_liftover.tsv"
WCPdata <- readr::read_tsv(file.path(input_dir, totalprotfile))
#

#get all WCP samples ids
WCP_KFid <- colnames(WCPdata)[4:ncol(WCPdata)]
#

#read phosphoproteomics data to data frame
phosprotfile <- "Hope_phosphosite_imputed_data_ischemia_removed_liftover.tsv"
PHOSdata <- readr::read_tsv(file.path(input_dir, phosprotfile))
#

#get all PHOS samples ids
PHOS_KFid <- colnames(PHOSdata)[9:ncol(PHOSdata)]
#

#read splicing events data to data frame
splicingfile <- "splicing_events-psi.tsv"
SPLICEdata <- readr::read_tsv(file.path(input_dir, splicingfile))
#

#get sample IDs from splicing file. these are the IDs used in the HGG splicing analysis
SPLICE_KFid <- colnames(SPLICEdata)[9:ncol(SPLICEdata)]
#

#function to draw violin plots of abundance distributions with selected gene abundance highlighted
MakeAbundanceViolinPlots <- function(select_gene, WCPdata, WCP_KFid){
  
  #create violin plot to show relative abundances of CKL1 across samples in which it is identified
  #tibble to store abundances and genes for violin plot
  AllSamples_AbdTotDataframe <- tibble("Gene"=character(),"Abd"=numeric(),"Condition"=character())
  gene_symbols<- "ApprovedGeneSymbol"
  #parse all samples to add to violin plot data frame
  for ( i in 1:length(WCP_KFid) ){
    
    #current sample ID
    sampleID<-WCP_KFid[i]
    
    #skip this sample if the selected gene is not identified
    if ( is.na(dplyr::filter(WCPdata, ApprovedGeneSymbol %in% select_gene) %>% select(sampleID)) ){
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
  #violin plot 
  file_name<-paste("SampleAbundanceDist", select_gene, "violin_NEW.pdf", sep="_")
  pdf(file.path(".", plot_dir, file_name), height=6, width=12)
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


#draw violin plots of abundance distributions with selected gene abundance highlighted
genes_to_plot<- SPLICEdata %>% dplyr::pull(gene)
#parse genes to plot
for ( i in 1:length(genes_to_plot) ){
  #gene to highlight in violin plots
  select_gene <- genes_to_plot[i]
  #skip if gene is not in WCP data
  if ( dim(dplyr::filter(WCPdata, ApprovedGeneSymbol %in% select_gene))[1]==0 ){
    next
  }
  #call function to draw violin plots of abundance distributions with selected gene abundance highlighted
  MakeAbundanceViolinPlots(select_gene, WCPdata, WCP_KFid)
}