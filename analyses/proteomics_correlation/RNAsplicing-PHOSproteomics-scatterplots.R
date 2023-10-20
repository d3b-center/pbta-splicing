

#load libraries
library(tidyverse)
library(openxlsx)

#Magrittr pipe
`%>%` <- dplyr::`%>%`

#set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
#directory containing input files
input_dir <- "./input"
#directory to store plots
plot_dir <- "./plots"
if(!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive=TRUE)
}
#directory containing module data
data_dir <- file.path(root_dir, "data")

#read histology file to data frame
histologyfile <- "Hope-GBM-histologies-base.tsv"
HISTdata <- readr::read_tsv(file.path(input_dir, histologyfile))
#

######read whole cell proteomics data to data frame
#####totalprotfile <- "Hope_proteome_imputed_data_liftover.tsv"
#####WCPdata <- readr::read_tsv(file.path(input_dir, totalprotfile))
######
#####
######get all WCP samples ids
#####WCP_KFid <- colnames(WCPdata)[4:ncol(WCPdata)]
######

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

#RNA abundance table
file_gene_counts = "gene-expression-rsem-tpm-collapsed.rds" 
#file_gene_counts = "gene-counts-rsem-expected_count-collapsed.rds" 
RNAdata_df <- readRDS(paste0(data_dir,"/", file_gene_counts)) 
RNAdata_df$gene<-row.names(RNAdata_df)
rownames(RNAdata_df)<-1:nrow(RNAdata_df)
RNAdata<-as_tibble(RNAdata_df)
RNAdata %>%
  relocate(gene) %>%
  head()
#

#make lookup table of splice KF ids and corresponding sample id
SpliceID_SampleID_lookup <- HISTdata %>%
  dplyr::filter(Kids_First_Biospecimen_ID%in%SPLICE_KFid) %>%
  dplyr::select(Kids_First_Biospecimen_ID,sample_id)

#make lookup table of phos KF ids and corresponding sample id
PhosID_SampleID_lookup <- HISTdata %>%
  dplyr::filter(Kids_First_Biospecimen_ID%in%PHOS_KFid) %>%
  dplyr::select(Kids_First_Biospecimen_ID,sample_id)



#function to draw scatter plots of psi values vs phos abundance for specific splicing events
MakePsiPhosScatterPlots <- function(select_gene, select_site, HISTdata, SPLICEdata, SPLICE_KFid, PHOSdata, PHOS_KFid, SpliceID_SampleID_lookup, PhosID_SampleID_lookup){

  
  #create scatter plot to compare splicing psi values to phosphorylation abundance
  #data frame to store psi and corresponding phos abundance values for each HOPE sample
  Psi_Phos_Data <- tibble("Psi"=numeric(),"Phos"=numeric())
  #parse each splice sample and find corresponding phos data
  for ( i in 1:length(SPLICE_KFid) ){
    
    #current splicing kids first ID
    currSpliceKFid <- SPLICE_KFid[i]
    
    #corresponding sample_id
    currSampleId <- dplyr::filter(SpliceID_SampleID_lookup, Kids_First_Biospecimen_ID %in% {{currSpliceKFid}}) %>% pull(sample_id)
    #if the sample_ID does not exist in the HOPE cohort move to next
    if ( length(currSampleId)==0 ){ next }
    
    #phosphorylation kids first ID corresponding to current sample_ID
    currPhosKFid <- dplyr::filter(PhosID_SampleID_lookup, sample_id %in% {{currSampleId}}) %>% pull(Kids_First_Biospecimen_ID)
    #if the kids first ID does not exist in the HOPE cohort move to next
    if ( length(currPhosKFid)==0 ){ next }
    
    #psi and phos data for current sample
    curr_psi <- dplyr::filter(SPLICEdata, gene %in% {{select_gene}}) %>% pull({{currSpliceKFid}})
    curr_phos <- dplyr::filter(PHOSdata, Site %in% {{select_site}}) %>% pull({{currPhosKFid}})
    
    #add psi and phos data to date frame
    Psi_Phos_Data <- Psi_Phos_Data %>% add_row(Psi={{curr_psi}},Phos={{curr_phos}})
    
  }
  #
  #scatter plot
  file_name<-paste("PsiPhosAbd", select_gene, "scatter.pdf", sep="_")
  pdf(file.path(".", plot_dir, file_name), height=6, width=6)
  Psi_Phos_plot<-ggplot(Psi_Phos_Data, aes(x=Psi, y=Phos)) +
    geom_point(size=2,color="darkorchid") +
    geom_smooth(method=lm, color="grey40") +
    theme_classic() + theme(axis.line = element_line(colour = 'black', linewidth = 1)) +
    scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(18.5, 20.25))
  print(Psi_Phos_plot)
  dev.off()
  
  return()
  
} #end MakePsiPhosScatterPlots function


#function to draw scatter plots of psi values vs phos abundance for specific splicing events
MakeRNAabdPhosScatterPlots <- function(select_gene, select_site, HISTdata, RNAdata, SPLICE_KFid, PHOSdata, PHOS_KFid, SpliceID_SampleID_lookup, PhosID_SampleID_lookup){
  
  
  #create scatter plot to compare RNA abundance values to phosphorylation abundance
  #data frame to store RNA abundance and corresponding phos abundance values for each HOPE sample
  RNA_Phos_Data <- tibble("RNA"=numeric(),"Phos"=numeric())
  #parse each splice sample and find corresponding phos data
  for ( i in 1:length(SPLICE_KFid) ){
    
    #current splicing kids first ID
    currSpliceKFid <- SPLICE_KFid[i]
    
    #corresponding sample_id
    currSampleId <- dplyr::filter(SpliceID_SampleID_lookup, Kids_First_Biospecimen_ID %in% {{currSpliceKFid}}) %>% pull(sample_id)
    #if the sample_ID does not exist in the HOPE cohort move to next
    if ( length(currSampleId)==0 ){ next }
    
    #phosphorylation kids first ID corresponding to current sample_ID
    currPhosKFid <- dplyr::filter(PhosID_SampleID_lookup, sample_id %in% {{currSampleId}}) %>% pull(Kids_First_Biospecimen_ID)
    #if the kids first ID does not exist in the HOPE cohort move to next
    if ( length(currPhosKFid)==0 ){ next }
    
    #psi and phos data for current sample
    curr_RNA <- dplyr::filter(RNAdata, gene %in% {{select_gene}}) %>% pull({{currSpliceKFid}})
    curr_phos <- dplyr::filter(PHOSdata, Site %in% {{select_site}}) %>% pull({{currPhosKFid}})
    
    #add psi and phos data to date frame
    RNA_Phos_Data <- RNA_Phos_Data %>% add_row(RNA={{curr_RNA}},Phos={{curr_phos}})
    
  }
  #
  #scatter plot
  file_name<-paste("RNAPhosAbd", select_gene, "scatter_TPM.pdf", sep="_")
  pdf(file.path(".", plot_dir, file_name), height=6, width=6)
  RNA_Phos_plot<-ggplot(RNA_Phos_Data, aes(x=RNA, y=Phos)) +
    geom_point(size=2,color="darkcyan") +
    geom_smooth(method=lm, color="grey40") +
    theme_classic() + theme(axis.line = element_line(colour = 'black', linewidth = 1)) +
    scale_x_continuous(limits = c(0, 200)) + scale_y_continuous(limits = c(18.5, 20.25))
  print(RNA_Phos_plot)
  dev.off()
  
  #
  
  return()
  
} #end MakeRNAabdPhosScatterPlots function

#compare RNA splicing and abundance and protein phosphorylation for selected gene and phospho site
select_gene<-"CLK1"
select_site<-"NP_001155879.1_182_200_1_1_S182"
#
#call function to draw scatter plots of psi values vs phos abundance for specific splicing events
MakePsiPhosScatterPlots(select_gene, select_site, HISTdata, SPLICEdata, SPLICE_KFid, PHOSdata, PHOS_KFid, SpliceID_SampleID_lookup, PhosID_SampleID_lookup)
#
#call function to draw scatter plots of RNA abundance values vs phos abundance for specific splicing events
MakeRNAabdPhosScatterPlots(select_gene, select_site, HISTdata, RNAdata, SPLICE_KFid, PHOSdata, PHOS_KFid, SpliceID_SampleID_lookup, PhosID_SampleID_lookup)

