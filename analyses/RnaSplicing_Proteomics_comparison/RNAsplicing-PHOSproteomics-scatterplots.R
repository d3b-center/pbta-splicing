

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
if(!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive=TRUE)
}


##read histology file to data frame
###histologyfile_orig <- "Hope-GBM-histologies-base.tsv"
###HISTdata <- readr::read_tsv(file.path(input_dir, histologyfile_orig))
histologyfile <- "histologies.tsv"
HISTdata <- readr::read_tsv(file.path(data_dir, histologyfile))
###head(HISTdata)
#


#read phosphoproteomics data to data frame
###phosprotfile <- "Hope_phosphosite_imputed_data_ischemia_removed_liftover.tsv"
###PHOSdata <- readr::read_tsv(file.path(input_dir, phosprotfile))
phosprotfile <- "hope-protein-imputed-phospho-expression-abundance.tsv.gz"
PHOSdata <- readr::read_tsv(file.path(data_dir, phosprotfile))
head(PHOSdata)
#

##get all PHOS samples ids
PHOS_KFid <- colnames(PHOSdata)[6:ncol(PHOSdata)]
print(PHOS_KFid)
#

##read splicing events data to data frame
#splicing events file prepared by naqvia with proteins of interest
splicingfile <- "splicing_events-psi.tsv"
SPLICEdata <- readr::read_tsv(file.path(input_dir, splicingfile))
print(SPLICEdata)
#



#get sample IDs from splicing file. these are the IDs used in the HGG splicing analysis
SPLICE_KFid <- colnames(SPLICEdata)[9:ncol(SPLICEdata)]
print(SPLICE_KFid)
#



#RNA abundance table. select table to compare based on RNA abundance as measured by TPM or EC
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
head(SpliceID_SampleID_lookup)
#
#make lookup table of phos KF ids and corresponding sample id
PhosID_SampleID_lookup <- HISTdata %>%
  dplyr::filter(Kids_First_Biospecimen_ID%in%PHOS_KFid) %>%
  dplyr::select(Kids_First_Biospecimen_ID,sample_id)
head(PhosID_SampleID_lookup)
#


#function to draw scatter plots of RNA psi or total abundance values vs phos abundance for specific splicing events
MakeRnaPhosScatterPlots <- function(select_gene, select_site, HISTdata, RNAdata, SPLICEdata, SPLICE_KFid, PHOSdata, PHOS_KFid, SpliceID_SampleID_lookup, PhosID_SampleID_lookup, plotfile, RNAvalue){

  #create scatter plot to compare splicing psi values to phosphorylation abundance
  #data frame to store psi and corresponding phos abundance values for each HOPE sample
  RNA_Phos_Data <- tibble("RNA"=numeric(),"Phos"=numeric())
  #parse each splice sample and find corresponding phos data
  for ( i in 1:length(SPLICE_KFid) ){
    
    #current splicing kids first ID
    currSpliceKFid <- SPLICE_KFid[i]
    
    #corresponding sample_id
    currSampleId <- dplyr::filter(SpliceID_SampleID_lookup, Kids_First_Biospecimen_ID %in% {{currSpliceKFid}}) %>% pull(sample_id)
    #if the sample_ID does not exist in the HOPE cohort move to next
    if ( length(currSampleId)>1 ){
      print(currSpliceKFid)
      print ("currSpliceKFid multiple in HOPE cohort")
    }
    if ( length(currSampleId)==0 ){
      next
    }
    
    #phosphorylation kids first ID corresponding to current sample_ID
    currPhosKFid <- dplyr::filter(PhosID_SampleID_lookup, sample_id %in% {{currSampleId}}) %>% pull(Kids_First_Biospecimen_ID)
    #if the kids first ID does not exist in the HOPE cohort move to next
    if ( length(currPhosKFid)>1 ){
      print(currPhosKFid)
      print("currPhosKFid multiple in phosphoproteoimcs file")
    }
    if ( length(currPhosKFid)==0 ){
      next
    }
    
    #psi and phos data for current sample
    curr_RNA <- dplyr::filter(RNAdata, gene %in% {{select_gene}}) %>% pull({{currSpliceKFid}})
    if ( length(curr_RNA)>1 ){ print("multiple RNA") } 
    curr_Psi <- dplyr::filter(SPLICEdata, gene %in% {{select_gene}}) %>% pull({{currSpliceKFid}})
    if ( length(curr_Psi)>1 ){ print("multiple Psi") }
    ###curr_phos <- dplyr::filter(PHOSdata, Site %in% {{select_site}}) %>% pull({{currPhosKFid}})
    curr_phos <- dplyr::filter(PHOSdata, GeneSymbol %in% {{select_gene}}, Site %in% {{select_site}}) %>% pull({{currPhosKFid}})
    if ( length(curr_phos)>1 ){ print("multiple phos") }
    ###print(curr_RNA)
    ###print(curr_Psi)
    ###print(curr_phos)
    #curr_phos <- 2^curr_phostmp
    
    #add psi and phos data to date frame
    #can filter by RNA and/or Psi values
    #change curr_RNA to filter based on TPM or EC value
    if ( curr_RNA>19 ){
      ifelse ( RNAvalue=="Abd",
               RNA_Phos_Data <- RNA_Phos_Data %>% add_row(RNA={{curr_RNA}},Phos={{curr_phos}}),
               ifelse(RNAvalue=="Psi",
                      RNA_Phos_Data <- RNA_Phos_Data %>% add_row(RNA={{curr_Psi}},Phos={{curr_phos}}),
                      stop("undef RNA value")
               )
      )
    }
    #

    
  }#end SPLICE_KFid for loop
  #
  
  #scatter plot
  pdf(file.path(".", plot_dir, plotfile), height=6, width=6)
  RNA_Phos_plot<-ggplot(RNA_Phos_Data, aes(x=RNA, y=Phos)) +
    geom_point(size=2,color="darkorchid") +
    geom_smooth(method=lm, color="grey40") +
    theme_classic() + theme(axis.line = element_line(colour = 'black', linewidth = 1))
  print(RNA_Phos_plot)
  dev.off()
  
  #+ scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(18.5, 20.25))
  
  return()
  
} #end MakeRnaPhosScatterPlots function


#compare RNA splicing and abundance and protein phosphorylation for selected gene and phospho site
select_gene<-"CLK1"
###select_site<-"NP_001155879.1_182_200_1_1_S182"
select_site<-"S182"

#head(PhosID_SampleID_lookup)
#testPhosKFid<-"BS_085TYY7Q_PHOS"
#test_curr_phos <- dplyr::filter(PHOSdata, GeneSymbol %in% {{select_gene}}, Site %in% {{select_site}}) %>% pull({{testPhosKFid}})
#length(test_curr_phos)

#
#call function to draw scatter plots of psi values vs phos abundance for specific splicing events
plotfile<-file_name<-paste("RnaPsiPhosAbundance_RnaTPMgt20", select_gene, "scatter_Jan18.pdf", sep="_")
MakeRnaPhosScatterPlots(select_gene, select_site, HISTdata, RNAdata, SPLICEdata, SPLICE_KFid, PHOSdata, PHOS_KFid, SpliceID_SampleID_lookup, PhosID_SampleID_lookup, plotfile, "Psi")
#
#call function to draw scatter plots of RNA abundance values vs phos abundance for specific splicing events
plotfile<-file_name<-paste("RnaTPMAbundancePhosAbundance_RnaTPMgt20", select_gene, "scatter_Jan18.pdf", sep="_")
MakeRnaPhosScatterPlots(select_gene, select_site, HISTdata, RNAdata, SPLICEdata, SPLICE_KFid, PHOSdata, PHOS_KFid, SpliceID_SampleID_lookup, PhosID_SampleID_lookup, plotfile, "Abd")
