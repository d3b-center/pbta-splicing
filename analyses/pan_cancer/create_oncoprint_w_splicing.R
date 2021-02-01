library(dplyr)
library(maftools)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

dataDir = "~/Desktop/AS-DMG/analyses/pan_cancer/results/"
 
file_fus <- "reformattted_fusions_for_oncoprint.total.onlywithWGS.brain-goi-list.txt"
fus_tab  = read.delim(paste0(dataDir, file_fus), sep = "\t", header=TRUE)

file_spl <- "reformatted_splicing_for_oncoprint.total.onlywithWGS.brain-goi-list.txt"
psi_tab  = read.delim(paste0(dataDir,file_spl), sep = "\t", header=TRUE)

file_maf <- "pbta-snv-consensus-mutation.maf.brain-goi-list.v2.txt"
maf_tab = read.delim(paste0(dataDir,file_maf), sep = "\t", header=TRUE)

clin_file = "~/Desktop/AS-DMG/data/pbta-histologies.tsv"
clin_tab = read.delim(paste0(clin_file), sep = "\t", header=TRUE)

clin_tab <- clin_tab %>% rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) 


##only required fields
maf_tab_req <- maf_tab %>%
  select(Tumor_Sample_Barcode,Hugo_Symbol,Variant_Classification, Variant_Type, Chromosome,Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2)
maf_obj <-read.maf(maf = maf_tab)


###merge with MAF
#file_pptc = "~/Desktop/test_maf.csv"
#pptc =  read.csv(file_pptc)
maf_fus <- bind_rows(maf_tab_req, fus_tab)
maf_fus_spl <- bind_rows(maf_fus, psi_tab)

subset(maf_fus_spl, Variant_Classification == "Splicing")



#maftools::oncoplot(maf = pptc)
maf = read.maf(maf = maf_fus_spl, clinicalData = clin_tab, vc_nonSyn = c(  "Frame_Shift_Del",
                                                                           "Frame_Shift_Ins",
                                                                           "Splice_Site",
                                                                           "Translation_Start_Site",
                                                                           "Nonsense_Mutation",
                                                                           "Nonstop_Mutation",
                                                                           "In_Frame_Del",
                                                                           "In_Frame_Ins",
                                                                           "Missense_Mutation", 
                                                                           "Fusion", 
                                                                           "Splicing"))


##summary
plotmafSummary(maf =maf)
oncoplot(maf = maf, clinicalFeatures = c("short_histology"), showTumorSampleBarcodes = F,
         drawRowBar = T,
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
         showTitle = F, logColBar = T)
         





