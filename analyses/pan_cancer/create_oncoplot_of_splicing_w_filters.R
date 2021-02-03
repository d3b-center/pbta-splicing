library(dplyr)
library(maftools)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

dataDir = "/Users/naqvia/Desktop/create_oncoprint/"
file_fus <- "reformattted_fusions_for_oncoprint.total.onlywithWGS.txt"
fus_tab  = read.delim(paste0(dataDir, file_fus), sep = "\t", header=TRUE)

file_spl <- "pan_cancer_splicing.thr10.report_all.formatted_for_oncoprint.v2.txt"
psi_tab  = read.delim(paste0(dataDir,file_spl), sep = "\t", header=TRUE)

file_maf <- "pbta-snv-consensus-mutation.maf.tsv"
maf_tab = read.delim(paste0(dataDir,file_maf), sep = "\t", header=TRUE)

clin_file = "~/Desktop/AS-DMG/data/pbta-histologies.tsv"
clin_tab = read.delim(paste0(clin_file), sep = "\t", header=TRUE)

clin_tab <- clin_tab %>% rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) 

## histology based filters
#all HGATs with WGS
sample_list <- read.delim("/Users/naqvia/Desktop/create_oncoprint/bsIDs.wgs.hgat_only.txt", sep = "\t",
                          header = F, as.is = T) 

#only DMG samples with WGS
sample_dmg_list <- read.delim("/Users/naqvia/Desktop/create_oncoprint/bsIDs.wgs.dmg_only.txt", sep = "\t",
                              header = F, as.is = T) 

## gene list filters
#kinase
gene_list_kinase <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/kinase_list.v2.txt", sep = "\t",
                        header = F, as.is = T)

#brain cancer genes
gene_list_brain <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/brain-goi-list-new.txt", sep = "\t",
                              header = F, as.is = T) 

#splicing factors
gene_list_sfs <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/splicing_factors.txt", sep = "\t",
                              header = F, as.is = T)

#SWI/SNF genes
gene_list_swisnf <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/swi_snf_genes.txt", sep = "\t",
                            header = F, as.is = T)

#top mis-spliced genes in DMG
gene_list_top_splicing_dmg_wgs <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/top10_dmg.wgs.splicing.txt", sep = "\t",
                                             header = F, as.is = T)  

#epigenetic related genes
gene_list_epi <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/epi_genes_list.grep.txt", sep = "\t",
                                             header = F, as.is = T) 

#tf related genes
gene_list_tf <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/tf_gene_list.txt", sep = "\t",
                                             header = F, as.is = T) 
## filter for only DMG samples
psi_tab_filtered_samples = psi_tab[psi_tab$Tumor_Sample_Barcode %in% sample_dmg_list$V1,]
fus_tab_filtered_samples = fus_tab[fus_tab$Tumor_Sample_Barcode %in% sample_dmg_list$V1,]
maf_tab_filtered_samples = maf_tab[maf_tab$Tumor_Sample_Barcode %in% sample_dmg_list$V1,]


#psi_tab_filtered_samples = psi_tab[psi_tab$Hugo_Symbol %in% c("SNHG14","OBSCN", "EPB41L1"),]
#fus_tab_filtered_samples = fus_tab[fus_tab$Hugo_Symbol %in% c("SNHG14","OBSCN", "EPB41L1"),]

## combine different types into one table for maf object
combine_tab     <- bind_rows(psi_tab_filtered_samples, maf_tab_filtered_samples,fus_tab_filtered_samples)



## filter for select genes from gene_list above
#maf_fus_filtered = maf_fus[maf_fus$Hugo_Symbol %in% gene_list$V1,]
#psi_tab_filtered = psi_tab[psi_tab$Hugo_Symbol %in% gene_list$V1,]
#fus_tab_filtered = fus_tab[fus_tab$Hugo_Symbol %in% gene_list$V1,]



#maf_fus_psi <- bind_rows(psi_tab_filtered, fus_tab_filtered)
#subset(maf_fus_psi, Hugo_Symbol == "SUPT5H")


## create and prepare maf object
maf = read.maf(maf = combine_tab, clinicalData = clin_tab, vc_nonSyn = c(  "Frame_Shift_Del",
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



col2hex(c('mediumseagreen', 'orange', 'green2', 'hotpink',
          'lightslateblue', 'midnightblue', 'gray33', 'turquoise3',
          'lightpink', 'lightsalmon', 'khaki1',
          'purple', 'black'))

colores = c("Missense_Mutation" = "#35978f", 
            "Nonsense_Mutation" = "#000000",
            "Frame_Shift_Del" = "#56B4E9", 
            "Frame_Shift_Ins" = "#FFBBFF", 
            "Splice_Site" = "#F0E442",
            "Translation_Start_Site" = "#191970",
            "Nonstop_Mutation" = "#545454",
            "In_Frame_Del" = "#CAE1FF",
            "In_Frame_Ins" = "#FFE4E1",
            "Stop_Codon_Ins" = "#CC79A7",
            "Start_Codon_Del" = "#56B4E9",
            "Fusion" = "#7B68EE",
            "Splicing" = "#f46d43")

mut.labels = c("Missense Mutation", 
               "Nonsense Mutation",
               "Frameshift Deletion", 
               "Frameshift Insertion", 
               "Splice Site Mutation",
               "Translation Start Site",
               "Nonstop Mutation",
               "Inframe Deletion",
               "Inframe Insertion",
               "Stop Codon Insertion",
               "Start Codon Deletion",
               "RNA Fusion",
               "Splicing")


## plot maf object
plotmafSummary(maf = maf)

# top mis-spliced genes (in DMG)
oncoplot(maf = maf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,genes=gene_list_top_splicing_dmg_wgs$V1, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = F, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# kinase-specific
oncoplot(maf = maf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,genes=gene_list_kinase$V1, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = F, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# brain cancer-specific
oncoplot(maf = maf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,genes=gene_list_brain$V1, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = F, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# splicing_factor-specific
oncoplot(maf = maf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,genes=gene_list_sfs$V1, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = F, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# swisnf-specific
oncoplot(maf = maf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,genes=gene_list_swisnf$V1, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = F, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= F)

# transcription_factor-specific
oncoplot(maf = maf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,genes=gene_list_tf$V1, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = F, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= F)

# epigenetic-specific
oncoplot(maf = maf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,genes=gene_list_epi$V1, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
         showTitle = T, colors = colores, logColBar = T, removeNonMutated= T)
