library(dplyr)
library(maftools)
library(vroom)
library(ggplot2)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

dataDir = "~/Desktop/AS-DMG/data/"
file_fus <- "reformattted_fusions_for_oncoprint.total.onlywithWGS.txt"
fus_tab  = read.delim(paste0(dataDir, file_fus), sep = "\t", header=TRUE)

file_spl <- "reformatted_splicing_for_oncoprint.total.onlywithWGS.thr10.hgat.all.v3_2.Rplot.txt"
psi_tab  = read.delim(paste0(dataDir,file_spl), sep = "\t", header=TRUE)

file_maf <- "pbta-merged-chop-method-consensus_somatic.maf"
#maf_tab = read.delim(paste0(dataDir,file_maf), sep = "\t", header=TRUE)

maf_tab <-vroom(paste0(dataDir,file_maf), delim = '\t')

#maf_tab %<>% 
#  mutate_each(funs(if(is.integer(.)) as.numeric(.) else .))

cn_tab = read.delim("~/Desktop/AS-DMG/data/pbta-consensus_seg_annotated_cn_autosomes.v2.tsv",sep="\t",header=T, as.is=T)

## reading it 
clin_file = "~/Desktop/AS-DMG/data/compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv"
clin_tab = read.delim(paste0(clin_file), sep = "\t", header=TRUE)

clin_tab <- clin_tab %>% rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) 

## histology based filters
#all HGATs with WGS
hgat_sample_list = "sample_names.hgat_w_both.txt"
sample_hgat_list <- read.delim(paste0(dataDir,hgat_sample_list), sep = "\t",
                          header = F, as.is = T) 

#only DMG samples with WGS
sample_dmg_list <- read.delim("/Users/naqvia/Desktop/create_oncoprint/bsIDs.wgs.dmg_only.txt", sep = "\t",
                              header = F, as.is = T) 

## gene list filters
#kinase
gene_list_kinase <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/kinase_list.txt", sep = "\t",
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
#psi_tab_filtered_samples = psi_tab[psi_tab$Tumor_Sample_Barcode %in% sample_dmg_list$V1,]
#fus_tab_filtered_samples = fus_tab[fus_tab$Tumor_Sample_Barcode %in% sample_dmg_list$V1,]
#maf_tab_filtered_samples = maf_tab[maf_tab$Tumor_Sample_Barcode %in% sample_dmg_list$V1,]

## filter for all HGAT samples
psi_tab_filtered_samples = psi_tab[psi_tab$Tumor_Sample_Barcode %in% sample_hgat_list$V1,]
fus_tab_filtered_samples = fus_tab[fus_tab$Tumor_Sample_Barcode %in% sample_hgat_list$V1,]
maf_tab_filtered_samples = maf_tab[maf_tab$Tumor_Sample_Barcode %in% sample_hgat_list$V1,]
cn_tab_filtered = cn_tab[cn_tab$Tumor_Sample_Barcode %in% sample_hgat_list$V1,]

## combine different types into one table for maf object
combine_tab     <- bind_rows(psi_tab_filtered_samples, maf_tab_filtered_samples, fus_tab_filtered_samples)
#combine_tab <- bind_rows(psi_tab, maf_tab,fus_tab)

## filter for select genes from gene_list above
#maf_fus_filtered = maf_fus[maf_fus$Hugo_Symbol %in% gene_list$V1,]
#psi_tab_filtered = psi_tab[psi_tab$Hugo_Symbol %in% gene_list$V1,]
#fus_tab_filtered = fus_tab[fus_tab$Hugo_Symbol %in% gene_list$V1,]

#maf_fus_psi <- bind_rows(psi_tab_filtered, fus_tab_filtered)
#subset(maf_fus_psi, Hugo_Symbol == "SUPT5H")


## create and prepare maf object
#maf = read.maf(maf = combine_tab, clinicalData = clin_tab, removeDuplicatedVariants = FALSE, vc_nonSyn = c("Frame_Shift_Del",
#                                                                           "Frame_Shift_Ins",
#                                                                           "Splice_Site",
#                                                                           "Translation_Start_Site",
#                                                                           "Nonsense_Mutation",
#                                                                           "Nonstop_Mutation",
#                                                                           "In_Frame_Del",
#                                                                           "In_Frame_Ins",
#                                                                           "Missense_Mutation",
#                                                                           "Fusion",
#                                                                           "Splicing"))

maf = read.maf(maf = combine_tab, clinicalData = clin_tab, removeDuplicatedVariants = FALSE, vc_nonSyn = c("Frame_Shift_Del",
                                                                                                             "Frame_Shift_Ins",
                                                                                                             "Splice_Site",
                                                                                                             "Translation_Start_Site",
                                                                                                             "Nonsense_Mutation",
                                                                                                             "Nonstop_Mutation",
                                                                                                             "In_Frame_Del",
                                                                                                             "In_Frame_Ins",
                                                                                                             "Missense_Mutation",
                                                                                                             "Fusion",
                                                                                                             "Splicing",
                                                                                                             "Multi_Hit"))

#maf2 = read.maf(maf = combine_tab, clinicalData = clin_tab, removeDuplicatedVariants = FALSE,
#                vc_nonSyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", 
#                              "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation",  
#                              "Stop_Codon_Ins", "Start_Codon_Del", "Fusion", "Multi_Hit", "Hom_Deletion",
#                              "Hem_Deletion", "Amplification", "Multi_Hit_Fusion", "Splicing"))


maf_total = read.maf(maf = combine_tab, clinicalData = clin_tab, removeDuplicatedVariants = FALSE, cnTable = cn_tab_filtered, 
                     vc_nonSyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", 
                                   "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation",  
                                   "Stop_Codon_Ins", "Start_Codon_Del", "Fusion", "Multi_Hit", "Del", "Amp",  "Splicing"))


# col2hex=c('mediumseagreen', 'orange', 'green2', 'hotpink',
#           'lightslateblue', 'midnightblue', 'gray33', 'turquoise3',
#           'lightpink','purple', 'black', 'brown')
# 
# 
# colores = c("Missense_Mutation" = "#35978f", 
#             "Nonsense_Mutation" = "#000000",
#             "Frame_Shift_Del" = "#56B4E9", 
#             "Frame_Shift_Ins" = "#FFBBFF", 
#             "Splice_Site" = "#F0E442",
#             "Translation_Start_Site" = "#191970",
#             "Nonstop_Mutation" = "#545454",
#             "In_Frame_Del" = "#CAE1FF",
#             "In_Frame_Ins" = "#FFE4E1",
#             "Fusion" = "#7B68EE",
#             "Splicing" = "#f46d43",
#             "Multi_Hit" = "#545454")
# 
# mut.labels = c("Missense Mutation", 
#                "Nonsense Mutation",
#                "Frameshift Deletion", 
#                "Frameshift Insertion", 
#                "Splice Site Mutation",
#                "Translation Start Site",
#                "Nonstop Mutation",
#                "Inframe Deletion",
#                "Inframe Insertion",
#                "Fusion",
#                "Splicing",
#                "Multi-hit")


#col2hex_2=c('mediumseagreen', 'orange', 'green2', 'hotpink',
#          'lightslateblue', 'midnightblue', 'gray33', 'turquoise3',
#          'lightpink', 'lightsalmon', 'khaki1',
#          'purple', 'black', 'dodgerblue2', 'red', 'cornflowerblue', 'lightcoral', 'red', 'mistyrose', 'lightsteelblue1', 'blue')

col2hex=c('mediumseagreen', 'orange', 'green2', 'hotpink',
          'lightslateblue', 'midnightblue', 'gray33', 'turquoise3',
          'lightpink','purple', 'black', 'red', 'darkred', 'brown')


colores = c("Missense_Mutation" = "#35978f", 
            "Nonsense_Mutation" = "#000000",
            "Frame_Shift_Del" = "#56B4E9", 
            "Frame_Shift_Ins" = "#FFBBFF", 
            "Splice_Site" = "#F0E442",
            "Translation_Start_Site" = "#191970",
            "Nonstop_Mutation" = "#545454",
            "In_Frame_Del" = "#CAE1FF",
            "In_Frame_Ins" = "#FFE4E1",
            "Fusion" = "#7B68EE",
            "Multi_Hit" = "#f46d43",
            "Del" = "#0072B2",
            "Amp" = "#c51b7d",
            "Splicing" = "#191970")

mut_labels = c("Missense Mutation", 
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
               "Multi-Hit",
               "Amp",
               "Del",
               "Splicing")

## plot maf object
plotmafSummary(maf = maf)

oncoplot(maf = maf_total,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T, top=30 , 
         annotationFontSize = 1, gene_mar = 3, barcode_mar = 3,
         sortByAnnotation = F, fontSize = .5, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# brain cancer-specific
sub_maf_brain = maftools::subsetMaf(maf = maf_total, genes = gene_list_brain$V1, mafObj = TRUE)
oncoplot(maf = sub_maf_brain,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,top = 20,
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T,  removeNonMutated = T)

# kinase-specific
sub_maf_kinase = maftools::subsetMaf(maf = maf_total, genes = gene_list_kinase$V1, mafObj = TRUE)
oncoplot(maf = sub_maf_kinase,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,top=20, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# splicing_factor-specific
sub_maf_sf = maftools::subsetMaf(maf = maf_total, genes = gene_list_sfs$V1, mafObj = TRUE)
oncoplot(maf = sub_maf_sf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,top = 20, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# swisnf-specific
sub_maf_swisnf = maftools::subsetMaf(maf = maf_total, genes = gene_list_swisnf$V1, mafObj = TRUE)
oncoplot(maf = sub_maf_swisnf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,top=20, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# transcription_factor-specific
sub_maf_tf = maftools::subsetMaf(maf = maf_total, genes = gene_list_tf$V1, mafObj = TRUE)
oncoplot(maf = sub_maf_tf,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,top=20, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
         showTitle = F, colors = colores, logColBar = T, removeNonMutated= T)

# epigenetic-specific
sub_maf_epi = maftools::subsetMaf(maf = maf_total, genes = gene_list_epi$V1, mafObj = TRUE)
oncoplot(maf = sub_maf_epi,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T,top=20, 
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = F, fontSize = 1, #legendFontSize = 2,
         showTitle = T, colors = colores, logColBar = T, removeNonMutated= T)

## custom genes
gene_custom = c("DPF2","NF1", "TP53", "CLK1", "MAST2", "FGFR1", "SMARCA4", "BAK1", "METTL26", "BIN1", "NTRK1","H3F3A" )
#sub_maf_cust = maftools::subsetMaf(maf = maf, genes = gene_custom, mafObj = TRUE)
oncoplot(maf = maf_total,clinicalFeatures = c("molecular_subtype"), showTumorSampleBarcodes = F, drawRowBar = T, genes =gene_custom,
         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
         sortByAnnotation = T, fontSize = 1, 
         showTitle = T, colors = colores, logColBar = T, removeNonMutated= T)

##look for occurrences of splicing
added_brain_list <- c(gene_list_brain$V1, "CLK1")
somaticInteractions(maf = maf_total , genes = added, pvalue = c(0.05, 0.1))
