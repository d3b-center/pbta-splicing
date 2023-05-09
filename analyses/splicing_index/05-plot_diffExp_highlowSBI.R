################################################################################
# 05-diffExp_highlowSBI.R
#
# generates volcanto plot of diff expression b/w low vs high splicing burden 
# tumors using output from generate_splicing_index_tab_using_tumors.pl
#
# written by Ammar Naqvi
#
# usage: Rscript 05-diffExp_highlowSBI.R
################################################################################

## libraries needed
suppressPackageStartupMessages({
  library("dplyr")
  library("EnhancedVolcano")
  library("DESeq2")
  })


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
file_volc_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_hgat_sbi.pdf")
file_volc_non_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_nonhgat_sbi.pdf")
file_tiff_volc_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_hgat_sbi.tiff")
file_tiff_volc_non_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_nonhgat_sbi.tiff")

sbi_coding_file  <- file.path(results_dir,"splicing_index.total.txt")
sbi_coding_df  <-  vroom(sbi_coding_file, comment = "#",delim="\t") %>% mutate(Kids_First_Biospecimen_ID=Sample) %>% filter(Histology=="HGG")

clin_file = file.path(data_dir, "histologies.tsv")
clin_df  <-  vroom(clin_file, comment = "#",delim="\t") 

sbi_ids_clin <- clin_df %>% inner_join(sbi_coding_df, by="Kids_First_Biospecimen_ID") 


## compute quantiles to define high vs low Exon 4 SBI tumors
quartiles_sbi <- quantile(sbi_ids_clin$SI, probs=c(.25, .75), na.rm = FALSE)
IQR_sbi <- IQR(sbi_ids_clin$SI)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

## subset tmb values and samples by high vs low SBI tumors
high_sbi_df <- dplyr::filter(sbi_ids_clin, SI > upper_sbi) %>% dplyr::mutate(SBI_level="high")
low_sbi_df  <- dplyr::filter(sbi_ids_clin, SI < lower_sbi) %>% dplyr::mutate(SBI_level="low")
high_vs_low_df <- rbind(low_sbi_df,high_sbi_df)

#count table for HGGs 

## get gene count table with  HGGs and stranded filter
stranded_samples_only <- high_vs_low_df %>% filter(RNA_library == "stranded")

file_gene_counts = "gene-counts-rsem-expected_count-collapsed.rds" 
count_data <- readRDS(paste0(data_dir,"/", file_gene_counts)) %>%  
              select(any_of(stranded_samples_only$Kids_First_Biospecimen_ID)) 

## HGAT differential gene expression analysis
# remove low expression genes
filtered.counts <- count_data[rowSums(count_data>=10) >= 69, ]
countTable <- filtered.counts

## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("Low",37), rep("High",32)))

condition = c(rep("Low",37), rep("High",32) )

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)
#res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

volc_hgat_plot <- EnhancedVolcano(res,
                  lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                  x = 'log2FoldChange',
                  y = 'padj',
                  ylim = c(0,30),
                  xlim = c(-3,3),
                  title = 'Low vs High SBI (HGGs)',
                  subtitle = NULL,
                  caption = NULL,
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 1,
                  labSize = 3,
                  typeConnectors = "closed",
                  #drawConnectors = TRUE,
                  widthConnectors = 0.15,
                  colConnectors = 'black')


# Save plot as PDF
pdf(file_volc_hgat_plot, 
    width = 8, height = 8)
volc_hgat_plot
dev.off()

# Save plot tiff version
tiff(file_tiff_volc_hgat_plot, height = 1800, width = 2400, res = 300)
print(volc_hgat_plot)
dev.off()

##non-HGGs
sbi_coding_file  <- file.path(results_dir,"splicing_index.total.txt")
sbi_coding_nonHGG_df  <-  vroom(sbi_coding_file, comment = "#",delim="\t") %>% mutate(Kids_First_Biospecimen_ID=Sample) %>% filter(Histology!="HGG")
sbi_ids_clin <- clin_df %>% inner_join(sbi_coding_nonHGG_df, by="Kids_First_Biospecimen_ID") 

## compute quantiles to define high vs low Exon 4 SBI tumors
quartiles_sbi <- quantile(sbi_ids_clin$SI, probs=c(.25, .75), na.rm = FALSE)
IQR_sbi <- IQR(sbi_ids_clin$SI)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

## subset tmb values and samples by high vs low SBI tumors
high_sbi_df <- dplyr::filter(sbi_ids_clin, SI > upper_sbi) %>% dplyr::mutate(SBI_level="high")
low_sbi_df  <- dplyr::filter(sbi_ids_clin, SI < lower_sbi) %>% dplyr::mutate(SBI_level="low")
high_vs_low_nonHGG_df <- rbind(low_sbi_df,high_sbi_df)

#count table for HGGs 
## get gene count table with  non-HGGs and stranded filter
stranded_samples_nonHGG_only <- high_vs_low_nonHGG_df %>% filter(RNA_library == "stranded")

file_gene_counts = "gene-counts-rsem-expected_count-collapsed.rds" 
count_data_nonHGG <- readRDS(paste0(data_dir,"/", file_gene_counts)) %>%  
  select(any_of(stranded_samples_nonHGG_only$Kids_First_Biospecimen_ID)) 

## non-HGG differential gene expression analysis
# remove low expression genes
filtered.counts <- count_data[rowSums(count_data_nonHGG>=5) >= 289, ]
countTable <- filtered.counts


## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("Low",145), rep("High",144)))

condition = c(rep("Low",145), rep("High",144) )

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)

res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

volc_non_hgat_plot <- EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,130),
                xlim = c(-3,3),
                title = 'Low vs High SBI (non HGGs)',
                subtitle = NULL,
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                #drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black')

# Save plot as PDF
pdf(file_volc_non_hgat_plot, 
    width = 8, height = 8)
volc_non_hgat_plot
dev.off()

# Save plot tiff version
tiff(file_tiff_volc_non_hgat_plot, height = 1800, width = 2400, res = 300)
print(volc_non_hgat_plot)
dev.off()
