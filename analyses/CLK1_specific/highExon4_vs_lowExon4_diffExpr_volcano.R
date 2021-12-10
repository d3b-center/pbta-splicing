library(dplyr)
library("sva")
library("EnhancedVolcano")
library("DESeq2")

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
#root_dir <- "/Users/naqvia/Desktop/pbta-splicing_git/pbta-splicing"
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## get data files and make table
input      = file.path(data_dir,"stranded_gene_counts_tab.tsv")
tab_rsem <- read.delim(input, header=TRUE, row.names=1)

head(tab_rsem)

tab_rsem = tab_rsem %>%
  select(c('gene_id', 'BS_Q13FQ8FV', 'BS_ZV1P6W9C','BS_WH8G4VFB','BS_NNPEC7W1','BS_PZVHMSYN', 
           'BS_DRY58DTF','BS_GXTFW99H','BS_E60JZ9Z3','BS_9CA93S6D'))


countTable <- tab_rsem[rowSums(tab_rsem>=10) > 9, ]
rownames(countTable) <- countTable$gene_id
countTable <- (countTable[-c(1)])

## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("High",5), rep("Low",4)))
condition = c(rep("High",5), rep("Low",4) )

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
                ylim = c(0,8),
                xlim = c(-4,4),
                title = 'High Exon 4  versus Low Exon 4 Tumors',
                pCutoff = 0.005,
                FCcutoff = 2,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                #drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black')

plot_file = file.path(plots_dir,"highExon4_vs_lowExon4.volano.pdf") 
ggsave(plot_file, width = 20, height = 15)

EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c("CYBB",
                              "NCF2",
                              "CIITA",
                              "IFI30",
                              "ADHFE1",
                              "PSMD9",
                              "PSMD13",
                              "TAP1",
                              "PSMA7",
                              "PSME3",
                              "PSMB1",
                              "PSMC3",
                              "PSMD1",
                              "PSME4",
                              "PSMB4",
                              "PSMB2",
                              "PSMD14",
                              "SEC13",
                              "PSMB8",
                              "PSMB9",
                              "PSMD8",
                              "PSMD5",
                              "PSMA6",
                              "PSMB3",
                              "PSMB5",
                              "PSMA2",
                              "PSMB6",
                              "PSMC1",
                              "PSMA5",
                              "PSMA4",
                              "PSME1",
                              "PSME2",
                              "SEC24D",
                              "PSMA3",
                              "SEC23A"),
                ylim = c(0,8),
                xlim = c(-4,4),
                title = 'High Exon 4  versus Low Exon 4 Tumors',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 2,
                labSize = 4,
                typeConnectors = "closed",
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'black')

file = "/Users/naqvia/Desktop/IR_proj/MYO1B-PSI_vs_MYO1B-Expr.txt"
tab = read.delim(file, sep="\t")
tab <- (tab[,2:3])

#convert NAs to 0
tab[is.na(tab)] <- 0

## remove columns that have stdv of 0
tab <- tab[, sapply(tab, function(x) { sd(x) != 0} )]
res <- cor((tab))

Q <- quantile(tab$breaks, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(tab$breaks)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

eliminated<- subset(tab, tab$breaks > (Q[1] - 1.5*iqr) & tab$breaks < (Q[2]+1.5*iqr))



ggscatter(tab, x="PSI", y="MYO1B", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          add.params = list(color = "blue",
                            fill = "lightgray"),
          ticks = TRUE,
          #xticks.by = .1, yticks.by = .1,
          xlab = "PSI", ylab = "HScore") + theme_Publication()
__DATA__
BS_Q13FQ8FV
BS_ZV1P6W9C
BS_WH8G4VFB
BS_NNPEC7W1
BS_PZVHMSYN

CLK1 low exon inclusion (group 2)
BS_XM1AHBDJ
BS_DRY58DTF
BS_GXTFW99H
BS_E60JZ9Z3
BS_9CA93S6D