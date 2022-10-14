# Splicing Factor Dysregulation

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is assess and identify splicing factor dysregulation

## Usage
### Perform differential gene expression analysis in HGGs and compute correlations of select splicing factors:
<br>**Run shell script to make tables and subsequent plots below**
```
./run_module.sh
```

*Input files:*
```
tpm_norm_vs_tumor.h3k28.SFs.v2.txt
tpm_norm_vs_tumor.hgg.SFs.txt
RSF11.rna_vs_prot.cptac.txt
SRSF11_psi.txt
SRSF11_splicing_vs_expr.txt
SRSF11rna_vs_phos.txt
RBM5rna_vs_phos.txt
```

## Folder content
* `run_module.sh` runs R scripts to perfom differential expression analysis and runs correlation scripts
* `01-volcano_plot_mRNA.R` Runs DESeq2 on count data from `input/tpm_norm_vs_tumor*` and generates volcano plots
* `02-SRSF11_plots.R` generates stacked barplot of PSI values from `SRSF11_psi.txt` and computes correlation and generates scatter plot of expr and PSI
* `03-plot_corr_SRSF11_RBM5.R` computes correlation and generates scatter plots of select gene expression, proteomics and phosphoprotein obtained from CPTAC and input files `input/*vs_phos.txt`
