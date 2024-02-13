# Splicing Factor Dysregulation

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to identify known splicing factors [(ensemble-annotated)](https://genome.cshlp.org/content/26/6/732.long#ref-15) that are dysregulated in high splicing burden HGGs, potentially explaining the global levels of splicing dysregulation across HGGs.

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `01-plot-diffExp_highlowSBI.R` performs differential gene expression on high vs low SBI tumors
* `02-plot_SFs_rna_vs_prot.R` generates heatmap of dysregulated splicing factors showing RNA and protein levels across brain tumors  

## Directory structure
```
.
├── 01-plot-diffExp_highlowSBI.R
├── 02-plot_SFs_rna_vs_prot.R
├── README.md
├── input
│   ├── CPTAC3-pbt_SF_family.xls
│   └── splicing_factors.txt
├── plots
│   ├── SF_RNA_vs_protein_levels_heatmap.pdf
│   ├── barplot_hgg_SFs.pdf
│   └── enhancedVolcano_hgg_sbi.pdf
├── results
│   └── diffSFs_sig_genes.txt
└── run_module.sh
```
