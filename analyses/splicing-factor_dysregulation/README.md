# Splicing Factor Dysregulation

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to identify known splicing factors [(ensemble-annotated)](https://genome.cshlp.org/content/26/6/732.long#ref-15) that are dysregulated in select HGGs, potentially explaining the global levels of splicing dysregulation across HGGs.

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `01-volcano_plot_mRNA.R` performs differential gene expression on splicing factors using DESeq2 and generates a volcano plot
* `02-identify_preprocess_SRSF11_psi.sh` formats and parses SRSF11 PSI values to use for plotting
* `03-SRSF11_plots.R` plots SRSF11 PSI values from above fo each sample and control
* `04-plot_SFs_rna_vs_prot.R` generates heatmap of dysregulated splicing factors showing RNA and protein levels across brain tumors  

## Directory structure
```
.
├── 01-volcano_plot_mRNA.R
├── 02-identify_preprocess_SRSF11_psi.sh
├── 03-SRSF11_plots.R
├── 04-plot_SFs_rna_vs_prot.R
├── README.md
├── input
│   ├── CPTAC3-pbt.xls
│   ├── rsem_counts.non_tumor.tsv
│   ├── se.srsf11.ctrl.incl.txt
│   ├── se.srsf11.incl.ctrl.tmp.txt
│   ├── se.srsf11.incl.tmp.txt
│   ├── se.srsf11.incl.txt
│   └── splicing_factors.txt
├── plots
│   ├── SRSF11_hgg_stacked.pdf
│   ├── enhancedVolcano_hggs_v_ctrl_SFs.pdf
│   └── heatmap_SFs_rna_prot.pdf
├── results
│   └── hggs_v_ctrl_SFs_sig_genes.txt
└── run_module.sh
```
