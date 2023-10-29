# CLK1 Splicing Correlations

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to correlate CLK1 Exon 4 splicing with splicing
burden, RNA expression and proteomics

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```
Input files:
```
./data/histologies.tsv
./data/splice-events-rmats.tsv.gz
./data/gene-counts-rsem-expected_count-collapsed.rds
./analyses/CLK1-splicing_correlations/input/CPTACT_SR_CLK1.xls
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-plot_highExon4_vs_lowExon4_and_SBI.R` correlate high vs low levels of overall splicing burden with CLK1 exon 4 inclusion levels
* `02-plot_splicing_vs_expr.R` correlate CLK1 exon 4 inclusion levels with RNA expression, including CLK1 and SRSF1
* `03-plot_CLK1-Ex4-splicing_vs_SRSF1-expr.R` correlate CLK1 RNA expression with RNA SRSF phospho proteomics across brain tumor types

## Directory structure
```
.
├── 01-plot_highExon4_vs_lowExon4_and_SBI.R
├── 02-plot_splicing_vs_expr.R
├── 03-plot_CLK1-Ex4-splicing_vs_SRSF1-expr.R
├── README.md
├── plots
│   ├── boxplot_high_vs_low_SBI.tiff
│   ├── heatmap_SR-phosp_CLK1-RNA_rna_prot.pdf
│   ├── CLK1-expr_vs_psi-midlineHGG.tiff
│   ├── CLK1-expr_vs_psi-totalHGG.tiff
│   ├── SRSF1-expr_vs_psi-midlineHGG.tiff
│   └── SRSF1-expr_vs_psi-totalHGG.tiff
└── run_module.sh
```
