# KNS42 Cell Line
This module investigates the dependency of CLK1 in pediatric HGG cell line KNS42 relative to all other brain tumor cell lines using DepMap portal data, and generates visualizations of KNS42 morpholino treated cells and functional validations (experimental) of CLK1 exon 4 splicing on cell proliferation.  

## Usage
To run analysis module from the command line:
```
bash run_module.sh
```

### Input
```
input/CLK1_CRISPR_depmap_score.csv
input/OmicsDefaultModelProfiles.csv
input/cell_prolif_res.tsv
input/qpcr-raw-triplicates.tsv
```

## Folder Content
```01-plot-and-process-depmap-data.R``` is script plotting CLK1 dependency scores and CLK1 transcript expression for brain tumor cell lines available in the DepMap Portal<br>
```02-plot_cell-proliferation-assay-res.R``` plots the results from the cell proliferation assay of treated vs ctrl/untreated cells<br>
```03-plot-qPCR-results.R``` visualizes qRT-PCR results

## Directory Structure
```
.
├── 01-plot-and-process-depmap-data.R
├── 02-plot_cell-proliferation-assay-res.R
├── 03-plot-qPCR-results.R
├── README.md
├── input
│   ├── 2023-03-22 162604_JLRmod.xls
│   ├── CLK1_CRISPR_depmap_score.csv
│   ├── OmicsDefaultModelProfiles.csv
│   ├── cell_prolif_res.tsv
│   └── qpcr-raw-triplicates.tsv
├── plots
│   ├── cell_prolif-line.pdf
│   ├── depmap_score_CLK1_vs_score_KNS42.pdf
│   ├── depmap_score_all_cell_lines.pdf
│   ├── depmap_score_cns_cell_lines.pdf
│   └── qPCR-morp.pdf
└── run_module.sh
```
