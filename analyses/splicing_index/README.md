# Splicing Index

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to compute the splicing index of each tumor (proportion of mis-spliced events)

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```
Input files (`data` folder):
```
histologies.tsv
splice-events-rmats.tsv.gz
independent-specimens.rnaseqpanel.primary.tsv
independent-specimens.rnaseqpanel.primary-plus.tsv
```

## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `01-generate_splicing-index_and_diff-events_table.pl` processes rMATs output and computes splicing burden index per sample stratified by splicing case
* `02-plot_splicing_burden_index.R` takes splicing index burden values and generates CDF plot categorized by histologies and denotes median values
* `03-identify_and_plot_histologies_by_SBI.R` identifies high vs low SBI tumors categorized by histologies
* `04-plot_total-splicing-cases.R` plots total splicing cases across samples  

## Directory structure
```
.
├── 01-generate_splicing-index_and_diff-events_table.pl
├── 02-plot_splicing_burden_index.R
├── 03-identify_and_plot_histologies_by_SBI.R
├── 04-plot_total-splicing-cases.R
├── README.md
├── plots
│   ├── hist_by_sbi_level_barplot.pdf
|   ├── sbi-plot-A3SS.pdf
│   ├── sbi-plot-A5SS.pdf
│   ├── sbi-plot-RI.pdf
│   ├── sbi-plot-SE.pdf
│   └── splice-types.pdf
├── results
│   ├── splicing_index.A3SS.txt
│   ├── splicing_index.A5SS.txt
│   ├── splicing_index.RI.txt
│   └── splicing_index.SE.txt
└── run_module.sh
```
