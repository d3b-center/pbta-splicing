# Splicing Index

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to compute the splicing index of each tumor (proportion of mis-spliced events)

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```
Input files (`data` folder, or TMB in `input` folder):
```
histologies.tsv
splice-events-rmats.tsv.gz
independent-specimens.rnaseqpanel.primary-plus.tsv
independent-specimens.rnaseqpanel.primary-plus.tsv
snv-mutation-tmb-coding.tsv
```

## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `01-generate_splicing-index_and_diff-events_table.pl` processes rMATs output and computes splicing burden index per sample stratified by splicing case
* `02-plot_splicing_burden_index.R` takes splicing index burden values and generates CDF plot categorized by histologies and denotes median values
* `03-identify_and_plot_histologies_by_SBI.R` identifies high vs low SBI tumors categorized by histologies
* `04-plot_total-splicing-cases.R` plots total splicing cases across samples  
* `05-plot-tmb-vs-sbi.R` plots splicing burden vs TMB by mutation status and cancer group 

## Directory structure
```.
├── 01-generate_splicing-index_and_diff-events_table.pl
├── 02-plot_splicing_burden_index.R
├── 03-identify_and_plot_histologies_by_SBI.R
├── 04-plot_total-splicing-cases.R
├── 05-plot-tmb-vs-sbi.R
├── README.md
├── input
│   └── snv-mutation-tmb-coding.tsv
├── plots
├── boxplot_sbi-tmb-by-cg.pdf
├── boxplot_sbi-tmb-by-mutation-status.pdf
├── corplot_sbi-tmb-by-cg.pdf
├── corplot_sbi-tmb.pdf
├── hist_by_sbi_level_barplot.pdf
├── sbi-plot-A3SS.pdf
├── sbi-plot-A5SS.pdf
├── sbi-plot-RI.pdf
├── sbi-plot-SE.pdf
└── splice-types.pdf
├── results
│   ├── histologies-plot-group.tsv
│   ├── splice_events.diff.A3SS.txt
│   ├── splice_events.diff.A5SS.txt
│   ├── splice_events.diff.RI.txt
│   ├── splice_events.diff.SE.txt
│   ├── splicing_index.A3SS.txt
│   ├── splicing_index.A5SS.txt
│   ├── splicing_index.RI.txt
│   └── splicing_index.SE.txt
└── run_module.sh
```
