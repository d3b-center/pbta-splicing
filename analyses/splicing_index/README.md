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
* `02-plot_splicing_burden_index.R` takes splicing index burden values and generates CDF plot

## Directory structure
```
.
├── 01-generate_splicing-index_and_diff-events_table.pl
├── 02-plot_splicing_burden_index.R
├── README.md
├── plots
│   ├── sbi-plot-A3SS.tiff
│   ├── sbi-plot-A5SS.tiff
│   ├── sbi-plot-RI.tiff
│   └── sbi-plot-SE.tiff
├── results
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
