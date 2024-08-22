# Cohort summary

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to summarize cohort in our downstream analyses.

## Usage
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```
Input files (`data` folder):
```
histologies.tsv
independent-specimens.rnaseqpanel.primary.tsv
```

## Folder content
* `01-generate-cohort-summary-circos-plot.R` creates circos plot of cohort highlighting histology, CNS region and reported gender

## Directory structure
```
.
├── 01-generate-cohort-summary-circos-plot.R
├── input
│   └── plot-mapping.tsv
├── plots
│   └── cohort_circos.pdf
└── results
    ├── histologies-plot-group.tsv
    └── plot_mapping.tsv
```
