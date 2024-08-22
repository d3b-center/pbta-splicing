# Histology-specific splicing

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify splicing signatures or unique splicing variants within each histology

## Usage
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
* `01-generate_hist_spec_events_tab.pl` processes rMATs output and identifies unique splicing events if it is in 2% of the histology-specific cohort
* `02-plot_histology-specific_splicing_events.R` takes result table frmo above and generates UpSetR plots for each splicing case
* `03-plot-histology-specific-norm-events.R` computes and plots the average number of unique events per histology

## Directory structure
```
.
├── 01-generate_hist_spec_events_tab.pl
├── 02-plot_histology-specific_splicing_events.R
├── 03-plot-histology-specific-norm-events.R
├── README.md
├── plots
│   ├── avg-uniq-hits.pdf
│   ├── upsetR_histology-specific.ei.pdf
│   └── upsetR_histology-specific.es.pdf
├── results
│   ├── recurrent_splice_events_by_histology.tsv
│   ├── unique_events-ei.tsv
│   └── unique_events-es.tsv
└── run_module.sh
```
