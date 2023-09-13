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

## Directory structure
```
.
├── 01-generate_hist_spec_events_tab.pl
├── 02-plot_histology-specific_splicing_events.R
├── plots
│   ├── upsetR_histology-specific.ei.pdf
│   ├── upsetR_histology-specific.ei.tiff
│   ├── upsetR_histology-specific.es.pdf
│   └── upsetR_histology-specific.es.tiff
├── results
│   ├── histology-specific_events.en.total.tsv
│   ├── histology-specific_events.es.total.tsv
│   └── splicing_events.hist-labeled_list.thr2freq.txt
└── run_module.sh
```
