# Splicing Overlapping Functional Sites

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify splicing events that result in loss/gain of functional sites (as defined by Uniprot)

## Usage
<br>**Run shell script to make tables and subsequent plots below**
```
./run_module.sh
```

## Input
```
data/histologies.tsv
data/rMATS_merged.comparison.tsv.gz
```

*input bed files obtained from unipro needed to  run*
```
unipDisulfBond.hg38col.bed
unipDomain.hg38.col.bed
unipLocSignal.hg38.col.bed
unipMod.hg38.col.bed
unipOther.hg38.col.bed
```

## Output
```
results/splicing_events.total.pos.intersectUnip.ggplot.txt
results/splicing_events.total.neg.intersectUnip.ggplot.txt
```

## Folder content
* `run_module.sh` takes the files from above and generates table with uniprot overlaps to be used for plotting
* `01-extract_recurrent_splicing_events_hgg.pl` processing output from rMATS with filters and constructs data table for all downstream analysis and output file to `results/splicing_events.total*`
* `02-run_bedtools_intersect.sh` runs bedtools intersect to find exon coordinates corresponding to Uniprot sites
* `03-format_for_ggplot.sh` formats and appends file into table for plotting
* `04-plot_splicing_across_functional_sites.R` generates ggplot violin plots of average dPSI per event identidied overlapping a functional site, outputting to `plots/*png`
* `05-plot-splice-patterns.R` generates plots of patterns observed categorized by functional sites

## Directory structure
```
.
├── 01-extract_recurrent_splicing_events_hgg.pl
├── 02-run_bedtools_intersect.sh
├── 03-format_for_ggplot.sh
├── 04-plot_splicing_across_functional_sites.R
├── 05-plot-splice-patterns.R
├── README.md
├── input
│   ├── unipDisulfBond.hg38col.bed
│   ├── unipDomain.hg38.col.bed
│   ├── unipLocSignal.hg38.col.bed
│   ├── unipMod.hg38.col.bed
│   └── unipOther.hg38.col.bed
├── plots
│   ├── dPSI_across_functional_sites_neg.pdf
│   ├── dPSI_across_functional_sites_pos.pdf
│   ├── flip_barplots.pdf
│   └── splicing_pattern_plot.pdf
├── results
│   ├── archive
│   ├── flip.tmp
│   ├── splicing_events.total.neg.bed
│   ├── splicing_events.total.neg.intersectUnip.ggplot.txt
│   ├── splicing_events.total.neg.tsv
│   ├── splicing_events.total.pos.bed
│   ├── splicing_events.total.pos.intersectUnip.ggplot.txt
│   └── splicing_events.total.pos.tsv
└── run_module.sh
```
