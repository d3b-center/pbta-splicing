# Splicing Overlapping Functional Sites

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify splicing events that result in loss/gain of functional sites (as defined by Uniprot) amongst identified differential splicing events

## Usage
### Make summary table of strong splicing events and relevant filtered tables that overlap functional sites:
<br>**Run shell script to make tables and subsequent plots below**
```
./run_module.sh
```

*Input files:*
```
data/histologies.tsv
data/splice-events-rmats.tsv.gz
```

*input bed files obtained from unipro needed to  run*
```
unipDisulfBond.hg38col.bed
unipDomain.hg38.col.bed
unipLocSignal.hg38.col.bed
unipMod.hg38.col.bed
unipOther.hg38.col.bed

```

*Output files:*
```
results/splice_events.diff.SE.HGG.txt
results/splicing_events.SE.total.HGG.pos.bed
results/splicing_events.SE.total.HGG.neg.bed
results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt
results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt

```

## Folder content
* `run_module.sh` takes the files from above and generates table with uniprot overlaps to be used for plotting
* `01-extract_recurrent_splicing_events_hgg.pl` processing output from rMATS with filters and constructs data table for all downstream analysis and output file to `results/splicing_events.total*`
* `02-run_bedtools_intersect.sh` runs bedtools intersect to find exon coordinates corresponding to Uniprot sites
* `03-format_for_ggplot.sh` formats and appends file into table for plotting
* `04-plot_splicing_across_functional_sites.R` generates ggplot violin plots of average dPSI per event identidied overlapping a functional site, outputting to `plots/*png`
* `05-plot-splice-patterns` generates plots for visualizing splicing event types into `plots` folder
* `06-generate-functional-site-tabs.pl` parses and formats result files from previous scripts into  table (tsv) for splice events hitting functional sites

## Directory structure
├── 01-extract_recurrent_splicing_events_hgg.pl
├── 02-run_bedtools_intersect.sh
├── 03-format_for_ggplot.sh
├── 04-plot_splicing_across_functional_sites.R
├── 05-plot-splice-patterns.R
├── 06-generate-functional-site-tabs.pl
├── README.md
├── input
│   ├── gene_lists.tsv
│   ├── unipDisulfBond.hg38col.bed
│   ├── unipDomain.hg38.col.bed
│   ├── unipLocSignal.hg38.col.bed
│   ├── unipMod.hg38.col.bed
│   └── unipOther.hg38.col.bed
├── plots
│   ├── dPSI_across_functional_sites.HGG.pdf
│   ├── dPSI_across_functional_sites_kinase.HGG.pdf
│   ├── kinases-ora-plot.pdf
│   └── splicing_pattern_plot.pdf
├── results
│   ├── kinases-functional_sites_included.txt
│   ├── kinases-functional_sites_skipped.txt
│   ├── splice-functional-sites-table.tsv
│   ├── splice_events.diff.SE.HGG.txt
│   ├── splicing_events.total.HGG.neg.intersectUnip.ggplot.txt
│   └── splicing_events.total.HGG.pos.intersectUnip.ggplot.txt
├── run_module.sh
