# Splicing Overlapping Functional Sites

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify splicing events that result in loss/gain of functional sites (as defined by Uniprot)

## Usage
### Make summary table of strong splicing events and relevant filtered tables that overlap functional sites:
<br>**Run shell script to make tables and subsequent plots below**
```
./run_module.sh
```

*Input files:*
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

*Output files:*
```
results/splicing_events.total.pos.intersectUnip.ggplot.txt
results/splicing_events.total.neg.intersectUnip.ggplot.txt
```

## Folder content
* `run_module.sh` takes the files from above and generates table with uniprot overlaps to be used for plotting
* `01-extract_recurrent_splicing_events_hgg.pl` processing output from rMATS with filters and constructs data table for all downstream analysis and output file to `results/splicing_events.total*`
* `01-run_bedtools_intersect.sh` runs bedtools intersect to find exon coordinates corresponding to Uniprot sites
* `02-format_for_ggplot.sh` formats and appends file into table for plotting
* `01-plot_splicing_across_functional_sites.R` generates ggplot violin plots of average dPSI per event identidied overlapping a functional site, outputting to `plots/*png`
* `02-plot-flip_mixed_events.R` generates plots for special flip splicing events into `plots` folder
