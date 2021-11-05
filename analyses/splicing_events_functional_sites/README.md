# Splicing Overlapping Functional Sites

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify splicing events that result in loss/gain of functional sites (as-defined by Uniprot)

## Usage
### Make summary table of strong splicing events and relevant filtered tables that overlap functional sites:
<br>**Run shell script to make final tables to be used for plotting below**
```
./run_module.sh
```

### Plot filtered tables that overlap functional sites using ggplot:
<br>**Run R script for plotting**
```
Rscript splicing_functional_sites.R
```

Input files:
```
data/v19_plus_20210311_pnoc.tsv
data/rMATS_merged.comparison.tsv
tables (*bed and *tsv) generated from `extract_recurrent_splicing_events.pl` run
unipDisulfBond.hg38col.bed
unipDomain.hg38.col.bed
unipLocSignal.hg38.col.bed
unipMod.hg38.col.bed
unipOther.hg38.col.bed
```

Output files:
```
results/splicing_events.total.neg.bed
results/splicing_events.total.pos.bed
results/splicing_events.total.neg.tsv
results/splicing_events.total.pos.tsv
results/splicing_events.total.neg.intersect*wo.txt
results/splicing_events.total.pos.intersect*wo.txt
```

![](plots/dPSI_distr_across_sites_positive_rec2.png)
<br>
![](plots/dPSI_distr_across_sites_negative_rec2.png)


## Folder content
* `run_module.sh` takes the files from above and generates table with uniprot overlaps to be used for plotting
* `extract_recurrent_splicing_events_hgg.pl` processing output from rMATS with filters and constructs data table for all downstream analysis and output file to `results/splicing_events.total*`
* `splicing_functional_sites.R` generates ggplot violin plots of average dPSI per event identidied overlapping a functional site, outputting to `plots/*png`
