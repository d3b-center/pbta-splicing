# Splicing index

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to compute the splicing index of each tumor (proportion of mis-spliced events)

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
./run_module.sh ../../data/v19_plus_20210311_pnoc_rna.tsv ../../data/merge_rMATS_splicing.SE.single.tsv
```
Input files:
```
../../data/v19_plus_20210311_pnoc_rna.tsv
../../data/merge_rMATS_splicing.SE.single.tsv
```
Output files:
```
results/splicing_index.wdPSI10_per_sample.txt
```

![](plots/SI_dpsi_thr10.png)
<br>


## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `splicing_index.R` takes SI matrix and plots CDF plots, outputting to `plots/*png`
