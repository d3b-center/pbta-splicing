# Splicing-based clustering

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to clustering based on PSI values

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
./run_module.sh
```
Input files:
```
data/pbta-histologies.RNA-Seq
```

Output files:
```
results/CC_groups_remDup.txt
```

![](plots/dPSI_distr_across_sites_positive_rec2.png)
<br>
![](plots/dPSI_distr_across_sites_negative_rec2.png)


## Folder content
* `run_module.sh` shell script to pre-process histology file and run analysis
* `consensus_clustering.R` takes PSI matrix and impliments consensus clustering, plots heatmap and stacked barplot of cluster members, outputting to `plots/*png`
