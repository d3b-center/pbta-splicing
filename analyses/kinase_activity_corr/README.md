## Clustered heatmap generation for kinase activity scores

Module author: Run Jin ([@runjin326](https://github.com/runjin326) 

This module contains two scrpits:
1) R notebook `kinase_activity_corr.Rmd` that explores different heatmap options using CLK1 PSI scores and either HGAT specific or pan cancer kinase activity scores
2) R script `kinase_activity_per_gene.R` that iterates through individual PSI scores of relevant RBP proteins and output heatmap and results

### Usage
The module does not have a run-script.sh but each script can be run individually as followed:
  
```
Rscript kinase_activity_corr.R
```

### Input files 
For `kinase_activity_corr.Rmd`: 
  - Input files (in the `input` folder of this analysis module) are: 
    ```
    supplementarytable_kinase.xlsx
    CLK1_PSI_grouping.txt
    CLK1_PSI.pan-cancer.groupings.txt
    ```
  - Data files (in the `data` folder of this repo): `pbta-histologies.tsv`

For `kinase_activity_per_gene.R`:
 - Input files (in the `input` folder of this analysis module) are: 
 1) files in `input/individual_files`
 2) `input/supplementarytable_kinase.xlsx`
 - Data files (in the `data` folder of this repo): `pbta-histologies.tsv`

### Output files 
For `kinase_activity_corr.Rmd`: 
  - Output file (in the `plots` folder of this analysis module):
    ```
    kinase_group_psi_hist_pan_all.pdf
    kinase_group_psi_hist.pdf
    kinase_group_psi_hist_pan_hist.pdf
    kinase_group_psi_hist_pan.pdf
    kinase_group_psi_hist_pan_sub.pdf
    ```
For `kinase_activity_per_gene.R`:
1) Heatmap output in `plots/individual_genes/sigificant_heatmap`
2) Correlation plot output in `plots/individual_genes/corrplot`
3) Correlation stats output in `results` 
All the files are named with the individual RBP gene of interest.
