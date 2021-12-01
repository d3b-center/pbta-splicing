## Clustered heatmap generation for kinase activity scores

This module contains R script that generate a heatmap of kinase activity scores (as published in @doi:10.1016/j.cell.2020.10.044) clustered by PSI scores, exon-4 groupings and molecular subtype

### Usage
The module does not have a run-script.sh but each script can be run individually as followed:
  
```
Rscript kinase_activity_corr.R
```

Input files (in the `input` folder of this analysis module):
```
supplementarytable_kinase.xlsx
CLK1_PSI_grouping.txt
```

Data files (in the `data` folder of this repo)
```
pbta-histologies.tsv
```

Output file (in the `plots` folder of this analysis module):
```
kinase_group_psi_hist.pdf
```

