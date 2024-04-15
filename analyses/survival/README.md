# Survival by splicing burden

Module authors: Ryan Corbett (@rjcorb), Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess survival in tumor histologies by overall splicing burden index and CLK1 exon 4 PSI

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run-survival-module.sh
```
Input files:
```
./data/histologies.tsv
./data/splice-events-rmats.tsv.gz
../splicing_index/results/splicing_index.SE.txt
```

## Folder content
* `run-survival-module.sh` shell script to run analysis
* `01-prepare-survival.Rmd` Merge splicing data with relevant histology data and define survival features
* `02-run-survival-SIgroup.Rmd` Assess histology-specific survival by splicing burden index (SBI) group (high vs low burden)
* `03-run-survival-SI.Rmd` Assess histology specific survival including SBI as a continuous variable
* `04-plot-survival.R` plot SBI survival models
* `05-survival-hgg-clk1-status.Rmd` Assess survival by CLK1 exon 4 PSI in HGG tumors


## Directory structure
```
├── 01-prepare-survival.Rmd
├── 01-prepare-survival.nb.html
├── 02-run-survival-SIgroup.Rmd
├── 02-run-survival-SIgroup.nb.html
├── 03-run-survival-SI.Rmd
├── 03-run-survival-SI.nb.html
├── 04-plot-survival.R
├── 05-survival-hgg-clk1-status.Rmd
├── 05-survival-hgg-clk1-status.nb.html
├── input
├── plots
│   ├── ATRT/
│   ├── CPG/
│   ├── DMG/
│   ├── EPN/
│   ├── GNG/
│   ├── HGG/
│   ├── LGG/
│   ├── MB/
│   ├── forest_DMG_EFS_add_subtype_clk1_status.pdf
│   ├── forest_DMG_OS_add_subtype_clk1_status.pdf
│   ├── forest_HGG_DMG_EFS_int_subtype_clk1_status.pdf
│   ├── forest_HGG_DMG_OS_add_subtype_clk1_status.pdf
│   ├── forest_HGG_DMG_OS_int_subtype_clk1_status.pdf
│   ├── forest_HGG_EFS_add_subtype_clk1_status.pdf
│   ├── forest_HGG_OS_add_subtype_clk1_status.pdf
│   ├── km_DMG_OS_EFS_CLK1_status.pdf
│   ├── km_HGG_DMG_OS_EFS_CLK1_status.pdf
│   └── km_HGG_OS_EFS_CLK1_status.pdf
├── results
│   ├── ATRT/
│   ├── CPG/
│   ├── DMG/
│   ├── EPN/
│   ├── GNG/
│   ├── HGG/
│   ├── LGG/
│   ├── MB/
│   ├── cox_DMG_EFS_additive_terms_subtype_clk1_status.RDS
│   ├── cox_DMG_OS_additive_terms_subtype_clk1_status.RDS
│   ├── cox_HGG_EFS_additive_terms_subtype_clk1_status.RDS
│   ├── cox_HGG_EFS_interaction_terms_subtype_clk1_status.RDS
│   ├── cox_HGG_OS_additive_terms_subtype_clk1_status.RDS
│   ├── cox_HGG_OS_interaction_terms_subtype_clk1_status.RDS
│   ├── logrank_DMG_EFS_CLK1_status.RDS
│   ├── logrank_DMG_OS_CLK1_status.RDS
│   ├── logrank_HGG_EFS_CLK1_status.RDS
│   ├── logrank_HGG_OS_CLK1_status.RDS
│   ├── splicing_indices_with_survival.tsv
│   └── subtypes-for-survival.tsv
├── run-survival-module.sh
└── util
    └── survival_models.R
```
