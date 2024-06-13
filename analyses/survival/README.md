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
* `06-survival_by_cluster.Rmd` Assess survival by splicing cluster assingment in all histologies, LGG, and HGG
* `07-run-survival-nf1-psi.Rmd` Generate survival models in histologies and subtypes by NF1-215 PSI
* `08-plot-survival-nf1-psi.Rmd` Plot survival models by NF1-215 PSI


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
├── 06-survival_by_cluster.Rmd
├── 06-survival_by_cluster.nb.html
├── 07-run-survival-nf1-psi.Rmd
├── 07-run-survival-nf1-psi.nb.html
├── 08-plot-survival-nf1-psi.Rmd
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
│   ├── forest_add_EFS_HGG_subtype_cluster_assignment.pdf
│   ├── forest_add_EFS_LGG_resection_subtype_cluster_assignment.pdf
│   ├── forest_add_EFS_resection_lgg_group_cluster_assignment.pdf
│   ├── forest_add_OS_HGG_subtype_cluster_assignment.pdf
│   ├── forest_add_OS_LGG_resection_subtype_cluster_assignment.pdf
│   ├── forest_add_OS_resection_lgg_group_cluster_assignment.pdf
│   ├── km_DMG_OS_EFS_CLK1_status.pdf
│   ├── km_EFS_LGG_cluster_assignment.pdf
│   ├── km_EFS_cluster_assignment.pdf
│   ├── km_HGG_DMG_OS_EFS_CLK1_status.pdf
│   ├── km_HGG_OS_EFS_CLK1_status.pdf
│   ├── km_OS_cluster_assignment.pdf
│   ├── km_hgg_EFS_cluster_assignment.pdf
│   ├── km_hgg_OS_cluster_assignment.pdf
│   ├── km_lgg_EFS_cluster_assignment.pdf
│   └── km_lgg_OS_cluster_assignment.pdf
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
│   ├── cox_EFS_additive_terms_resection_lgg_group_cluster.RDS
│   ├── cox_HGG_EFS_additive_terms_subtype_clk1_status.RDS
│   ├── cox_HGG_EFS_interaction_terms_subtype_clk1_status.RDS
│   ├── cox_HGG_OS_additive_terms_subtype_clk1_status.RDS
│   ├── cox_HGG_OS_interaction_terms_subtype_clk1_status.RDS
│   ├── cox_OS_additive_terms_resection_lgg_group_cluster.RDS
│   ├── cox_hgg_EFS_additive_terms_subtype_cluster.RDS
│   ├── cox_hgg_OS_additive_terms_resection_subtype_cluster.RDS
│   ├── cox_hgg_OS_additive_terms_subtype_cluster.RDS
│   ├── cox_lgg_EFS_additive_terms_resection_subtype_cluster.RDS
│   ├── cox_lgg_OS_additive_terms_resection_subtype_cluster.RDS
│   ├── logrankEFS_cluster_assignment.RDS
│   ├── logrank_DMG_EFS_CLK1_status.RDS
│   ├── logrank_DMG_OS_CLK1_status.RDS
│   ├── logrank_HGG_EFS_CLK1_status.RDS
│   ├── logrank_HGG_OS_CLK1_status.RDS
│   ├── logrank_OS_cluster_assignment.RDS
│   ├── logrank_hgg_EFS_cluster_assignment.RDS
│   ├── logrank_hgg_OS_cluster_assignment.RDS
│   ├── logrank_lgg_EFS_cluster_assignment.RDS
│   ├── logrank_lgg_OS_cluster_assignment.RDS
│   ├── splicing_indices_with_survival.tsv
│   └── subtypes-for-survival.tsv
├── run-survival-module.sh
└── util
    └── survival_models.R
```
