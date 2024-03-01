# Write manuscript tables

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is write and generate suppl tables for manuscript

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
Rscript --vanilla make-suppl-tables.R
```

## Folder content
* `make-suppl-tables.R` script to generate suppl tables in xls format for manuscript

## Directory structure
```
.
├── input
│   └── CNS_primary_site_match.json
├── make-suppl-tables.R
└── output
    ├── TableS1-histologies.xlsx
    ├── TableS2-histology-specific-splice-events.xlsx
    ├── TableS3-DeSeq2-sbi-SFs.xlsx
    ├── TableS4-functional-sites.xlsx
    └── TableS5-CLK1-ex4-splicing-impact-morpholino.xlsx
```
