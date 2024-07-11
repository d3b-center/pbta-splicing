# Oncoprint

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess CLK1 exon 4 PSI levels in the context of mutations

## Usage
### How to Run:
```
<<<<<<< HEAD
Rscript --vanilla 01-oncoprint.R
=======
bash run-oncoprint.sh
>>>>>>> main
```

#### Input:
* MAF - T/N: ```./data/snv-consensus-plus-hotspots.maf.tsv.gz``` <br>
* MAF - T only: ```./data/snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz``` <br>
* clinical info: ```./data/histologies.tsv```  <br>
* genes of interest: ```input/oncoprint-goi-lists-OpenPedCan-gencode-v39.csv``` <br>
* independent specimen file (RNA): ```./data/independent-specimens.rnaseqpanel.primary.tsv```  <br>
* mutation colors: ```input/mutation-colors.R``` <br>
* splicing events: ```./data/splice-events-rmats.tsv.gz```  <br>
* tumor mutation burden: ```./input/snv-mutation-tmb-coding.tsv``` <br>

#### Output:
```plots/oncoprint.pdf``` <br>

## Scripts
* `01-oncoprint.R` generates oncoprint with mutation frequencies with CLK1 exon 4 PSI, gender, molecular subtype, CNS region and mutation status information across pediatric HGGs, as well as enrichment of CLK1 high/low tumors by gene alteration 
<<<<<<< HEAD
=======
* `02-oncoprint-SFs.R` generates oncoprint with SF mutation frequencies with CLK1 exon 4 PSI, gender, molecular subtype, CNS region and mutation status information across pediatric HGGs, as well as enrichment of CLK1 high/low tumors by gene alteration 
>>>>>>> main

## Directory Structure
```
.
├── 01-oncoprint.R
├── README.md
├── input
│   ├── mutation-colors.R
│   ├── oncoprint-goi-lists-OpenPedCan-gencode-v39.csv
│   └── snv-mutation-tmb-coding.tsv
├── plots
│   └── oncoprint.pdf
└── results
    └── clk1_high_low_mutation_counts.tsv
```
