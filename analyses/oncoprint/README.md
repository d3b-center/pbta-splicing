# Oncoprint

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess CLK1 exon 4 PSI levels in the context of 

## Usage
### How to Run:
```
Rscript --vanilla 01-co-occurence-interactions.R
```

#### Input:
MAF: ```./data/snv-consensus-plus-hotspots.maf.tsv.gz``` <br>
clinical info: ```./data/histologies.tsv```  <br>
genes of interest: ```input/oncoprint-goi-lists-OpenPedCan-gencode-v39.csv``` <br>
independent specimen file (RNA): ```./data/independent-specimens.rnaseqpanel.primary.tsv```  <br>
independent specimen file (DNA): ```./data/independent-specimens.wgswxspanel.primary.prefer.wgs.tsv```  <br>
mutation colors: ```input/mutation-colors.R``` <br>
splicing events: ```./data/splice-events-rmats.tsv.gz```  <br>
tumor only maf: ```./data/snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz```  <br>
tumor mutation burden `<- ```./data/snv-mutation-tmb-coding.tsv``` <br>

#### Output:
```plots/oncoprint.pdf``` <br>

## Scripts
* `01-oncoprint.R` generates oncoprint with mutation frequencies with CLK1 exon 4 PSI, gender, molecular subtype, CNS region and mutation status information across pediatric HGGs 

## Directory Structure
```
.
├── 01-oncoprint.R
├── README.md
├── input
│   ├── mutation-colors.R
│   └── oncoprint-goi-lists-OpenPedCan-gencode-v39.csv
├── plots
│   └── oncoprint.pdf
├── results
└── util
```
