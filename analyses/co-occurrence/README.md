# Co-occurrence

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify co-occurring mutations in differentially spliced genes in HGG tumors.

## Usage
### How to Run:
```
Rscript --vanilla 01-co-occurence-interactions.R
```

#### Input:
maf: ```./data/snv-consensus-plus-hotspots.maf.tsv.gz``` <br>
clinical info: ```./data/histologies.tsv```)  <br>
independent specimen file (RNA): ```./data/independent-specimens.rnaseqpanel.primary.tsv```  <br>
independent specimen file (DNA): ```./data/independent-specimens.wgswxspanel.primary.prefer.wgs.tsv```  <br>
differential splicing events: ```./analyses/splicing_events_functional_sites/results/splice_events.diff.SE.HGG.txt```  <br>

#### Output:
```results/co-occurence-sign.tsv``` <br>

## Scripts
* `01-co-occurence-interactions.R` identifies genes that have signficant co-occuring variants in differentially spliced genes in HGG tumors

## Directory Structure
```
.
├── 01-co-occurrence-interactions.R
├── README.md
├── plots
└── results
    └── co-occurence-sign.tsv
```
