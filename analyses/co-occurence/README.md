# Co-occurrence 

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify co-occurring mutations in differentially spliced genes in HGG tumors.

## Usage
### Make summary table of strong splicing events and relevant filtered tables that overlap functional sites:
<br>**Run shell script to make tables and subsequent plots below**
```
Rscript --vanilla 01-co-occurence-interactions.R
```

*Input files:*
maf: ```./data/snv-consensus-plus-hotspots.maf.tsv.gz```
clinical info: ```./data/histologies.tsv```)
independent specimen file (RNA): ```./data/independent-specimens.rnaseqpanel.primary.tsv```
independent specimen file (DNA): ```./data/independent-specimens.wgswxspanel.primary.prefer.wgs.tsv```
differential splicing events: ```./analyses/splicing_events_functional_sites/results/splice_events.diff.SE.HGG.txt```

*Output files:*
```results/co-occurence-sign.tsv```

## Folder content
* `01-co-occurence-interactions.R` identifies genes that have signficant co-occuring variants in differentially spliced genes in HGG tumors
