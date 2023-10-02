## Usage

To run analysis module from the command line:

```
Rscript 01-plot-and-process-depmap-data.R
```
## Summary
This module investigates the dependency of CLK1 in pediatric HGG cell line KNS42 relative to all other brain tumor cell lines. It produces two plots using DepMap portal showing 1. CLK1 dependency score highlighting KNS42 relative to all other brain tumor cell lines, and 2. Sina plot comparing the distributions of CRISPR scores betwenn cell lines with high vs low CLK1 Exon 4-containing transcript expression.

```01-plot-and-process-depmap-data.R``` is script plotting CLK1 dependency scores and CLK1 transcript expression for brain tumor cell lines available in the DepMap Portal

## Results
We confirm that CLK1 is essential to KNS42. Furthermore, high ex4 transcript expression is significantly more essential than low exon 4 transcript expression.

## Directory Structure
```
├── 01-plot-and-process-depmap-data.R
├── input
│   ├── CLK1_CRISPR_depmap_score.csv
│   └── OmicsDefaultModelProfiles.csv
└── plots
    ├── depmap_score_CLK1_vs_score_KNS42.tiff
    └── depmap_score_cell-lines.tiff
```
