# long-read-CLK1-validation

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to visualize long (Oxford Nanopore Technologies) 
along with illumina short RNA-seq reads. This module helps and shows
full-length transcript (with and without Exon 4) validations of CLK1.

## Usage
### script to run analysis
<br>**Run R script to make stacked barplots of inclsuion/skipping isoforms**
```
Rscript --vanilla 01-plot_ont_vs_short_CLK1-Ex4_psi.R
```
Input files:
```
ont_vs_rmats.ggplot.tsv
```

## Folder content
* `01-plot_ont_vs_short_CLK1-Ex4_psi.R` plot CLK1 Exon 4 included/skipped isoforms categorized by RNA-seq sequencing strategy (long vs short sequencing)

## Directory structure
```
.
├── 01-plot_ont_vs_short_CLK1-Ex4_psi.R
├── input
│   └── ont_vs_rmats.ggplot.tsv
└── plots
    └── isoform_stackedbarplot.tiff
```
