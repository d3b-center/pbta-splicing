# long-read-CLK1-validation

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to visualize long (Oxford Nanopore Technologies)
along with illumina short RNA-seq reads. This module helps and shows
full-length transcript (with and without Exon 4) validations of CLK1.

## Usage
### script to run analysis
<br>**Run bash script to make stacked barplots of inclsuion/skipping isoforms**
```
bash run_module.sh
```
Input files:
```
7316_1763.CLK1.aln
7316_1769.CLK1.aln
KNS42.CLK1.aln
```

## Folder content
* `00-pre-process-aln.sh` processes and re-formats *aln files for plotting script
* `01-plot_ont_vs_short_CLK1-Ex4_psi.R` plot CLK1 Exon 4 included/skipped isoforms categorized by RNA-seq sequencing strategy (long vs short sequencing)

## Directory structure
```
.
├── 00-pre-process-aln.sh
├── 01-plot_ont_vs_short_CLK1-Ex4_psi.R
├── README.md
├── input
│   ├── 7316_1763.CLK1.aln
│   ├── 7316_1763.CLK1.processed.txt
│   ├── 7316_1769.CLK1.aln
│   ├── 7316_1769.CLK1.processed.txt
│   ├── KNS42.CLK1.aln
│   ├── KNS42.CLK1.processed.txt
│   └── ont_vs_rmats.ggplot.tsv
├── plots
│   └── isoform-stacked-barplot.pdf
└── run_module.sh
```
