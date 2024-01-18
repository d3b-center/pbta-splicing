# Total Protein Abundance Distributions

Module authors: Joseph Dybas (@JosephDybas)

The purpose of this module is to compare the total protein distributions for each of the samples in the HOPE cohort.

## Usage
### Scripts and files to run analysis
<br>**Run R scripts to make plots below**

Input files (`./input` directory)
```
splicing_events-psi.tsv: File containing differential splicing events and Psi values for each sample
hope-protein-imputed-prot-expression-abundance.tsv.gz: File containing total proteome abundance from HOPE cohort
gbm-protein-imputed-prot-expression-abundance.tsv.gz: File containing total proteome abundance from GBM cohort
```
## Folder content
* `WCPproteomics-abundance-volcanoplots.R` Creates violin plots of the distributions of protein abundances in each sample with abundance of specific protein highlighted

### Directory structure
```
.
├── WCPproteomics-abundance-violinplots.R
├── README.md
├── input
│   ├── splicing_events-psi.tsv
├── plots
│   ├── SampleAbundanceDist_CLK1_violin.pdf
│   ├── SampleAbundanceDist_YAP1_violin.pdf
│   ├── SampleAbundanceDist_CHRD_violin.pdf
│   ├── SampleAbundanceDist_ATRX_violin.pdf
```

