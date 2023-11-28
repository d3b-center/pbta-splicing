# Total Protein Abundance Distributions

Module authors: Joseph Dybas (@JosephDybas)

The purpose of this module is to compare the total protein distributions for each of the samples in the HOPE cohort.

## Usage
### Scripts and files to run analysis
<br>**Run R scripts to make plots below**

Input files (`./input` directory)
```
splicing_events-psi.tsv: File containing differential splicing events and Psi values for each sample
Hope_proteome_imputed_data_liftover.tsv: File containing total proteome abundance from Hope and GBM cohorts
Hope_phosphosite_imputed_data_ischemia_removed_liftover.tsv: File containing phospho protein abundance from Hope and GBM cohorts
Hope-GBM-histologies-base.tsv: Histologies file with relevant IDs
```
## Folder content
* `WCPproteomics-abundance-volcanoplots.R` Creates violin plots of the distributions of protein abundances in each sample with abundance of specific protein highlighted

### Directory structure
```
.
├── WCPproteomics-abundance-volcanoplots.R
├── RNAsplicing-PHOSproteomics-scatterplots.R
├── README.md
├── input
│   ├── splicing_events-psi.tsv
│   ├── Hope_proteome_imputed_data_liftover.tsv
│   ├── Hope_phosphosite_imputed_data_ischemia_removed_liftover.tsv
│   └── Hope-GBM-histologies-base.tsv
├── plots
│   ├── SampleAbundanceDist_CLK1_violin.pdf
│   ├── SampleAbundanceDist_YAP1_violin.pdf
│   ├── SampleAbundanceDist_CHRD_violin.pdf
│   ├── SampleAbundanceDist_ATRX_violin.pdf
```

