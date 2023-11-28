# Proteomics_correlation

Module authors: Joseph Dybas (@JosephDybas)

The purpose of this module is to compare the differential splicing in HGGs total and phospho proteomics abundances.

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
* `WCPproteomics-abundance-volcanoplots.R` creates volcano plots of the distributions of protein abundances in identified samples with abundance of specific gene highlighted
* `RNAsplicing-PHOSproteomics-scatterplots.R` Creates scatter plot of differential splicing Psi value and corresponding phosphorylation abundance for specified gene, splicing event, and phosphorylation event

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
│   ├── PsiPhosAbd_CLK1_scatter.pdf
```

