# Total Protein Abundance Distributions

Module authors: Joseph Dybas (@JosephDybas)

The purpose of this module is to compare the total and phospho protein data with the RNA differential splicing data

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
* `WCPproteomics-abundance-volcanoplots.R` Creates violin plots of the distributions of total protein abundances in each sample with abundance of specific protein highlighted
* `RNAsplicing-PHOSproteomics-scatterplots.R` Creates scatter plots comparing RNA splicing (Psi) or abundance values (TPM or EC) with phosphorylation event abundance for specified splicing event

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
│   ├── SampleAbundanceDist_CLK1_violin.pdf - violin plots of total protein abundance distributions with CLK1 highlighted
│   ├── SampleAbundanceDist_YAP1_violin.pdf - violin plots of total protein abundance distributions with YAP1 highlighted
│   ├── SampleAbundanceDist_CHRD_violin.pdf - violin plots of total protein abundance distributions with CCRD highlighted
│   ├── SampleAbundanceDist_ATRX_violin.pdf - violin plots of total protein abundance distributions with ATRX highlighted
│   ├── RnaPsiPhosAbundance_RnaTPMgt1_CLK1_scatter.pdf - scatterplot comparing CLK1 Psi and phosphorylation abundance for splicing event in exon 4 on S182
│   ├── RnaTPMAbundancePhosAbundance_RnaTPMgt1_CLK1_scatter.pdf - scatterplot comparing CLK1 TPM and phosphorylation abundance for splicing event in exon 4 on S182
│   ├── RnaPsiPhosAbundance_RnaTPMgt20_CLK1_scatter.pdf - scatterplot comparing CLK1 Psi and phosphorylation abundance for splicing event in exon 4 on S182, filtered to include only TPM >= 20
│   ├── RnaTPMAbundancePhosAbundance_RnaTPMgt20_CLK1_scatter.pdf - scatterplot comparing CLK1 TPM and phosphorylation abundance for splicing event in exon 4 on S182, filtered to include only TMP >= 20
```

