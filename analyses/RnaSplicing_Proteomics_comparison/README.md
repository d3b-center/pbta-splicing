# Total Protein Abundance Distributions

Module authors: Joseph Dybas (@JosephDybas)

The purpose of this module is to compare the total and phospho protein data with the RNA differential splicing data

## Usage
### Scripts and files to run analysis
<br>**Run R scripts to make plots below**

Input files (`./input` directory)
```
splicing_events-psi.tsv: File containing differential splicing events and Psi values for each sample
Hope-GBM-histologies-base.tsv: Histologies file with relevant IDs
```
## Folder content
* `RNAsplicing-PHOSproteomics-scatterplots.R` Creates scatter plots comparing RNA splicing (Psi) or abundance values (TPM or EC) with phosphorylation event abundance for specified splicing event

### Directory structure
```
.
├── RNAsplicing-PHOSproteomics-scatterplots.R
├── README.md
├── input
│   ├── splicing_events-psi.tsv
│   └── Hope-GBM-histologies-base.tsv
├── plots
│   ├── RnaPsiPhosAbundance_RnaTPMgt1_CLK1_scatter.pdf - scatterplot comparing CLK1 Psi and phosphorylation abundance for splicing event in exon 4 on S182
│   ├── RnaTPMAbundancePhosAbundance_RnaTPMgt1_CLK1_scatter.pdf - scatterplot comparing CLK1 TPM and phosphorylation abundance for splicing event in exon 4 on S182
│   ├── RnaPsiPhosAbundance_RnaTPMgt20_CLK1_scatter.pdf - scatterplot comparing CLK1 Psi and phosphorylation abundance for splicing event in exon 4 on S182, filtered to include only TPM >= 20
│   ├── RnaTPMAbundancePhosAbundance_RnaTPMgt20_CLK1_scatter.pdf - scatterplot comparing CLK1 TPM and phosphorylation abundance for splicing event in exon 4 on S182, filtered to include only TMP >= 20
```

