# OpenPBTA Oncoprint

## Work in progress

This module is a work in progress.
The code was written to be updated as consensus, filtered, and/or prioritized data becomes available.

**The plots included as PNGs should be regarded as proof of concept, rather than for interpretation.**



## Usage

To run this module from the command line as intended, use:

Preprocessing and create matrix of inclusion/psi levels
```
./create_matrix_of_PSI.pl pbta-histologies.tsv filtered_samples_files.txt
```
input:  histology file, file paths for rMATS output
output: matrix of psi/inclusion levels for each splice change in each sample

Run consensus clustering method and compute clusters
```
Rscript consensus_clustering.R
```
input:  tables generated from `create_matrix_of_PSI.pl` run
output: clustering plots

Genrate splicing index plot, similar to TMB

```
./create_matrix_of_PSI.pl pbta-histologies.tsv filtered_samples_files.txt
/generate_splicing_index_tab.pl pbta-histologies.tsv filtered_samples_files.txt
Rscript consensus_clustering.R
```
input:  histology file, file paths for rMATS output, results/splicing_index.wdPSI.txt
output: matrix of psi/inclusion levels for each splice change in each sample

## Folder content

* `00-map-to-sample_id.R` prepares MAF, focal CN (the output of the `focal-cn-file-preparation` module), and standardized fusion files for use with `01-plot-oncoprint.R`.
  * The `Tumor_Sample_Barcode` column in the output corresponds to the `sample_id` column in the histologies file
  * We remove ambiguous `sample_id` -- i.e., where there are more than two tumor biospecimens that map to the same sample identifier.
  * Filtering via an [independent specimen file](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples) is optional, but highly recommended.
* `01-plot-oncoprint.R` takes the files from above and optionally a gene list and creates an oncoprint.
