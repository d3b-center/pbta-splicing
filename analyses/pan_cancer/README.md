# Pan Cancer

### Work in progress

This module is a work in progress.

## Usage
### Consensus clustering:
<br>Preprocessing and create matrix of inclusion/psi levels
```
./create_matrix_of_PSI.pl pbta-histologies.tsv filtered_samples_files.txt
```

Input files:
```
data/pbta-histologies.tsv
data/filtered_samples_files.txt
```

Output files:
```
results/pan_cancer_splicing.thr10.report_all.txt
```

Run consensus clustering method and compute clusters
```
Rscript consensus_clustering.R
```

Input files:
```
tables generated from `create_matrix_of_PSI.pl` run
```

Output files:
```
plots/*png
```

### Compute splicing index clustering:
<br>Genrate splicing index plot, similar to TMB

```
/generate_splicing_index_tab.pl pbta-histologies.RNA-Seq.initial.tsv filtered_samples_files.txt
```

Input files:
```
data/pbta-histologies.tsv
data/filtered_samples_files.txt
```

Output files:
```
results/splicing_index.wdPSI.txt
```

```
Rscript splicing_index.R
```

Input files:
```
results/splicing_index.wdPSI.txt
```

Output files:
```
plots/splicing_index_cdf_current.pdf
```

### Oncoplot with mutations, fusions, and splicing for HGATs:
<br>Preprocessing and create matrix of inclusion/psi levels
```
./create_matrix_of_PSI_HGATs.pl pbta-histologies.RNA-Seq.initial.tsv filtered_samples_files.txt
```

Input files:
```
data/pbta-histologies.tsv
data/filtered_samples_files.txt
```

Output files:
```
results/hgat.diffsplicing.psi.txt
```

```
Rscript create_oncoplot_of_splicing_w_filters.R
```

Input files:
```
results/hgat.diffsplicing.psi.txt
```

Output files:
```
plots/oncoplot*pdf
```

## Folder content

* `00-map-to-sample_id.R` prepares MAF, focal CN (the output of the `focal-cn-file-preparation` module), and standardized fusion files for use with `01-plot-oncoprint.R`.
  * The `Tumor_Sample_Barcode` column in the output corresponds to the `sample_id` column in the histologies file
  * We remove ambiguous `sample_id` -- i.e., where there are more than two tumor biospecimens that map to the same sample identifier.
  * Filtering via an [independent specimen file](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples) is optional, but highly recommended.
* `01-plot-oncoprint.R` takes the files from above and optionally a gene list and creates an oncoprint.
