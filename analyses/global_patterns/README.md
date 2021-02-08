# Global patterns in aberrant splicing

Module authors: Ammar Naqvi (@naqvia)

Generate tables and plots to assess global patterns of aberrant splicing, including splicing types, splicing index, and splicing heterogeneity.


## Usage
### Quantify splicing cases:
<br>**Quantify splicing cases (SE, MXE, RI, A533, A3SS)**
```
./quant_splicing.pl <input>
```

### Compute splicing index clustering:
<br>**Genrate splicing index plot for DMGs, similar to TMB**

```
/generate_splicing_index_DMG_tab.pl <pbta-histologies.RNA-Seq.initial.tsv> <filtered_samples_files.txt>
```
```
Rscript splicing_index.R
```

Input files:
```
data/pbta-histologies.tsv
data/filtered_samples_files.txt
results/splicing_index_DMG.wdPSI.txt
```

Output files:
```
plots/splicing_index_cdf_DMG.pdf
```

![](plots/splicing_index_cdf_current.png)

### Compute aberrant splicing heterogeneity:
<br>**Genrate barplot showing aberrant splicing heterogeneity**

```
/generate_splicing_index_DMG_tab.pl <pbta-histologies.RNA-Seq.initial.tsv> <filtered_samples_files.txt>
```
```
Rscript splicing_index.R
```

Input files:
```
data/pbta-histologies.tsv
data/filtered_samples_files.txt
results/splicing_index_DMG.wdPSI.txt
```

Output files:
```
plots/splicing_index_cdf_DMG.pdf
```

![](plots/splicing_index_cdf_current.png)



## Folder content

* `quant_splicing.pl` quantifies signficant splicing cases from rMATS output `results/splicing_cases.txt`
* `create_matrix_of_PSI_HGATs.pl` prepares, filters, and constructs data table with healthy comparisons of dPSI for downstream , outputs file to `results/hgat.diffsplicing.psi.txt`
* `Rscript splicing_index.R` takes the files from above and generates CDF plot in `plots/splicing_index.png`
