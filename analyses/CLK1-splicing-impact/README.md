# CLK1 Exon 4 Splicing Impact

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify differential splicing and expression between CLK1 vs CLK1 morphilino treated KNS42 cell lines. The morpholino treated cell lines forces CLK1 exon 4 incluson, thereby increasing CLK1 protein levels. The results of this analysis will allow us to identify downstream impact on splicing and gene expression of CLK1 splicing dysregulation or CLK1 targets. The overall goal will help us determine the impact of CLK1 aberrant splicing on tumor progression and oncogenesis.



## Usage
### Make summary table of strong splicing events and relevant filtered tables that overlap functional sites:
<br>**Run script to make tables and subsequent plots below**
```
Rscript --vanilla 01-diffExpr-ctrl_vs_morph.R
```

*Input files:*
```
data/ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv
input/hgnc_complete_set.txt
input/genelistreference.txt
input/oncoprint-goi-lists-OpenPedCan-gencode-v39.csv
input_dir/RBP_known.txt
input_dir/kinase_known.txt
input_dir/epi_known.txt
```

*Output:*
```
results/ctrl_vs_treated.de.tsv
results/ctrl_vs_treated.de.regl.tsv
plots/ctrl_vs_clk1-morp_volcano.pdf
```

## Folder content
* `01-diffExpr-ctrl_vs_morph.R` performs differential expression analysis on ctrl vs treated cells using DESeq2. It also subsets based on genes of interests (transcription factors, kinases, etc)
