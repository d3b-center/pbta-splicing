# CLK1 Exon 4 Splicing Impact

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify differential splicing and expression between CLK1 vs CLK1 morphilino treated KNS42 cell lines. The morpholino treated cell lines forces CLK1 exon 4 incluson, thereby increasing CLK1 protein levels. The results of this analysis will allow us to identify downstream impact on splicing and gene expression of CLK1 splicing dysregulation or CLK1 targets. The overall goal will help us determine the impact of CLK1 aberrant splicing on tumor progression and oncogenesis.



## Usage
### Make summary table of strong splicing events and relevant filtered tables that overlap functional sites:
<br>**Run script**
```
bash run_module.sh
```

*Input files:*
```
data/ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv
input/hgnc_complete_set.txt
input/genelistreference.txt
input/oncoprint-goi-lists-OpenPedCan-gencode-v39.csv
input_dir/RBP_known.txt
input_dir/epi_known.txt
```

*Output:*
```
results/ctrl_vs_treated.de.tsv
results/ctrl_vs_treated.de.regl.tsv
plots/CLK1_morph_gsea_ridgeplot.tiff
plots/CLK1_status_gsea_dotplot.tiff
plots/CLK1_targets_ora_dotplot.tiff
plots/ctrl_vs_clk1-morp_volcano.pdf
plots/dPSI_distr_ei.pdf
plots/dPSI_distr_es.pdf
plots/gene-fam-DE-plot.pdf
```

## Folder content
* `01-diffExpr-ctrl_vs_morph.R` performs differential expression analysis on ctrl vs treated cells using DESeq2. It also subsets based on genes of interests (transcription factors, kinases, etc)
* `02-gsea-analysis.R` performs gsea on differentially expressed genes from previous script
* `03-plot_diff-splice-events.R` plots PSI distributions of differential splice events
* `04-ora-analysis.R` performs an over-representation analysis on above splice events

## Directory structure
```
.
├── 01-diffExpr-ctrl_vs_morph.R
├── 02-gsea-analysis.R
├── 03-plot_diff-splice-events.R
├── 04-ora-analysis.R
├── CLK1_splicing_impact
├── README.md
├── input
│   ├── RBP_known.txt
│   ├── epi_known.txt
│   ├── genelistreference.txt
│   ├── hgnc_complete_set.txt
│   ├── morpholno.merged.rmats.tsv
│   └── oncoprint-goi-lists-OpenPedCan-gencode-v39.csv
├── plots
│   ├── CLK1_morph_gsea_ridgeplot.tiff
│   ├── CLK1_status_gsea_dotplot.tiff
│   ├── CLK1_targets_ora_dotplot.tiff
│   ├── ctrl_vs_clk1-morp_volcano.pdf
│   ├── dPSI_distr_ei.pdf
│   ├── dPSI_distr_es.pdf
│   └── gene-fam-DE-plot.pdf
├── results
│   ├── clk1_morph_gsea_results.csv
│   ├── ctrl_vs_treated.de.formatted.tsv
│   ├── ctrl_vs_treated.de.regl.tsv
│   └── ctrl_vs_treated.de.tsv
└── run_module.sh
```
