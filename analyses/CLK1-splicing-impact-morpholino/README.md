# CLK1 Exon 4 Splicing Impact

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to identify differential splicing and expression between CLK1 vs CLK1 morphilino treated KNS42 cell lines. The morpholino treated cell lines forces CLK1 exon 4 incluson, thereby increasing CLK1 protein levels. The results of this analysis will allow us to identify downstream impact on splicing and gene expression of CLK1 splicing dysregulation or CLK1 targets. The overall goal will help us determine the impact of CLK1 aberrant splicing on tumor progression and oncogenesis.

## Usage
### Run script
```
bash run_module.sh
```

### Input
```
data/ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv
../../splicing_events_functional_sites/input/unip*bed

input
├── RBP_known.txt
├── base_excision_repair.txt
├── dna_repair_all.txt
├── epi_known.txt
├── homologous_recombination.txt
├── mismatch_repair.txt
├── morpholno.merged.rmats.tsv
├── nonhomologous_end_joining.txt
├── nucleotide_excision_repair.txt
└── oncoprint-goi-lists-OpenPedCan-gencode-v39.csv

```

### Output
```
results/ctrl_vs_morpho.rsem.genes.collapsed.rds
results/ctrl_vs_treated.de.formatted.tsv
results/ctrl_vs_treated.de.tsv
results/expr_collapsed_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
results/expr_collapsed_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
results/expr_collapsed_clk1_ctrl_morpho_kegg_gsva_scores.tsv
results/expr_collapsed_de_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
results/expr_collapsed_de_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
results/expr_collapsed_de_clk1_ctrl_morpho_kegg_gsva_scores.tsv
results/gsva_score_diff_dna_repair.tsv
results/gsva_score_diff_dna_repair_de.tsv
results/gsva_score_diff_hallmark.tsv
results/gsva_score_diff_hallmark_de.tsv
results/gsva_score_diff_kegg.tsv
results/gsva_score_diff_kegg_de.tsv
results/splicing_events.morpho.intersectUnip.ggplot.txt
plots/CLK1_targets_ora_dotplot.func-sites.pdf
plots/CLK1_targets_ora_dotplot.pdf
plots/ctrl_vs_clk1-morp_volcano.pdf
plots/dPSI-distr-func-goi.pdf
plots/dPSI-distr-func.pdf
plots/dPSI_distr.pdf
plots/gene-fam-DE-plot.pdf
plots/gsva_heatmap_dna_repair.pdf
plots/gsva_heatmap_dna_repair_de.pdf
plots/gsva_heatmap_hallmark.pdf
plots/gsva_heatmap_hallmark_de.pdf
plots/gsva_heatmap_kegg.pdf
plots/gsva_heatmap_kegg_de.pdf
```

## Folder content
* `00-get-splice-transcripts.R` identifies the transcripts which are being spliced, using the rMATs morpholino results file
* `01-diffExpr-ctrl_vs_morph.R` performs differential expression analysis on ctrl vs treated cells using DESeq2. It also subsets based on genes of interests (transcription factors, kinases, etc)
* `02-plot_diff-splice-events.R` plots PSI distributions of differential splice events
* `03-run_bedtools_intersect-morpho.sh` formats rMATs output into bed format, and runs bedtools to intersect with uniprot db files
* `04-plot_diff-func-splice-events.R` plots PSI distributions of differential splice events categorized by functionnal site and then subsets cancer genes
* `05-ora-analysis.R` performs an over-representation analysis on above splice events that hit functional sites
* `06-conduct-gsva-analysis.R` performs GSVA on HALLMARK, KEGG, and DNA repair pathways from https://pubmed.ncbi.nlm.nih.gov/29617664/
* `07-run-gsva-comparisons.Rmd` performs CLK1 morpholino vs non targeting morpholino cell line treatment comparisons of GSVA scores HALLMARK, KEGG, and DNA repair pathways from https://pubmed.ncbi.nlm.nih.gov/29617664/.
* `08-intersection-dex-des.R` intersects dysregulated genes to assess overlap between genes that are both differentially spliced and expressed. Generates Venn diagram and performs ORA of overlapping genes. 
* `09-plot_total-splicing-cases.R` plots number of splicing events per type per treatment.

## Directory structure
```
.
├── 01-diffExpr-ctrl_vs_morph.R
├── 02-plot_diff-splice-events.R
├── 03-run_bedtools_intersect-morpho.sh
├── 04-plot_diff-func-splice-events.R
├── 05-ora-analysis.R
├── 06-conduct-gsva-analysis.R
├── 07-run-gsva-comparisons.Rmd
├── 07-run-gsva-comparisons.html
├── 08-intersection-dex-des.R
├── 09-plot_total-splicing-cases.R 
├── README.md
├── input
│   ├── RBP_known.txt
│   ├── base_excision_repair.txt
│   ├── dna_repair_all.txt
│   ├── epi_known.txt
│   ├── homologous_recombination.txt
│   ├── mismatch_repair.txt
│   ├── morpholno.merged.rmats.tsv
│   ├── nonhomologous_end_joining.txt
│   ├── nucleotide_excision_repair.txt
│   └── oncoprint-goi-lists-OpenPedCan-gencode-v39.csv
├── plots
│   ├── CLK1_targets_ora_dotplot.func-sites.pdf
│   ├── CLK1_targets_ora_dotplot.pdf
│   ├── ctrl_vs_clk1-morp_volcano.pdf
│   ├── dPSI-distr-func-goi.pdf
│   ├── dPSI-distr-func.pdf
│   ├── dPSI_distr.pdf
│   ├── gene-fam-DE-plot.pdf
│   ├── gsva_heatmap_dna_repair.pdf
│   ├── gsva_heatmap_dna_repair_de.pdf
│   ├── gsva_heatmap_hallmark.pdf
│   ├── gsva_heatmap_hallmark_de.pdf
│   ├── gsva_heatmap_kegg.pdf
│   ├── gsva_heatmap_kegg_de.pdf
│   └── splice-types.pdf
├── results
│   ├── ctrl_vs_morpho.rsem.genes.collapsed.rds
│   ├── ctrl_vs_treated.de.formatted.tsv
│   ├── ctrl_vs_treated.de.tsv
│   ├── expr_collapsed_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
│   ├── expr_collapsed_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
│   ├── expr_collapsed_clk1_ctrl_morpho_kegg_gsva_scores.tsv
│   ├── expr_collapsed_de_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
│   ├── expr_collapsed_de_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
│   ├── expr_collapsed_de_clk1_ctrl_morpho_kegg_gsva_scores.tsv
│   ├── gsva_score_diff_dna_repair.tsv
│   ├── gsva_score_diff_dna_repair_de.tsv
│   ├── gsva_score_diff_hallmark.tsv
│   ├── gsva_score_diff_hallmark_de.tsv
│   ├── gsva_score_diff_kegg.tsv
│   ├── gsva_score_diff_kegg_de.tsv
│   └── splicing_events.morpho.intersectUnip.ggplot.txt
└── run_module.sh
```
