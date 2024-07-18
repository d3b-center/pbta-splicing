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

## Folder content
* `00-get-splice-transcripts.R` identifies the transcripts which are being spliced, using the rMATs morpholino results file
* `01-diffExpr-ctrl_vs_morph.R` performs differential expression analysis on ctrl vs treated cells using DESeq2. It also subsets based on genes of interests (transcription factors, kinases, etc)
* `02-plot_diff-splice-events.R` plots PSI distributions of differential splice events
* `03-run_bedtools_intersect-morpho.sh` formats rMATs output into bed format, and runs bedtools to intersect with uniprot db files
* `04-plot_diff-func-splice-events.R` plots PSI distributions of differential splice events categorized by functionnal site and then subsets cancer genes
* `05-ora-analysis.R` performs an over-representation analysis on above splice events that hit functional sites
* `06-conduct-gsva-analysis.R` performs GSVA on CLK1 morpholino and control morpholino DE and DS events. DS events are delineated by all functional events or all functional events in onco/tsgs. This module uses HALLMARK, KEGG, and DNA repair pathways from https://pubmed.ncbi.nlm.nih.gov/29617664/.
* `07-run-gsva-comparisons.Rmd` performs CLK1 morpholino vs non targeting morpholino cell line treatment comparisons of GSVA scores HALLMARK, KEGG, and DNA repair pathways from https://pubmed.ncbi.nlm.nih.gov/29617664/.
* `08-intersection-dex-des.R` intersects dysregulated genes to assess overlap between genes that are both differentially spliced and expressed. Generates Venn diagram and performs ORA of overlapping genes.
* `09-plot_total-splicing-cases.R` plots number of splicing events per type per treatment.
* `10-crispr-screen-intersection.R` plots the intersection of CLK1 targets and genes that are essential in HGGs (via CRISPR scores obtained from [CCMA](https://data.mendeley.com/datasets/rnfs539pfw/3), doi: 10.17632/rnfs539pfw.3)

## Directory structure
```
.
├── 00-get-splice-transcripts.R
├── 01-diffExpr-ctrl_vs_morph.R
├── 02-plot_diff-splice-events.R
├── 03-run_bedtools_intersect-morpho.sh
├── 04-plot_diff-func-splice-events.R
├── 05-ora-analysis.R
├── 06-conduct-gsva-analysis.R
├── 07-run-gsva-comparisons.Rmd
├── 07-run-gsva-comparisons.html
├── 08-intersection-dex-des.R
├── 09-crispr-screen-intersection.R
├── 09-plot_total-splicing-cases.R
├── 10-crispr-screen-intersection.R
├── README.md
├── input
│   ├── CCMA_crispr_genedependency_042024.csv
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
│   ├── CLK1_ds-dex-targets_ora_dotplot-func.pdf
│   ├── CLK1_ds-dex-targets_ora_dotplot.pdf
│   ├── CLK1_targets_ora_dotplot.func-sites.pdf
│   ├── CLK1_targets_ora_dotplot.pdf
│   ├── clk1-targets-crispr-cbtn-lines.pdf
│   ├── clk1-targets-crispr_cbtn_lines-sign.pdf
│   ├── ctrl_vs_clk1-morp_volcano.pdf
│   ├── dPSI-distr-func-goi.pdf
│   ├── dPSI-distr-func.pdf
│   ├── dPSI_distr.pdf
│   ├── des-dex-venn-func.pdf
│   ├── des-dex-venn.pdf
│   ├── gene-fam-DE-plot.pdf
│   ├── gsva_heatmap_dna_repair.pdf
│   ├── gsva_heatmap_dna_repair_de.pdf
│   ├── gsva_heatmap_dna_repair_sp.pdf
│   ├── gsva_heatmap_dna_repair_sp_onc.pdf
│   ├── gsva_heatmap_hallmark.pdf
│   ├── gsva_heatmap_hallmark_de.pdf
│   ├── gsva_heatmap_hallmark_sp.pdf
│   ├── gsva_heatmap_hallmark_sp_onc.pdf
│   ├── gsva_heatmap_kegg.pdf
│   ├── gsva_heatmap_kegg_de.pdf
│   ├── gsva_heatmap_kegg_sp.pdf
│   ├── gsva_heatmap_kegg_sp_onc.pdf
│   └── splice-types.pdf
├── results
│   ├── clk1-targets-crispr-cbtn-lines.txt
│   ├── common_genes_de_ds_functional.txt
│   ├── ctrl_vs_morpho.rsem.genes.collapsed.rds
│   ├── ctrl_vs_treated.de.formatted.tsv
│   ├── ctrl_vs_treated.de.tsv
│   ├── differential_splice_by_goi_category.tsv
│   ├── expr_collapsed_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
│   ├── expr_collapsed_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
│   ├── expr_collapsed_clk1_ctrl_morpho_kegg_gsva_scores.tsv
│   ├── expr_collapsed_de_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
│   ├── expr_collapsed_de_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
│   ├── expr_collapsed_de_clk1_ctrl_morpho_kegg_gsva_scores.tsv
│   ├── expr_splice_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
│   ├── expr_splice_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
│   ├── expr_splice_clk1_ctrl_morpho_kegg_gsva_scores.tsv
│   ├── expr_splice_onco_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv
│   ├── expr_splice_onco_clk1_ctrl_morpho_hallmark_gsva_scores.tsv
│   ├── expr_splice_onco_clk1_ctrl_morpho_kegg_gsva_scores.tsv
│   ├── gene_categories.tsv
│   ├── gsva_score_diff_dna_repair.tsv
│   ├── gsva_score_diff_dna_repair_de.tsv
│   ├── gsva_score_diff_dna_repair_sp.tsv
│   ├── gsva_score_diff_dna_repair_sp_onc.tsv
│   ├── gsva_score_diff_hallmark.tsv
│   ├── gsva_score_diff_hallmark_de.tsv
│   ├── gsva_score_diff_hallmark_sp.tsv
│   ├── gsva_score_diff_hallmark_sp_onc.tsv
│   ├── gsva_score_diff_kegg.tsv
│   ├── gsva_score_diff_kegg_de.tsv
│   ├── gsva_score_diff_kegg_sp.tsv
│   ├── gsva_score_diff_kegg_sp_onc.tsv
│   ├── morpholino_splice_NF1_0.1_psi_diff_transcripts.tsv
│   ├── splice-events-significant.tsv
│   ├── splicing_events.morpho.A3SS.intersectUnip.ggplot.txt
│   ├── splicing_events.morpho.A5SS.intersectUnip.ggplot.txt
│   ├── splicing_events.morpho.RI.intersectUnip.ggplot.txt
│   ├── splicing_events.morpho.SE.intersectUnip.ggplot.txt
│   └── splicing_events.morpho.intersectUnip.ggplot.txt
├── run_module.sh
```
