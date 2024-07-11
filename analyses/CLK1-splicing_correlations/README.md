# CLK1 Splicing Correlations

Module authors: Ammar Naqvi (@naqvia), Jo Lynne Rokita (@jharenza)

The purpose of this module is to correlate CLK1 Exon 4 splicing with splicing
burden, RNA expression and proteomics

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run_module.sh
```
Input files:
```
./data/histologies.tsv
./data/splice-events-rmats.tsv.gz
./data/clk1-splice-events-rmats.tsv.gz
./data/gene-counts-rsem-expected_count-collapsed.rds
./analyses/CLK1-splicing_correlations/input/CPTAC3-pbt.xls
../splicing_index/results/splicing_index.SE.txt
../splicing-factor_dysregulation/input/splicing_factors.txt
./data/cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
./data/cptac-protein-imputed-prot-expression-abundance.tsv.gz
./data/hope-protein-imputed-phospho-expression-abundance.tsv.gz
./data/hope-protein-imputed-prot-expression-abundance.tsv.gz
```

## Folder content
* `run_module.sh` shell script to run analysis
* `00-subset-CLK1-from-rMATs.R` subset rMATs output for CLK1
* `01-plot_highExon4_vs_lowExon4_and_SBI.R` correlate high vs low levels of overall splicing burden with CLK1 exon 4 inclusion levels
* `02-plot_splicing_vs_expr.R` correlate CLK1 exon 4 inclusion levels with RNA expression, including CLK1, SRSF1, SRSF2, SRSF10
* `03-plot_CLK1-Ex4-splicing_vs_SRSF1-expr.R` correlate CLK1 RNA expression with RNA SRSFs, phospho, and total proteomics across brain tumor types
* `04-CLK1_PSI_plots.R` generates stacked barplot of relative inclusion/skipping of exon 4 across midline HGG tumors with differential splicing in CLK1.
* `05-CLK-SRSF-expr-correlations.R` generates scatter plots of CLK and SRSF transcript abundance vs. CLK1 exon 4 PSI and transcript abundance
* `06-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R` generates summary heatmap of spearman correlation coefficients between CLK1 PSI and transcript abundance against CLK and SRSF RNA, total protein, and phosphoprotein expression
* `07-plot-clk1ex4-hgg-normals.R` plots CLK1 exon 4 transcript expression for HGGs, GTEX normals, fetal brain, and TGEN pediatric normals.
<<<<<<< HEAD
=======
* `08-CLK1-impact-NF1-splicing.Rmd` plots CLK1 and NF1 psi, transcript expression, and gene expression correlations for each cancer group and as a whole.
* `09-clk1-nf1-protein-correlations.R` plots NF1 total protein correlations and NF1 phosphosite correlations to all data from 08.
* `10-clk1-nf1-single-sample-heatmap.R` plots single DMG sample heatmap of CLK1/NF1 RNA, splice, protein z-scores.
>>>>>>> main

## Directory structure
```
.
├── 00-subset-CLK1-from-rMATs.R
├── 01-plot_highExon4_vs_lowExon4_and_SBI.R
├── 02-plot_splicing_vs_expr.R
├── 03-plot_SR-phosp_vs_CLK1-RNA.R
├── 04-CLK1_PSI_plots.R
├── 05-CLK-SRSF-expr-correlations.R
├── 06-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R
├── 07-plot-clk1ex4-hgg-normals.R
<<<<<<< HEAD
=======
├── 08-CLK1-impact-NF1-splicing.Rmd
├── 08-CLK1-impact-NF1-splicing.html
├── 09-clk1-nf1-protein-correlations.R
├── 10-clk1-nf1-single-sample-heatmap.R
>>>>>>> main
├── README.md
├── input
│   └── CPTAC3-pbt.xls
├── plots
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK1-201_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
<<<<<<< HEAD
=======
│   ├── CLK1-NF1-all-cors.pdf
│   ├── CLK1-NF1-cors-Atypical Teratoid Rhabdoid Tumor.pdf
│   ├── CLK1-NF1-cors-Choroid plexus tumor.pdf
│   ├── CLK1-NF1-cors-Craniopharyngioma.pdf
│   ├── CLK1-NF1-cors-DIPG or DMG.pdf
│   ├── CLK1-NF1-cors-Ependymoma.pdf
│   ├── CLK1-NF1-cors-Germ cell tumor.pdf
│   ├── CLK1-NF1-cors-High-grade glioma.pdf
│   ├── CLK1-NF1-cors-Low-grade glioma.pdf
│   ├── CLK1-NF1-cors-Medulloblastoma.pdf
│   ├── CLK1-NF1-cors-Meningioma.pdf
│   ├── CLK1-NF1-cors-Mesenchymal tumor.pdf
│   ├── CLK1-NF1-cors-Mixed neuronal-glial tumor.pdf
│   ├── CLK1-NF1-cors-Neurofibroma plexiform.pdf
│   ├── CLK1-NF1-cors-Non-neoplastic tumor.pdf
│   ├── CLK1-NF1-cors-Other CNS embryonal tumor.pdf
│   ├── CLK1-NF1-cors-Other high-grade glioma.pdf
│   ├── CLK1-NF1-cors-Other tumor.pdf
│   ├── CLK1-NF1-cors-Schwannoma.pdf
│   ├── CLK1-NF1-psi-expr-correlation-heatmap.pdf
│   ├── CLK1-NF1-subset-cors.pdf
│   ├── CLK1-exon4-NF1-psi-all-cns-boxplot.pdf
│   ├── CLK1-psi-all-cns-boxplot.pdf
>>>>>>> main
│   ├── CLK1-psi-expr-correlation-heatmap.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_201_exp_DMG.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_201_exp_HGG.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_exp_DMG.pdf
│   ├── CLK1_SRSF_phospho_vs_CLK1_exp_HGG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_201_exp_DMG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_201_exp_HGG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_exp_DMG.pdf
│   ├── CLK1_SRSF_prot_vs_CLK1_exp_HGG.pdf
│   ├── CLK1_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── CLK1_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_exp_all_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_exp_midline_hgg.pdf
│   ├── CLK1_exp_vs_SRSF_exp_other_hgg.pdf
<<<<<<< HEAD
=======
│   ├── CLK1ex4-hgg-gtex.pdf
>>>>>>> main
│   ├── CLK2_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK2_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK2_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK3_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK3_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK3_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK4_exp_vs_SRSF_SRPK_exp_all_hgg.pdf
│   ├── CLK4_exp_vs_SRSF_SRPK_exp_midline_hgg.pdf
│   ├── CLK4_exp_vs_SRSF_SRPK_exp_other_hgg.pdf
│   ├── CLK_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── CLK_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── CLK_exp_vs_CLK1_psi_other_hgg.pdf
<<<<<<< HEAD
│   ├── CLK1ex4-hgg-normals.pdf
=======
│   ├── NF1 pS2796-phos-cor-CLK1-201 Exon4 PSI-DMG.pdf
│   ├── NF1 pS2796-phos-cor-CLK1-201 Exon4 PSI-HGG.pdf
│   ├── NF1 pS2796-phos-cor-CLK1-201 Log2 TPM-DMG.pdf
│   ├── NF1 pS2796-phos-cor-CLK1-201 Log2 TPM-HGG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-202 Exon23a PSI-DMG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-202 Exon23a PSI-HGG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-202 Log2 TPM-DMG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-202 Log2 TPM-HGG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-215 Log2 TPM-DMG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-215 Log2 TPM-HGG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-215 PSI-DMG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-215 PSI-HGG.pdf
│   ├── NF1 pS2796-phos-cor-NF1-DMG.pdf
│   ├── NF1 pS2796-phos-cor-Total CLK1 Log2 TPM-DMG.pdf
│   ├── NF1 pS2796-phos-cor-Total CLK1 Log2 TPM-HGG.pdf
│   ├── NF1 pS2796-phos-cor-Total NF1 Log2 TPM-DMG.pdf
│   ├── NF1 pS2796-phos-cor-Total NF1 Log2 TPM-HGG.pdf
│   ├── NF1 pS864-phos-cor-CLK1-201 Exon4 PSI-DMG.pdf
│   ├── NF1 pS864-phos-cor-CLK1-201 Exon4 PSI-HGG.pdf
│   ├── NF1 pS864-phos-cor-CLK1-201 Log2 TPM-DMG.pdf
│   ├── NF1 pS864-phos-cor-CLK1-201 Log2 TPM-HGG.pdf
│   ├── NF1 pS864-phos-cor-NF1-202 Exon23a PSI-DMG.pdf
│   ├── NF1 pS864-phos-cor-NF1-202 Exon23a PSI-HGG.pdf
│   ├── NF1 pS864-phos-cor-NF1-202 Log2 TPM-DMG.pdf
│   ├── NF1 pS864-phos-cor-NF1-202 Log2 TPM-HGG.pdf
│   ├── NF1 pS864-phos-cor-NF1-215 Log2 TPM-DMG.pdf
│   ├── NF1 pS864-phos-cor-NF1-215 Log2 TPM-HGG.pdf
│   ├── NF1 pS864-phos-cor-NF1-215 PSI-DMG.pdf
│   ├── NF1 pS864-phos-cor-NF1-215 PSI-HGG.pdf
│   ├── NF1 pS864-phos-cor-Total CLK1 Log2 TPM-DMG.pdf
│   ├── NF1 pS864-phos-cor-Total CLK1 Log2 TPM-HGG.pdf
│   ├── NF1 pS864-phos-cor-Total NF1 Log2 TPM-DMG.pdf
│   ├── NF1 pS864-phos-cor-Total NF1 Log2 TPM-HGG.pdf
│   ├── NF1-protein-cor-CLK1-201 Exon4 PSI-DMG.pdf
│   ├── NF1-protein-cor-CLK1-201 Exon4 PSI-HGG.pdf
│   ├── NF1-protein-cor-CLK1-201 Log2 TPM-DMG.pdf
│   ├── NF1-protein-cor-CLK1-201 Log2 TPM-HGG.pdf
│   ├── NF1-protein-cor-CLK1-201-DMG.pdf
│   ├── NF1-protein-cor-NF1-202 Exon23a PSI-DMG.pdf
│   ├── NF1-protein-cor-NF1-202 Exon23a PSI-HGG.pdf
│   ├── NF1-protein-cor-NF1-202 Log2 TPM-DMG.pdf
│   ├── NF1-protein-cor-NF1-202 Log2 TPM-HGG.pdf
│   ├── NF1-protein-cor-NF1-215 Log2 TPM-DMG.pdf
│   ├── NF1-protein-cor-NF1-215 Log2 TPM-HGG.pdf
│   ├── NF1-protein-cor-NF1-215 PSI-DMG.pdf
│   ├── NF1-protein-cor-NF1-215 PSI-HGG.pdf
│   ├── NF1-protein-cor-NF1-DMG.pdf
│   ├── NF1-protein-cor-NF1-HGG.pdf
│   ├── NF1-protein-cor-Total CLK1 Log2 TPM-DMG.pdf
│   ├── NF1-protein-cor-Total CLK1 Log2 TPM-HGG.pdf
│   ├── NF1-protein-cor-Total NF1 Log2 TPM-DMG.pdf
│   ├── NF1-protein-cor-Total NF1 Log2 TPM-HGG.pdf
>>>>>>> main
│   ├── SRPK_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRPK_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRPK_exp_vs_CLK1_psi_other_hgg.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF_exp_vs_CLK1_psi_other_hgg.pdf
│   ├── SR_phos_CLK1_exp_heatmap.pdf
│   ├── all_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_exome_capture.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_poly-A stranded.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_poly-A.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_stranded.pdf
<<<<<<< HEAD
=======
│   ├── clk1-nf1-single-sample-exp-protein-heatmap-dmg.pdf
>>>>>>> main
│   ├── dmg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf
│   ├── dmg_SBI_high_vs_low_CLK1_exome_capture.pdf
│   ├── dmg_SBI_high_vs_low_CLK1_poly-A stranded.pdf
│   ├── dmg_SBI_high_vs_low_CLK1_poly-A.pdf
│   ├── dmg_SBI_high_vs_low_CLK1_stranded.pdf
│   ├── other_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf
│   ├── other_hgg_SBI_high_vs_low_CLK1_exome_capture.pdf
│   ├── other_hgg_SBI_high_vs_low_CLK1_poly-A stranded.pdf
│   ├── other_hgg_SBI_high_vs_low_CLK1_poly-A.pdf
│   └── other_hgg_SBI_high_vs_low_CLK1_stranded.pdf
├── results
│   ├── all_hgg-mean_clk1_psi.txt
│   ├── clk1-exon4-psi-hgg.tsv
<<<<<<< HEAD
│   ├── clk1-splice-events-rmats.tsv
│   ├── dmg-mean_clk1_psi.txt
│   ├── hgg-dmg-clk-srsf-expression-phosphorylation.tsv
│   ├── mean_clk1_psi.txt
=======
│   ├── clk1-nf1-psi-exp-df.rds
│   ├── clk1-splice-events-rmats.tsv
│   ├── dmg-mean_clk1_psi.txt
│   ├── hgg-dmg-clk-nf1-expression-phosphorylation.tsv
│   ├── hgg-dmg-clk-srsf-expression-phosphorylation.tsv
│   ├── mean_clk1_psi.txt
│   ├── nf1-splice-events-rmats.tsv
>>>>>>> main
│   └── other_hgg-mean_clk1_psi.txt
├── run_module.sh
└── util
    └── function-create-scatter-plot.R
<<<<<<< HEAD
```
=======
```
>>>>>>> main
