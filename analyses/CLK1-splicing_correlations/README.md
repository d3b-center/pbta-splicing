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
./data/gene-counts-rsem-expected_count-collapsed.rds
./analyses/CLK1-splicing_correlations/input/CPTAC3-pbt.xls
../splicing_index/results/splicing_index.SE.txt
../splicing-factor_dysregulation/input/splicing_factors.txt
```

## Folder content
* `run_module.sh` shell script to run analysis
* `00-subset-CLK1-from-rMATs.R` subset rMATs output for CLK1
* `01-plot_highExon4_vs_lowExon4_and_SBI.R` correlate high vs low levels of overall splicing burden with CLK1 exon 4 inclusion levels
* `02-plot_splicing_vs_expr.R` correlate CLK1 exon 4 inclusion levels with RNA expression, including CLK1, SRSF1, SRSF2, SRSF10
* `03-plot_CLK1-Ex4-splicing_vs_SRSF1-expr.R` correlate CLK1 RNA expression with RNA SRSFs, phospho, and total proteomics across brain tumor types
* `04-CLK1_PSI_plots.R` generates stacked barplot of relative inclusion/skipping of exon 4 across midline HGG tumors with differential splicing in CLK1.
* `04-plot-diffExp_highlowSBI.R` perform differential gene expression on high vs low SBI HGG tumors
* `05-CLK1-NF1-NMD.Rmd` correlate and ivnestigate NF1 NMD transcripts splicing found to be differentially spliced in morpholino in our patient samples


## Directory structure
```
.
├── 00-subset-CLK1-from-rMATs.R
├── 01-plot_highExon4_vs_lowExon4_and_SBI.R
├── 02-plot_splicing_vs_expr.R
├── 03-plot_SR-phosp_vs_CLK1-RNA.R
├── 04-CLK1_PSI_plots.R
├── README.md
├── input
│   └── CPTAC3-pbt.xls
├── plots
│   ├── CLK1_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── CLK1_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF10_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF1_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_all_hgg.pdf
│   ├── SRSF2_exp_vs_CLK1_psi_midline_hgg.pdf
│   ├── SR_phos_CLK1_exp_heatmap.pdf
│   ├── all_hgg_CLK1_exon4_inclusion_fraction_hgg_stacked.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_exome_capture.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_poly-A stranded.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_poly-A.pdf
│   ├── all_hgg_SBI_high_vs_low_CLK1_stranded.pdf
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
│   ├── dmg-mean_clk1_psi.txt
│   ├── mean_clk1_psi.txt
│   └── other_hgg-mean_clk1_psi.txt
├── run_module.sh
└── util
    └── function-create-scatter-plot.R
```
